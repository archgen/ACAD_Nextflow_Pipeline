#!/usr/bin/env nextflow

// Show help message
if (params.help) {    
    helpMessage()
    exit 0
}

// Check validity of input TSV - custom function defined at bottom of script
checkTSV(params.sampleSheet, 8)

// Check and import fastq files - custom function efined at bottom of script
// This will make a list object, where each element is a file and each file has values for each column
// of the TSV file
fq = getFastq(params.sampleSheet)

// Split fastq files into separate input channels - each process needs its own channel (can't reuse same channel)
( fq, ch_fastp ) = fq.into(2)

// Parsing 'genome.config' into a key-value list (really a dict)
referenceMap = defineReferenceMap()

println referenceMap
/*
STEP 1: QC the data with FastP
*/

process runFastp {

    // This is what the job name will appear as in the terminal 
    tag { sampleName + ' - Fastp' }

    // Directory that results will be copied to - can change to 'move' if you want to clean up the working dir
    publishDir "${params.outDir}/FastP/${sampleName}_${experimentName}", mode: 'copy'

    // Specifying the fastq channel as a set (each set value is a per-sample column value from TSV)
    input:
    set experimentName, 
    sampleName, 
    libraryName, 
    unitName, 
    platformName, 
    runName, 
    file(reads) from ch_fastp

    // Output to copy to the publishDir
    // Anything not specified here will be created but remain in the working directory (not copied)
    output:
    file "${sampleName}_${experimentName}_${libraryName}_${runName}_fastp.{html,json}"
    file "${sampleName}_${experimentName}_${libraryName}_${runName}_fastp_{R1,R2}.fastq.gz"
    set experimentName, 
    sampleName, 
    libraryName, 
    unitName, 
    platformName, 
    runName,
    file("${sampleName}_${experimentName}_${libraryName}_${runName}_collapsed_fastp.fastq.gz") into result_FastP

    // Code to execute
    script:
    """
    fastp \
    --thread ${task.cpus} \
    -g -x -y -p -V -c \
    --in1 ${reads[0]} \
    --in2 ${reads[1]} \
    -R '${sampleName} – ${experimentName}' \
    -h ${sampleName}_${experimentName}_${libraryName}_${runName}_fastp.html \
    -j ${sampleName}_${experimentName}_${libraryName}_${runName}_fastp.json \
    --merge \
    --merged_out ${sampleName}_${experimentName}_${libraryName}_${runName}_collapsed_fastp.fastq.gz \
    --out1 ${sampleName}_${experimentName}_${libraryName}_${runName}_fastp_R1.fastq.gz \
    --out2 ${sampleName}_${experimentName}_${libraryName}_${runName}_fastp_R2.fastq.gz
    """
}

/*
STEP 2: Align the data 
*/

process runBwaAln {

    tag { sampleName + ' - BWA-aln' }

    publishDir "${params.outDir}/BWA-aln/${sampleName}_${experimentName}", mode: 'copy'

    input:
    set experimentName, 
    sampleName, 
    libraryName, 
    unitName, 
    platformName, 
    runName, 
    file(fastqCollapsed) from result_FastP
    val ref_fasta_basename from params.genome
    file ref_fasta from referenceMap.genomeFile
    file ref_idx from Channel.fromPath( [ params.genomeFile, '.*' ].join() ).collect() // Copies index files to wd

    output:
    set experimentName, 
        sampleName, 
        libraryName, 
        unitName, 
        platformName, 
        runName, 
        file("${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted.bam"), 
        file("${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted.bam.bai") into results_bwa

    // Loading Phoenix modules - separate with colon
    module 'BWA/0.7.15-foss-2017a:SAMtools/1.9-foss-2016b'

    script:
    """
    bwa aln \
    -t 8 \
    ${ref_fasta} \
    ${fastqCollapsed} \
    -n 0.04 \
    -l 1024 \
    -k 2 \
    -f ${sampleName}_${experimentName}_${libraryName}_${runName}.sai

    bwa samse \
    -r "@RG\tID:${runName}\tPL:${platformName}\tPU:${unitName}\tSM:${sampleName}" \
    ${ref_fasta} \
    ${sampleName}_${experimentName}_${libraryName}_${runName}.sai \
    ${fastqCollapsed} | \
    samblaster | \
    samtools sort \
        -@ ${task.cpus} \
        -O BAM  \
        -o ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted.bam \
        -

    samtools index \
    ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted.bam
"""
}

// Split channel into 3 for IdxStats, PreSeq and DeDup
( input_SamToolsStats, input_preSeq, input_DeDup ) = results_bwa.into(3)

/*
STEP 3: De-duplication
*/

process runDeDup {

    tag { sampleName + ' - deDup' }

    publishDir "${params.outDir}/DeDup/${sampleName}_${experimentName}", mode: 'copy'

    input:
    set experimentName,
    sampleName,
    libraryName,
    unitName,
    platformName,
    runName, 
    file(bam),
    file(bamBai) from input_DeDup
    val deDupJar from params.deDup
    val java from params.java
    val ref_fasta_basename from params.genome
    
    output:
    set experimentName,
        sampleName,
        libraryName,
        unitName,
        platformName,
        runName,
        file("${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_rmdup_sorted.bam"),
        file("${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_rmdup_sorted.bam.bai") into results_DeDup
        file("${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted.hist")
        file("${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted.dedup.json")
        file("${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted.log")
    // Loading Phoenix modules - separate with colon
    module 'Java/1.8.0_121:SAMtools/1.9-foss-2016b'

    script:
    """
    ${java} -jar ${deDupJar} -i ${bam} -m -u -o . 
    samtools sort \
        -@ ${task.cpus} \
        -O BAM \
        -o ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_rmdup_sorted.bam \
        ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup.bam

    samtools index \
    ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_rmdup_sorted.bam
"""
}

// Split channel into 2 for MitoMD and IndelRealignment
( input_MitoRetrieve, input_indelRealign ) = results_DeDup.into(2)

/*
STEP 4a: Mito-Retrieve
*/

process runMitoRetrieve {

    tag { sampleName + ' - mitoRetrieve' }

    publishDir "${params.outDir}/mitoGenomes/${sampleName}_${experimentName}", mode: 'copy'

    input:
    set experimentName,
    sampleName,
    libraryName,
    unitName,
    platformName,
    runName,
    file(bam),
    file(bamBai) from input_MitoRetrieve
    val ref_fasta_basename from params.genome
    file ref_fasta from referenceMap.genomeFile
    file ref_idx from Channel.fromPath( [ params.genomeFile, '.*' ].join() ).collect() // Copies index files to wd

    output:
    set experimentName,
        sampleName,
        libraryName,
        unitName,
        platformName,
        runName,
        file("${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_MT.bam"),
        file("${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_MT.bam.bai") into results_MitoRetrieve

    // Loading Phoenix modules - separate with colon
    module 'SAMtools/1.9-foss-2016b'

    script:
    """
    samtools view -b ${bam} MT > ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_MT.bam 

    samtools index \
    ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_MT.bam
"""
}

( input_mitoRSRS ) = results_MitoRetrieve.into(1)

/*
STEP 4b: Remap mito to RSRS 
*/

process runMitoBwaAlnRSRS {

    tag { sampleName + ' - MitoBWA-aln-RSRS' }

    publishDir "${params.outDir}/mitoGenomes", mode: 'copy'

    input:
    set experimentName, 
    sampleName, 
    libraryName, 
    unitName, 
    platformName, 
    runName, 
    file(bam),
    file(bamBai) from input_mitoRSRS
    val ref_fasta_basename from params.genome
    file ref_fasta from referenceMap.genomeFile
    file ref_idx from Channel.fromPath( [ params.genomeFile, '.*' ].join() ).collect() // Copies index files to wd
    

    output:
    set experimentName, 
        sampleName, 
        libraryName, 
        unitName, 
        platformName, 
        runName, 
        file("${sampleName}_${experimentName}_${libraryName}_${runName}_RSRS_collapsed_bwaALN_sorted_rmdup_MT_MD.bam"), 
        file("${sampleName}_${experimentName}_${libraryName}_${runName}_RSRS_collapsed_bwaALN_sorted_rmdup_MT_MD.bam.bai") into results_mitoRSRS

    // Loading Phoenix modules - separate with colon
    module 'BWA/0.7.15-foss-2017a:SAMtools/1.9-foss-2016b'

    script:
    """
    bwa aln \
    -t 8 \
    /home/a1222106/fastdir/a1222106/Refs/RSRS_mtDNA/REFHumanRSRS.fasta \
    -b ${bam} \
    -n 0.04 \
    -l 1024 \
    -k 2 \
    -f ${sampleName}_${experimentName}_${libraryName}_${runName}.sai

    bwa samse \
    -r "@RG\tID:${runName}\tPL:${platformName}\tPU:${unitName}\tSM:${sampleName}" \
    /home/a1222106/fastdir/a1222106/Refs/RSRS_mtDNA/REFHumanRSRS.fasta \
    ${sampleName}_${experimentName}_${libraryName}_${runName}.sai \
    ${bam} \
    | \
    samtools sort \
        -@ ${task.cpus} \
        -O BAM  \
        -o ${sampleName}_${experimentName}_${libraryName}_${runName}_RSRS_collapsed_bwaALN_sorted_rmdup_MT.bam \
        -

        samtools index \
    ${sampleName}_${experimentName}_${libraryName}_${runName}_RSRS_collapsed_bwaALN_sorted_rmdup_MT.bam

    samtools calmd -b ${sampleName}_${experimentName}_${libraryName}_${runName}_RSRS_collapsed_bwaALN_sorted_rmdup_MT.bam /home/a1222106/fastdir/a1222106/Refs/RSRS_mtDNA/REFHumanRSRS.fasta > ${sampleName}_${experimentName}_${libraryName}_${runName}_RSRS_collapsed_bwaALN_sorted_rmdup_MT_MD.bam
    samtools index \
    ${sampleName}_${experimentName}_${libraryName}_${runName}_RSRS_collapsed_bwaALN_sorted_rmdup_MT_MD.bam 
    
    rm -r ${params.outDir}/mitoGenomes/${sampleName}_${experimentName}
"""
}

/*
STEP 5: Statistics on the BAM file
*/

process runSamtoolsIdxStats {
    tag { sampleName + ' - Samtools stats'}

    publishDir "${params.outDir}/Samtools-Stats/${sampleName}_${experimentName}", mode: 'copy'

    input:
    set experimentName, 
        sampleName, 
        libraryName, 
        unitName, 
        platformName, 
        runName, 
        file(bam), 
        file(bamBai) from input_SamToolsStats
    val ref_fasta_basename from params.genome

    output:
    file ("${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted.idxStats") // Not needed for downstream == no channel
    file ("${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted.samStats") // Not needed for downstream == no channel
    file ("${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted.covHist.txt") // as above
    //file '*.png'  

    // Loading Phoenix modules - separate with colon
    module 'gnuplot/5.0.3-foss-2016b'

    script:
    """
    samtools idxstats ${bam} > ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted.idxStats
    samtools stats ${bam} > ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted.samStats
    samtools coverage ${bam} -m > ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted.covHist.txt
    """
}

/*
STEP 6: Preseq on the BAM file   
*/

process runPreseq {

    tag { sampleName + ' - Preseq'}

    publishDir "${params.outDir}/PreSeq/${sampleName}_${experimentName}", mode: 'copy'

    input:
    set experimentName, 
        sampleName, 
        libraryName, 
        unitName, 
        platformName, 
        runName,
        file(bam),
        file(bamBai) from input_preSeq
    val ref_fasta_basename from params.genome

    output:
    file '*.txt' // Doesn't need its own output channel --- .{ComplexityCurve,YieldCurve,CoverageCurve}

    script:
    """
    preseq c_curve -B -o ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted.ComplexityCurve.txt ${bam}
    preseq lc_extrap -B -o ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted.YieldCurve.txt ${bam}
    """
}

/*
STEP 7: Indel Realignment
*/

process runIndelRealignment {

    // This is what the job name will appear as in the terminal
    tag { sampleName + ' - IndelRealignment'}

    // Directory where the results will be copied to 
    publishDir "${params.outDir}/IndelRealign/${sampleName}_${experimentName}", mode: 'copy'

    // List and definition of all input variables for the IndelRealignment process
    input:
    set experimentName, 
    sampleName, 
    libraryName, 
    unitName, 
    platformName, 
    runName, 
    file(bam), 
    file(bamBai) from input_indelRealign
    val ref_fasta_basename from params.genome
    file ref_fasta from referenceMap.genomeFile
    file ref_idx from referenceMap.genomeIndex //Channel.fromPath( [ params.genomeFile, '.*' ].join() ).collect() // Copies index files to wd
    file ref_dict from referenceMap.genomeDict
    val gatk3_jar from params.gatk3
    val java from params.java
   
    output:
    set experimentName, 
    sampleName, 
    libraryName,
    unitName,
    platformName, 
    runName, 
    file("${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal.bam"),
    file("${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal.bam.bai") into results_IndelReal  
   
    module 'GATK/4.0.0.0-Java-1.8.0_121:picard/2.2.4-Java-1.8.0_71'

    script:
    """ 
    ${java} -Xmx8G -jar ${gatk3_jar} -T RealignerTargetCreator -R ${ref_fasta} --num_threads ${task.cpus} --mismatchFraction 0.30 --maxIntervalSize 650 --allow_potentially_misencoded_quality_scores -I ${bam} -o ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal.intervals

    ${java} -Xmx8G -jar ${gatk3_jar} -T IndelRealigner -R ${ref_fasta} -model USE_READS -compress 0 --filter_bases_not_stored --allow_potentially_misencoded_quality_scores -I ${bam} -targetIntervals ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal.intervals -o ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal.bam

    cp -v ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal.bai ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal.bam.bai 
"""
}

// Split channel into 2 for MapDamage and TrimBam 
( input_MapDamage, input_TrimBam ) = results_IndelReal.into(2)

/*
STEP 8: MapDamage
*/

process runMapDamage {

    tag { sampleName + ' - MapDamage'}

    publishDir "${params.outDir}/MapDamage/${sampleName}_${experimentName}", mode: 'copy'

    input:
    set experimentName, 
    sampleName, 
    libraryName, 
    unitName, 
    platformName, 
    runName, 
    file(bam), 
    file(bamBai) from input_MapDamage 
    val ref_fasta_basename from params.genome
    file ref_fasta from referenceMap.genomeFile
    file ref_idx from referenceMap.genomeIndex
    publishDir = publishDir
 
    output:
    set experimentName,
    sampleName,
    libraryName,
    unitName,
    platformName,
    runName,
    file("${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal_Fragmisincorporation_plot.pdf") // Not needed for downstream == no channel
    file("${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal_Length_plot.pdf") 
    file("${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal_misincorporation.txt")
    file("${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal_dnacomp.txt")
    file("${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal_dnacomp_genome.csv")  
    file("${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal_5pCtoT_freq.txt")
    file("${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal_3pGtoA_freq.txt")
    file("${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal_lgdistribution.txt")

    // Installing modules

    script:
    """ 
    mapDamage -i ${bam} -r ${ref_fasta} -m 25 -d ./ --title="${sampleName}_${libraryName}_${experimentName}_${runName}" 
    
   cp -v Fragmisincorporation_plot.pdf ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal_Fragmisincorporation_plot.pdf
   cp -v Length_plot.pdf ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal_Length_plot.pdf
   cp -v misincorporation.txt ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal_misincorporation.txt
   cp -v dnacomp.txt ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal_dnacomp.txt
   cp -v dnacomp_genome.csv ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal_dnacomp_genome.csv
   cp -v 5pCtoT_freq.txt ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal_5pCtoT_freq.txt
   cp -v 3pGtoA_freq.txt ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal_3pGtoA_freq.txt
   cp -v lgdistribution.txt ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal_lgdistribution.txt
"""
}

/*
STEP 9: TrimBam
*/

process runTrimBam {

    tag { sampleName + ' - TrimBam'}

    publishDir "${params.outDir}/TrimBam", mode: 'copy'

    input:
    set experimentName,
    sampleName,
    libraryName,
    unitName,
    platformName,
    runName,
    file(bam),
    file(bamBai) from input_TrimBam 
    val ref_fasta_basename from params.genome
    file ref_fasta from referenceMap.genomeFile
    file ref_idx from referenceMap.genomeIndex

    output:
    set experimentName,
    sampleName,
    libraryName,
    unitName,
    platformName,
    runName,
    file("${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal.trim3_2ends.bam"),
    file("${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal.trim3_2ends.bam.bai") into results_TrimBam
    
    // Installing modules
    module 'SAMtools/1.9-foss-2016b'

    script:
    """
    bam trimBam \
    ${bam} \
    ${sampleName}_${experimentName}_${libraryName}_${runName}_tmp.bam \
    3 \
    --ignoreStrand \
    --clip

    samtools sort \
        -@ ${task.cpus} \
        -O BAM  \
        -o ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal.trim3_2ends.bam \
        ${sampleName}_${experimentName}_${libraryName}_${runName}_tmp.bam 

    samtools index \
    ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal.trim3_2ends.bam
    """
}


// Split channel into 2 for QualiMap and BedTools
( input_QualiMap, input_BedTools ) = results_TrimBam.into(2)

/*
STEP 10: QualiMap
*/

process runQualimap {

    tag { sampleName + ' - QualiMap'}

    publishDir "${params.outDir}/QualiMap/${sampleName}_${experimentName}", mode: 'copy'
    
    input:
    set experimentName,
    sampleName,
    libraryName,
    unitName,
    platformName,
    runName,
    file(bam),
    file(bamBai) from input_QualiMap
    val ref_fasta_basename from params.genome
    file ref_fasta from referenceMap.genomeFile
    file ref_idx from referenceMap.genomeIndex
    java_mem = 8G
    publishDir = publishDir
 
    output:
    file "${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal.trim3_2ends_qualimapReport.html" // Not needed for downstream == no channel
    file "${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal.trim3_2ends_genome_results.txt" // Not needed for downstream == no channel

    // Installing modules
    
    script:
    """
    qualimap bamqc -bam ${bam} -gff /home/a1222106/fastdir/Refs/qualimap_1240k_positions_XY.bed -c -nt ${task.cpus} --skip-duplicated --skip-dup-mode 0 --java-mem-size=12G -outdir ./
    cp -v *qualimapReport.html ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal.trim3_2ends_qualimapReport.html
    cp -v *genome_results.txt ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal.trim3_2ends_genome_results.txt
    """
} 

/*
STEP 11: BedTools
*/

process runBedTools {

    tag { sampleName + ' - BedTools'}

    publishDir "${params.outDir}/BedTools/${sampleName}_${experimentName}", mode: 'copy'

    input:
    set experimentName,
    sampleName,
    libraryName,
    unitName,
    platformName,
    runName,
    file(bam),
    file(bamBai) from input_BedTools
    val ref_fasta_basename from params.genome
    file ref_fasta from referenceMap.genomeFile
    file ref_idx from referenceMap.genomeIndex
    publishDir = publishDir

    output:
    file "*"

    script:
    """
    bedtools coverage -a /home/a1222106/fastdir/Refs/qualimap_1240k_positions_XY.bed -b ${bam} > ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal.trim3_2ends.BedTools.tsv
    bedtools coverage -a /home/a1222106/fastdir/Refs/qualimap_1240k_positions_XY.bed -b ${bam} -d > ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal.trim3_2ends.depth.tsv
    bedtools coverage -a /home/a1222106/fastdir/Refs/qualimap_1240k_positions_XY.bed -b ${bam} -hist > ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal.trim3_2ends.hist.tsv
    bedtools coverage -a /home/a1222106/fastdir/Refs/qualimap_1240k_positions_XY.bed -b ${bam} -counts > ${sampleName}_${experimentName}_${libraryName}_${runName}_${ref_fasta_basename}_collapsed_bwaALN_sorted_rmdup_IndelReal.trim3_2ends.counts.tsv
"""
}


/* 
Introspection
    - Welcome to comment the below out if you don't want
    - https://www.nextflow.io/docs/latest/metadata.html
*/
workflow.onComplete {
    println """Pipeline execution summary
    ---------------------------
    Completed at : ${workflow.complete}
    Duration     : ${workflow.duration}
    Success      : ${workflow.success}
    Work Dir     : ${workflow.workDir}
    Exit status  : ${workflow.exitStatus}
    Error report : ${workflow.errorReport ?: '-'}
    """
}

/*
Pipeline functions
*/

// Help message
def helpMessage() {
    
    log.info pipelineHeader()
    log.info "\n"
    log.info """

    Usage:

    Typical command to run the pipeline:

    \$ nextflow run fusion/main.nf -profile slurm --outDir ./out_dir --sampleSheet ./sampleSheet.csv --genome GRCh37 -w ./custom_wd

    Arguments:
        -profile <local/slurm>                      Running locally on HPC (can take values local/slurm)
        -w </path/to/working/dir>                   Path to specified working directory (default: ./work)
        --outDir </path/to/output/dir>              Output directory path
        --sampelSheet </path/to/samplesheet.csv>    File path to sample sheet containing samples to run
        --genome <GRCh37>                           Genome version to use

    """
    .stripIndent()
}


// Header for help page
def pipelineHeader() {

    // Log colors ANSI codes
    c_reset = "\033[0m";
    c_dim = "\033[2m";
    c_green = "\033[0;32m";
    c_lblue = "\033[1;34m";
    c_yellow = "\033[0;33m";
    c_lred = "\033[1;31m";
    c_red = "\033[0;31m";
    c_purple = "\033[0;35m";

    return """
-${c_dim}------------------------------------------------------------------------------------------${c_reset}-

                    ${c_green}         8888888b.  888b    888        d8888${c_reset} 
                    ${c_green}         888  "Y88b 8888b   888       d88888${c_reset}
                    ${c_yellow}         888    888 88888b  888      d88P888${c_reset}
                    ${c_lred} 8888b.  888    888 888Y88b 888     d88P 888${c_reset} 
                    ${c_red}    "88b 888    888 888 Y88b888    d88P  888${c_reset} 
                    ${c_purple}.d888888 888    888 888  Y88888   d88P   888${c_reset} 
                    ${c_lblue}888  888 888  .d88P 888   Y8888  d8888888888${c_reset} 
                    ${c_lblue}"Y888888 8888888P"  888    Y888 d88P     888${c_reset} 
                                       
                                                            
                                  A Nextflow Pipeline                           

-${c_dim}------------------------------------------------------------------------------------------${c_reset}-

    """.stripIndent()

}

// Check sample sheet validity
def checkTSV(sampleSheet, nCol) {
    // Check sample sheet exists and isn't empty
    ss = file(sampleSheet)
    if (!ss.exists()) exit 1, "Sample sheet file doesn't exist: ${sampleSheet}"
    if (ss.isEmpty()) exit 1, "Sample sheet file is empty: ${sampleSheet}"

    // Check number of columns in TSV
    Channel
        .fromPath(sampleSheet)
        .splitCsv(header: true, sep: '\t')
        .subscribe {val ->

            // Check that the sample sheet has correct no. fields
            if (val.size() != nCol) exit 1, "TSV row is incorrect - incorrect column number: ${val}"
            
            // Check TSV fields are complete (no empty fields)
            val.each { entry ->
                if(!entry.value) exit 1, "TSV value missing - Key: ${entry.key} Value: ${entry.value}"
            }
        }
    
    return true
    
}

// Import fastq files using samplesheet
def getFastq(sampleSheet) {
    // Read lines of sampleSheet
    ssFile = file(sampleSheet)
    fqChannel = Channel
                    .from(ssFile)
                    .splitCsv(header: true, sep: '\t')
                    .map { row ->

                        // Check if reads files exist
                        if (!file(row.fastqRead1).exists()) exit 1, "Reads file doesn't exist: ${row.fastqRead1}"
                        if (!file(row.fastqRead2).exists()) exit 1, "Reads file doesn't exist: ${row.fastqRead2}"

                        // Check if file is empty
                        if (file(row.fastqRead1).isEmpty()) exit 1, "Reads file is empty: ${row.fastqRead1}"
                        if (file(row.fastqRead2).isEmpty()) exit 1, "Reads file is empty: ${row.fastqRead2}"

                        // Returning list object
                        lst = [ row.experimentName, 
                                row.sampleName, 
                                row.libraryName, 
                                row.unitName,
                                row.platformName,
                                row.runName,
                                tuple( file(row.fastqRead1), file(row.fastqRead2) ) ]
                        return lst
                    }
}

// This extracts the required information from the genomes.config file
def checkParamReturnFile(item) {
    params."${item}" = params.genomes[params.genome]."${item}"

    f = file(params."${item}")
    if (!f.exists()) exit 1, "File in genome config does not exist: ${f}"
    if (f.isEmpty()) exit 1, "File in genome config file is empty: ${f}"

    return f
}

// This calls the function above to extract information from the 
// genomes.config file and assign it to a variable
def defineReferenceMap() {
  // Check genome version is in genome configuration file
  if (!(params.genome in params.genomes)) exit 1, "Genome ${params.genome} not found in configuration"
  
  // Return genome map 
  lst = [
    // FASTA genome reference
    'genomeFile'       : checkParamReturnFile("genomeFile"),
    // genome .fai file
    'genomeIndex'      : checkParamReturnFile("genomeIndex"),
    // Genome Dictionary
    'genomeDict'       : checkParamReturnFile("genomeDict"),
  ]

  return lst
}
