#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process SPLASH_MERGE_FASTQ {

    tag "${sample_id}"
    label 'process_medium'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.fastq.gz"), emit: merged_reads

    script:
    """
    pigz -dc ${reads} | pigz -c > merged_${sample_id}.fastq.gz
    """

}

process SPLASH_PEAR {

    tag "${sample_id}"
    label 'process_medium'
    container 'quay.io/biocontainers/pear:0.9.6--hb1d24b7_13' // run pear in docker biocontainer
    publishDir "${params.outdir}/logs/${sample_id}/", mode: "copy", pattern: "*.log"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.assembled.fastq.gz"), emit: fastq
    path "*.log"                                      , emit: log

    script:
    """
    pear --forward-fastq ${reads[0]} --reverse-fastq ${reads[1]} --output ${sample_id} --threads $task.cpus > ${sample_id}_pear.log

    gzip ${sample_id}.assembled.fastq
    """

}

process SPLASH_FASTP_DEDUPLICATION {

    tag "${sample_id}"
    label 'process_medium'    
    container 'quay.io/biocontainers/fastp:0.24.2--heae3180_0' // run fastp in docker biocontainer
    publishDir "${params.outdir}/logs/${sample_id}/", mode: "copy", pattern: "*.log"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.fastq.gz"), emit: fastq
    path "*.log"                            , emit: log

    script:
    """
    fastp --thread $task.cpus --dedup --disable_adapter_trimming --in1 ${reads} --out1 deduplicated_${sample_id}.fastq.gz 2> ${sample_id}_fastp_dedup.log
    """

}

process SPLASH_TRUNCATE_FASTQ_SEQID {
    // removes characters following the space from sequence identifier lines:
    // @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<sample number>
    //                                                                       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ removes this
    // needed to prevent errors during IDENTIFY_HYBRIDS step (pblat alignments omit these characters)

    tag "${sample_id}"
    label 'process_medium'

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.trunc_ids.fastq.gz"), emit: fastq

    script:
    """
    pigz -dc ${reads} | sed '/^@/s/ .*//' | pigz > ${sample_id}.trunc_ids.fastq.gz
    """
}

process SPLASH_CUTADAPT {

    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/logs/${sample_id}/", mode: "copy", pattern: "*.log"

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.trimmed.fastq.gz"), emit: fastq
        tuple val(sample_id), path("*.cutadapt.log"), emit: log

    script:
    args = " -j ${task.cpus}"
    args += " -a " + params.adapter 
    args += " -q " + params.min_quality
    args += " --minimum-length " + params.min_readlength
    args += " -o ${sample_id}.trimmed.fastq.gz"
    args += " --cut " + params.five_prime_trim_length

    """
    cutadapt $args $reads > ${sample_id}.cutadapt.log
    """
}

process SPLASH_TRUNCATE_FASTA_SEQID {
    // removes characters following the space from sequence identifier lines (similar to SPLASH_TRUNCATE_FASTQ_SEQID)
    // needed so that sequence names are correctly assigned in downstream processing (viz. SPLASH_MAKE_PSEUDO_TRACKS and ANALYSE_STRUCTURES:CHUNK_SEQUENCES)

    tag "${ref_name}"
    label 'process_low'
    
    input:
        path(fasta)
    
    output:
        path("*_strippedID.fa")
    
    script:
    ref_name = fasta.getSimpleName()
    """
    sed '/^>/s/ .*//' ${fasta} >  ${ref_name}_strippedID.fa
    """
}

process SPLASH_INDEX_FASTA {

    tag "${ref_name}"
    label 'process_low'
    
    input:
        path(fasta)

    output:
        path("*.fai"), emit: fai

    script:
    ref_name = fasta.getSimpleName()
    """
    samtools faidx ${fasta} -o ${fasta}.fai
    """
}

process SPLASH_MAKE_PSEUDO_TRACKS {
    // generate a (pseudo) transcript annotation (gtf) file 
    // used by tosca for converting from transcript to genomic coordinates

    tag "${ref_name}"
    label 'process_low'

    input:
        path(fasta)

    output:
        path("*.gtf"), emit: gtf

    script:
    ref_name = fasta.getSimpleName()
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(Biostrings))
    suppressPackageStartupMessages(library(GenomicRanges))
    suppressPackageStartupMessages(library(rtracklayer))

    fasta_sequences = readDNAStringSet("$fasta")

    transcriptome.gr = GRanges(
        seqnames = names(fasta_sequences), # Assumes sequence identifier lines have no spaces (i.e. no metadata, just name)
        ranges = IRanges(start = rep(1, length(fasta_sequences)), width = width(fasta_sequences)),  # Entire sequence is the range
        strand = rep("+", length(fasta_sequences)),  # RNA on positive strand
        fasta_id = names(fasta_sequences)  # Store the FASTA sequence identifier as fasta_id
    )

    export.gff2(transcriptome.gr, "${ref_name}.gtf")

    """
    
}