#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process SPLASH_FASTP_DEDUPLICATION {

    tag "${sample_id}"
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
    publishDir "${params.outdir}/logs/${sample_id}/", mode: "copy", pattern: "*.log"

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.trimmed.fastq.gz"), emit: fastq
        path("*.cutadapt.log"), emit: log

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

process SPLASH_INDEX_FASTA {

    tag "${ref_name}"
    
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

    stripped_seqnames = sapply(strsplit(names(fasta_sequences), " "), function(x) x[1]) # Assign first item of sequence identifier as name

    transcriptome.gr = GRanges(
        seqnames = stripped_seqnames,
        ranges = IRanges(start = rep(1, length(fasta_sequences)), width = width(fasta_sequences)),  # Entire sequence is the range
        strand = rep("+", length(fasta_sequences)),  # RNA on positive strand
        fasta_id = stripped_seqnames  # Store the FASTA sequence name as fasta_id
    )

    export.gff2(transcriptome.gr, "${ref_name}.gtf")

    """
    
}