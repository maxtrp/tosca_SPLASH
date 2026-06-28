#!/usr/bin/env nextflow

/*
========================================================================================
                                    amchakra/tosca
========================================================================================
hiCLIP/proximity ligation analysis pipeline.
 #### Homepage / Documentation
 https://github.com/amchakra/tosca
----------------------------------------------------------------------------------------
*/

// Define DSL2
nextflow.enable.dsl=2

// Processes
include { hiclipheader } from './modules/utils.nf'
include { METADATA } from './workflows/metadata.nf'
include { CUTADAPT } from './modules/cutadapt.nf'
include { PREMAP } from './workflows/premap.nf'
include { GET_HYBRIDS } from './workflows/gethybrids.nf'
include { TRACK_READ_FATE } from './modules/logging.nf'
include { GET_NON_HYBRIDS } from './modules/getnonhybrids.nf'
include { PROCESS_HYBRIDS } from './workflows/processhybrids.nf'
include { GET_VISUALISATIONS } from './workflows/getvisualisations.nf'
include { GET_ATLAS } from './workflows/getatlas.nf'
include { MAKE_REPORT } from './workflows/makereport.nf'

// SPLASH-specific processes
include { SPLASH_FASTP_DEDUPLICATION } from './modules/splashspecific.nf'
include { SPLASH_TRUNCATE_FASTQ_SEQID } from './modules/splashspecific.nf'
include { SPLASH_CUTADAPT } from './modules/splashspecific.nf'
include { SPLASH_INDEX_FASTA as SPLASH_INDEX_FASTA_GENOME } from './modules/splashspecific.nf'
include { SPLASH_INDEX_FASTA as SPLASH_INDEX_FASTA_TRANSCRIPT } from './modules/splashspecific.nf'
include { SPLASH_MAKE_PSEUDO_TRACKS } from './modules/splashspecific.nf'
include { SPLASH_TRUNCATE_FASTA_SEQID } from './modules/splashspecific.nf'
include { SPLASH_MERGE_FASTQ } from './modules/splashspecific.nf'
include { SPLASH_PEAR } from './modules/splashspecific.nf'

// Genome variables
// if(params.org && params.genomesdir) {

//     params.genome_fai = params.genomes[ params.org ].genome_fai
//     params.transcript_gtf = params.genomes[ params.org ].transcript_gtf
//     params.regions_gtf = params.genomes[ params.org ].regions_gtf

// //     params.genome_fai = params.genomes[ params.org ].genome_fai
// //     params.transcript_fa = params.genomes[ params.org ].transcript_fa
// //     params.transcript_fai = params.genomes[ params.org ].transcript_fai
// //     params.transcript_gtf = params.genomes[ params.org ].transcript_gtf
// //     params.star_genome = params.genomes[ params.org ].star_genome
// //     params.regions_gtf = params.genomes[ params.org ].regions_gtf

//     if(!params.genome_fai) { exit 1, "--genome_fai is not specified." }
//     if(!params.transcript_gtf) { exit 1, "--transcript_gtf is not specified." }
//     if(!params.regions_gtf) { exit 1, "--regions_gtf is not specified." }

// //     if(!params.genome_fai) { exit 1, "--genome_fai is not specified." } 
// //     if(!params.transcript_fa) { exit 1, "--transcript_fa is not specified." } 
// //     if(!params.transcript_fai) { exit 1, "--transcript_fai is not specified." } 
// //     if(!params.transcript_gtf) { exit 1, "--transcript_gtf is not specified." } 
// //     if(!params.regions_gtf) { exit 1, "--regions_gtf is not specified." } 

// // }


// // Create channels for static files
// ch_genome_fai = Channel.fromPath(params.genome_fai, checkIfExists: true)
// ch_transcript_gtf = Channel.fromPath(params.transcript_gtf, checkIfExists: true)
ch_regions_gtf = Channel.fromPath(params.regions_gtf, checkIfExists: true)

// If not making atlas
if(!params.atlas) {

    if(params.org && params.genomesdir) {

        params.transcript_fa = params.genomes[ params.org ].transcript_fa
        params.transcript_fai = params.genomes[ params.org ].transcript_fai
        params.star_genome = params.genomes[ params.org ].star_genome

    } else {

        if(!params.transcript_fa) { exit 1, "--transcript_fa is not specified." }
        // if(!params.transcript_fai) { exit 1, "--transcript_fai is not specified." }

    }


    // Create channels for static files
    ch_transcript_fa = Channel.fromPath(params.transcript_fa, checkIfExists: true)
    // ch_transcript_fai = Channel.fromPath(params.transcript_fai, checkIfExists: true)

    // Channels for optional inputs
    if(!params.skip_premap) {
        ch_star_genome = Channel.fromPath(params.star_genome, checkIfExists: true)
    } else {
        ch_star_genome = Channel.empty()
    }

    if(params.goi) {
        ch_goi = Channel.fromPath(params.goi, checkIfExists: true)
    } else {
        ch_goi = Channel.empty()
    }

}

// Channel for MultiQC config
ch_multiqc_config = Channel.fromPath(params.multiqc_config, checkIfExists: true)

// Show header
log.info hiclipheader()

def summary = [:]
summary['Output directory'] = params.outdir
summary['Trace directory'] = params.tracedir
if(workflow.repository) summary['Pipeline repository'] = workflow.repository
if(workflow.revision) summary['Pipeline revision'] = workflow.revision
summary['Pipeline directory'] = workflow.projectDir
summary['Working dir'] = workflow.workDir
summary['Run name'] = workflow.runName
summary['Profile'] = workflow.profile
if(workflow.container) summary['Container'] = workflow.container
if(params.keep_intermediates) summary['Keep intermediates'] = params.keep_intermediates
if(!params.keep_cache) summary['Keep cache'] = params.keep_cache
log.info summary.collect { k,v -> "${k.padRight(25)}: $v" }.join("\n")
log.info "-\033[2m---------------------------------------------------------------\033[0m-"

def settings = [:]
settings['Organism'] = params.org
if(params.skip_qc) { settings['Skip QC'] = params.skip_qc }
// if(params.skip_atlas) { settings['Skip atlas generation'] = params.skip_atlas }
if(params.skip_premap) { settings['Skip premapping'] = params.skip_premap }
settings['Adapter sequence'] = params.adapter
settings['Minimum read quality'] = params.min_quality
settings['Minimum read length'] = params.min_readlength
settings['FASTQ split size'] = params.split_size
settings['Minimum e-value'] = params.evalue
settings['Maximum hits/read'] = params.maxhits
settings['Deduplication method'] = params.dedup_method
if(params.dedup_method != 'none') settings['UMI separator'] = params.umi_separator
if(params.slurm) settings['Use SLURM'] = params.slurm
settings['Clustering chunk number'] = params.chunk_number
settings['Clustering sample size'] = params.sample_size
settings['Clustering overlap'] = params.percent_overlap
settings['Analyse structures'] = params.analyse_structures
if(params.analyse_structures) settings['Analyse clusters only'] = params.clusters_only
if(params.analyse_structures) settings['Analyse shuffled energies'] = params.shuffled_energies

if(params.goi) { settings['Genes of interest'] = params.goi }
if(params.goi) { settings['Bin size for contact maps'] = params.bin_size }
if(params.goi) { settings['Breaks for arcs'] = params.breaks }
log.info settings.collect { k,v -> "${k.padRight(25)}: $v" }.join("\n")
log.info "-----------------------------------------------------------------"

// Pipeline
workflow {

    if(params.merge_fastq && params.assemble_pe_reads){ 
        exit 1, "Cannot set both --merge_fastq and --assemble_pe_reads to true. Only enable one. "
        }

    if(params.atlas) {

        Channel.fromPath(params.input)
               .map { path -> [ params.atlas, path ] }
               .groupTuple(by: 0)
               .set { ch_all_hybrids }
            //    .view()
            // [atlas, [hybrids_1.tsv.gz, hybrids_2.tsv.gz, hybrids_3.tsv.gz, ...]

        GET_ATLAS(ch_all_hybrids, ch_transcript_gtf, ch_regions_gtf, ch_genome_fai)

    } else {

        /*
        PREPARE INPUTS
        */
        METADATA(params.input) // Get fastq paths
        // CUTADAPT(METADATA.out) // Trim adapters

        if(params.merge_fastq) {
            // Concatenate reads if sample is from multiple fastq files (e.g. one fastq per lane of flowcell) 
            merged_reads_ch = SPLASH_MERGE_FASTQ(METADATA.out)
            SPLASH_FASTP_DEDUPLICATION(merged_reads_ch) // Remove PCR duplicates using fastp
        } else if(params.assemble_pe_reads) {
            assembled_reads_ch = SPLASH_PEAR(METADATA.out) // Assemble paired-end reads into a single fastq
            SPLASH_FASTP_DEDUPLICATION(assembled_reads_ch.fastq) // Remove PCR duplicates using fastp
        } else {
            SPLASH_FASTP_DEDUPLICATION(METADATA.out) // Remove PCR duplicates using fastp
        }        

        SPLASH_TRUNCATE_FASTQ_SEQID(SPLASH_FASTP_DEDUPLICATION.out.fastq) // Modify seq ids
        SPLASH_CUTADAPT(SPLASH_TRUNCATE_FASTQ_SEQID.out.fastq) // Trim adapters

        ch_transcript_fa = SPLASH_TRUNCATE_FASTA_SEQID(params.transcript_fa) // Modify seq ids
        ch_genome_fai = SPLASH_INDEX_FASTA_GENOME(ch_transcript_fa) // same as ch_transcript_fai
        ch_transcript_fai = SPLASH_INDEX_FASTA_TRANSCRIPT(ch_transcript_fa) // same as ch_genome_fai
        ch_transcript_gtf = SPLASH_MAKE_PSEUDO_TRACKS(ch_transcript_fa) // needed for converting transcript coords to genomic coords 
                                                                // (which are the same when aligning to a transcriptome)
        /*
        IDENTIFY HYBRIDS
        */
        if (!params.skip_premap) {
            // PREMAP(CUTADAPT.out.fastq, ch_star_genome)
            PREMAP(SPLASH_CUTADAPT.out.fastq, ch_star_genome)
            ch_for_hybrids = PREMAP.out.fastq
            ch_premap_log = PREMAP.out.logs
        } else {
            // ch_for_hybrids = CUTADAPT.out.fastq
            ch_for_hybrids = SPLASH_CUTADAPT.out.fastq
            ch_premap_log = Channel.empty()
        }

        GET_HYBRIDS(ch_for_hybrids, ch_transcript_fa) // Identify hybrids

        /*
        TRACK READ FATE
        */
        // ch_for_read_fate = CUTADAPT.out.log
        ch_for_read_fate = SPLASH_CUTADAPT.out.log
        .join(ch_premap_log, by: 0, remainder: true)
        .join(GET_HYBRIDS.out.logs, by: 0)
        .map { tuple ->
            def sample_id = tuple[0]
            def logs = tuple[1..-1].findAll { it != null }
            [sample_id, logs]
        }
            // .view { "Channel for TRACK_READ_FATE: $it" }

        TRACK_READ_FATE(ch_for_read_fate)

        /*
        IDENTIFY NON-HYBRIDS
        */
        if (params.merge_fastq) {
            GET_NON_HYBRIDS(GET_HYBRIDS.out.hybrids.join(merged_reads_ch))
        } else if(params.assemble_pe_reads) {
            GET_NON_HYBRIDS(GET_HYBRIDS.out.hybrids.join(assembled_reads_ch.fastq))
        } else {
            GET_NON_HYBRIDS(GET_HYBRIDS.out.hybrids.join(METADATA.out))
        }

        /*
        PROCESS HYBRIDS
        */
        PROCESS_HYBRIDS(GET_HYBRIDS.out.hybrids, ch_transcript_fa, ch_transcript_gtf, ch_regions_gtf)

        // /*
        // GET ATLAS
        // */
        // if(!params.skip_atlas) {
        //     GET_ATLAS(PROCESS_HYBRIDS.out.hybrids, ch_transcript_gtf, ch_regions_gtf, ch_genome_fai)
        // }

        /*
        GET VISUALISATIONS
        */
        GET_VISUALISATIONS(PROCESS_HYBRIDS.out.hybrids, PROCESS_HYBRIDS.out.clusters, ch_genome_fai, ch_transcript_fai, ch_goi)

        /*
        MAKE REPORT
        */
        if(!params.skip_qc) {
            // ch_input_logs = params.skip_premap ? Channel.of([]) : PREMAP.out.logs.collect()
            // ch_input_logs = params.skip_premap ? CUTADAPT.out.log.collect { it[1] } : PREMAP.out.logs.collect { it[1] }
            ch_input_logs = params.skip_premap ? SPLASH_CUTADAPT.out.log.collect { it[1] } : PREMAP.out.logs.collect { it[1] }

            MAKE_REPORT(
                ch_input_logs,
                GET_HYBRIDS.out.logs.collect { it[1..-1].flatten() },
                GET_HYBRIDS.out.raw_hybrids.collect { it[1] },
                PROCESS_HYBRIDS.out.hybrids.collect { it[1] },
                PROCESS_HYBRIDS.out.clusters.collect { it[1] },
                ch_multiqc_config
            )
        }

    }

}

workflow.onComplete {

    if (workflow.success) {
        log.info "-\033[0;34m[Tosca]\033[0;32m Pipeline completed successfully\033[0m-\n"
    } else {
        log.info "-\033[0;34m[Tosca]\033[1;91m Pipeline completed with errors\033[0m\n"
    }

}