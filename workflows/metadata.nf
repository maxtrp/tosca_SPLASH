#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

workflow METADATA {
    take: csv
    main:

        if(params.interleave_pe_reads) {
            Channel
                .fromPath( csv )
                .splitCsv(header:true)
                .map { row -> [ row.sample, [file(row.fastq1, checkIfExists: true), file(row.fastq2, checkIfExists: true)] ]  }
                .set { data }
        } else{
            Channel
                .fromPath( csv )
                .splitCsv(header:true)
                .map { row -> [ row.sample, file(row.fastq, checkIfExists: true) ]  }
                .set { data }            
        }

    emit:
        data
}