### Addidional parameters

The tweaked version of the pipeline has a few extra parameters required to process SPLASH data:

- `--five_prime_trim_length`: The number of bases to trim off the beginning of reads - default `3`. Required for removing the template-switching oligo added during [SMARTer smRNA-Seq](https://www.takarabio.com/documents/User%20Manual/SMARTer%20smRNA/SMARTer%20smRNA-Seq%20Kit%20for%20Illumina%20User%20Manual.pdf) library prep.
- `--merge_fastq`: Concatenates fastq files matching a file path containing a wildcard in the samplesheet (e.g. `sample_L00*_R1.fastq.gz`) into a single fastq - default `false`. It is useful for consolidating a sample split across multiple lanes of a sequening run. This param can't be used with paired-end data. 
- `--assemble_pe_reads`: Assembles paired-end reads into a single fastq file using [PEAR](https://doi.org/10.1093/bioinformatics/btt593) - default `false`. A slightly modified samplesheet structure is required to pass paired-end data to the pipeline:
```
sample,fastq1,fastq2
sample1,/path/to/file1_R1.fastq.gz,/path/to/file1_R2.fastq.gz
sample2,/path/to/file2_R1.fastq.gz,/path/to/file2_R2.fastq.gz
sample3,/path/to/file3_R1.fastq.gz,/path/to/file3_R2.fastq.gz
```

### Example usage

Sample containing a single RNA of interest (single-end data):

```
nextflow run main.nf \
    -profile SPLASH_singRNA \
    --input /path/to/samplesheet.csv \
    --goi /path/to/goi.txt \
    --transcript_fa /path/to/reference.fa \
    --outdir results/ \
```

Sample containing multiple RNAs of interest (single-end data):

```
nextflow run main.nf \
    -profile SPLASH_multiRNA \
    --input /path/to/samplesheet.csv \
    --goi /path/to/goi.txt \
    --transcript_fa /path/to/reference.fa \
    --outdir results/ \
```

Sample containing multiple RNAs of interest (paired-end data):

```
nextflow run main.nf \
    --assemble_pe_reads \
    -profile SPLASH_multiRNA \
    --input /path/to/samplesheet.csv \
    --goi /path/to/goi.txt \
    --transcript_fa /path/to/reference.fa \
    --outdir results/ \
```