#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.samples = "${params.samples ?: 'samples.tsv'}"
params.outdir  = "${params.outdir ?: 'results'}"
params.salmon_index = "${params.salmon_index ?: '/path/to/salmon/index'}"

// channel: read sample sheet
samples_ch = Channel.fromPath(params.samples)
  .splitCsv(sep:'\t', header:true)
  .map { row -> tuple(row.sample_id, file(row.fastq_1), file(row.fastq_2)) }

// include modules
include { FASTQC }        from './modules/fastqc'
include { TRIMGALORE }    from './modules/trim_galore'
include { SALMON_QUANT }  from './modules/salmon'
include { MULTIQC }       from './modules/multiqc'
include { STAR_ALIGN }    from './modules/star'
include { FEATURECOUNTS } from './modules/featurecounts'
include { MERGE_COUNTS }  from './modules/merge_counts'
include { DESEQ2 }        from './modules/deseq2'

// workflow definition
workflow {
  // FASTQC on raw reads
  raw_qc = FASTQC(samples_ch)

  // trimming
  trimmed_ch = TRIMGALORE(samples_ch)

  if (params.aligner == 'salmon') {
    quant_ch = SALMON_QUANT(trimmed_ch)
    MULTIQC(raw_qc.mix(trimmed_ch.map{ it[1..2] }).mix(quant_ch))
  }
  else if (params.aligner == 'star') {
    bam_ch = STAR_ALIGN(trimmed_ch)
    counts_ch = FEATURECOUNTS(bam_ch)
    
    merged_counts = MERGE_COUNTS(counts_ch.collect())

    // run DESeq2 using merged counts + design file
    deseq_results = DESEQ2( merged_counts.mix( Channel.fromPath(params.design) ) )

    MULTIQC(raw_qc.mix(trimmed_ch.map{ it[1..2] }).mix(bam_ch).mix(counts_ch).mix(merged_counts).mix(deseq_results))

  }

  // collect QC + quant outputs
  MULTIQC(raw_qc.mix(trimmed_ch.map{ it[1..2] }).mix(quant_ch))
}
