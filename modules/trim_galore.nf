process TRIMGALORE {
  tag { sample_id }
  publishDir "${params.outdir}/trimmed", mode: 'copy'

  input:
    tuple val(sample_id), path(read1), path(read2)

  output:
    tuple val(sample_id), path("${sample_id}_R1_trimmed.fq.gz"), path("${sample_id}_R2_trimmed.fq.gz")

  script:
  """
  trim_galore --paired --gzip -o . ${read1} ${read2}

  mv ${sample_id}_R1_val_1.fq.gz ${sample_id}_R1_trimmed.fq.gz || true
  mv ${sample_id}_R2_val_2.fq.gz ${sample_id}_R2_trimmed.fq.gz || true
  """
}
