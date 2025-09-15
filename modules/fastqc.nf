process FASTQC {
  tag { sample_id }
  publishDir "${params.outdir}/fastqc", mode: 'copy'

  input:
    tuple val(sample_id), path(read1), path(read2)

  output:
    path("${sample_id}_fastqc.zip")

  script:
  """
  fastqc -o . ${read1} ${read2}
  """
}
