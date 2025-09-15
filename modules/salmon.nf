process SALMON_QUANT {
  tag { sample_id }
  publishDir "${params.outdir}/salmon", mode: 'copy'

  input:
    tuple val(sample_id), path(read1), path(read2)

  output:
    path("${sample_id}_quant")

  script:
  """
  salmon quant -i ${params.salmon_index} -l A \
    -1 ${read1} -2 ${read2} \
    -p ${task.cpus} -o ${sample_id}_quant
  """
}
