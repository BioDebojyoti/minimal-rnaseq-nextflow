process FEATURECOUNTS {
  tag { sample_id }
  publishDir "${params.outdir}/counts", mode: 'copy'

  input:
    tuple val(sample_id), path(bam)

  output:
    path("${sample_id}_counts.txt")

  script:
  """
  featureCounts -T ${task.cpus} -a ${params.gtf} -o ${sample_id}_counts.txt ${bam}
  """
}
