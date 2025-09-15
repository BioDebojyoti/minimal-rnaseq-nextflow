process SALMON_INDEX {
  tag "Salmon index"

  publishDir "${params.outdir}/reference/salmon", mode: 'copy'

  input:
    path fasta
    path gtf

  output:
    path "salmon_index"

  script:
  """
  salmon index \
    -t ${fasta} \
    -i salmon_index \
    -g ${gtf}
  """
}
