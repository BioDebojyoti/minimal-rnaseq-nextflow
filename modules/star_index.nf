process STAR_INDEX {
  tag "STAR index"

  publishDir "${params.outdir}/reference/star", mode: 'copy'

  input:
    path fasta
    path gtf

  output:
    path "star_index"

  script:
  """
  mkdir star_index
  STAR --runThreadN ${task.cpus} \
       --runMode genomeGenerate \
       --genomeDir star_index \
       --genomeFastaFiles ${fasta} \
       --sjdbGTFfile ${gtf} \
       --sjdbOverhang 99
  """
}
