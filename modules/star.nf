process STAR_ALIGN {
  conda 'envs/star.yml'
  
  tag { sample_id }
  publishDir "${params.outdir}/star", mode: 'copy'

  input:
    tuple val(sample_id), path(read1), path(read2)
    path star_index  

  output:
    tuple val(sample_id), path("${sample_id}.bam")

  script:
  """
  STAR --genomeDir ${star_index} \
       --readFilesIn ${read1} ${read2} \
       --readFilesCommand zcat \
       --runThreadN ${task.cpus} \
       --outSAMtype BAM SortedByCoordinate \
       --outFileNamePrefix ${sample_id}.

  mv ${sample_id}.Aligned.sortedByCoord.out.bam ${sample_id}.bam
  """
}
