process MULTIQC {
  publishDir "${params.outdir}/multiqc", mode: 'copy'

  input:
    path(qc_dirs)

  output:
    path("multiqc_report.html")

  script:
  """
  multiqc ${qc_dirs} -o .
  """
}
