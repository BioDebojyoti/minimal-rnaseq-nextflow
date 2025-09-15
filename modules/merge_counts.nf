process MERGE_COUNTS {
  conda 'envs/deseq2.yml'
  
  publishDir "${params.outdir}/counts", mode: 'copy'

  input:
    path count_files

  output:
    path("gene_counts_matrix.txt")

  script:
  """
  Rscript - <<'EOF'
  suppressMessages({
    library(DESeq2)
  })

  files <- Sys.glob("${count_files}/*.txt")
  dflist <- lapply(files, function(f) {
    df <- read.table(f, header=TRUE, comment.char="#", stringsAsFactors=FALSE)
    df <- df[, c("Geneid","Length","Counts")]
    colnames(df)[3] <- gsub("_counts.txt\$", '', basename(f))
    df
  })
  merged <- Reduce(function(x,y) merge(x,y,by=c("Geneid","Length")), dflist)
  write.table(merged, file="gene_counts_matrix.txt", sep="\t", quote=FALSE, row.names=FALSE)
EOF
  """
}
