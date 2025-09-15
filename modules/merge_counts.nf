process MERGE_COUNTS {
  publishDir "${params.outdir}/counts", mode: 'copy'

  input:
    path(count_files)

  output:
    path("gene_counts_matrix.txt")

  script:
  """
  # featureCounts outputs tables with headers + comments starting with '#'
  # Use paste/join in R or Python, but here’s a simple awk-based approach.

  # Extract the header from the first file
  awk 'BEGIN{OFS="\\t"} !/^#/ {if(NR==1){header="Geneid"}; next} END{print header}' $(ls ${count_files} | head -n1)

  # Merge counts using Rscript for robustness
  Rscript - <<'EOF'
  files <- Sys.glob("*.txt")
  dflist <- lapply(files, function(f) {
    df <- read.table(f, header=TRUE, comment.char="#", stringsAsFactors=FALSE)
    df <- df[, c("Geneid", "Length", "Counts")]
    colnames(df)[3] <- gsub("_counts.txt$", "", basename(f))
    df
  })
  merged <- Reduce(function(x,y) merge(x,y,by=c("Geneid","Length")), dflist)
  write.table(merged, file="gene_counts_matrix.txt", sep="\\t", quote=FALSE, row.names=FALSE)
EOF
  """
}
