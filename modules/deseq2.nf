process DESEQ2 {
  publishDir "${params.outdir}/deseq2", mode: 'copy'

  input:
    path(count_matrix)
    path(design)

  output:
    path("deseq2_results.tsv")
    path("deseq2_plots.pdf")

  script:
  """
  Rscript - <<'EOF'
  suppressMessages({
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
  })

  # read counts and metadata
  counts <- read.table("${count_matrix}", header=TRUE, sep="\\t", row.names=1, check.names=FALSE)
  coldata <- read.table("${design}", header=TRUE, sep="\\t", row.names=1)

  # filter counts to numeric matrix
  count_mat <- as.matrix(counts[ , !(colnames(counts) %in% c("Length")) ])

  # align metadata
  coldata <- coldata[colnames(count_mat), , drop=FALSE]

  # build DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData=count_mat,
                                colData=coldata,
                                design=~condition)

  # run DESeq2
  dds <- DESeq(dds)
  res <- results(dds)

  # write results
  write.table(as.data.frame(res), file="deseq2_results.tsv", sep="\\t", quote=FALSE)

  # QC / plots
  pdf("deseq2_plots.pdf")
  plotMA(res, main="DESeq2 MA-plot")
  vsd <- vst(dds, blind=TRUE)
  plotPCA(vsd, intgroup="condition")
  dev.off()
EOF
  """
}
