

suppressPackageStartupMessages({
  library(DESeq2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})


count_file <- "results/counts/gene_counts.txt"
metadata_file <- "docs/samples.tsv"
output_dir <- "results/deseq2"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


fc <- read.delim(
  count_file,
  comment.char = "#",
  check.names = FALSE
)

samples <- read.delim(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Extract raw gene counts (featureCounts columns 7 onward)
count_matrix <- as.matrix(fc[, 7:ncol(fc)])

rownames(count_matrix) <- fc$Geneid

# Convert BAM filenames into sample IDs
colnames(count_matrix) <- sub(
  "\\.sorted\\.bam$",
  "",
  basename(colnames(count_matrix))
)


stopifnot(
  ncol(fc) == 12,
  nrow(samples) == 6,
  !anyDuplicated(samples$sample_id),
  !anyDuplicated(rownames(count_matrix)),
  identical(colnames(count_matrix), samples$sample_id),
  !anyNA(samples$condition),
  setequal(
    unique(samples$condition),
    c("Control", "OSM_overexpression")
  ),
  all(is.finite(count_matrix)),
  all(count_matrix >= 0),
  all(count_matrix == floor(count_matrix))
)

storage.mode(count_matrix) <- "integer"

# Control reference
samples$condition <- factor(
  samples$condition,
  levels = c("Control", "OSM_overexpression")
)

rownames(samples) <- samples$sample_id


keep <- rowSums(count_matrix) >= 10

filtered_counts <- count_matrix[keep, ]

cat("Genes before filtering:", nrow(count_matrix), "\n")
cat("Genes retained:", nrow(filtered_counts), "\n")
cat("Genes removed:", sum(!keep), "\n")


dds <- DESeqDataSetFromMatrix(
  countData = filtered_counts,
  colData = samples,
  design = ~ condition
)


dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c(
    "condition",
    "OSM_overexpression",
    "Control"
  ),
  alpha = 0.05
)

res_ordered <- res[order(res$padj), ]



ensembl_ids <- rownames(res_ordered)

ensembl_keys <- sub(
  "\\.[0-9]+$",
  "",
  ensembl_ids
)


gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_keys,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "asNA"
)

entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_keys,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "asNA"
)


res_df <- as.data.frame(res_ordered)

res_df$ENSEMBL <- ensembl_ids
res_df$SYMBOL <- unname(gene_symbols[ensembl_keys])
res_df$ENTREZID <- unname(entrez_ids[ensembl_keys])


res_df <- res_df[
  ,
  c(
    "ENSEMBL",
    "SYMBOL",
    "ENTREZID",
    setdiff(
      colnames(res_df),
      c("ENSEMBL", "SYMBOL", "ENTREZID")
    )
  )
]


write.csv(
  res_df,
  file = file.path(
    output_dir,
    "DESeq2_OSM_vs_Control.csv"
  ),
  row.names = FALSE,
  na = ""
)


significant_degs <- res_df[
  !is.na(res_df$padj) & res_df$padj < 0.05,
]

write.csv(
  significant_degs,
  file = file.path(
    output_dir,
    "significant_DEGs_padj05.csv"
  ),
  row.names = FALSE,
  na = ""
)

# padj + fold-change cutoff

significant_degs_lfc1 <- significant_degs[
  abs(significant_degs$log2FoldChange) >= 1,
]

write.csv(
  significant_degs_lfc1,
  file = file.path(
    output_dir,
    "significant_DEGs_padj05_lfc1.csv"
  ),
  row.names = FALSE,
  na = ""
)


normalized_counts <- counts(
  dds,
  normalized = TRUE
)

write.csv(
  normalized_counts,
  file = file.path(
    output_dir,
    "normalized_counts.csv"
  ),
  row.names = TRUE
)


saveRDS(
  dds,
  file = file.path(
    output_dir,
    "DESeq2_dataset.rds"
  )
)

writeLines(
  capture.output(summary(res)),
  file.path(output_dir, "DESeq2_summary.txt")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo.txt")
)


cat("\n===== DESeq2 ANALYSIS COMPLETE =====\n")
cat("Comparison: OSM_overexpression vs Control\n")
cat("Design: ~ condition\n")
cat("Reference condition: Control\n")
cat("Genes tested:", nrow(res), "\n")

cat(
  "Significant DEGs (adjusted p-value < 0.05):",
  sum(res$padj < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "Significant DEGs with |log2FC| >= 1:",
  nrow(significant_degs_lfc1),
  "\n"
)

cat(
  "Genes with mapped symbols:",
  sum(!is.na(res_df$SYMBOL)),
  "\n"
)

cat("Annotation database version:",
    as.character(packageVersion("org.Hs.eg.db")), "\n")

cat("Results directory:", output_dir, "\n")
