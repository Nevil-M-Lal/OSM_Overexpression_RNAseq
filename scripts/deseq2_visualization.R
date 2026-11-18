

suppressPackageStartupMessages({
  library(DESeq2)
})


dds <- readRDS("results/deseq2/DESeq2_dataset.rds")

output_dir <- "results/deseq2/plots"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(
  "condition" %in% colnames(colData(dds)),
  ncol(dds) >= 3
)


res <- results(
  dds,
  contrast = c("condition", "OSM_overexpression", "Control")
)

res_df <- as.data.frame(res)

# Load annotated results

annotated_results <- read.csv(
  "results/deseq2/DESeq2_OSM_vs_Control.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(
  all(c("ENSEMBL", "SYMBOL") %in% colnames(annotated_results))
)


vsd <- varianceStabilizingTransformation(
  dds,
  blind = FALSE
)

vst_counts <- assay(vsd)

sample_info <- as.data.frame(colData(dds))

stopifnot(
  identical(colnames(vst_counts), rownames(sample_info)),
  all(is.finite(vst_counts))
)

sample_colors <- ifelse(
  sample_info$condition == "Control",
  "steelblue",
  "firebrick"
)


gene_variances <- apply(vst_counts, 1, var)

top_variable_genes <- order(
  gene_variances,
  decreasing = TRUE
)[1:min(500, length(gene_variances))]

pca <- prcomp(
  t(vst_counts[top_variable_genes, , drop = FALSE]),
  center = TRUE,
  scale. = FALSE
)

percent_var <- round(
  100 * pca$sdev^2 / sum(pca$sdev^2),
  1
)

png(
  file.path(output_dir, "PCA.png"),
  width = 2000,
  height = 1400,
  res = 200
)

par(mar = c(5, 4, 4, 10), xpd = NA)

plot(
  pca$x[, 1],
  pca$x[, 2],
  col = sample_colors,
  pch = 19,
  cex = 1.5,
  xlab = paste0("PC1 (", percent_var[1], "%)"),
  ylab = paste0("PC2 (", percent_var[2], "%)"),
  main = "PCA of RNA-seq Samples"
)

text(
  pca$x[, 1],
  pca$x[, 2],
  labels = rownames(sample_info),
  pos = 3,
  cex = 0.8
)

legend(
  "topright",
  inset = c(-0.45, 0),
  legend = c("Control", "OSM overexpression"),
  col = c("steelblue", "firebrick"),
  pch = 19,
  bty = "n"
)

dev.off()


sample_cor <- cor(
  vst_counts,
  method = "pearson"
)

png(
  file.path(output_dir, "sample_correlation_heatmap.png"),
  width = 1600,
  height = 1400,
  res = 200
)

heatmap(
  sample_cor,
  scale = "none",
  symm = TRUE,
  margins = c(9, 9),
  main = "Sample-to-Sample Correlation"
)

dev.off()


png(
  file.path(output_dir, "MA_plot.png"),
  width = 1600,
  height = 1200,
  res = 200
)

plotMA(
  res,
  alpha = 0.05,
  main = "MA Plot: OSM vs Control"
)

dev.off()


valid <- !is.na(res_df$padj) &
         !is.na(res_df$log2FoldChange) &
         is.finite(res_df$log2FoldChange)

volcano_data <- res_df[valid, ]

neg_log10_padj <- -log10(
  pmax(volcano_data$padj, .Machine$double.xmin)
)

is_upregulated <- volcano_data$padj < 0.05 &
                  volcano_data$log2FoldChange >= 1

is_downregulated <- volcano_data$padj < 0.05 &
                    volcano_data$log2FoldChange <= -1

volcano_colors <- rep("grey60", nrow(volcano_data))

volcano_colors[is_upregulated] <- "firebrick"
volcano_colors[is_downregulated] <- "steelblue"

png(
  file.path(output_dir, "volcano_plot.png"),
  width = 1800,
  height = 1400,
  res = 200
)

plot(
  volcano_data$log2FoldChange,
  neg_log10_padj,
  col = volcano_colors,
  pch = 19,
  cex = 0.5,
  xlab = "log2 Fold Change",
  ylab = "-log10 Adjusted P-value",
  main = "Volcano Plot: OSM vs Control"
)

abline(
  v = c(-1, 1),
  h = -log10(0.05),
  lty = 2,
  col = "black"
)

legend(
  "topright",
  legend = c(
    "Upregulated",
    "Downregulated",
    "Other genes"
  ),
  col = c("firebrick", "steelblue", "grey60"),
  pch = 19
)

dev.off()


significant_genes <- res_df[
  !is.na(res_df$padj) & res_df$padj < 0.05,
]

significant_genes <- significant_genes[
  order(significant_genes$padj),
]

top_genes <- head(
  rownames(significant_genes),
  30
)

deg_matrix <- vst_counts[
  top_genes,
  ,
  drop = FALSE
]

# Drop zero-variance genes

gene_sd <- apply(deg_matrix, 1, sd)

deg_matrix <- deg_matrix[
  is.finite(gene_sd) & gene_sd > 0,
  ,
  drop = FALSE
]

stopifnot(nrow(deg_matrix) > 0)

# Z-scores

deg_zscore <- t(scale(t(deg_matrix)))

# Match gene symbols

annotation_index <- match(
  rownames(deg_zscore),
  annotated_results$ENSEMBL
)

gene_labels <- annotated_results$SYMBOL[
  annotation_index
]

# Ensembl fallback

missing_symbols <- is.na(gene_labels) |
                   gene_labels == ""

gene_labels[missing_symbols] <- rownames(deg_zscore)[
  missing_symbols
]

# Handle duplicate symbols

duplicated_symbols <- duplicated(gene_labels) |
                      duplicated(gene_labels, fromLast = TRUE)

gene_labels[duplicated_symbols] <- paste0(
  gene_labels[duplicated_symbols],
  " (",
  rownames(deg_zscore)[duplicated_symbols],
  ")"
)

png(
  file.path(output_dir, "DEG_heatmap.png"),
  width = 1800,
  height = 1800,
  res = 200
)

heatmap(
  deg_zscore,
  scale = "none",
  labRow = gene_labels,
  ColSideColors = sample_colors,
  margins = c(9, 12),
  main = "Top Significant DEGs"
)

dev.off()


writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo.txt")
)


cat("\n===== VISUALISATION COMPLETE =====\n")

cat("PCA: PCA.png\n")
cat("Correlation heatmap: sample_correlation_heatmap.png\n")
cat("MA plot: MA_plot.png\n")
cat("Volcano plot: volcano_plot.png\n")
cat("DEG heatmap: DEG_heatmap.png\n")

cat("\nOutput directory:", output_dir, "\n")
