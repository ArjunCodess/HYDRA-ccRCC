source("analysis/00_config.R")
source("analysis/functions/io.R")
source("analysis/functions/plotting.R")

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

se <- read_required_rds(FILES$tcga_se)
counts <- SummarizedExperiment::assay(se, "unstranded")
coldata <- as.data.frame(SummarizedExperiment::colData(se)) |>
  tibble::rownames_to_column("sample_barcode") |>
  mutate(condition = factor(shortLetterCode, levels = c("NT", "TP")))

sample_summary <- coldata |>
  count(sample_type, shortLetterCode, name = "n_samples")
write_csv_atomic(sample_summary, file.path(DIRS$tables, "tcga_kirc_sample_summary.csv"))

keep <- rowSums(counts >= THRESHOLDS$min_count) >= THRESHOLDS$min_samples
dds <- DESeqDataSetFromMatrix(
  countData = counts[keep, ],
  colData = coldata,
  design = ~ condition
)

vst_mat <- assay(vst(dds, blind = TRUE))
write_rds_atomic(vst_mat, FILES$tcga_vst)

pca <- prcomp(t(vst_mat), scale. = FALSE)
pca_df <- as.data.frame(pca$x[, 1:2]) |>
  tibble::rownames_to_column("sample_barcode") |>
  left_join(coldata, by = "sample_barcode")

variance <- pca$sdev^2 / sum(pca$sdev^2)
p <- ggplot(pca_df, aes(PC1, PC2, color = sample_type)) +
  geom_point(size = 2, alpha = 0.85) +
  labs(
    title = "TCGA-KIRC RNA-seq PCA",
    x = paste0("PC1 (", round(variance[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(variance[2] * 100, 1), "%)")
  ) +
  theme_hydra()

ggsave(file.path(DIRS$figures, "tcga_kirc_pca.png"), p, width = 7, height = 5, dpi = 300)

message("QC complete. Saved VST matrix and PCA figure.")
