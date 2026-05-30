source("analysis/00_config.R")
source("analysis/functions/io.R")

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(DESeq2)
  library(dplyr)
  library(readr)
})

se <- read_required_rds(FILES$tcga_se)
counts <- SummarizedExperiment::assay(se, "unstranded")
coldata <- as.data.frame(SummarizedExperiment::colData(se)) |>
  tibble::rownames_to_column("sample_barcode") |>
  mutate(condition = factor(shortLetterCode, levels = c("NT", "TP")))

keep <- rowSums(counts >= THRESHOLDS$min_count) >= THRESHOLDS$min_samples

dds <- DESeqDataSetFromMatrix(
  countData = counts[keep, ],
  colData = coldata,
  design = ~ condition
)

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "TP", "NT"))
res_df <- as.data.frame(res) |>
  tibble::rownames_to_column("gene_id") |>
  arrange(padj) |>
  mutate(
    significant = !is.na(padj) &
      padj < THRESHOLDS$deg_fdr &
      abs(log2FoldChange) >= THRESHOLDS$deg_abs_log2fc
  )

write_csv_atomic(res_df, FILES$tcga_deg)
write_rds_atomic(dds, file.path(DIRS$processed, "tcga_kirc_deseq2_dds.rds"))

message("DESeq2 complete. Significant DEG count: ", sum(res_df$significant, na.rm = TRUE))
