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
coef_name <- "condition_TP_vs_NT"
if (!coef_name %in% resultsNames(dds)) {
  stop(
    "Expected DESeq2 coefficient is unavailable: ", coef_name,
    ". Available coefficients: ", paste(resultsNames(dds), collapse = ", ")
  )
}
res_apeglm <- lfcShrink(dds, coef = coef_name, type = "apeglm")
apeglm_lfc <- setNames(res_apeglm$log2FoldChange, rownames(res_apeglm))
apeglm_lfc_se <- setNames(res_apeglm$lfcSE, rownames(res_apeglm))

res_df <- as.data.frame(res) |>
  tibble::rownames_to_column("gene_id") |>
  arrange(padj) |>
  mutate(
    log2FoldChange_apeglm = unname(apeglm_lfc[gene_id]),
    lfcSE_apeglm = unname(apeglm_lfc_se[gene_id]),
    significant = !is.na(padj) &
      padj < THRESHOLDS$deg_fdr &
      abs(log2FoldChange) >= THRESHOLDS$deg_abs_log2fc
  )

write_csv_atomic(res_df, FILES$tcga_deg)
write_rds_atomic(dds, file.path(DIRS$processed, "tcga_kirc_deseq2_dds.rds"))

message("DESeq2 complete. Significant DEG count: ", sum(res_df$significant, na.rm = TRUE))
