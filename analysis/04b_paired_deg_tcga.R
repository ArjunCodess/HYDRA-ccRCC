source("analysis/00_config.R")
source("analysis/functions/io.R")
source("analysis/functions/patient_samples.R")

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(DESeq2)
  library(dplyr)
  library(readr)
})

se <- read_required_rds(FILES$tcga_se)
counts <- SummarizedExperiment::assay(se, "unstranded")
meta <- as.data.frame(SummarizedExperiment::colData(se)) |>
  tibble::rownames_to_column("sample_barcode") |>
  mutate(patient_barcode = substr(sample_barcode, 1, 12),
         condition = factor(shortLetterCode, levels = c("NT", "TP")))
meta <- select_tcga_patient_samples(meta, counts)
paired <- meta |>
  count(patient_barcode) |>
  filter(n == 2L) |>
  pull(patient_barcode)
meta <- meta |> filter(patient_barcode %in% paired)
stopifnot(!anyDuplicated(paste(meta$patient_barcode, meta$condition)),
          nrow(meta) == 2L * length(paired), length(paired) > 0L)
counts <- counts[, meta$sample_barcode, drop = FALSE]
meta <- as.data.frame(meta)
rownames(meta) <- meta$sample_barcode
meta$patient_barcode <- factor(meta$patient_barcode)
keep <- rowSums(counts >= THRESHOLDS$min_count) >= THRESHOLDS$min_samples
dds <- DESeqDataSetFromMatrix(counts[keep, ], meta, design = ~ patient_barcode + condition)
dds <- DESeq(dds)
result <- as.data.frame(results(dds, contrast = c("condition", "TP", "NT"))) |>
  tibble::rownames_to_column("gene_id") |>
  mutate(significant = !is.na(padj) & padj < THRESHOLDS$deg_fdr &
           abs(log2FoldChange) >= THRESHOLDS$deg_abs_log2fc)
write_csv_atomic(result, file.path(DIRS$tables, "tcga_kirc_paired_deseq2_tumor_vs_normal.csv"))
write_csv_atomic(tibble(metric = c("paired_patients", "paired_samples", "tested_genes", "significant_genes"),
                        value = c(length(paired), nrow(meta), nrow(result), sum(result$significant))),
                 file.path(DIRS$tables, "tcga_kirc_paired_deg_summary.csv"))
message("Paired TCGA DE complete: ", length(paired), " patients; ", sum(result$significant), " significant genes.")
