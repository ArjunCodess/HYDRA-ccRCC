# GSE53757 pair audit and unpaired limma sensitivity.
# Does not replace the primary adjacent-row paired table.

source("analysis/00_config.R")
source("analysis/functions/io.R")
source("analysis/05_deg_geo.R")

suppressPackageStartupMessages({
  library(Biobase)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

series <- readRDS(file.path(DIRS$processed, "gse53757_series_matrix.rds"))
pdata <- Biobase::pData(series[[1]])
titles <- as.character(pdata$title)
if (length(titles) %% 2L != 0L) stop("GSE53757 sample count is not even.")

token_of <- function(title) {
  stripped <- gsub("\\s+", "", title)
  sub("[TNtn]$", "", stripped)
}
stage_of <- function(source_name) {
  hit <- str_match(source_name, "Stage\\s+([1-4])")[, 2]
  ifelse(is.na(hit), NA_character_, paste0("Stage ", hit))
}

left <- seq(1L, length(titles), by = 2L)
audit <- tibble(
  pair_index = seq_along(left),
  left_title = titles[left],
  right_title = titles[left + 1L],
  left_source = as.character(pdata$source_name_ch1[left]),
  right_source = as.character(pdata$source_name_ch1[left + 1L]),
  left_tissue = as.character(pdata$`tissue:ch1`[left]),
  right_tissue = as.character(pdata$`tissue:ch1`[left + 1L])
) |>
  mutate(
    left_token = token_of(left_title),
    right_token = token_of(right_title),
    same_title_token = left_token == right_token,
    left_stage = stage_of(left_source),
    right_stage = stage_of(right_source),
    same_stage = !is.na(left_stage) & left_stage == right_stage,
    opposite_tissue = left_tissue != right_tissue
  )
write_csv_atomic(audit, file.path(DIRS$tables, "gse53757_pair_audit.csv"))
write_csv_atomic(tibble(
  metric = c("adjacent_pairs", "opposite_tissue_pairs", "same_title_token_pairs",
             "same_stage_pairs", "shared_patient_field"),
  value = c(nrow(audit), sum(audit$opposite_tissue), sum(audit$same_title_token),
            sum(audit$same_stage), 0)
), file.path(DIRS$tables, "gse53757_pair_audit_summary.csv"))

message("Running unpaired GSE53757 limma.")
run_geo_limma("GSE53757", paired = FALSE)

paired <- read_csv(file.path(DIRS$tables, "gse53757_limma_tumor_vs_normal.csv"), show_col_types = FALSE) |>
  select(symbol, paired_log2fc = log2FoldChange, paired_pvalue = pvalue, paired_padj = padj)
unpaired <- read_csv(file.path(DIRS$tables, "gse53757_unpaired_limma_tumor_vs_normal.csv"), show_col_types = FALSE) |>
  select(symbol, unpaired_log2fc = log2FoldChange, unpaired_pvalue = pvalue, unpaired_padj = padj)
gse40435 <- read_csv(file.path(DIRS$tables, "gse40435_limma_tumor_vs_normal.csv"), show_col_types = FALSE) |>
  select(symbol, gse40435_log2fc = log2FoldChange, gse40435_pvalue = pvalue)
repro <- read_csv(file.path(DIRS$tables, "reproducible_deg_tcga_gse40435_gse53757.csv"), show_col_types = FALSE)
high <- read_csv(file.path(DIRS$tables, "high_confidence_candidate_genes.csv"), show_col_types = FALSE)

alt <- repro |>
  select(symbol, tcga_significant, tcga_direction, gse40435_direction, gse40435_nominal) |>
  inner_join(unpaired, by = "symbol") |>
  mutate(
    unpaired_direction = sign(unpaired_log2fc),
    same_direction_unpaired = !is.na(unpaired_direction) & unpaired_direction == tcga_direction,
    unpaired_nominal = unpaired_pvalue < 0.05,
    nominal_support = (gse40435_nominal & gse40435_direction == tcga_direction) |
      (same_direction_unpaired & unpaired_nominal),
    reproducible_unpaired_gate = tcga_significant &
      gse40435_direction == tcga_direction &
      same_direction_unpaired &
      nominal_support
  )

high_rows <- high |>
  select(symbol) |>
  left_join(paired, by = "symbol") |>
  left_join(unpaired, by = "symbol") |>
  left_join(alt |> select(symbol, reproducible_unpaired_gate), by = "symbol") |>
  mutate(
    same_sign = sign(paired_log2fc) == sign(unpaired_log2fc),
    paired_nominal = paired_pvalue < 0.05,
    unpaired_nominal = unpaired_pvalue < 0.05,
    paired_effect = abs(paired_log2fc) >= THRESHOLDS$min_geo_abs_log2fc,
    unpaired_effect = abs(unpaired_log2fc) >= THRESHOLDS$min_geo_abs_log2fc
  )
write_csv_atomic(high_rows, file.path(DIRS$tables, "gse53757_unpaired_sensitivity.csv"))

paired_repro <- repro$symbol[which(repro$reproducible_deg)]
unpaired_repro <- alt$symbol[which(alt$reproducible_unpaired_gate)]
shared <- intersect(paired_repro, unpaired_repro)
jaccard <- length(shared) / length(union(paired_repro, unpaired_repro))
write_csv_atomic(tibble(
  metric = c("high_confidence_genes", "high_confidence_still_reproducible",
             "high_confidence_same_sign", "paired_reproducible_deg",
             "unpaired_gate_reproducible_deg", "reproducible_jaccard"),
  value = c(nrow(high_rows), sum(high_rows$reproducible_unpaired_gate, na.rm = TRUE),
            sum(high_rows$same_sign, na.rm = TRUE), length(paired_repro),
            length(unpaired_repro), jaccard)
), file.path(DIRS$tables, "gse53757_unpaired_summary.csv"))
message("Unpaired sensitivity written. High-confidence genes still reproducible: ",
        sum(high_rows$reproducible_unpaired_gate, na.rm = TRUE), " of ", nrow(high_rows))
