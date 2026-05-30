source("analysis/00_config.R")

suppressPackageStartupMessages({
  library(readr)
})

required_files <- c(
  file.path(DIRS$tables, "tcga_kirc_sample_summary.csv"),
  file.path(DIRS$tables, "gse40435_sample_summary.csv"),
  file.path(DIRS$tables, "gse53757_sample_summary.csv"),
  FILES$tcga_deg,
  file.path(DIRS$tables, "gse40435_limma_tumor_vs_normal.csv"),
  file.path(DIRS$tables, "gse53757_limma_tumor_vs_normal.csv"),
  file.path(DIRS$tables, "reproducible_deg_tcga_gse40435_gse53757.csv"),
  FILES$tcga_survival,
  FILES$tcga_enrichment,
  file.path(DIRS$tables, "candidate_gene_evidence_table.csv"),
  file.path(DIRS$tables, "strict_candidate_genes.csv"),
  file.path(DIRS$tables, "high_confidence_candidate_genes.csv"),
  file.path(DIRS$tables, "high_confidence_candidate_pathway_summary.csv"),
  file.path(DIRS$tables, "candidate_summary.csv"),
  file.path(DIRS$figures, "tcga_kirc_pca.png"),
  file.path(DIRS$figures, "tcga_kirc_volcano.png"),
  file.path(DIRS$figures, "tcga_kirc_discordance.png")
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required pipeline outputs: ", paste(missing_files, collapse = ", "))
}

if (dir.exists(file.path("results", "logs"))) {
  stop("results/logs exists, but the pipeline must not create a logs directory.")
}

summary <- read_csv(file.path(DIRS$tables, "candidate_summary.csv"), show_col_types = FALSE)
required_metrics <- c(
  "reproducible_deg",
  "main_stage_grade_complete_prognostic",
  "ph_pass",
  "sensitivity_pass",
  "strict_candidate",
  "high_confidence_candidate"
)

missing_metrics <- setdiff(required_metrics, summary$metric)
if (length(missing_metrics) > 0) {
  stop("candidate_summary.csv is missing metrics: ", paste(missing_metrics, collapse = ", "))
}

values <- setNames(summary$value, summary$metric)
if (values[["high_confidence_candidate"]] > values[["strict_candidate"]]) {
  stop("High-confidence candidate count exceeds strict candidate count.")
}
if (values[["strict_candidate"]] > values[["sensitivity_pass"]]) {
  stop("Strict candidate count exceeds sensitivity-pass count.")
}

message("Output validation complete.")
