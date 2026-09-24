# Tables and the cross-cohort forest from results that already exist.

source("analysis/00_config.R")
source("analysis/functions/io.R")
source("analysis/functions/plotting.R")
source("analysis/functions/tcga_metadata.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
})

predictions <- read_csv(file.path(DIRS$tables, "nested_cv_predictions.csv"), show_col_types = FALSE)
fold_events <- predictions |>
  group_by(repeat_id, fold) |>
  summarise(n_test = n(), test_events = sum(os_event), .groups = "drop")
write_csv_atomic(fold_events, file.path(DIRS$tables, "nested_cv_fold_events.csv"))
write_csv_atomic(tibble(
  metric = c("test_events_min", "test_events_median", "test_events_max"),
  value = c(min(fold_events$test_events), median(fold_events$test_events), max(fold_events$test_events))
), file.path(DIRS$tables, "nested_cv_fold_event_summary.csv"))

clinical <- read_csv(FILES$tcga_clinical, show_col_types = FALSE) |>
  transmute(stage = normalize_stage(ajcc_pathologic_stage),
            grade = normalize_grade(tumor_grade))
write_csv_atomic(bind_rows(
  clinical |> count(stage, name = "n") |> mutate(factor = "stage", level = stage) |> select(factor, level, n),
  clinical |> count(grade, name = "n") |> mutate(factor = "grade", level = grade) |> select(factor, level, n)
), file.path(DIRS$tables, "clinical_factor_counts.csv"))

gse <- read_csv(file.path(DIRS$tables, "external_survival_gse29609.csv"), show_col_types = FALSE) |>
  filter(external_present) |>
  transmute(symbol, cohort = "GSE29609", log_hr = external_log_hr,
            ci_low = log(external_hr_ci_low), ci_high = log(external_hr_ci_high))
em <- read_csv(file.path(DIRS$tables, "external_survival_emtab1980.csv"), show_col_types = FALSE) |>
  filter(external_present) |>
  transmute(symbol, cohort = "E-MTAB-1980", log_hr = external_log_hr,
            ci_low = log(external_hr_ci_low), ci_high = log(external_hr_ci_high))
tcga <- read_csv(file.path(DIRS$tables, "high_confidence_candidate_genes.csv"), show_col_types = FALSE) |>
  transmute(symbol, cohort = "TCGA-KIRC", log_hr = main_log_hr,
            ci_low = log(main_hr_ci_low), ci_high = log(main_hr_ci_high))
forest <- bind_rows(tcga, gse, em) |>
  mutate(cohort = factor(cohort, levels = c("TCGA-KIRC", "GSE29609", "E-MTAB-1980")),
         label = if_else(symbol %in% c("DDC", "TCIRG1"), paste0(symbol, " *"), symbol))
gene_order <- tcga |> arrange(log_hr) |> pull(symbol)
forest$label <- factor(forest$label, levels = rev(if_else(gene_order %in% c("DDC", "TCIRG1"),
                                                          paste0(gene_order, " *"), gene_order)))
write_csv_atomic(forest |> mutate(label = as.character(label)),
                 file.path(DIRS$tables, "external_log_hr_forest.csv"))
p <- ggplot(forest, aes(log_hr, label, color = cohort)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y",
                width = 0, linewidth = 0.4, position = position_dodge(width = 0.6)) +
  geom_point(position = position_dodge(width = 0.6), size = 1.8) +
  scale_color_manual(values = c("TCGA-KIRC" = "#205D72", "GSE29609" = "#B85C38", "E-MTAB-1980" = "#3E6B48")) +
  labs(x = "log hazard ratio", y = NULL, color = NULL,
       title = "High-confidence genes across cohorts",
       subtitle = "Stars mark DDC and TCIRG1. Intervals are model confidence intervals.") +
  theme_hydra()
ggsave(file.path(DIRS$figures, "external_log_hr_forest.png"), p, width = 8, height = 7, dpi = 300)

sample <- read_csv(file.path(DIRS$tables, "tcga_kirc_sample_summary.csv"), show_col_types = FALSE)
paired <- read_csv(file.path(DIRS$tables, "tcga_kirc_paired_deg_summary.csv"), show_col_types = FALSE)
gse_s <- read_csv(file.path(DIRS$tables, "external_survival_gse29609_summary.csv"), show_col_types = FALSE)
em_s <- read_csv(file.path(DIRS$tables, "external_survival_emtab1980_summary.csv"), show_col_types = FALSE)
tx <- read_csv(file.path(DIRS$tables, "tracerx_multiregion_study_summary.csv"), show_col_types = FALSE)
cm <- read_csv(file.path(DIRS$tables, "checkmate025_study_summary.csv"), show_col_types = FALSE)
geo_n <- function(accession) {
  x <- read_csv(file.path(DIRS$tables, paste0(tolower(accession), "_sample_summary.csv")), show_col_types = FALSE)
  unique(x$n_patients)
}
value_of <- function(data, metric) data$value[data$metric == metric][1]
cohort <- tibble(
  dataset = c("TCGA-KIRC tumors", "TCGA-KIRC normals", "TCGA-KIRC pairs",
              "GSE40435", "GSE53757", "GSE29609", "E-MTAB-1980", "TRACERx Renal", "CheckMate 025"),
  role = c("discovery expression and survival", "discovery normals", "paired tumor-normal sensitivity",
           "paired expression replication", "paired expression replication",
           "previously inspected survival check", "previously inspected survival check",
           "region sampling sensitivity", "treatment-interaction check"),
  patients = c(sample$n_samples[sample$sample_type == "Primary Tumor"],
               sample$n_samples[sample$sample_type == "Solid Tissue Normal"],
               value_of(paired, "paired_patients"),
               geo_n("GSE40435"), geo_n("GSE53757"),
               value_of(gse_s, "gse29609_samples"), value_of(em_s, "emtab1980_samples"),
               value_of(tx, "tracerx_survival_linked_patients"),
               value_of(cm, "checkmate025_rna_linked_patients")),
  samples_or_regions = c(sample$n_samples[sample$sample_type == "Primary Tumor"],
                         sample$n_samples[sample$sample_type == "Solid Tissue Normal"],
                         value_of(paired, "paired_samples"),
                         geo_n("GSE40435") * 2, geo_n("GSE53757") * 2,
                         value_of(gse_s, "gse29609_samples"), value_of(em_s, "emtab1980_samples"),
                         value_of(tx, "tracerx_primary_regions"),
                         value_of(cm, "checkmate025_rna_linked_patients")),
  events = c(NA, NA, NA, NA, NA,
             value_of(gse_s, "gse29609_events"), value_of(em_s, "emtab1980_events"),
             value_of(tx, "tracerx_survival_events"),
             value_of(cm, "os_events")),
  platform = c("TCGA STAR counts", "TCGA STAR counts", "TCGA STAR counts",
               "GPL570", "GPL570", "GEO expression array", "Illumina HumanHT-12",
               "TRACERx TPM", "CheckMate 025 RNA"),
  previously_inspected = c("no", "no", "no", "no", "no", "yes", "yes", "yes", "yes"),
  script = c("analysis/04_deg_tcga.R", "analysis/03_qc_tcga.R", "analysis/04b_paired_deg_tcga.R",
             "analysis/05_deg_geo.R", "analysis/05_deg_geo.R", "analysis/13_external_survival_gse29609.R",
             "analysis/14_external_survival_emtab1980.R", "analysis/20_tracerx_multiregion_transportability.R",
             "analysis/21_checkmate025_treatment_interaction.R")
)
write_csv_atomic(cohort, file.path(DIRS$tables, "cohort_dictionary.csv"))
message("Evidence display tables and forest written.")
