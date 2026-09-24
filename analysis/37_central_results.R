# Central displays from saved results. Does not refit DESeq2 or the nested CV.
# Ridge boundaries are recorded from the implementation and from the saved folds.

source("analysis/00_config.R")
source("analysis/functions/io.R")
source("analysis/functions/plotting.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(tidyr)
})

tab <- function(name) read_csv(file.path(DIRS$tables, name), show_col_types = FALSE)
value_of <- function(data, metric) {
  hit <- data$value[data$metric == metric]
  if (length(hit) != 1L || !is.finite(hit)) stop("Missing metric: ", metric)
  hit
}

genes <- tab("high_confidence_candidate_genes.csv")
if (nrow(genes) != 23L) stop("Expected 23 high-confidence genes.")

geo_a <- tab("gse40435_limma_tumor_vs_normal.csv") |>
  transmute(symbol, gse40435_pvalue = pvalue, gse40435_fdr = padj) |>
  distinct(symbol, .keep_all = TRUE)
geo_b <- tab("gse53757_limma_tumor_vs_normal.csv") |>
  transmute(symbol, gse53757_pvalue = pvalue, gse53757_fdr = padj) |>
  distinct(symbol, .keep_all = TRUE)
gse <- tab("external_survival_gse29609.csv") |>
  transmute(symbol, gse29609_present = external_present, gse29609_log_hr = external_log_hr,
            gse29609_p = external_p_value, gse29609_fdr = external_fdr,
            gse29609_same_direction = external_same_direction)
em <- tab("external_survival_emtab1980.csv") |>
  transmute(symbol, emtab_present = external_present, emtab_log_hr = external_log_hr,
            emtab_p = external_p_value, emtab_fdr = external_fdr,
            emtab_same_direction = external_same_direction)
composition <- tab("candidate_clinical_composition_sensitivity.csv") |>
  transmute(symbol, marker_log_hr = composition_adjusted_log_hr,
            marker_fdr = composition_adjusted_fdr)
purity <- tab("candidate_direct_tumor_purity_sensitivity.csv") |>
  transmute(symbol, purity_log_hr = gene_log_hr, purity_fdr = gene_fdr)
frequency <- tab("nested_gene_selection_frequency.csv") |>
  transmute(symbol, nested_folds = folds)
tracerx <- tab("tracerx_one_region_cox_summary.csv") |>
  filter(scenario == "fixed_subset_regions") |>
  transmute(symbol, tracerx_fixed_direction = same_tcga_direction_fraction,
            tracerx_fixed_median_log_hr = median_log_hr)
checkmate <- tab("checkmate025_candidate_treatment_interactions.csv") |>
  filter(endpoint == "OS", model == "age_sex_mskcc_adjusted") |>
  transmute(symbol, os_interaction_log_hr = interaction_log_hr,
            os_interaction_fdr = interaction_fdr)
any_interaction <- tab("checkmate025_candidate_treatment_interactions.csv") |>
  group_by(symbol) |>
  summarise(any_interaction_fdr_lt_0.05 = any(interaction_fdr < 0.05), .groups = "drop")
aliquot <- tab("aliquot_sensitivity_summary.csv")
retained_all <- Reduce(intersect, strsplit(aliquot$retained, ";"))

evidence <- genes |>
  left_join(geo_a, by = "symbol") |>
  left_join(geo_b, by = "symbol") |>
  left_join(gse, by = "symbol") |>
  left_join(em, by = "symbol") |>
  left_join(composition, by = "symbol") |>
  left_join(purity, by = "symbol") |>
  left_join(frequency, by = "symbol") |>
  left_join(tracerx, by = "symbol") |>
  left_join(checkmate, by = "symbol") |>
  left_join(any_interaction, by = "symbol") |>
  mutate(
    nested_folds = replace_na(nested_folds, 0L),
    gse40435_nominal = gse40435_pvalue < 0.05 & sign(gse40435_log2fc) == sign(tcga_log2fc),
    gse53757_nominal = gse53757_pvalue < 0.05 & sign(gse53757_log2fc) == sign(tcga_log2fc),
    marker_fdr_support = marker_fdr < 0.05,
    purity_fdr_support = purity_fdr < 0.05,
    retained_under_every_aliquot_rule = symbol %in% retained_all,
    evidence_profile = case_when(
      symbol == "DDC" ~ "external reversal",
      symbol == "CRYL1" ~ "transportable association on one inspected cohort",
      symbol == "ACADM" ~ "concordant metabolic association",
      !gse29609_same_direction & gse29609_fdr < 0.05 ~ "external reversal",
      !marker_fdr_support ~ "composition-sensitive association",
      gse29609_same_direction & emtab_same_direction ~ "directionally concordant association",
      TRUE ~ "discovery association with mixed external direction"
    )
  ) |>
  select(symbol, evidence_score, tcga_log2fc, main_log_hr, main_hr, main_fdr,
         gse40435_log2fc, gse40435_pvalue, gse40435_nominal,
         gse53757_log2fc, gse53757_pvalue, gse53757_nominal,
         gse29609_present, gse29609_log_hr, gse29609_fdr, gse29609_same_direction,
         emtab_present, emtab_log_hr, emtab_fdr, emtab_same_direction,
         marker_fdr, marker_fdr_support, purity_fdr, purity_fdr_support,
         nested_folds, tracerx_fixed_direction, tracerx_fixed_median_log_hr,
         os_interaction_log_hr, os_interaction_fdr, any_interaction_fdr_lt_0.05,
         retained_under_every_aliquot_rule, evidence_profile)
if (anyNA(evidence$gse40435_pvalue) || anyNA(evidence$os_interaction_fdr)) {
  stop("Candidate evidence matrix is missing a joined column.")
}
if (sum(evidence$any_interaction_fdr_lt_0.05) != 0L) {
  stop("An interaction FDR below 0.05 contradicts the locked CheckMate result.")
}
write_csv_atomic(evidence, file.path(DIRS$tables, "candidate_evidence_matrix.csv"))

ridge_folds <- tab("nested_benchmark_folds.csv") |> filter(strategy == "ridge_eligible")
ridge_summary <- tab("nested_benchmark_summary.csv") |> filter(strategy == "ridge_eligible")
nested <- tab("nested_cv_summary.csv")
hydra <- tab("nested_benchmark_summary.csv") |> filter(strategy == "hydra")
unpaired <- tab("gse53757_unpaired_summary.csv")
clearcode <- tab("published_signature_summary.csv") |>
  filter(cohort == "TCGA-KIRC", arm == "clinical_plus_clearcode34")
if (nrow(ridge_folds) != 50L || any(ridge_folds$fallback_to_clinical) || any(!is.finite(ridge_folds$lambda))) {
  stop("Ridge arm is missing a finite training lambda.")
}
if (abs(hydra$mean_delta_c - nested$mean_delta_c) > 1e-12) {
  stop("HYDRA benchmark increment does not match the nested CV increment.")
}
spec <- tibble(
  item = c(
    "model", "alpha", "lambda_rule", "inner_criterion", "ties", "standardize",
    "clinical_columns", "gene_columns", "expression", "size_factors",
    "inner_folds", "inner_seed", "outer_folds", "geo_gate", "what_this_does_not_do"
  ),
  definition = c(
    "glmnet Cox, clinical covariates unpenalized, genes ridge-penalized",
    "0",
    "lambda.min",
    "partial-likelihood deviance inside the training fold",
    "efron",
    "glmnet standardizes every training column and applies those moments to the test rows",
    "age, sex, stage, and grade; penalty factor 0; factor levels from the complete-case cohort; values from the training rows only",
    "training-fold reproducible DEGs with nonzero training standard deviation; penalty factor 1; none are removed by the penalty",
    "log2(count / size factor + 1), not the full-data variance-stabilizing transform",
    "training DESeq2 size factors; test size factors use the training geometric-mean reference",
    "five event-stratified folds of the training patients only",
    "RESAMPLING$seed + 31000 + repeat_id * 10 + fold, drawn before cv.glmnet",
    "the saved ten-by-five patient splits; test outcomes are not an input",
    "fixed full GEO tables; GEO cohorts are not re-split",
    "not the prespecified one-gene acceptance test; not a sparse gene selection; patient bootstrap does not refit lambda"
  )
)
write_csv_atomic(spec, file.path(DIRS$tables, "ridge_specification.csv"))
write_csv_atomic(ridge_folds |>
  transmute(repeat_id, fold, lambda, n_genes = n_model_genes,
            train_repro_genes, fallback_to_clinical),
  file.path(DIRS$tables, "ridge_fold_audit.csv"))
write_csv_atomic(tibble(
  metric = c("folds", "fallback_folds", "lambda_min", "lambda_median", "lambda_max",
             "mean_n_genes", "min_n_genes", "max_n_genes",
             "mean_c", "mean_delta_c", "delta_c_ci_low", "delta_c_ci_high",
             "mean_brier3", "mean_delta_brier3", "delta_brier_ci_low", "delta_brier_ci_high",
             "mean_calibration_slope"),
  value = c(nrow(ridge_folds), sum(ridge_folds$fallback_to_clinical),
            min(ridge_folds$lambda), median(ridge_folds$lambda), max(ridge_folds$lambda),
            mean(ridge_folds$n_model_genes), min(ridge_folds$n_model_genes), max(ridge_folds$n_model_genes),
            ridge_summary$mean_c, ridge_summary$mean_delta_c,
            ridge_summary$patient_bootstrap_ci_low, ridge_summary$patient_bootstrap_ci_high,
            ridge_summary$mean_brier3, ridge_summary$mean_delta_brier3,
            ridge_summary$patient_bootstrap_brier_ci_low, ridge_summary$patient_bootstrap_brier_ci_high,
            ridge_summary$mean_calibration_slope)
), file.path(DIRS$tables, "ridge_multigene_summary.csv"))

checks <- tibble(
  check = c("high_confidence_n", "hydra_matches_nested", "ridge_no_fallback",
            "ridge_delta", "unpaired_keeps_23", "unpaired_jaccard", "clearcode_genes",
            "no_interaction_fdr"),
  passed = c(
    nrow(evidence) == 23L,
    abs(hydra$mean_delta_c - nested$mean_delta_c) < 1e-12,
    all(!ridge_folds$fallback_to_clinical),
    abs(ridge_summary$mean_delta_c - 0.025688259109311752) < 1e-12,
    value_of(unpaired, "high_confidence_still_reproducible") == 23,
    value_of(unpaired, "reproducible_jaccard") > 0.99,
    clearcode$genes_used == 33,
    sum(evidence$any_interaction_fdr_lt_0.05) == 0L
  )
)
if (any(!checks$passed)) stop("Numerical consistency failed: ", paste(checks$check[!checks$passed], collapse = ", "))
write_csv_atomic(checks, file.path(DIRS$tables, "numerical_consistency.csv"))

arithmetic <- tab("tcga_aliquot_arithmetic.csv")
funnel <- tab("candidate_summary.csv")
spine <- tibble(
  step = c("Tumor aliquots", "Patients", "Reproducible DEGs", "Strict candidates", "High-confidence"),
  n = c(value_of(arithmetic, "tumor_aliquots"), value_of(arithmetic, "tumor_patients"),
        value_of(funnel, "reproducible_deg"), value_of(funnel, "strict_candidate"),
        value_of(funnel, "high_confidence_candidate")),
  note = c("8 extra aliquots", "one tumor each", "both GEO cohorts", "survival gates", "one-gene rule")
)
spine$step <- factor(spine$step, levels = rev(spine$step))
p_spine <- ggplot(spine, aes(n, step)) +
  geom_col(fill = "#205D72", width = 0.72) +
  geom_text(aes(label = format(n, big.mark = ",")), hjust = -0.08, size = 3.4) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(x = NULL, y = NULL, title = "What the funnel keeps",
       subtitle = "Counts are gates, not a claim that the last gate predicts.") +
  theme_hydra(base_size = 11) +
  theme(panel.grid.major.y = element_blank())

bench <- tab("nested_benchmark_summary.csv") |>
  filter(strategy != "clinical") |>
  transmute(arm = recode(strategy,
                         hydra = "HYDRA gene",
                         survival_only = "Survival only",
                         de_only = "DE only",
                         matched_control = "Matched control",
                         ridge_eligible = "Ridge on reproducible DEGs"),
            delta = mean_delta_c, lo = patient_bootstrap_ci_low, hi = patient_bootstrap_ci_high)
bench <- bind_rows(bench, tibble(
  arm = "ClearCode34",
  delta = clearcode$delta_c, lo = clearcode$delta_c_ci_low, hi = clearcode$delta_c_ci_high
))
arm_levels <- c("HYDRA gene", "Survival only", "DE only", "Matched control",
                "ClearCode34", "Ridge on reproducible DEGs")
bench$arm <- factor(bench$arm, levels = rev(arm_levels))
p_bench <- ggplot(bench, aes(delta, arm)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_vline(xintercept = 0.01, linetype = "dashed", linewidth = 0.4, color = "#9b3d2a") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0, linewidth = 0.5) +
  geom_point(size = 2.4, color = "#205D72") +
  labs(x = "Held-out concordance minus clinical", y = NULL,
       title = "Same patients, same splits",
       subtitle = "Dashed line is the prespecified +0.01 rule. Intervals are conditional on the saved predictions.") +
  theme_hydra(base_size = 11)

pairing <- tab("gse53757_unpaired_sensitivity.csv")
p_pair <- ggplot(pairing, aes(paired_log2fc, unpaired_log2fc)) +
  geom_abline(linewidth = 0.3, color = "grey40") +
  geom_point(color = "#205D72", size = 2.2) +
  labs(x = "Paired log2 fold change", y = "Unpaired log2 fold change",
       title = "GSE53757 pairing sensitivity",
       subtitle = "All 23 high-confidence genes. Jaccard overlap of the reproducible lists is 0.998.") +
  theme_hydra(base_size = 11)
ggsave(file.path(DIRS$figures, "gse53757_pairing_sensitivity.png"), p_pair,
       width = 6.2, height = 5.4, dpi = 300)

forest <- tab("external_log_hr_forest.csv")
forest$reversal <- forest$symbol %in% c("DDC", "TCIRG1")
forest$label <- factor(forest$label, levels = unique(forest$label[order(match(forest$symbol, genes$symbol[order(genes$main_log_hr)]))]))
# Rebuild labels in TCGA log-HR order from the saved forest file.
tcga_order <- forest |> filter(cohort == "TCGA-KIRC") |> arrange(log_hr) |> pull(symbol)
forest <- forest |>
  mutate(label = if_else(symbol %in% c("DDC", "TCIRG1"), paste0(symbol, "  reversal"), symbol),
         label = factor(label, levels = rev(if_else(tcga_order %in% c("DDC", "TCIRG1"),
                                                    paste0(tcga_order, "  reversal"), tcga_order))))
p_forest <- ggplot(forest, aes(log_hr, label, color = cohort)) +
  geom_rect(data = data.frame(label = factor(c("DDC  reversal", "TCIRG1  reversal"),
                                              levels = levels(forest$label))),
            aes(ymin = as.numeric(label) - 0.45, ymax = as.numeric(label) + 0.45),
            xmin = -Inf, xmax = Inf, inherit.aes = FALSE, fill = "#F4D6C6", alpha = 0.85) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y",
                width = 0, linewidth = 0.4, position = position_dodge(width = 0.55)) +
  geom_point(position = position_dodge(width = 0.55), size = 1.9) +
  scale_color_manual(values = c("TCGA-KIRC" = "#205D72", "GSE29609" = "#B85C38", "E-MTAB-1980" = "#3E6B48")) +
  labs(x = "log hazard ratio", y = NULL, color = NULL,
       title = "External survival does not copy TCGA",
       subtitle = "Shaded rows are DDC and TCIRG1. GSE29609 has 17 deaths; E-MTAB-1980 has 23.") +
  theme_hydra()
ggsave(file.path(DIRS$figures, "external_log_hr_forest.png"), p_forest,
       width = 8.2, height = 7.2, dpi = 300)

png(file.path(DIRS$figures, "master_funnel_benchmark.png"), width = 11.4, height = 5.6, units = "in", res = 300)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2, widths = grid::unit(c(1, 1.25), "null"))))
print(p_spine, vp = grid::viewport(layout.pos.col = 1))
print(p_bench, vp = grid::viewport(layout.pos.col = 2))
dev.off()

message("Central evidence table, ridge audit, and figures written.")
