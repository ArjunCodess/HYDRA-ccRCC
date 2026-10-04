source("analysis/00_config.R")
suppressPackageStartupMessages({ library(readr); library(dplyr); library(jsonlite) })
tab <- function(name) read_csv(file.path(DIRS$tables, paste0(name, ".csv")), show_col_types = FALSE)
ledger <- tab("candidate_ledger")
all_counts <- tab("ledger_funnel_counts")
funnel <- filter(all_counts, role == "serial_discovery_selection")
repro <- tab("reproducible_deg_tcga_gse40435_gse53757")
candidate <- tab("candidate_gene_evidence_table")
cox <- tab("tcga_kirc_cox_models")
stopifnot(!anyDuplicated(ledger$gene_id), !anyDuplicated(ledger$symbol),
  setequal(ledger$gene_id, repro$tcga_gene_id), setequal(ledger$symbol, repro$symbol))
aligned <- candidate[match(ledger$gene_id, candidate$tcga_gene_id), ]
for (pair in list(c("log_hr", "main_log_hr"), c("log_hr_ci_low", "main_log_hr_ci_low"),
    c("log_hr_ci_high", "main_log_hr_ci_high"), c("p_value", "main_p_value"),
    c("adjusted_p", "main_fdr"), c("ph_p_value", "main_ph_p_value"), c("n", "main_n"), c("events", "main_events"))) {
  stopifnot(isTRUE(all.equal(ledger[[pair[1]]], aligned[[pair[2]]], tolerance = 1e-12)))
}
# Derive all gates from estimates, rather than accepting copied candidate flags.
expected <- candidate |>
  mutate(discovery_de = !is.na(tcga_padj) & tcga_padj < THRESHOLDS$deg_fdr & abs(tcga_log2fc) >= THRESHOLDS$deg_abs_log2fc,
    reproducible_de = discovery_de & same_direction_count >= 2 & nominal_support_count >= 1,
    adjusted_survival = reproducible_de & !is.na(main_fdr) & main_fdr < THRESHOLDS$strict_survival_fdr,
    covariate_sensitivity = adjusted_survival &
      sign(stage_complete_log_hr) == sign(main_log_hr) & sign(grade_complete_log_hr) == sign(main_log_hr) &
      stage_complete_p_value < 0.05 & grade_complete_p_value < 0.05,
    strict_rebuilt = covariate_sensitivity & abs(main_log_hr) >= THRESHOLDS$min_abs_log_hr &
      abs(gse40435_log2fc) >= THRESHOLDS$min_geo_abs_log2fc & abs(gse53757_log2fc) >= THRESHOLDS$min_geo_abs_log2fc,
    priority_rebuilt = strict_rebuilt & main_fdr < THRESHOLDS$high_confidence_survival_fdr & abs(main_log_hr) >= THRESHOLDS$high_confidence_abs_log_hr)
stages <- c("discovery_de", "reproducible_de", "adjusted_survival", "covariate_sensitivity", "strict_candidate", "priority_shortlist")
expected_cols <- c(stages[1:4], "strict_rebuilt", "priority_rebuilt")
for (i in seq_along(stages)) {
  genes <- expected$symbol[expected[[expected_cols[i]]] %in% TRUE]
  stopifnot(setequal(genes, ledger$symbol[ledger[[stages[i]]] %in% TRUE]),
    funnel$count[funnel$stage == stages[i]] == length(genes))
}
stopifnot(funnel$count[funnel$stage == "mapped_tcga"] == nrow(ledger),
  all(diff(funnel$count) <= 0), all(funnel$removed_from_previous[-1] == -diff(funnel$count)),
  setequal(ledger$symbol[ledger$priority_shortlist], tab("high_confidence_candidate_genes")$symbol),
  setequal(ledger$symbol[ledger$strict_candidate], tab("strict_candidate_genes")$symbol))
for (model in unique(cox$model_type)) {
  x <- filter(cox, model_type == model)
  stopifnot(isTRUE(all.equal(x$fdr, p.adjust(x$p_value, "BH"), tolerance = 1e-12)))
}
for (spec in list(c("external_survival_gse29609", "external_p_value", "external_fdr"),
    c("external_survival_emtab1980", "external_p_value", "external_fdr"),
    c("external_survival_emtab1980", "external_adjusted_p_value", "external_adjusted_fdr"),
    c("candidate_direct_tumor_purity_sensitivity", "gene_p_value", "gene_fdr"),
    c("candidate_clinical_composition_sensitivity", "composition_adjusted_p_value", "composition_adjusted_fdr"))) {
  x <- tab(spec[1])
  stopifnot(isTRUE(all.equal(x[[spec[3]]], p.adjust(x[[spec[2]]], "BH"), tolerance = 1e-12)))
}
legacy_funnel <- tab("evidence_funnel")
stopifnot(identical(as.numeric(legacy_funnel$count), as.numeric(funnel$count[-1])))
summary <- tab("candidate_summary")
for (pair in list(c("reproducible_deg", "reproducible_de"), c("main_stage_grade_complete_prognostic", "adjusted_survival"),
    c("sensitivity_pass", "covariate_sensitivity"), c("strict_candidate", "strict_candidate"), c("high_confidence_candidate", "priority_shortlist"))) {
  stopifnot(summary$value[summary$metric == pair[1]] == funnel$count[funnel$stage == pair[2]])
}
for (field in grep("_json$", names(ledger), value = TRUE)) stopifnot(all(vapply(ledger[[field]], jsonlite::validate, logical(1))))
new_sources <- c(harmonized_cohort_models_json = "limitations_harmonized_cox",
  heterogeneity_json = "limitations_heterogeneity",
  missing_covariate_scenarios_json = "limitations_missing_covariate_scenarios",
  source_site_sensitivity_json = "limitations_source_site_cox",
  plate_de_sensitivity_json = "limitations_plate_de_candidates")
for (field in names(new_sources)) {
  source_rows <- tab(new_sources[[field]])
  stopifnot(field %in% names(ledger))
  for (gene_symbol in unique(source_rows$symbol)) {
    stored <- fromJSON(ledger[[field]][ledger$symbol == gene_symbol])
    expected_rows <- as.data.frame(source_rows[source_rows$symbol == gene_symbol, ])
    stopifnot(isTRUE(all.equal(stored, expected_rows, check.attributes = FALSE, tolerance = 1e-12)))
  }
  missing_rows <- !ledger$symbol %in% source_rows$symbol
  stopifnot(all(ledger[[field]][missing_rows] == '{"status":"not_assessed"}'))
}
stopifnot(all(ledger$effect_direction[is.na(ledger$log_hr)] == "not_assessed"),
  all(ledger$ph_check[is.na(ledger$ph_p_value)] == "not_assessed"))

# Existing macros are audited even if no optional declarations have been added.
macro_map <- tibble(macro = c("MappedDeg", "ReproDeg", "MainSurvival", "SensitivityPass", "StrictGenes", "HighGenes"), stage = stages)
if (file.exists("paper/ledger_claims.csv")) macro_map <- bind_rows(macro_map, read_csv("paper/ledger_claims.csv", show_col_types = FALSE)) |> distinct()
macros <- readLines("paper/results_macros.tex", warn = FALSE)
for (i in seq_len(nrow(macro_map))) {
  prefix <- paste0("\\newcommand{\\", macro_map$macro[i], "}{")
  lines <- macros[startsWith(macros, prefix)]
  stopifnot(length(lines) == 1)
  value <- gsub(",", "", substr(lines, nchar(prefix) + 1, nchar(lines) - 1))
  stopifnot(as.numeric(value) == funnel$count[funnel$stage == macro_map$stage[i]])
}
# Optional declarations: each named-gene assertion must name a ledger stage.
if (file.exists("paper/ledger_gene_claims.csv")) {
  claims <- read_csv("paper/ledger_gene_claims.csv", show_col_types = FALSE)
  stopifnot(all(c("claim", "symbol", "stage") %in% names(claims)), all(claims$stage %in% stages))
  for (i in seq_len(nrow(claims))) stopifnot(claims$symbol[i] %in% ledger$symbol[ledger[[claims$stage[i]]] %in% TRUE])
}
procedures <- tab("funnel_external_summary")
review_macros <- readLines("paper/review_macros.tex", warn = FALSE)
assert_review_macro <- function(name, expected) {
  prefix <- paste0("\\newcommand{\\", name, "}{")
  lines <- review_macros[startsWith(review_macros, prefix)]
  stopifnot(length(lines) == 1L)
  value <- as.numeric(substr(lines, nchar(prefix) + 1, nchar(lines) - 1))
  stopifnot(value == expected)
  if (name %in% all_counts$stage) stopifnot(all_counts$count[all_counts$stage == name] == expected)
}
specs <- list(c("ReviewEmDeMapped", "de_only", "E-MTAB-1980", "present"),
  c("ReviewEmDeSame", "de_only", "E-MTAB-1980", "same_direction"),
  c("ReviewEmCompleteMapped", "complete_rule", "E-MTAB-1980", "present"),
  c("ReviewEmCompleteSame", "complete_rule", "E-MTAB-1980", "same_direction"),
  c("ReviewGseDeMapped", "de_only", "GSE29609", "present"),
  c("ReviewGseDeSame", "de_only", "GSE29609", "same_direction"),
  c("ReviewGseCompleteMapped", "complete_rule", "GSE29609", "present"),
  c("ReviewGseCompleteSame", "complete_rule", "GSE29609", "same_direction"),
  c("ReviewEmDeFdr", "de_only", "E-MTAB-1980", "same_direction_fdr"),
  c("ReviewEmCompleteFdr", "complete_rule", "E-MTAB-1980", "same_direction_fdr"))
for (spec in specs) {
  row <- procedures[procedures$rule == spec[2] & procedures$cohort == spec[3], ]
  stopifnot(nrow(row) == 1L)
  assert_review_macro(spec[1], row[[spec[4]]])
}
em <- tab("external_survival_emtab1980"); gse <- tab("external_survival_gse29609")
marker <- tab("candidate_clinical_composition_sensitivity")
purity <- tab("candidate_direct_tumor_purity_sensitivity")
bootstrap <- tab("candidate_cox_bootstrap_summary")
assert_review_macro("ReviewMarkerRetained", sum(marker$composition_adjusted_fdr < 0.05, na.rm = TRUE))
assert_review_macro("ReviewPhFlags", sum(ledger$priority_shortlist & ledger$ph_p_value < 0.05, na.rm = TRUE))
assert_review_macro("ReviewEmMapped", sum(em$external_present))
assert_review_macro("ReviewEmSame", sum(em$external_present & em$external_same_direction, na.rm = TRUE))
assert_review_macro("ReviewEmFdr", sum(em$external_strict_support, na.rm = TRUE))
assert_review_macro("ReviewGseMapped", sum(gse$external_present))
assert_review_macro("ReviewGseSame", sum(gse$external_present & gse$external_same_direction, na.rm = TRUE))
assert_review_macro("ReviewBootstrapIntervals", sum(bootstrap$ci_excludes_zero, na.rm = TRUE))
assert_review_macro("ReviewPurityRetained", sum(purity$gene_fdr < 0.05, na.rm = TRUE))
assert_review_macro("ReviewDiscoveryAttritionPct", round(100 * (1 - sum(ledger$priority_shortlist) / sum(ledger$discovery_de)), 2))
assert_review_macro("ReviewMappedAttritionPct", round(100 * (1 - sum(ledger$priority_shortlist) / nrow(ledger)), 2))
stopifnot(all(c("ACADM", "DDC", "TCIRG1") %in% ledger$symbol[ledger$priority_shortlist]))
reversed <- gse$symbol[gse$external_present & !gse$external_same_direction & gse$external_fdr < 0.05]
stopifnot(setequal(reversed, c("DDC", "TCIRG1")))
marker_loss <- marker$symbol[marker$composition_adjusted_fdr >= 0.05]
stopifnot(setequal(marker_loss, c("DDC", "KL", "GJB1", "CLCN5", "PODXL", "HIBCH")))
acadm <- em[em$symbol == "ACADM", ]
stopifnot(nrow(acadm) == 1L, acadm$external_strict_support)
table_lines <- readLines("paper/review_candidate_table.tex", warn = FALSE)
named_rows <- sub(" .*", "", table_lines[grepl("^[A-Z0-9]+ & ", table_lines)])
stopifnot(length(named_rows) == sum(ledger$priority_shortlist), setequal(named_rows, ledger$symbol[ledger$priority_shortlist]))
main_text <- paste(readLines("paper/main.tex", warn = FALSE), collapse = "\n")
supplement_text <- paste(readLines("paper/supplement.tex", warn = FALSE), collapse = "\n")
benchmark <- tab("nested_benchmark_summary")
ridge <- benchmark[benchmark$strategy == "ridge_eligible", ]
control <- benchmark[benchmark$strategy == "matched_control", ]
if (grepl("Its Brier difference interval spans zero", main_text, fixed = TRUE)) {
  stopifnot(ridge$patient_bootstrap_brier_ci_low <= 0, ridge$patient_bootstrap_brier_ci_high >= 0)
}
if (grepl("The interval sits above zero and remains well below 0.01", supplement_text, fixed = TRUE)) {
  stopifnot(control$patient_bootstrap_ci_low > 0, control$patient_bootstrap_ci_high < .01)
}
if (grepl("Both left the mean three-year Brier score worse", supplement_text, fixed = TRUE)) {
  stopifnot(all(benchmark$mean_delta_brier3[benchmark$strategy %in% c("survival_only", "de_only")] > 0))
}
if (grepl("the latter was farther from one", supplement_text, fixed = TRUE)) {
  slopes <- benchmark$mean_calibration_slope[match(c("clinical", "hydra"), benchmark$strategy)]
  stopifnot(abs(slopes[2] - 1) > abs(slopes[1] - 1))
}
message("PASS: candidate ledger, reconstructed selection gates, BH families, source gene sets, JSON evidence and manuscript funnel macros.")
