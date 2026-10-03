source("analysis/00_config.R")
source("analysis/functions/io.R")
suppressPackageStartupMessages({ library(dplyr); library(readr); library(jsonlite) })

# This script joins recorded results. It does not fit new models or add gates.
tab <- function(name) read_csv(file.path(DIRS$tables, paste0(name, ".csv")), show_col_types = FALSE)
repro <- tab("reproducible_deg_tcga_gse40435_gse53757")
candidate <- tab("candidate_gene_evidence_table")
cox <- tab("tcga_kirc_cox_models")
stopifnot(!anyDuplicated(repro$symbol), !anyDuplicated(repro$tcga_gene_id))
stopifnot(setequal(repro$tcga_gene_id, candidate$tcga_gene_id))

# Recheck the Cox multiplicity family separately for each recorded model.
for (model in unique(cox$model_type)) {
  rows <- cox$model_type == model
  stopifnot(isTRUE(all.equal(cox$fdr[rows], p.adjust(cox$p_value[rows], "BH"), tolerance = 1e-12)))
}

# Each JSON field preserves the complete source row(s), including unavailable
# estimates as null. A missing assessment is a labeled object rather than [].
pack <- function(data, key = "symbol") {
  pieces <- split(data, data[[key]])
  vapply(pieces, function(x) as.character(toJSON(x, dataframe = "rows", na = "null", digits = NA)), character(1))
}
attach_json <- function(ledger, name, data, key = "symbol", ledger_key = "symbol") {
  lookup <- pack(data, key)
  value <- unname(lookup[ledger[[ledger_key]]])
  value[is.na(value)] <- '{"status":"not_assessed"}'
  ledger[[name]] <- value
  ledger
}
ledger <- candidate |>
  transmute(gene_id = tcga_gene_id, symbol,
    cohort = "TCGA-KIRC", model = "stage_grade_complete Cox; standardized tumor expression",
    discovery_log2fc = tcga_log2fc, discovery_adjusted_p = tcga_padj,
    discovery_de = tcga_significant, reproducible_de = reproducible_deg,
    adjusted_survival = reproducible_deg & prognostic,
    covariate_sensitivity = reproducible_deg & prognostic & sensitivity_pass,
    strict_candidate, priority_shortlist = high_confidence_candidate,
    effect_direction = case_when(is.na(main_log_hr) ~ "not_assessed", main_log_hr > 0 ~ "higher_hazard", main_log_hr < 0 ~ "lower_hazard", TRUE ~ "zero"),
    log_hr = main_log_hr, log_hr_ci_low = main_log_hr_ci_low, log_hr_ci_high = main_log_hr_ci_high,
    hr = main_hr, hr_ci_low = main_hr_ci_low, hr_ci_high = main_hr_ci_high,
    p_value = main_p_value, adjusted_p = main_fdr,
    multiplicity_adjustment = if_else(is.na(main_p_value), "not_assessed", "BH within recorded Cox model over reproducible DE genes"),
    n = main_n, events = main_events, ph_p_value = main_ph_p_value,
    ph_check = case_when(is.na(main_ph_p_value) ~ "not_assessed", main_ph_p_value < THRESHOLDS$ph_min_p ~ "nominal_diagnostic_flag", TRUE ~ "diagnostic_not_flagged"),
    final_status = case_when(high_confidence_candidate ~ "priority_shortlist_requires_external_interpretation",
      strict_candidate ~ "strict_discovery_candidate", !tcga_significant ~ "not_discovery_de",
      !reproducible_deg ~ "not_reproducible_de", !prognostic ~ "not_adjusted_survival",
      !sensitivity_pass ~ "not_covariate_sensitivity", TRUE ~ "not_strict_effect_thresholds")) |>
  arrange(symbol)
ledger <- attach_json(ledger, "de_evidence_json", repro)
ledger <- attach_json(ledger, "cox_models_json", cox, "gene_id", "gene_id")
sources <- c(external_gse29609_json = "external_survival_gse29609", external_emtab1980_json = "external_survival_emtab1980",
  bootstrap_json = "candidate_cox_bootstrap_summary", conditional_prediction_json = "candidate_cv_clinical_increment",
  purity_sensitivity_json = "candidate_direct_tumor_purity_sensitivity", composition_sensitivity_json = "candidate_clinical_composition_sensitivity",
  cell_source_json = "hpa_candidate_cell_source_summary", evidence_profile_json = "candidate_evidence_matrix")
for (name in names(sources)) ledger <- attach_json(ledger, name, tab(sources[[name]]))
external <- tab("external_survival_emtab1980")
ledger$external_replication <- ifelse(ledger$symbol %in% external$symbol[external$external_strict_support %in% TRUE],
  "E-MTAB-1980 same-direction FDR support; see cohort-specific evidence", ifelse(ledger$priority_shortlist, "mixed_or_unavailable; see cohort-specific evidence", "not_assessed"))
ledger$purity_assessment <- ifelse(ledger$purity_sensitivity_json == '{"status":"not_assessed"}', "not_assessed", "assessed; see purity_sensitivity_json")
ledger$cell_source_assessment <- ifelse(ledger$cell_source_json == '{"status":"not_assessed"}', "not_assessed", "normal-tissue HPA reference; see cell_source_json")
write_csv_atomic(ledger, file.path(DIRS$tables, "candidate_ledger.csv"))

stages <- c("mapped_tcga", "discovery_de", "reproducible_de", "adjusted_survival", "covariate_sensitivity", "strict_candidate", "priority_shortlist")
counts <- c(nrow(ledger), vapply(stages[-1], function(x) sum(ledger[[x]], na.rm = TRUE), numeric(1)))
funnel <- tibble(stage = stages, count = counts, previous_count = c(NA, head(counts, -1)),
  removed_from_previous = c(NA, -diff(counts)), role = "serial_discovery_selection")
write_csv_atomic(funnel, file.path(DIRS$tables, "ledger_funnel_counts.csv"))

# Procedure comparisons use the independently generated matched-list table.
procedures <- tab("funnel_external_summary")
procedure_value <- function(rule, cohort, field) {
  row <- procedures[procedures$rule == rule & procedures$cohort == cohort, ]
  stopifnot(nrow(row) == 1L)
  row[[field]]
}
em <- tab("external_survival_emtab1980")
gse <- tab("external_survival_gse29609")
marker <- tab("candidate_clinical_composition_sensitivity")
purity <- tab("candidate_direct_tumor_purity_sensitivity")
bootstrap <- tab("candidate_cox_bootstrap_summary")
priority <- ledger[ledger$priority_shortlist, ]
review_values <- c(
  ReviewEmDeMapped = procedure_value("de_only", "E-MTAB-1980", "present"),
  ReviewEmDeSame = procedure_value("de_only", "E-MTAB-1980", "same_direction"),
  ReviewEmCompleteMapped = procedure_value("complete_rule", "E-MTAB-1980", "present"),
  ReviewEmCompleteSame = procedure_value("complete_rule", "E-MTAB-1980", "same_direction"),
  ReviewGseDeMapped = procedure_value("de_only", "GSE29609", "present"),
  ReviewGseDeSame = procedure_value("de_only", "GSE29609", "same_direction"),
  ReviewGseCompleteMapped = procedure_value("complete_rule", "GSE29609", "present"),
  ReviewGseCompleteSame = procedure_value("complete_rule", "GSE29609", "same_direction"),
  ReviewEmDeFdr = procedure_value("de_only", "E-MTAB-1980", "same_direction_fdr"),
  ReviewEmCompleteFdr = procedure_value("complete_rule", "E-MTAB-1980", "same_direction_fdr"),
  ReviewMarkerRetained = sum(marker$composition_adjusted_fdr < 0.05, na.rm = TRUE),
  ReviewPhFlags = sum(priority$ph_check == "nominal_diagnostic_flag"),
  ReviewEmMapped = sum(em$external_present),
  ReviewEmSame = sum(em$external_present & em$external_same_direction, na.rm = TRUE),
  ReviewEmFdr = sum(em$external_strict_support, na.rm = TRUE),
  ReviewGseMapped = sum(gse$external_present),
  ReviewGseSame = sum(gse$external_present & gse$external_same_direction, na.rm = TRUE),
  ReviewBootstrapIntervals = sum(bootstrap$ci_excludes_zero, na.rm = TRUE),
  ReviewPurityRetained = sum(purity$gene_fdr < 0.05, na.rm = TRUE))
annotations <- tibble(stage = names(review_values), count = unname(review_values), previous_count = NA_real_,
  removed_from_previous = NA_real_, role = "annotation_not_serial_gate")
write_csv_atomic(bind_rows(funnel, annotations), file.path(DIRS$tables, "ledger_funnel_counts.csv"))
macro <- function(name, value) paste0("\\newcommand{\\", name, "}{", value, "}")
review_macros <- c("% Generated by analysis/39_candidate_ledger.R; do not edit.",
  vapply(seq_along(review_values), function(i) macro(names(review_values)[i], review_values[i]), character(1)),
  macro("ReviewDiscoveryAttritionPct", sprintf("%.2f", 100 * (1 - nrow(priority) / sum(ledger$discovery_de)))),
  macro("ReviewMappedAttritionPct", sprintf("%.2f", 100 * (1 - nrow(priority) / nrow(ledger)))))
writeLines(review_macros, "paper/review_macros.tex")

# A compact complete shortlist table is regenerated from the ledger's evidence.
fmt <- function(x, digits = 2) ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
qfmt <- function(x) ifelse(is.na(x), "NA", formatC(x, digits = 2, format = "g"))
rows <- vapply(seq_len(nrow(priority)), function(i) {
  gene <- priority$symbol[i]
  e <- em[match(gene, em$symbol), ]; g <- gse[match(gene, gse$symbol), ]; m <- marker[match(gene, marker$symbol), ]
  paste0(paste(c(gene, paste0(fmt(priority$hr[i]), " [", fmt(priority$hr_ci_low[i]), ", ", fmt(priority$hr_ci_high[i]), "]"),
    qfmt(priority$adjusted_p[i]), fmt(e$external_hr), qfmt(e$external_fdr), fmt(g$external_hr),
    qfmt(g$external_fdr), qfmt(m$composition_adjusted_fdr), qfmt(priority$ph_p_value[i])), collapse = " & "), " \\\\")
}, character(1))
writeLines(c("% Generated by analysis/39_candidate_ledger.R; do not edit.",
  "\\begin{table}[htbp]", "\\centering\\scriptsize",
  "\\caption{Complete discovery shortlist. TCGA hazard ratios and 95\\% intervals are adjusted for age, sex, stage and grade. External hazard ratios are unadjusted. The marker column reports the composition-adjusted TCGA FDR; PH is the primary TCGA proportional-hazards diagnostic p-value. NA denotes unavailable evidence. None of these columns establishes clinical utility.}",
  "\\label{tab:review-candidates}", "\\resizebox{\\linewidth}{!}{\\begin{tabular}{lrrrrrrrr}", "\\toprule",
  "Gene & TCGA HR [95\\% CI] & TCGA FDR & E-MTAB HR & E-MTAB FDR & GSE29609 HR & GSE29609 FDR & Marker FDR & PH $p$ \\\\",
  "\\midrule", rows, "\\bottomrule", "\\end{tabular}}", "\\end{table}"), "paper/review_candidate_table.tex")
message("Candidate ledger: ", nrow(ledger), " mapped genes; ", tail(counts, 1), " priority genes. External, bootstrap, prediction, purity and cell source are annotations, not serial gates.")
