source("analysis/00_config.R")
suppressPackageStartupMessages({library(readr); library(dplyr)})
tab <- function(name) read_csv(file.path(DIRS$tables, paste0(name, ".csv")), show_col_types = FALSE)
priority <- tab("high_confidence_candidate_genes")
h <- tab("limitations_harmonized_cox")
stopifnot(nrow(h) == 198L, all(h$status == "fit"), all(is.finite(h$log_hr)),
          all(h$ci_low <= h$hr & h$ci_high >= h$hr), all(h$symbol %in% priority$symbol))
for (cohort in unique(h$cohort)) for (model in unique(h$model)) {
  rows <- h |> filter(.data$cohort == .env$cohort, .data$model == .env$model)
  stopifnot(isTRUE(all.equal(rows$fdr, p.adjust(rows$p_value, "BH", n = 23), tolerance = 1e-12)))
}
# Independently reproduce the recorded unadjusted external effects; these
# sensitivities must not silently change probe collapsing, sample alignment or coding.
for (spec in list(c("GSE29609", "external_survival_gse29609"), c("E-MTAB-1980", "external_survival_emtab1980"))) {
  original <- tab(spec[2])
  compared <- h |> filter(cohort == spec[1], model == "unadjusted") |>
    inner_join(original |> select(symbol, external_log_hr), by = "symbol")
  stopifnot(max(abs(compared$log_hr - compared$external_log_hr)) < 1e-8)
}
m <- tab("limitations_missing_covariate_scenarios")
s <- tab("limitations_source_site_cox")
stopifnot(nrow(m) == 184L, all(m$n == 529L), all(m$events == 173L), all(m$status == "fit"),
          nrow(s) == 23L, all(s$n == 517L), all(s$events == 170L), all(s$status == "fit"))
strata <- tab("limitations_source_site_strata")
stopifnot(sum(strata$n) == 517L, sum(strata$events) == 170L, nrow(strata) == unique(s$site_strata))
q <- tab("limitations_heterogeneity")
for (symbol in unique(q$symbol)) for (model in unique(q$model)) {
  rows <- h |> filter(.data$symbol == .env$symbol, .data$model == .env$model)
  weights <- 1/rows$se^2
  expected <- sum(weights * (rows$log_hr - sum(weights*rows$log_hr)/sum(weights))^2)
  result <- q |> filter(.data$symbol == .env$symbol, .data$model == .env$model)
  stopifnot(nrow(result) == 1L, abs(expected - result$q) < 1e-10)
}
overlap <- tab("limitations_accession_overlap")
stopifnot(nrow(overlap) == 10L, all(overlap$shared_accession_ids >= 0L),
          all(grepl("unverified", overlap$patient_independence)))
macro_lines <- readLines("paper/limitations_macros.tex")
assert_macro <- function(name, value) stopifnot(paste0("\\newcommand{\\", name, "}{", value, "}") %in% macro_lines)
assert_macro("LimitEmAgeGrade", sum(h$cohort == "E-MTAB-1980" & h$model == "age_grade" & h$fdr < .05 & h$same_primary_direction))
assert_macro("LimitEmAgeGradeT", sum(h$cohort == "E-MTAB-1980" & h$model == "age_grade_t" & h$fdr < .05 & h$same_primary_direction))
ms <- tab("limitations_missing_covariate_summary")
assert_macro("LimitMissingRetained", sum(ms$all_same_direction & ms$all_fdr_supported))
assert_macro("LimitSiteRetained", sum(s$fdr < .05 & s$same_primary_direction))
reversals <- h |> filter(cohort == "GSE29609", model != "unadjusted", symbol %in% c("DDC", "TCIRG1"))
stopifnot(nrow(reversals) == 4L, all(!reversals$same_primary_direction), all(reversals$fdr >= .05))
acadm_adjusted <- h |> filter(cohort == "E-MTAB-1980", model == "age_grade_t", symbol == "ACADM")
stopifnot(nrow(acadm_adjusted) == 1L, acadm_adjusted$same_primary_direction, acadm_adjusted$fdr >= .05)
audit <- tab("normalization_primary_correction_audit")
stopifnot(nrow(audit) == 50L, all(audit$original_prediction_max_error < 1e-9),
          all(audit$training_expression_max_change == 0), any(audit$max_gene_lp_change > 0))
message("PASS: harmonized estimates reproduce original cohorts; BH, complete-case populations, Q, scenario coverage, macros and normalization provenance.")
