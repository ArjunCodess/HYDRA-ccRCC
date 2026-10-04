# Post hoc threat-specific sensitivities; never redefine the discovery shortlist.
source("analysis/00_config.R")
source("analysis/functions/io.R")
source("analysis/functions/tcga_metadata.R")
source("analysis/functions/patient_samples.R")
suppressPackageStartupMessages({
  library(Biobase); library(readxl); library(AnnotationDbi); library(org.Hs.eg.db)
  library(readr); library(dplyr); library(tibble); library(survival)
})
set.seed(RESAMPLING$seed + 46L)
priority <- read_csv(file.path(DIRS$tables, "high_confidence_candidate_genes.csv"), show_col_types = FALSE)
stopifnot(nrow(priority) == 23L, !anyDuplicated(priority$symbol))
audit <- read_selected_tcga_coldata(FILES$tcga_coldata, FILES$tcga_counts) |>
  filter(sample_type == "Primary Tumor")
stopifnot(nrow(audit) == 533L, !anyDuplicated(audit$patient_barcode))
clinical_raw <- read_csv(FILES$tcga_clinical, show_col_types = FALSE)
tcga <- audit |> select(sample = sample_barcode, patient = patient_barcode) |>
  inner_join(clinical_raw, by = c("patient" = "submitter_id")) |>
  transmute(sample, patient, os_time, os_event,
    age = as.numeric(age_at_diagnosis) / 365.25, sex = factor(gender),
    stage = factor(normalize_stage(ajcc_pathologic_stage), levels = paste("Stage", c("I", "II", "III", "IV"))),
    grade_high = case_when(normalize_grade(tumor_grade) %in% c("G1", "G2") ~ 0,
                          normalize_grade(tumor_grade) %in% c("G3", "G4") ~ 1, TRUE ~ NA_real_),
    t_high = case_when(grepl("^T[12]", ajcc_pathologic_t) ~ 0,
                       grepl("^T[34]", ajcc_pathologic_t) ~ 1, TRUE ~ NA_real_)) |>
  filter(is.finite(os_time), os_time > 0, os_event %in% 0:1)
stopifnot(nrow(tcga) == 529L,
          sum(complete.cases(tcga[, c("age", "sex", "stage", "grade_high")])) == 517L)
vst <- read_required_rds(FILES$tcga_vst)
tcga_expr <- vst[priority$tcga_gene_id, tcga$sample, drop = FALSE]
rownames(tcga_expr) <- priority$symbol
rm(vst)

gse <- read_required_rds(file.path(DIRS$processed, "gse29609_series_matrix.rds"))
pd <- pData(gse); features <- fData(gse)
symbol_column <- names(features)[grepl("gene symbol", names(features), ignore.case = TRUE)][1]
stopifnot(!is.na(symbol_column))
symbols <- trimws(sub(" /// .*| // .*|;.*|,.*", "", as.character(features[[symbol_column]])))
gse_expr <- as.data.frame(exprs(gse)) |> mutate(symbol = symbols) |>
  filter(symbol %in% priority$symbol) |> group_by(symbol) |>
  summarise(across(where(is.numeric), mean), .groups = "drop") |> column_to_rownames("symbol") |> as.matrix()
gse_clin <- tibble(sample = rownames(pd), patient = rownames(pd),
  os_time = as.numeric(pd[["survival time:ch1"]]), os_event = as.numeric(pd[["death (1=yes, 0=no):ch1"]]),
  age = as.numeric(pd[["age at diagnosis (y):ch1"]]),
  grade_high = as.numeric(pd[["fuhrman grade:ch1"]]) >= 3,
  t_high = as.numeric(pd[["t (tnm stage):ch1"]]) >= 3) |>
  filter(is.finite(os_time), os_time > 0, os_event %in% 0:1)

em_path <- file.path(DIRS$raw, "arrayexpress", "E-MTAB-1980")
em_raw <- read_excel(file.path(em_path, "supplementary_table_1.xlsx"), skip = 1)
em_clin <- em_raw |> transmute(sample = as.character(.data[["sample ID"]]), patient = sample,
  age = as.numeric(Age), sex = factor(toupper(Sex)),
  os_time = as.numeric(.data[["observation period (month)"]]),
  os_event = as.integer(tolower(outcome) == "dead"),
  grade_high = as.numeric(.data[["Fuhrman grade"]]) >= 3,
  stage_text = as.character(.data[["Stage at diagnosis"]]),
  t_high = case_when(grepl("^pT[12]", stage_text) ~ 0,
                     grepl("^pT[34]", stage_text) ~ 1, TRUE ~ NA_real_),
  molecular_subtype = .data[["gene expression profile"]]) |>
  filter(molecular_subtype %in% c("ccA", "ccB"), is.finite(os_time), os_time > 0, os_event %in% 0:1)
em_matrix <- read_tsv(file.path(em_path, "ccRCC_exp_log_quantile_normalized.txt"), show_col_types = FALSE, progress = FALSE)
em_matrix$refseq <- sub("\\.\\d+$", "", em_matrix$SystematicName)
mapping <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(em_matrix$refseq[!is.na(em_matrix$refseq)]),
                                keytype = "REFSEQ", columns = "SYMBOL") |>
  filter(SYMBOL %in% priority$symbol) |> distinct(REFSEQ, SYMBOL)
em_expr <- em_matrix |> inner_join(mapping, by = c("refseq" = "REFSEQ")) |>
  select(SYMBOL, all_of(em_clin$sample)) |> mutate(across(-SYMBOL, as.numeric)) |>
  group_by(SYMBOL) |> summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop") |>
  column_to_rownames("SYMBOL") |> as.matrix()
cohorts <- list("TCGA-KIRC" = list(clinical = tcga, expr = tcga_expr),
                "GSE29609" = list(clinical = gse_clin, expr = gse_expr),
                "E-MTAB-1980" = list(clinical = em_clin, expr = em_expr))
availability <- bind_rows(lapply(names(cohorts), function(cohort) {
  dat <- cohorts[[cohort]]$clinical
  bind_rows(lapply(c("age", "sex", "grade_high", "t_high"), function(variable) {
    available <- variable %in% names(dat)
    tibble(cohort, variable, n = nrow(dat), available,
           missing = if (available) sum(is.na(dat[[variable]])) else nrow(dat))
  }))
}))
write_csv_atomic(availability, file.path(DIRS$tables, "limitations_covariate_availability.csv"))

fit_record <- function(dat, formula, symbol, cohort, model) {
  vars <- all.vars(formula)
  dat <- dat[complete.cases(dat[, vars, drop = FALSE]), , drop = FALSE]
  if (anyDuplicated(dat$patient)) stop("Duplicate patients in sensitivity model.")
  warning_text <- character()
  fit <- withCallingHandlers(tryCatch(coxph(formula, data = dat), error = identity),
    warning = function(w) {warning_text <<- c(warning_text, conditionMessage(w)); invokeRestart("muffleWarning")})
  if (inherits(fit, "error")) return(tibble(symbol, cohort, model, n = nrow(dat), events = sum(dat$os_event),
    log_hr = NA_real_, se = NA_real_, hr = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
    p_value = NA_real_, ph_p_value = NA_real_,
    status = conditionMessage(fit), warning = paste(warning_text, collapse = "; ")))
  beta <- unname(coef(fit)["expr"]); se <- sqrt(vcov(fit)["expr", "expr"])
  ph <- tryCatch(cox.zph(fit)$table["expr", "p"], error = function(e) NA_real_)
  tibble(symbol, cohort, model, n = fit$n, events = fit$nevent,
         log_hr = beta, se, hr = exp(beta), ci_low = exp(beta - 1.96*se), ci_high = exp(beta + 1.96*se),
         p_value = 2*pnorm(-abs(beta/se)), ph_p_value = ph, status = ifelse(is.finite(beta+se), "fit", "nonfinite"),
         warning = paste(unique(warning_text), collapse = "; "))
}
harmonized <- bind_rows(lapply(names(cohorts), function(cohort) {
  dat <- cohorts[[cohort]]$clinical; expr <- cohorts[[cohort]]$expr
  bind_rows(lapply(intersect(priority$symbol, rownames(expr)), function(symbol) {
    # Standardization before model-specific exclusions matches the existing association convention.
    dat$expr <- as.numeric(scale(expr[symbol, dat$sample]))
    bind_rows(fit_record(dat, Surv(os_time, os_event) ~ expr, symbol, cohort, "unadjusted"),
              fit_record(dat, Surv(os_time, os_event) ~ expr + age + grade_high, symbol, cohort, "age_grade"),
              fit_record(dat, Surv(os_time, os_event) ~ expr + age + grade_high + t_high, symbol, cohort, "age_grade_t"))
  }))
})) |> group_by(cohort, model) |> mutate(fdr = p.adjust(p_value, "BH", n = nrow(priority))) |> ungroup() |>
  left_join(priority |> select(symbol, main_log_hr), by = "symbol") |>
  mutate(same_primary_direction = sign(log_hr) == sign(main_log_hr))
write_csv_atomic(harmonized, file.path(DIRS$tables, "limitations_harmonized_cox.csv"))
heterogeneity <- harmonized |> filter(status == "fit", is.finite(se), se > 0) |>
  group_by(symbol, model) |> group_modify(function(dat, key) {
    w <- 1/dat$se^2; center <- sum(w*dat$log_hr)/sum(w)
    q <- sum(w*(dat$log_hr-center)^2); df <- nrow(dat)-1L
    tibble(cohorts = nrow(dat), q, df, p_value = if (df > 0) pchisq(q, df, lower.tail = FALSE) else NA_real_,
           i2 = if (q > 0 && df > 0) max(0, (q-df)/q) else 0,
           scope = "heterogeneity of cohort-standardized coefficients; post-selection exploratory, not a pooled biological effect")
  }) |> ungroup() |> group_by(model) |> mutate(fdr = p.adjust(p_value, "BH", n = nrow(priority))) |> ungroup()
write_csv_atomic(heterogeneity, file.path(DIRS$tables, "limitations_heterogeneity.csv"))

# Eight explicit completion scenarios assess covariate exclusion. They are not
# multiple imputation, mathematical worst-case bounds, or a missingness model.
scenarios <- expand.grid(grade_value = 0:1, stage_value = c("Stage I", "Stage IV"), age_quantile = c(.1, .9))
missing_results <- bind_rows(lapply(seq_len(nrow(scenarios)), function(i) {
  scenario <- scenarios[i, ]; dat <- tcga
  dat$age[is.na(dat$age)] <- quantile(dat$age, scenario$age_quantile, na.rm = TRUE)
  dat$grade_high[is.na(dat$grade_high)] <- scenario$grade_value
  dat$stage[is.na(dat$stage)] <- scenario$stage_value
  bind_rows(lapply(priority$symbol, function(symbol) {
    dat$expr <- as.numeric(scale(tcga_expr[symbol, dat$sample]))
    fit_record(dat, Surv(os_time, os_event) ~ expr + age + sex + stage + grade_high,
               symbol, "TCGA-KIRC", sprintf("completion_%02d", i))
  })) |> mutate(grade_completion = scenario$grade_value, stage_completion = scenario$stage_value,
                age_quantile = scenario$age_quantile)
})) |> group_by(model) |> mutate(fdr = p.adjust(p_value, "BH")) |> ungroup() |>
  left_join(priority |> select(symbol, main_log_hr), by = "symbol") |>
  mutate(same_primary_direction = sign(log_hr) == sign(main_log_hr))
write_csv_atomic(missing_results, file.path(DIRS$tables, "limitations_missing_covariate_scenarios.csv"))
missing_summary <- missing_results |> group_by(symbol) |>
  summarise(scenarios = n(), all_same_direction = all(same_primary_direction),
    all_fdr_supported = all(fdr < .05), min_log_hr = min(log_hr), max_log_hr = max(log_hr),
    max_fdr = max(fdr), .groups = "drop")
write_csv_atomic(missing_summary, file.path(DIRS$tables, "limitations_missing_covariate_summary.csv"))

# A collection-site sensitivity addresses site-associated clinical/expression
# differences. Barcode sites are not technical batches or measured confounders.
site_data <- tcga |>
  filter(complete.cases(age, sex, stage, grade_high)) |>
  mutate(source_site = substr(patient, 6L, 7L))
site_counts <- site_data |> count(source_site, name = "n")
large_sites <- site_counts$source_site[site_counts$n >= 10L]
site_data$site_group <- ifelse(site_data$source_site %in% large_sites, site_data$source_site, "small_sites")
site_result <- bind_rows(lapply(priority$symbol, function(symbol) {
  dat <- site_data
  dat$expr <- as.numeric(scale(tcga_expr[symbol, dat$sample]))
  fit_record(dat, Surv(os_time, os_event) ~ expr + age + sex + stage + grade_high + strata(site_group),
             symbol, "TCGA-KIRC", "source_site_stratified")
})) |> mutate(fdr = p.adjust(p_value, "BH", n = nrow(priority))) |>
  left_join(priority |> select(symbol, main_log_hr), by = "symbol") |>
  mutate(same_primary_direction = sign(log_hr) == sign(main_log_hr),
         site_strata = length(unique(site_data$site_group)),
         scope = "post-selection collection-site sensitivity; not technical batch correction")
write_csv_atomic(site_result, file.path(DIRS$tables, "limitations_source_site_cox.csv"))
write_csv_atomic(site_data |> group_by(site_group) |> summarise(n = n(), events = sum(os_event), .groups = "drop"),
                 file.path(DIRS$tables, "limitations_source_site_strata.csv"))

# Shared accession IDs can detect documented reuse. Disjoint ID namespaces do
# not establish patient independence, so both kinds of evidence are separate.
geo_sets <- lapply(c("gse40435", "gse53757", "gse29609"), function(id) {
  x <- read_required_rds(file.path(DIRS$processed, paste0(id, "_series_matrix.rds")))
  if (is.list(x)) x <- x[[1L]]
  data.frame(cohort = toupper(id), sample = rownames(pData(x)), stringsAsFactors = FALSE)
})
sample_ids <- bind_rows(geo_sets, tibble(cohort = "TCGA-KIRC", sample = audit$sample_barcode),
                        tibble(cohort = "E-MTAB-1980", sample = em_clin$sample))
pairs <- combn(unique(sample_ids$cohort), 2, simplify = FALSE)
overlap <- bind_rows(lapply(pairs, function(pair) {
  a <- sample_ids$sample[sample_ids$cohort == pair[1]]; b <- sample_ids$sample[sample_ids$cohort == pair[2]]
  shared <- intersect(a,b)
  tibble(cohort_a = pair[1], cohort_b = pair[2], shared_accession_ids = length(shared),
         shared_ids = paste(shared, collapse = ";"),
         patient_independence = "unverified; accession-ID screen cannot exclude different identifiers for the same patient")
}))
write_csv_atomic(overlap, file.path(DIRS$tables, "limitations_accession_overlap.csv"))
stopifnot(all(harmonized$status == "fit"), all(missing_results$status == "fit"), all(site_result$status == "fit"))
support <- function(cohort, model) sum(harmonized$cohort == cohort & harmonized$model == model &
  harmonized$fdr < .05 & harmonized$same_primary_direction, na.rm = TRUE)
macro <- function(name, value) paste0("\\newcommand{\\", name, "}{", value, "}")
writeLines(c("% Generated by analysis/46_limitations_sensitivity.R.",
  macro("LimitEmAgeGrade", support("E-MTAB-1980", "age_grade")),
  macro("LimitEmAgeGradeT", support("E-MTAB-1980", "age_grade_t")),
  macro("LimitMissingRetained", sum(missing_summary$all_same_direction & missing_summary$all_fdr_supported)),
  macro("LimitSiteRetained", sum(site_result$fdr < .05 & site_result$same_primary_direction)),
  macro("LimitSiteStrata", length(unique(site_data$site_group)))), "paper/limitations_macros.tex")
summary_rows <- harmonized |> group_by(cohort, model) |>
  summarise(mapped = n(), n = min(n), events = min(events),
    same_fdr = sum(fdr < .05 & same_primary_direction), reversed_fdr = sum(fdr < .05 & !same_primary_direction),
    ph_flags = sum(ph_p_value < .05, na.rm = TRUE), .groups = "drop")
table_rows <- apply(as.data.frame(summary_rows), 1, function(row) paste0(paste(row, collapse = " & "), " \\\\"))
table_rows <- gsub("age_grade_t", "age + grade + T", table_rows, fixed = TRUE)
table_rows <- gsub("age_grade", "age + grade", table_rows, fixed = TRUE)
writeLines(c("% Generated by analysis/46_limitations_sensitivity.R.",
  "\\begin{table}[htbp]\\centering\\small\\setlength{\\tabcolsep}{3pt}",
  "\\caption{Exploratory harmonized adjustment of the fixed shortlist. Each cohort/model has a 23-hypothesis BH family, including unmapped genes. Same and opposite count FDR-supported effects relative to the primary TCGA direction; N counts patients. PH flags are nominal expression diagnostics, not exclusion gates. Grade and T are collapsed to low/high categories. These post-selection comparisons do not establish transportability.}",
  "\\label{tab:harmonized}\\begin{tabular}{llrrrrrr}\\toprule",
  "Cohort & Model & Mapped & N & Events & Same & Opposite & PH \\\\",
  "\\midrule", table_rows, "\\bottomrule\\end{tabular}\\end{table}"), "paper/limitations_table.tex")
writeLines(capture.output(sessionInfo()), "environment/limitations_sensitivity_sessionInfo.txt")
message("Threat-specific harmonized, heterogeneity, missing-covariate and accession-overlap analyses complete.")
