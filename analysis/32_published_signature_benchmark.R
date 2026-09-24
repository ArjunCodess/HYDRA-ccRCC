# Published-signature comparison. Does not replace the one-gene acceptance test.
# ClearCode34 symbols and subtypes are Table 1 of Brooks et al., Eur Urol 2014.
# ccA genes are signed -1 and ccB genes +1 so a higher score tracks the poor-prognosis subtype.

source("analysis/00_config.R")
source("analysis/functions/io.R")
source("analysis/functions/tcga_metadata.R")
source("analysis/functions/patient_samples.R")
source("analysis/functions/nested_benchmark.R")

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(dplyr)
  library(readr)
  library(readxl)
  library(stringr)
  library(tibble)
  library(survival)
})

published <- read_csv("analysis/clearcode34_genes.csv", show_col_types = FALSE)
stopifnot(nrow(published) == 34L, !anyDuplicated(published$published_symbol))

map_symbols <- function(published) {
  alias <- AnnotationDbi::select(org.Hs.eg.db, keys = published$published_symbol,
                                 keytype = "ALIAS", columns = c("ALIAS", "SYMBOL")) |>
    filter(!is.na(SYMBOL), SYMBOL != ALIAS) |>
    distinct(ALIAS, SYMBOL)
  if (anyDuplicated(alias$ALIAS)) stop("A ClearCode34 alias maps to more than one current symbol.")
  published |>
    left_join(alias, by = c("published_symbol" = "ALIAS")) |>
    mutate(mapped_symbol = if_else(!is.na(SYMBOL), SYMBOL, published_symbol)) |>
    select(-SYMBOL)
}

score_block <- function(pred, lp_col, risk_col, horizon) {
  ci <- function(lp) unname(concordance(Surv(pred$os_time, pred$os_event) ~ lp, reverse = TRUE)$concordance)
  cens <- survfit(Surv(os_time, 1L - os_event) ~ 1, data = pred)
  censor_surv <- function(t) {
    i <- findInterval(t, cens$time)
    ifelse(i == 0L, 1, cens$surv[pmax(i, 1L)])
  }
  target <- as.numeric(pred$os_time <= horizon & pred$os_event == 1L)
  weight <- ifelse(target == 1, 1 / pmax(censor_surv(pmax(pred$os_time - 1, 0)), 0.01),
                   ifelse(pred$os_time > horizon, 1 / pmax(censor_surv(horizon), 0.01), 0))
  brier <- function(risk) mean(weight * (target - risk)^2)
  slope <- function(lp) tryCatch(unname(coef(coxph(Surv(pred$os_time, pred$os_event) ~ lp, data = pred))[1]),
                                 error = function(e) NA_real_)
  c(c_index = ci(pred[[lp_col]]), brier = brier(pred[[risk_col]]), slope = slope(pred[[lp_col]]))
}

risk_from_fit <- function(fit, newdata, horizon) {
  lp <- as.numeric(predict(fit, newdata = newdata, type = "lp", reference = "zero"))
  bh <- basehaz(fit, centered = FALSE)
  h <- if (any(bh$time <= horizon)) tail(bh$hazard[bh$time <= horizon], 1) else 0
  list(lp = lp, risk = 1 - exp(-h * exp(lp)))
}

bootstrap_delta <- function(pred, lp_col, risk_col, clinical_lp, clinical_risk, horizon, draws = 1000L) {
  ids <- unique(pred$patient_barcode)
  deltas <- vapply(seq_len(draws), function(i) {
    sampled <- sample(ids, length(ids), replace = TRUE)
    weights <- as.integer(table(factor(sampled, levels = ids)))
    x <- pred[rep(seq_len(nrow(pred)), weights[match(pred$patient_barcode, ids)]), , drop = FALSE]
    sig <- score_block(x, lp_col, risk_col, horizon)
    clin <- score_block(x, clinical_lp, clinical_risk, horizon)
    unname(sig["c_index"] - clin["c_index"])
  }, numeric(1))
  quantile(deltas, c(0.025, 0.975))
}

mapped <- map_symbols(published)
write_csv_atomic(mapped, file.path(DIRS$tables, "clearcode34_symbol_map.csv"))
if (sum(!is.na(mapped$mapped_symbol)) < 30L) {
  stop("Fewer than 30 ClearCode34 genes mapped. Refusing to score a partial list.")
}

se <- read_required_rds(FILES$tcga_se)
counts <- assay(se, "unstranded")
meta <- as.data.frame(colData(se)) |>
  tibble::rownames_to_column("sample_barcode") |>
  mutate(patient_barcode = substr(sample_barcode, 1, 12),
         condition = factor(shortLetterCode, levels = c("NT", "TP")))
meta <- select_tcga_patient_samples(meta, counts)
clinical <- read_csv(FILES$tcga_clinical, show_col_types = FALSE) |>
  transmute(patient_barcode = submitter_id, os_time, os_event,
            age = suppressWarnings(as.numeric(age_at_diagnosis)) / 365.25,
            sex = factor(gender), stage = factor(normalize_stage(ajcc_pathologic_stage)),
            grade = factor(case_when(normalize_grade(tumor_grade) %in% c("G1", "G2") ~ "Low grade",
                                     normalize_grade(tumor_grade) %in% c("G3", "G4") ~ "High grade",
                                     TRUE ~ NA_character_)))
patients <- meta |>
  filter(condition == "TP") |>
  inner_join(clinical, by = "patient_barcode") |>
  filter(is.finite(os_time), os_time > 0, os_event %in% c(0L, 1L),
         is.finite(age), !is.na(sex), !is.na(stage), !is.na(grade)) |>
  arrange(patient_barcode)
symbols <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(sub("\\..*$", "", rownames(counts))),
                                 keytype = "ENSEMBL", columns = "SYMBOL") |>
  filter(!is.na(SYMBOL)) |> arrange(ENSEMBL, SYMBOL) |> distinct(ENSEMBL, .keep_all = TRUE)
gene_rows <- tibble(gene_id = rownames(counts), ensembl = sub("\\..*$", "", rownames(counts))) |>
  left_join(symbols, by = c("ensembl" = "ENSEMBL")) |>
  filter(SYMBOL %in% mapped$mapped_symbol) |>
  group_by(SYMBOL) |> slice(1) |> ungroup()
present <- mapped |> filter(mapped_symbol %in% gene_rows$SYMBOL)
expr <- log2(sweep(counts[gene_rows$gene_id, patients$sample_barcode, drop = FALSE], 2,
                   colSums(counts[, patients$sample_barcode, drop = FALSE]), "/") * 1e6 + 1)
rownames(expr) <- gene_rows$SYMBOL[match(rownames(expr), gene_rows$gene_id)]

signed_score <- function(train_idx, idx, genes) {
  used <- genes |> filter(mapped_symbol %in% rownames(expr))
  train_expr <- expr[used$mapped_symbol, train_idx, drop = FALSE]
  center <- rowMeans(train_expr)
  scale <- apply(train_expr, 1, sd)
  scale[scale == 0 | !is.finite(scale)] <- NA_real_
  z <- sweep(sweep(expr[used$mapped_symbol, idx, drop = FALSE], 1, center, "-"), 1, scale, "/")
  colMeans(z * used$sign, na.rm = TRUE)
}

n_repeats <- 10L
fold_rows <- list()
pred_rows <- list()
for (repeat_id in seq_len(n_repeats)) {
  set.seed(RESAMPLING$seed + 22L + repeat_id)
  fold_id <- event_stratified_folds(patients$os_event, 5L)
  for (fold in seq_len(5L)) {
    train_idx <- which(fold_id != fold)
    test_idx <- which(fold_id == fold)
    train <- patients[train_idx, ]
    test <- patients[test_idx, ]
    train$score <- signed_score(train_idx, train_idx, present)
    test$score <- signed_score(train_idx, test_idx, present)
    clinical_fit <- coxph(Surv(os_time, os_event) ~ age + sex + stage + grade, data = train)
    score_fit <- coxph(Surv(os_time, os_event) ~ age + sex + stage + grade + score, data = train)
    c_hat <- risk_from_fit(clinical_fit, test, 1095)
    s_hat <- risk_from_fit(score_fit, test, 1095)
    pred_rows[[length(pred_rows) + 1L]] <- tibble(
      patient_barcode = test$patient_barcode, os_time = test$os_time, os_event = test$os_event,
      repeat_id = repeat_id, fold = fold, clinical_lp = c_hat$lp, score_lp = s_hat$lp,
      clinical_risk3 = c_hat$risk, score_risk3 = s_hat$risk)
    fold_rows[[length(fold_rows) + 1L]] <- tibble(repeat_id, fold, n_test = nrow(test),
                                                  test_events = sum(test$os_event))
  }
  message("ClearCode34 repeat ", repeat_id)
}
pred <- bind_rows(pred_rows)
repeat_metrics <- pred |>
  group_by(repeat_id) |>
  group_modify(function(x, key) {
    clin <- score_block(x, "clinical_lp", "clinical_risk3", 1095)
    sig <- score_block(x, "score_lp", "score_risk3", 1095)
    tibble(clinical_c = unname(clin["c_index"]), signature_c = unname(sig["c_index"]),
           delta_c = unname(sig["c_index"] - clin["c_index"]),
           clinical_brier3 = unname(clin["brier"]), signature_brier3 = unname(sig["brier"]),
           delta_brier3 = unname(sig["brier"] - clin["brier"]),
           clinical_slope = unname(clin["slope"]), signature_slope = unname(sig["slope"]))
  }) |>
  ungroup()
set.seed(RESAMPLING$seed + 32L)
ci <- bootstrap_delta(pred, "score_lp", "score_risk3", "clinical_lp", "clinical_risk3", 1095)
tcga_row <- tibble(
  cohort = "TCGA-KIRC", arm = "clinical_plus_clearcode34", evaluation = "nested",
  n = n_distinct(pred$patient_barcode), events = NA_real_,
  c_index = mean(repeat_metrics$signature_c),
  delta_c = mean(repeat_metrics$delta_c),
  delta_c_ci_low = ci[1], delta_c_ci_high = ci[2],
  brier = mean(repeat_metrics$signature_brier3),
  delta_brier3 = mean(repeat_metrics$delta_brier3),
  calibration_slope = mean(repeat_metrics$signature_slope),
  genes_used = sum(!is.na(present$mapped_symbol))
)

# E-MTAB-1980 in-sample comparison on the previously inspected cohort.
em_expr_path <- file.path(DIRS$raw, "arrayexpress", "E-MTAB-1980", "ccRCC_exp_log_quantile_normalized.txt")
em_clin_path <- file.path(DIRS$raw, "arrayexpress", "E-MTAB-1980", "supplementary_table_1.xlsx")
em_raw <- readxl::read_excel(em_clin_path, sheet = 1, skip = 1)
em <- em_raw |>
  transmute(patient_barcode = as.character(.data[["sample ID"]]),
            age = suppressWarnings(as.numeric(.data[["Age"]])),
            subtype = as.character(.data[["gene expression profile"]]),
            os_event = as.integer(str_to_lower(as.character(.data[["outcome"]])) == "dead"),
            os_time = suppressWarnings(as.numeric(.data[["observation period (month)"]])),
            t_high = str_detect(as.character(.data[["Stage at diagnosis"]]), "^pT(?:3|4)"),
            metastatic = str_detect(as.character(.data[["Stage at diagnosis"]]), "M1")) |>
  filter(subtype %in% c("ccA", "ccB"), is.finite(os_time), os_time > 0, !is.na(os_event))
em_mat <- read_tsv(em_expr_path, show_col_types = FALSE, progress = FALSE, name_repair = "minimal")
refseq <- sub("\\.\\d+$", "", as.character(em_mat$SystematicName))
ref_map <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(refseq[refseq != ""]),
                                 keytype = "REFSEQ", columns = "SYMBOL") |>
  filter(!is.na(SYMBOL)) |> distinct(REFSEQ, .keep_all = TRUE)
sample_columns <- intersect(em$patient_barcode, names(em_mat))
em_gene <- em_mat |>
  mutate(refseq = sub("\\.\\d+$", "", as.character(SystematicName))) |>
  inner_join(ref_map, by = c("refseq" = "REFSEQ")) |>
  select(SYMBOL, all_of(sample_columns)) |>
  mutate(across(all_of(sample_columns), as.numeric)) |>
  group_by(SYMBOL) |>
  summarise(across(all_of(sample_columns), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
gene_vector <- function(symbol) {
  row <- em_gene |> filter(SYMBOL == symbol)
  if (nrow(row) != 1L) return(rep(NA_real_, nrow(em)))
  vals <- as.numeric(row[1, sample_columns])
  names(vals) <- sample_columns
  as.numeric(scale(vals[em$patient_barcode]))
}
high <- read_csv(file.path(DIRS$tables, "high_confidence_candidate_genes.csv"), show_col_types = FALSE) |>
  arrange(desc(evidence_score))
mapped_expr <- lapply(high$symbol, gene_vector)
mapped_ok <- vapply(mapped_expr, function(x) any(is.finite(x)), logical(1))
if (!any(mapped_ok)) stop("No high-confidence gene is present in E-MTAB-1980.")
top_symbol <- high$symbol[which(mapped_ok)[1]]
panel <- high[mapped_ok, ]
panel$expr <- mapped_expr[mapped_ok]
panel_score <- colMeans(do.call(rbind, lapply(seq_len(nrow(panel)), function(i) {
  sign(panel$main_log_hr[i]) * panel$expr[[i]]
})), na.rm = TRUE)
em$subtype_ccb <- as.integer(em$subtype == "ccB")
em$top_expr <- gene_vector(top_symbol)
em$panel_score <- panel_score
em_complete <- em |> filter(is.finite(age), !is.na(t_high), !is.na(metastatic),
                            is.finite(top_expr), is.finite(panel_score))
if (nrow(em_complete) < 50L) stop("E-MTAB-1980 complete cases collapsed to ", nrow(em_complete), ".")
fit_arm <- function(formula) {
  fit <- coxph(formula, data = em_complete)
  hat <- risk_from_fit(fit, em_complete, 36)
  em_complete |>
    mutate(lp = hat$lp, risk = hat$risk) |>
    select(patient_barcode, os_time, os_event, lp, risk)
}
arms <- list(
  clinical_limited = fit_arm(Surv(os_time, os_event) ~ age + t_high + metastatic),
  subtype = fit_arm(Surv(os_time, os_event) ~ subtype_ccb),
  clinical_plus_subtype = fit_arm(Surv(os_time, os_event) ~ age + t_high + metastatic + subtype_ccb),
  clinical_plus_top_gene = fit_arm(Surv(os_time, os_event) ~ age + t_high + metastatic + top_expr),
  clinical_plus_signed_panel = fit_arm(Surv(os_time, os_event) ~ age + t_high + metastatic + panel_score),
  unadjusted_top_gene = fit_arm(Surv(os_time, os_event) ~ top_expr),
  unadjusted_signed_panel = fit_arm(Surv(os_time, os_event) ~ panel_score)
)
clinical_pred <- arms$clinical_limited
em_rows <- bind_rows(lapply(names(arms), function(arm) {
  pred_arm <- arms[[arm]]
  point <- score_block(pred_arm, "lp", "risk", 36)
  base <- if (arm == "clinical_limited") point else score_block(clinical_pred, "lp", "risk", 36)
  set.seed(RESAMPLING$seed + 32L)
  if (arm == "clinical_limited") {
    ci <- c(NA_real_, NA_real_)
    delta <- 0
  } else {
    joined <- pred_arm |>
      transmute(patient_barcode, os_time, os_event, score_lp = lp, score_risk = risk) |>
      inner_join(clinical_pred |> transmute(patient_barcode, clinical_lp = lp, clinical_risk = risk),
                 by = "patient_barcode")
    ci <- bootstrap_delta(joined, "score_lp", "score_risk", "clinical_lp", "clinical_risk", 36)
    delta <- unname(point["c_index"] - base["c_index"])
  }
  tibble(cohort = "E-MTAB-1980", arm = arm, evaluation = "in_sample_previously_inspected",
         n = nrow(pred_arm), events = sum(pred_arm$os_event),
         c_index = unname(point["c_index"]), delta_c = delta,
         delta_c_ci_low = ci[1], delta_c_ci_high = ci[2],
         brier = unname(point["brier"]), delta_brier3 = unname(point["brier"] - base["brier"]),
         calibration_slope = unname(point["slope"]),
         genes_used = if (arm == "clinical_plus_signed_panel") nrow(panel) else NA_real_)
}))
summary <- bind_rows(tcga_row, em_rows)
write_csv_atomic(summary, file.path(DIRS$tables, "published_signature_summary.csv"))
write_csv_atomic(repeat_metrics, file.path(DIRS$tables, "clearcode34_repeat_metrics.csv"))
write_csv_atomic(pred, file.path(DIRS$tables, "clearcode34_predictions.csv"))
message("Published-signature benchmark written. Top full-data gene: ", top_symbol,
        ". ClearCode34 genes mapped: ", sum(!is.na(present$mapped_symbol)))
