source("analysis/00_config.R")
source("analysis/functions/io.R")
source("analysis/functions/tcga_metadata.R")
source("analysis/functions/patient_samples.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(dplyr)
  library(readr)
  library(survival)
})

set.seed(RESAMPLING$seed + 22L)
n_repeats <- as.integer(Sys.getenv("HYDRA_NESTED_REPEATS", "10"))
n_null <- as.integer(Sys.getenv("HYDRA_NULL_SIMS", "200"))
stopifnot(n_repeats >= 1L, n_null >= 0L)
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
stopifnot(!anyDuplicated(patients$patient_barcode), nrow(patients) > 400L)

symbols <- AnnotationDbi::select(org.Hs.eg.db,
  keys = unique(sub("\\..*$", "", rownames(counts))),
  keytype = "ENSEMBL", columns = "SYMBOL") |>
  filter(!is.na(SYMBOL)) |> arrange(ENSEMBL, SYMBOL) |>
  distinct(ENSEMBL, .keep_all = TRUE)
geo_a <- read_csv(file.path(DIRS$tables, "gse40435_limma_tumor_vs_normal.csv"),
                  show_col_types = FALSE) |>
  select(symbol, gse40435_log2fc = log2FoldChange, gse40435_pvalue = pvalue)
geo_b <- read_csv(file.path(DIRS$tables, "gse53757_limma_tumor_vs_normal.csv"),
                  show_col_types = FALSE) |>
  select(symbol, gse53757_log2fc = log2FoldChange, gse53757_pvalue = pvalue)
geo <- tibble(tcga_gene_id = rownames(counts),
              ensembl = sub("\\..*$", "", rownames(counts))) |>
  left_join(symbols, by = c("ensembl" = "ENSEMBL")) |>
  rename(symbol = SYMBOL) |>
  inner_join(geo_a, by = "symbol") |>
  inner_join(geo_b, by = "symbol")
tumor_counts <- counts[, patients$sample_barcode, drop = FALSE]
all_normal <- meta |> filter(condition == "NT")

folds <- function(event, k) {
  out <- integer(length(event))
  for (e in c(0L, 1L)) {
    idx <- which(event == e)
    out[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  out
}

make_fold_data <- function(train_idx, test_idx) {
  train_patient <- patients$patient_barcode[train_idx]
  normals <- all_normal |> filter(!patient_barcode %in% patients$patient_barcode[test_idx])
  train_meta <- bind_rows(patients[train_idx, names(meta)], normals) |>
    mutate(condition = factor(condition, levels = c("NT", "TP")))
  train_raw <- counts[, train_meta$sample_barcode, drop = FALSE]
  keep <- rowSums(train_raw >= THRESHOLDS$min_count) >= THRESHOLDS$min_samples
  train_raw <- train_raw[keep, , drop = FALSE]
  train_meta <- as.data.frame(train_meta)
  rownames(train_meta) <- train_meta$sample_barcode
  dds <- DESeqDataSetFromMatrix(train_raw, train_meta, design = ~ condition)
  # Match the primary DESeq2 fit. If its optional Cook outlier replacement
  # fails in a fold, retry without replacement and record that deviation.
  replacement_fallback <- FALSE
  dds <- tryCatch(DESeq(dds, quiet = TRUE), error = function(e) {
    if (!grepl("is.finite(partial)", conditionMessage(e), fixed = TRUE)) stop(e)
    message("DESeq2 outlier replacement failed: ", conditionMessage(e))
    replacement_fallback <<- TRUE
    DESeq(dds, quiet = TRUE, minReplicatesForReplace = Inf)
  })
  de <- as.data.frame(results(dds, contrast = c("condition", "TP", "NT")))
  de$gene_id <- rownames(de)
  de <- de |>
    inner_join(geo, by = c("gene_id" = "tcga_gene_id")) |>
    group_by(symbol) |>
    arrange(padj, .by_group = TRUE) |>
    slice_head(n = 1L) |>
    ungroup() |>
    filter(!is.na(padj), padj < THRESHOLDS$deg_fdr,
           abs(log2FoldChange) >= THRESHOLDS$deg_abs_log2fc,
           sign(log2FoldChange) == sign(gse40435_log2fc),
           sign(log2FoldChange) == sign(gse53757_log2fc),
           (gse40435_pvalue < 0.05 | gse53757_pvalue < 0.05))
  if (nrow(de) == 0L) return(list(de = de, expr = matrix(numeric(), 0, nrow(patients)),
                                 replacement_fallback = replacement_fallback))

  # Estimate library size using training samples only, then apply the frozen
  # geometric-mean reference to held-out tumors.
  train_sf <- sizeFactors(dds)
  geo_mean <- exp(rowMeans(log(train_raw)))
  test_raw <- tumor_counts[rownames(train_raw), test_idx, drop = FALSE]
  test_sf <- DESeq2::estimateSizeFactorsForMatrix(test_raw, geoMeans = geo_mean)
  all_sf <- c(setNames(train_sf, colnames(train_raw)),
              setNames(test_sf, colnames(test_raw)))
  chosen <- intersect(de$gene_id, rownames(train_raw))
  expr <- log2(sweep(counts[chosen, patients$sample_barcode, drop = FALSE],
                     2, all_sf[patients$sample_barcode], "/") + 1)
  stopifnot(all(is.finite(expr)))
  list(de = de, expr = expr, replacement_fallback = replacement_fallback)
}

gene_fit <- function(dat, expr, covariates) {
  dat$expr <- expr
  tryCatch(suppressWarnings(coxph(
    as.formula(paste("Surv(os_time, os_event) ~", paste(c("expr", covariates), collapse = "+"))),
    data = dat
  )), error = function(e) NULL)
}

select_gene <- function(fd, train_idx, times, events) {
  de <- fd$de
  if (nrow(de) == 0L) return(NA_character_)
  dat <- patients[train_idx, ]
  dat$os_time <- times[train_idx]
  dat$os_event <- events[train_idx]
  main_model <- t(vapply(de$gene_id, function(id) {
    x <- as.numeric(fd$expr[id, train_idx])
    if (!is.finite(sd(x)) || sd(x) == 0) return(c(NA_real_, NA_real_))
    x <- as.numeric(scale(x))
    fit <- gene_fit(dat, x, c("age", "sex", "stage", "grade"))
    if (is.null(fit)) return(c(NA_real_, NA_real_))
    c(unname(coef(fit)["expr"]), summary(fit)$coefficients["expr", "Pr(>|z|)"])
  }, numeric(2)))
  colnames(main_model) <- c("main_beta", "main_p")
  de <- bind_cols(de, as.data.frame(main_model))
  de$main_fdr <- p.adjust(de$main_p, method = "BH")
  de <- de |>
    filter(is.finite(main_fdr), main_fdr < THRESHOLDS$high_confidence_survival_fdr,
           abs(main_beta) >= THRESHOLDS$high_confidence_abs_log_hr,
           abs(gse40435_log2fc) >= THRESHOLDS$min_geo_abs_log2fc,
           abs(gse53757_log2fc) >= THRESHOLDS$min_geo_abs_log2fc)
  if (nrow(de) == 0L) return(NA_character_)
  sensitivity <- t(vapply(seq_len(nrow(de)), function(i) {
    x <- as.numeric(scale(fd$expr[de$gene_id[i], train_idx]))
    fits <- lapply(list(c("age", "sex", "stage"), c("age", "sex", "grade")),
                   function(vars) gene_fit(dat, x, vars))
    if (any(vapply(fits, is.null, logical(1)))) return(rep(NA_real_, 4L))
    c(vapply(fits, function(f) unname(coef(f)["expr"]), numeric(1)),
      vapply(fits, function(f) summary(f)$coefficients["expr", "Pr(>|z|)"], numeric(1)))
  }, numeric(4)))
  colnames(sensitivity) <- c("stage_beta", "grade_beta", "stage_p", "grade_p")
  de <- bind_cols(de, as.data.frame(sensitivity)) |>
    filter(sign(stage_beta) == sign(main_beta), sign(grade_beta) == sign(main_beta),
           stage_p < 0.05, grade_p < 0.05) |>
    mutate(evidence_score = -log10(pmax(main_fdr, .Machine$double.xmin)) +
             abs(main_beta) + pmin(abs(log2FoldChange), 5) / 5 +
             pmin(abs(gse40435_log2fc), 3) / 3 +
             pmin(abs(gse53757_log2fc), 3) / 3) |>
    arrange(desc(evidence_score), main_fdr, gene_id)
  if (nrow(de) == 0L) NA_character_ else de$gene_id[1]
}

predict_fold <- function(fd, train_idx, test_idx, selected, times, events) {
  train <- patients[train_idx, ]
  test <- patients[test_idx, ]
  train$os_time <- times[train_idx]
  train$os_event <- events[train_idx]
  clinical_fit <- coxph(Surv(os_time, os_event) ~ age + sex + stage + grade, data = train)
  gene_model <- clinical_fit
  if (!is.na(selected)) {
    train_mean <- mean(fd$expr[selected, train_idx])
    train_sd <- sd(fd$expr[selected, train_idx])
    train$expr <- (fd$expr[selected, train_idx] - train_mean) / train_sd
    test$expr <- (fd$expr[selected, test_idx] - train_mean) / train_sd
    gene_model <- gene_fit(train, train$expr, c("age", "sex", "stage", "grade"))
    if (is.null(gene_model)) gene_model <- clinical_fit
  }
  risk <- function(fit, newdata) {
    lp <- as.numeric(predict(fit, newdata = newdata, type = "lp", reference = "zero"))
    bh <- basehaz(fit, centered = FALSE)
    h3 <- if (any(bh$time <= 1095)) tail(bh$hazard[bh$time <= 1095], 1) else 0
    list(lp = lp, risk3 = 1 - exp(-h3 * exp(lp)))
  }
  c_pred <- risk(clinical_fit, test)
  g_pred <- risk(gene_model, test)
  tibble(patient_barcode = test$patient_barcode, os_time = times[test_idx],
         os_event = events[test_idx], clinical_lp = c_pred$lp,
         gene_lp = g_pred$lp, clinical_risk3 = c_pred$risk3,
         gene_risk3 = g_pred$risk3, selected_gene = selected)
}

score_repeat <- function(pred) {
  ci <- function(lp) unname(concordance(Surv(pred$os_time, pred$os_event) ~ lp,
                                        reverse = TRUE)$concordance)
  # KM censoring weights; patients censored before three years contribute zero.
  cens <- survfit(Surv(os_time, 1L - os_event) ~ 1, data = pred)
  censor_surv <- function(t) {
    i <- findInterval(t, cens$time)
    ifelse(i == 0L, 1, cens$surv[pmax(i, 1L)])
  }
  horizon <- 1095
  target <- as.numeric(pred$os_time <= horizon & pred$os_event == 1L)
  weight <- ifelse(target == 1, 1 / pmax(censor_surv(pmax(pred$os_time - 1, 0)), 0.01),
                   ifelse(pred$os_time > horizon, 1 / pmax(censor_surv(horizon), 0.01), 0))
  brier <- function(risk) mean(weight * (target - risk)^2)
  slope <- function(lp) tryCatch(unname(coef(coxph(Surv(pred$os_time, pred$os_event) ~ lp))[1]),
                                  error = function(e) NA_real_)
  tibble(clinical_c = ci(pred$clinical_lp), gene_c = ci(pred$gene_lp),
         delta_c = gene_c - clinical_c,
         clinical_brier3 = brier(pred$clinical_risk3),
         gene_brier3 = brier(pred$gene_risk3),
         delta_brier3 = gene_brier3 - clinical_brier3,
         clinical_calibration_slope = slope(pred$clinical_lp),
         gene_calibration_slope = slope(pred$gene_lp))
}

predictions <- vector("list", n_repeats)
fold_summary <- list()
first_fold_cache <- vector("list", 5L)
for (repeat_id in seq_len(n_repeats)) {
  set.seed(RESAMPLING$seed + 22L + repeat_id)
  fold_id <- folds(patients$os_event, 5L)
  fold_predictions <- vector("list", 5L)
  for (fold in seq_len(5L)) {
    train_idx <- which(fold_id != fold)
    test_idx <- which(fold_id == fold)
    message("Nested selection: repeat ", repeat_id, "/", n_repeats, " fold ", fold, "/5")
    fd <- make_fold_data(train_idx, test_idx)
    selected <- select_gene(fd, train_idx, patients$os_time, patients$os_event)
    fold_predictions[[fold]] <- predict_fold(fd, train_idx, test_idx, selected,
                                              patients$os_time, patients$os_event) |>
      mutate(repeat_id = repeat_id, fold = fold)
    fold_summary[[length(fold_summary) + 1L]] <- tibble(repeat_id, fold,
      n_train = length(train_idx), n_test = length(test_idx),
      train_de_genes = nrow(fd$de), selected_gene = selected,
      outlier_replacement_fallback = fd$replacement_fallback,
      no_gene_selected = is.na(selected))
    if (repeat_id == 1L) first_fold_cache[[fold]] <- list(fd = fd, train = train_idx, test = test_idx)
    gc()
  }
  predictions[[repeat_id]] <- bind_rows(fold_predictions)
}
pred <- bind_rows(predictions)
scores <- pred |> group_by(repeat_id) |> group_modify(~ score_repeat(.x)) |> ungroup()

set.seed(RESAMPLING$seed + 23L)
patient_ids <- unique(pred$patient_barcode)
boot <- bind_rows(lapply(seq_len(1000L), function(i) {
  sampled <- sample(patient_ids, length(patient_ids), replace = TRUE)
  x <- pred |> filter(patient_barcode %in% sampled)
  # Replication weights preserve duplicated patients in the bootstrap draw.
  weights <- as.integer(table(factor(sampled, levels = patient_ids)))
  x <- x[rep(seq_len(nrow(x)), weights[match(x$patient_barcode, patient_ids)]), ]
  tibble(draw = i, delta_c = mean(x |> group_by(repeat_id) |>
                                   group_modify(~ score_repeat(.x)) |> pull(delta_c)))
}))
summary <- tibble(n_patients = nrow(patients), repeats = n_repeats, folds = 5L,
  no_gene_folds = sum(vapply(fold_summary, function(x) x$no_gene_selected, logical(1))),
  outlier_replacement_fallbacks = sum(vapply(fold_summary,
    function(x) x$outlier_replacement_fallback, logical(1))),
  mean_delta_c = mean(scores$delta_c),
  patient_bootstrap_ci_low = quantile(boot$delta_c, 0.025),
  patient_bootstrap_ci_high = quantile(boot$delta_c, 0.975),
  mean_delta_brier3 = mean(scores$delta_brier3),
  cv_acceptance = mean_delta_c >= 0.01 & quantile(boot$delta_c, 0.025) > 0 &
    mean(scores$delta_brier3) <= 0)
write_csv_atomic(pred, file.path(DIRS$tables, "nested_cv_predictions.csv"))
write_csv_atomic(bind_rows(fold_summary), file.path(DIRS$tables, "nested_cv_folds.csv"))
write_csv_atomic(scores, file.path(DIRS$tables, "nested_cv_repeat_metrics.csv"))
write_csv_atomic(boot, file.path(DIRS$tables, "nested_cv_patient_bootstrap.csv"))
write_csv_atomic(summary, file.path(DIRS$tables, "nested_cv_summary.csv"))

# The null generates survival under a fitted clinical-only model and resamples
# censoring times. The same fold-level DE and survival selection is then rerun.
baseline <- coxph(Surv(os_time, os_event) ~ age + sex + stage + grade, data = patients)
null_lp <- as.numeric(predict(baseline, newdata = patients, type = "lp", reference = "zero"))
base_hazard <- tail(basehaz(baseline, centered = FALSE)$hazard, 1)
censor_times <- patients$os_time[patients$os_event == 0L]
null_rows <- vector("list", n_null)
set.seed(RESAMPLING$seed + 24L)
for (sim in seq_len(n_null)) {
  simulated_event_time <- -log(runif(nrow(patients))) / (base_hazard * exp(null_lp)) *
    max(patients$os_time)
  simulated_censor <- sample(censor_times, nrow(patients), replace = TRUE)
  sim_time <- pmin(simulated_event_time, simulated_censor)
  sim_event <- as.integer(simulated_event_time <= simulated_censor)
  # One held-out fold per simulation cycles across the five outer folds. The
  # complete training selection is repeated, while 200 simulations remain
  # computationally feasible on the cached public inputs.
  fold <- ((sim - 1L) %% 5L) + 1L
  f <- first_fold_cache[[fold]]
  selected <- select_gene(f$fd, f$train, sim_time, sim_event)
  sim_pred <- predict_fold(f$fd, f$train, f$test, selected, sim_time, sim_event)
  null_rows[[sim]] <- score_repeat(sim_pred) |>
    mutate(simulation = sim, fold = fold, selected_gene = selected)
  if (sim %% 10L == 0L) message("Clinical-only null simulation ", sim, "/", n_null)
}
null_table <- bind_rows(null_rows)
if (n_null > 0L) write_csv_atomic(null_table, file.path(DIRS$tables, "nested_cv_clinical_null.csv"))
message("Nested CV and clinical-only null complete.")
