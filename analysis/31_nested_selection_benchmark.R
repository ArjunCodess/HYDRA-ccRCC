# Nested comparison of the HYDRA rule with simpler selectors.
# Same patients, seeds, and event-stratified folds as analysis/22_nested_cv.R.
# Does not replace the prespecified concordance criterion.
#
# Full run:
#   Rscript analysis/31_nested_selection_benchmark.R
# One repeat or fold, for a smoke test:
#   HYDRA_REPEAT_IDS=1 HYDRA_FOLD_IDS=1 Rscript analysis/31_nested_selection_benchmark.R

source("analysis/00_config.R")
source("analysis/functions/io.R")
source("analysis/functions/tcga_metadata.R")
source("analysis/functions/patient_samples.R")
source("analysis/functions/nested_benchmark.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(survival)
  library(glmnet)
})

say <- function(...) {
  message(...)
  flush(stderr())
}

n_repeats <- as.integer(Sys.getenv("HYDRA_BENCHMARK_REPEATS", "10"))
screen_n <- as.integer(Sys.getenv("HYDRA_SURVIVAL_SCREEN", "2000"))
stopifnot(n_repeats >= 1L, screen_n >= 50L)
worker_spec <- Sys.getenv("HYDRA_REPEAT_IDS", "")
worker_mode <- nzchar(worker_spec)
repeat_ids <- if (worker_mode) {
  as.integer(strsplit(worker_spec, ",", fixed = TRUE)[[1]])
} else {
  seq_len(n_repeats)
}
fold_spec <- Sys.getenv("HYDRA_FOLD_IDS", "")
fold_ids <- if (nzchar(fold_spec)) {
  as.integer(strsplit(fold_spec, ",", fixed = TRUE)[[1]])
} else {
  seq_len(5L)
}
if (anyNA(repeat_ids) || anyDuplicated(repeat_ids) || any(!repeat_ids %in% seq_len(n_repeats))) {
  stop("Invalid HYDRA_REPEAT_IDS.")
}
if (anyNA(fold_ids) || anyDuplicated(fold_ids) || any(!fold_ids %in% seq_len(5L))) {
  stop("Invalid HYDRA_FOLD_IDS.")
}
if (!identical(fold_ids, seq_len(max(fold_ids)))) {
  stop("HYDRA_FOLD_IDS must be a prefix of 1,2,3,4,5 so the DESeq2 random stream matches the primary nested CV.")
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
primary_folds <- read_csv(file.path(DIRS$tables, "nested_cv_folds.csv"), show_col_types = FALSE)
symbol_of <- function(gene_id) {
  if (length(gene_id) != 1L || is.na(gene_id)) return(NA_character_)
  hit <- symbols$SYMBOL[match(sub("\\..*$", "", gene_id), symbols$ENSEMBL)]
  if (length(hit) == 0L || is.na(hit)) NA_character_ else unname(hit)
}

make_fold_data <- function(train_idx, test_idx) {
  normals <- all_normal |> filter(!patient_barcode %in% patients$patient_barcode[test_idx])
  train_meta <- bind_rows(patients[train_idx, names(meta)], normals) |>
    mutate(condition = factor(condition, levels = c("NT", "TP")))
  train_raw <- counts[, train_meta$sample_barcode, drop = FALSE]
  keep <- rowSums(train_raw >= THRESHOLDS$min_count) >= THRESHOLDS$min_samples
  train_raw <- train_raw[keep, , drop = FALSE]
  train_meta <- as.data.frame(train_meta)
  rownames(train_meta) <- train_meta$sample_barcode
  dds <- DESeqDataSetFromMatrix(train_raw, train_meta, design = ~ condition)
  replacement_fallback <- FALSE
  fitted <- NULL
  for (attempt in seq_len(3L)) {
    fit <- tryCatch(
      if (attempt == 1L) DESeq(dds, quiet = TRUE)
      else DESeq(dds, quiet = TRUE, minReplicatesForReplace = Inf),
      error = identity
    )
    if (!inherits(fit, "error")) {
      fitted <- fit
      break
    }
    if (!identical(conditionMessage(fit),
                   "default method not implemented for type 'expression'") ||
        !identical(deparse(conditionCall(fit)), "is.finite(partial)")) stop(fit)
    replacement_fallback <- TRUE
    say("DESeq2 trimmed-mean failure in fold, attempt ", attempt,
        "; retrying without outlier replacement.")
    gc()
  }
  if (is.null(fitted)) stop("DESeq2 failed after three fold-fit attempts.")
  dds <- fitted
  raw_de <- as.data.frame(results(dds, contrast = c("condition", "TP", "NT")))
  raw_de$gene_id <- rownames(raw_de)
  mapped <- tibble(gene_id = raw_de$gene_id,
                   ensembl = sub("\\..*$", "", raw_de$gene_id)) |>
    left_join(symbols, by = c("ensembl" = "ENSEMBL")) |>
    rename(symbol = SYMBOL)
  de_tcga <- raw_de |>
    left_join(mapped |> select(gene_id, symbol), by = "gene_id") |>
    mutate(group_key = if_else(is.na(symbol) | symbol == "", gene_id, symbol)) |>
    group_by(group_key) |>
    arrange(padj, gene_id, .by_group = TRUE) |>
    slice_head(n = 1L) |>
    ungroup()
  de <- raw_de |>
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

  train_sf <- sizeFactors(dds)
  geo_mean <- exp(rowMeans(log(train_raw)))
  test_raw <- tumor_counts[rownames(train_raw), test_idx, drop = FALSE]
  test_sf <- DESeq2::estimateSizeFactorsForMatrix(test_raw, geoMeans = geo_mean)
  all_sf <- c(setNames(train_sf, colnames(train_raw)),
              setNames(test_sf, colnames(test_raw)))
  expr_all <- log2(sweep(counts[rownames(train_raw), patients$sample_barcode, drop = FALSE],
                         2, all_sf[patients$sample_barcode], "/") + 1)
  stopifnot(all(is.finite(expr_all)),
            identical(colnames(expr_all), patients$sample_barcode))
  list(de = de, de_tcga = de_tcga, expr_all = expr_all,
       replacement_fallback = replacement_fallback)
}

same_gene <- function(a, b) {
  if (is.na(a) && is.na(b)) TRUE else identical(a, b)
}

clinical_matrix <- function(data) {
  mm <- stats::model.matrix(~ age + sex + stage + grade, data = data)
  mm[, colnames(mm) != "(Intercept)", drop = FALSE]
}

predict_gene <- function(clinical_fit, expr_all, gene_id, train_idx, test_idx, train, test) {
  if (is.na(gene_id)) return(list(risk = cox_risk(clinical_fit, test), fallback = TRUE))
  x_train <- as.numeric(expr_all[gene_id, train_idx])
  if (!is.finite(stats::sd(x_train)) || stats::sd(x_train) == 0) {
    return(list(risk = cox_risk(clinical_fit, test), fallback = TRUE))
  }
  scaled <- scale(x_train)
  train$expr <- as.numeric(scaled)
  test$expr <- as.numeric((as.numeric(expr_all[gene_id, test_idx]) -
                             attr(scaled, "scaled:center")) / attr(scaled, "scaled:scale"))
  fit <- gene_cox(train, train$expr, c("age", "sex", "stage", "grade"))
  if (is.null(fit)) return(list(risk = cox_risk(clinical_fit, test), fallback = TRUE))
  list(risk = cox_risk(fit, test), fallback = FALSE)
}

checkpoint_dir <- file.path(DIRS$processed, "nested_benchmark_checkpoints")
dir.create(checkpoint_dir, showWarnings = FALSE, recursive = TRUE)
checkpoint_inputs <- c(
  "analysis/31_nested_selection_benchmark.R",
  "analysis/functions/nested_benchmark.R",
  "analysis/00_config.R",
  "analysis/functions/patient_samples.R",
  FILES$tcga_se,
  FILES$tcga_clinical,
  file.path(DIRS$tables, "gse40435_limma_tumor_vs_normal.csv"),
  file.path(DIRS$tables, "gse53757_limma_tumor_vs_normal.csv")
)
spec <- list(files = unname(tools::md5sum(checkpoint_inputs)),
             screen_n = screen_n, seed = RESAMPLING$seed,
             ridge = "glmnet-5.0-efron-alpha0-lambda.min")
checkpoint_path <- function(repeat_id, fold) {
  file.path(checkpoint_dir, sprintf("repeat_%02d_fold_%d.rds", repeat_id, fold))
}

run_fold <- function(repeat_id, fold, fold_id) {
  path <- checkpoint_path(repeat_id, fold)
  if (file.exists(path)) {
    saved <- readRDS(path)
    if (identical(saved$spec, spec)) {
      assign(".Random.seed", saved$rng_after_deseq, envir = .GlobalEnv)
      say("Reused benchmark checkpoint: repeat ", repeat_id, " fold ", fold)
      return(saved)
    }
  }
  started <- proc.time()[["elapsed"]]
  train_idx <- which(fold_id != fold)
  test_idx <- which(fold_id == fold)
  say("Benchmark fold: repeat ", repeat_id, "/", n_repeats, " fold ", fold, "/5")
  fd <- make_fold_data(train_idx, test_idx)
  rng_after_deseq <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  train <- patients[train_idx, ]
  test <- patients[test_idx, ]
  train$os_time <- patients$os_time[train_idx]
  train$os_event <- patients$os_event[train_idx]
  expr_train <- fd$expr_all[, train_idx, drop = FALSE]
  clinical_fit <- coxph(Surv(os_time, os_event) ~ age + sex + stage + grade, data = train)
  primary_gene <- primary_folds$selected_gene[
    primary_folds$repeat_id == repeat_id & primary_folds$fold == fold
  ]
  if (length(primary_gene) != 1L) stop("Primary nested-CV fold gene is missing.")
  hydra <- select_hydra_gene(fd$de, expr_train, train, THRESHOLDS)
  if (!same_gene(hydra$gene_id, primary_gene)) {
    stop("HYDRA gene ", hydra$gene_id, " does not match primary nested CV gene ",
         primary_gene, " for repeat ", repeat_id, " fold ", fold, ".")
  }
  say("HYDRA gene matches primary nested CV: ",
      ifelse(is.na(hydra$gene_id), "none", hydra$gene_id))

  resid <- residuals(clinical_fit, type = "martingale")
  abs_cor <- abs(residual_correlation(expr_train, resid))
  survival_arm <- select_survival_only_gene(expr_train, train, abs_cor, THRESHOLDS, screen_n)
  de_arm <- select_de_only_gene(fd$de_tcga, THRESHOLDS$deg_fdr, THRESHOLDS$deg_abs_log2fc)
  if (!is.na(de_arm$gene_id) && !de_arm$gene_id %in% rownames(fd$expr_all)) {
    stop("DE-only gene is missing from the normalized matrix.")
  }
  target <- if (!is.na(hydra$gene_id)) {
    mean(as.numeric(expr_train[hydra$gene_id, ]))
  } else if (nrow(fd$de) > 0L) {
    stats::median(rowMeans(expr_train[fd$de$gene_id, , drop = FALSE]))
  } else {
    NA_real_
  }
  control_id <- select_matched_control(
    rowMeans(expr_train), abs_cor, target,
    excluded = unique(c(hydra$gene_id, survival_arm$gene_id, de_arm$gene_id))
  )

  clin_mm <- clinical_matrix(patients)
  train_sd <- apply(clin_mm[train_idx, , drop = FALSE], 2, stats::sd)
  if (any(!is.finite(train_sd) | train_sd == 0)) {
    stop("Training fold is missing a clinical factor level.")
  }
  ridge <- NULL
  ridge_fallback <- FALSE
  repro_ids <- intersect(fd$de$gene_id, rownames(expr_train))
  if (length(repro_ids) >= 1L) {
    gene_train <- t(expr_train[repro_ids, , drop = FALSE])
    gene_sd <- apply(gene_train, 2, stats::sd)
    keep_gene <- is.finite(gene_sd) & gene_sd > 0
    repro_ids <- repro_ids[keep_gene]
  }
  if (length(repro_ids) >= 1L) {
    x_train <- cbind(clin_mm[train_idx, , drop = FALSE],
                     t(expr_train[repro_ids, , drop = FALSE]))
    x_test <- cbind(clin_mm[test_idx, , drop = FALSE],
                    t(fd$expr_all[repro_ids, test_idx, drop = FALSE]))
    colnames(x_train) <- colnames(x_test) <- c(
      paste0("c", seq_len(ncol(clin_mm))),
      paste0("g", seq_along(repro_ids))
    )
    ridge <- tryCatch(with_random_seed(
      RESAMPLING$seed + 31000L + repeat_id * 10L + fold,
      fit_ridge_cox(
        x_train = x_train, x_test = x_test,
        time = train$os_time, event = train$os_event,
        foldid = event_stratified_folds(train$os_event, 5L),
        n_unpenalized = ncol(clin_mm)
      )
    ), error = function(e) e)
    if (inherits(ridge, "error")) {
      say("Ridge arm failed: ", conditionMessage(ridge))
      ridge <- NULL
      ridge_fallback <- TRUE
    }
  } else {
    ridge_fallback <- TRUE
  }

  clinical_risk <- cox_risk(clinical_fit, test)
  finish_arm <- function(fitted, gene, n_pass, n_model, rank = NA_integer_, lambda = NA_real_) {
    list(risk = fitted$risk, gene = gene, n_pass = n_pass, n_model = n_model,
         fallback = isTRUE(fitted$fallback) || is.na(gene),
         rank = rank, lambda = lambda)
  }
  arms <- list(
    clinical = list(risk = clinical_risk, gene = NA_character_, n_pass = 0L,
                    n_model = 0L, fallback = FALSE, rank = NA_integer_, lambda = NA_real_),
    hydra = finish_arm(
      predict_gene(clinical_fit, fd$expr_all, hydra$gene_id, train_idx, test_idx, train, test),
      hydra$gene_id, hydra$n_pass, as.integer(!is.na(hydra$gene_id))
    ),
    survival_only = finish_arm(
      predict_gene(clinical_fit, fd$expr_all, survival_arm$gene_id, train_idx, test_idx, train, test),
      survival_arm$gene_id, survival_arm$n_pass, as.integer(!is.na(survival_arm$gene_id)),
      rank = survival_arm$screen_rank
    ),
    de_only = finish_arm(
      predict_gene(clinical_fit, fd$expr_all, de_arm$gene_id, train_idx, test_idx, train, test),
      de_arm$gene_id, de_arm$n_pass, as.integer(!is.na(de_arm$gene_id))
    ),
    ridge_eligible = if (is.null(ridge)) {
      list(risk = clinical_risk, gene = NA_character_, n_pass = 0L, n_model = 0L,
           fallback = TRUE, rank = NA_integer_, lambda = NA_real_)
    } else {
      list(risk = list(lp = ridge$lp_test, risk3 = ridge$risk_test),
           gene = NA_character_, n_pass = ridge$n_genes, n_model = ridge$n_genes,
           fallback = FALSE, rank = NA_integer_, lambda = ridge$lambda)
    },
    matched_control = finish_arm(
      predict_gene(clinical_fit, fd$expr_all, control_id, train_idx, test_idx, train, test),
      control_id, as.integer(!is.na(control_id)), as.integer(!is.na(control_id))
    )
  )
  if (ridge_fallback) arms$ridge_eligible$fallback <- TRUE

  elapsed <- proc.time()[["elapsed"]] - started
  predictions <- bind_rows(lapply(names(arms), function(strategy_name) {
    arm <- arms[[strategy_name]]
    tibble(patient_barcode = test$patient_barcode,
           os_time = patients$os_time[test_idx],
           os_event = patients$os_event[test_idx],
           lp = arm$risk$lp,
           risk3 = arm$risk$risk3,
           repeat_id = repeat_id, fold = fold, strategy = strategy_name,
           selected_gene = arm$gene)
  }))
  fold_rows <- bind_rows(lapply(names(arms), function(strategy_name) {
    arm <- arms[[strategy_name]]
    tibble(repeat_id = repeat_id, fold = fold, strategy = strategy_name,
           selected_gene = arm$gene,
           selected_symbol = symbol_of(arm$gene),
           n_pass = arm$n_pass,
           n_model_genes = arm$n_model,
           fallback_to_clinical = arm$fallback,
           screen_rank = arm$rank,
           lambda = arm$lambda,
           outlier_replacement_fallback = fd$replacement_fallback,
           train_universe = nrow(fd$expr_all),
           train_repro_genes = nrow(fd$de),
           primary_selected_gene = primary_gene,
           primary_match = if (strategy_name == "hydra") TRUE else NA,
           elapsed_seconds = elapsed)
  }))
  saved <- list(spec = spec, rng_after_deseq = rng_after_deseq,
                predictions = predictions, fold_rows = fold_rows)
  write_rds_atomic(saved, path)
  rm(fd, expr_train)
  gc()
  saved
}

# The repeat seed is consumed only by the outer fold assignment, matching
# analysis/22_nested_cv.R. Later random draws restore the previous seed.
for (repeat_id in repeat_ids) {
  set.seed(RESAMPLING$seed + 22L + repeat_id)
  fold_id <- event_stratified_folds(patients$os_event, 5L)
  for (fold in fold_ids) run_fold(repeat_id, fold, fold_id)
}

expected <- expand.grid(fold = seq_len(5L), repeat_id = seq_len(n_repeats))
expected$present <- file.exists(mapply(checkpoint_path, expected$repeat_id, expected$fold))
if (worker_mode || !all(expected$present)) {
  say("Benchmark checkpoints saved. ", sum(expected$present), " of ",
      nrow(expected), " folds are complete.")
  quit(save = "no", status = 0L)
}

pieces <- lapply(seq_len(nrow(expected)), function(i) {
  saved <- readRDS(checkpoint_path(expected$repeat_id[i], expected$fold[i]))
  if (!identical(saved$spec, spec)) stop("Checkpoint spec changed while aggregating.")
  saved
})
predictions <- bind_rows(lapply(pieces, function(x) x$predictions))
fold_rows <- bind_rows(lapply(pieces, function(x) x$fold_rows))
if (any(fold_rows$strategy == "hydra" & !fold_rows$primary_match)) {
  stop("Aggregated HYDRA genes do not match the primary nested CV.")
}
if (anyDuplicated(predictions |> select(repeat_id, strategy, patient_barcode))) {
  stop("A patient appears twice inside one repeat and strategy.")
}

score_group <- function(df) {
  scored <- df |>
    group_by(repeat_id, strategy) |>
    group_modify(~ {
      s <- score_one(.x$os_time, .x$os_event, .x$lp, .x$risk3)
      tibble(c = s$c, brier3 = s$brier3, calibration_slope = s$calibration_slope)
    }) |>
    ungroup()
  clinical_scores <- scored |>
    filter(strategy == "clinical") |>
    select(repeat_id, clinical_c = c, clinical_brier3 = brier3,
           clinical_calibration_slope = calibration_slope)
  scored |>
    inner_join(clinical_scores, by = "repeat_id") |>
    mutate(delta_c = c - clinical_c, delta_brier3 = brier3 - clinical_brier3)
}
repeat_metrics <- score_group(predictions)
if (any(repeat_metrics$strategy == "clinical" & abs(repeat_metrics$delta_c) > 1e-8)) {
  stop("Clinical arm concordance increment is not zero.")
}

patient_ids <- patients$patient_barcode
strategies <- setdiff(unique(predictions$strategy), "clinical")
boot <- with_random_seed(RESAMPLING$seed + 31L, {
  bind_rows(lapply(seq_len(1000L), function(draw) {
    sampled <- sample(patient_ids, length(patient_ids), replace = TRUE)
    weights <- as.integer(table(factor(sampled, levels = patient_ids)))
    expand <- function(strategy) {
      x <- predictions |> filter(.data$strategy == .env$strategy)
      w <- weights[match(x$patient_barcode, patient_ids)]
      x[rep(seq_len(nrow(x)), w), , drop = FALSE]
    }
    clinical_x <- expand("clinical")
    bind_rows(lapply(strategies, function(strategy) {
      scored <- score_group(bind_rows(clinical_x |> mutate(strategy = "clinical"),
                                      expand(strategy)))
      tibble(draw = draw, strategy = strategy,
             delta_c = mean(scored$delta_c[scored$strategy == strategy]),
             delta_brier3 = mean(scored$delta_brier3[scored$strategy == strategy]))
    }))
  }))
})

summary <- repeat_metrics |>
  group_by(strategy) |>
  summarise(n_patients = n_distinct(predictions$patient_barcode),
            repeats = n_repeats,
            folds = 5L,
            mean_c = mean(c),
            mean_brier3 = mean(brier3),
            mean_calibration_slope = mean(calibration_slope),
            mean_delta_c = mean(delta_c),
            mean_delta_brier3 = mean(delta_brier3),
            .groups = "drop")
boot_summary <- boot |>
  group_by(strategy) |>
  summarise(patient_bootstrap_ci_low = quantile(delta_c, 0.025),
            patient_bootstrap_ci_high = quantile(delta_c, 0.975),
            patient_bootstrap_brier_ci_low = quantile(delta_brier3, 0.025),
            patient_bootstrap_brier_ci_high = quantile(delta_brier3, 0.975),
            .groups = "drop")
gene_summary <- fold_rows |>
  group_by(strategy) |>
  summarise(folds_with_gene = sum(!is.na(selected_gene)),
            distinct_genes = n_distinct(selected_gene[!is.na(selected_gene)]),
            fallback_folds = sum(fallback_to_clinical),
            .groups = "drop")
summary <- summary |>
  left_join(boot_summary, by = "strategy") |>
  left_join(gene_summary, by = "strategy") |>
  mutate(patient_bootstrap_ci_low = if_else(strategy == "clinical", 0, patient_bootstrap_ci_low),
         patient_bootstrap_ci_high = if_else(strategy == "clinical", 0, patient_bootstrap_ci_high),
         survival_screen_n = screen_n,
         primary_hydra_mismatches = 0L)

save_csv <- function(data, path) {
  if (file.exists(path)) unlink(path)
  write_csv_atomic(data, path)
}
save_csv(predictions, file.path(DIRS$tables, "nested_benchmark_predictions.csv"))
save_csv(fold_rows, file.path(DIRS$tables, "nested_benchmark_folds.csv"))
save_csv(repeat_metrics, file.path(DIRS$tables, "nested_benchmark_repeat_metrics.csv"))
save_csv(boot, file.path(DIRS$tables, "nested_benchmark_patient_bootstrap.csv"))
save_csv(summary, file.path(DIRS$tables, "nested_benchmark_summary.csv"))

plot_df <- repeat_metrics |>
  filter(strategy != "clinical") |>
  mutate(strategy = factor(strategy, levels = c(
    "hydra", "survival_only", "de_only", "ridge_eligible", "matched_control"
  ), labels = c(
    "HYDRA", "Survival only", "DE only", "Ridge on eligible DEGs", "Matched control"
  )))
figure <- ggplot(plot_df, aes(strategy, delta_c)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = 0.01, color = "#9b3d2a", linetype = "dashed", linewidth = 0.5) +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 1),
             color = "#216b75", size = 2) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.45, linewidth = 0.4,
               color = "#173f46", fatten = 0) +
  labs(x = NULL, y = "Held-out concordance minus clinical",
       title = "Nested selection benchmark",
       subtitle = "Same outer splits. Dashed line: prespecified 0.01 increment. Points are repeats.") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))
ggsave(file.path(DIRS$figures, "nested_selection_benchmark.png"), figure,
       width = 8.4, height = 4.8, dpi = 180)
say("Nested selection benchmark complete.")
