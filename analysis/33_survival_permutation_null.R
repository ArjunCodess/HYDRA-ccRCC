# Permutation null for the high-confidence survival gate.
# Reuses the repeat-1 nested DESeq2 fits. Does not replace the exponential null.
# Seed: RESAMPLING$seed + 33.

source("analysis/00_config.R")
source("analysis/functions/io.R")
source("analysis/functions/tcga_metadata.R")
source("analysis/functions/patient_samples.R")

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(dplyr)
  library(readr)
  library(tibble)
  library(survival)
})

n_perm <- as.integer(Sys.getenv("HYDRA_PERMUTATIONS", "200"))
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

gene_fit <- function(dat, expr, covariates) {
  dat$expr <- expr
  tryCatch(suppressWarnings(coxph(
    as.formula(paste("Surv(os_time, os_event) ~", paste(c("expr", covariates), collapse = "+"))),
    data = dat
  )), error = function(e) NULL)
}

select_gene <- function(fd, train_idx, times, events) {
  de <- fd$de
  if (nrow(de) == 0L) return(tibble(gene_id = NA_character_, main_fdr = NA_real_, main_beta = NA_real_))
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
  if (nrow(de) == 0L) return(tibble(gene_id = NA_character_, main_fdr = NA_real_, main_beta = NA_real_))
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
  if (nrow(de) == 0L) {
    tibble(gene_id = NA_character_, main_fdr = NA_real_, main_beta = NA_real_)
  } else {
    de |> slice(1) |> select(gene_id, main_fdr, main_beta)
  }
}

folds <- lapply(seq_len(5L), function(fold) {
  path <- sprintf("data/processed/nested_cv_fold_checkpoints/repeat_01_fold_%d.rds", fold)
  saved <- readRDS(path)
  if (is.null(saved$cache$fd)) stop("Repeat-1 fold cache is missing: ", path)
  saved
})
observed <- select_gene(folds[[1]]$cache$fd, folds[[1]]$cache$train,
                        patients$os_time, patients$os_event)
if (!identical(observed$gene_id, folds[[1]]$summary$selected_gene)) {
  stop("Permutation-null selector does not reproduce repeat-1 fold 1. Got ",
       observed$gene_id, " expected ", folds[[1]]$summary$selected_gene)
}
message("Selector check matched ", observed$gene_id)

set.seed(RESAMPLING$seed + 33L)
rows <- vector("list", n_perm)
for (i in seq_len(n_perm)) {
  fold <- ((i - 1L) %% 5L) + 1L
  saved <- folds[[fold]]
  train <- saved$cache$train
  order_train <- sample(train)
  times <- patients$os_time
  events <- patients$os_event
  times[train] <- patients$os_time[order_train]
  events[train] <- patients$os_event[order_train]
  picked <- select_gene(saved$cache$fd, train, times, events)
  rows[[i]] <- tibble(simulation = i, fold = fold, gene_id = picked$gene_id,
                      main_fdr = picked$main_fdr, main_beta = picked$main_beta,
                      gene_selected = !is.na(picked$gene_id))
  if (i %% 10L == 0L) message("Permutation ", i, " of ", n_perm)
}
result <- bind_rows(rows)
write_csv_atomic(result, file.path(DIRS$tables, "survival_permutation_null.csv"))
selected <- result |> filter(gene_selected)
write_csv_atomic(tibble(
  metric = c("simulations", "folds_cycled", "gene_selected", "seed"),
  value = c(nrow(result), 5, sum(result$gene_selected), RESAMPLING$seed + 33L)
), file.path(DIRS$tables, "survival_permutation_null_summary.csv"))
message("Permutation null selected a gene in ", sum(result$gene_selected), " of ", n_perm, " simulations.")
