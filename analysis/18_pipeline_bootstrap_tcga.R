source("analysis/00_config.R")
source("analysis/functions/io.R")
source("analysis/functions/tcga_metadata.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(survival)
  library(broom)
})

vst_mat <- read_required_rds(FILES$tcga_vst)
coldata <- read_csv(FILES$tcga_coldata, show_col_types = FALSE)
clinical <- read_csv(FILES$tcga_clinical, show_col_types = FALSE)
candidates <- read_csv(
  file.path(DIRS$tables, "candidate_gene_evidence_table.csv"),
  show_col_types = FALSE
)

clinical_surv <- clinical |>
  transmute(
    patient_barcode = submitter_id,
    os_time,
    os_event,
    age = suppressWarnings(as.numeric(age_at_diagnosis)) / 365.25,
    sex = factor(gender),
    stage = factor(normalize_stage(ajcc_pathologic_stage)),
    grade = factor(case_when(
      normalize_grade(tumor_grade) %in% c("G1", "G2") ~ "Low grade",
      normalize_grade(tumor_grade) %in% c("G3", "G4") ~ "High grade",
      TRUE ~ NA_character_
    ))
  )

sample_data <- coldata |>
  filter(sample_type == "Primary Tumor") |>
  mutate(patient_barcode = tcga_patient_barcode(sample_barcode)) |>
  inner_join(clinical_surv, by = "patient_barcode") |>
  filter(
    !is.na(os_time),
    os_time > 0,
    !is.na(os_event),
    !is.na(age),
    !is.na(sex),
    !is.na(stage),
    !is.na(grade),
    sample_barcode %in% colnames(vst_mat)
  ) |>
  distinct(patient_barcode, .keep_all = TRUE) |>
  mutate(across(c(sex, stage, grade), droplevels))

gene_table <- candidates |>
  filter(reproducible_deg, tcga_gene_id %in% rownames(vst_mat)) |>
  select(
    symbol,
    tcga_gene_id,
    gse40435_log2fc,
    gse53757_log2fc,
    full_main_log_hr = main_log_hr,
    frozen_high_confidence = high_confidence_candidate
  )

fit_expr <- function(gene_id, indices, covars, ph = FALSE) {
  barcodes <- sample_data$sample_barcode[indices]
  raw_expr <- as.numeric(vst_mat[gene_id, barcodes])
  expr_sd <- sd(raw_expr)
  if (!is.finite(expr_sd) || expr_sd <= 0) return(NULL)

  dat <- sample_data[indices, ] |>
    mutate(expr = (raw_expr - mean(raw_expr)) / expr_sd) |>
    select(os_time, os_event, expr, all_of(covars)) |>
    mutate(across(any_of(c("sex", "stage", "grade")), droplevels))

  factors <- intersect(c("sex", "stage", "grade"), covars)
  if (nrow(dat) < 100 || sum(dat$os_event) < 25 ||
      any(vapply(dat[factors], nlevels, integer(1)) < 2)) return(NULL)

  form <- as.formula(paste("Surv(os_time, os_event) ~ expr +", paste(covars, collapse = " + ")))
  fit <- tryCatch(coxph(form, data = dat), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  term <- tidy(fit) |> filter(term == "expr")
  if (nrow(term) != 1 || !is.finite(term$estimate) || !is.finite(term$p.value)) return(NULL)

  ph_p <- NA_real_
  if (ph) {
    zph <- tryCatch(cox.zph(fit), error = function(e) NULL)
    if (!is.null(zph) && "expr" %in% rownames(zph$table)) ph_p <- zph$table["expr", "p"]
  }

  tibble(
    log_hr = term$estimate,
    p_value = term$p.value,
    ph_p_value = ph_p,
    n = fit$n,
    events = fit$nevent
  )
}

fit_oob_increment <- function(gene_id, inbag_indices, oob_indices) {
  train <- sample_data[inbag_indices, ]
  test <- sample_data[oob_indices, ]
  train_raw <- as.numeric(vst_mat[gene_id, train$sample_barcode])
  test_raw <- as.numeric(vst_mat[gene_id, test$sample_barcode])
  train_sd <- sd(train_raw)
  if (!is.finite(train_sd) || train_sd <= 0 || nrow(test) < 50 || sum(test$os_event) < 15) {
    return(NULL)
  }

  train$expr <- (train_raw - mean(train_raw)) / train_sd
  test$expr <- (test_raw - mean(train_raw)) / train_sd

  clinical_fit <- tryCatch(
    coxph(Surv(os_time, os_event) ~ age + sex + stage + grade, data = train, x = TRUE),
    error = function(e) NULL
  )
  gene_fit <- tryCatch(
    coxph(Surv(os_time, os_event) ~ expr + age + sex + stage + grade, data = train, x = TRUE),
    error = function(e) NULL
  )
  if (is.null(clinical_fit) || is.null(gene_fit)) return(NULL)

  clinical_lp <- tryCatch(
    predict(clinical_fit, newdata = test, type = "lp", reference = "zero"),
    error = function(e) rep(NA_real_, nrow(test))
  )
  gene_lp <- tryCatch(
    predict(gene_fit, newdata = test, type = "lp", reference = "zero"),
    error = function(e) rep(NA_real_, nrow(test))
  )
  complete <- is.finite(clinical_lp) & is.finite(gene_lp)
  if (sum(complete) < 50 || sum(test$os_event[complete]) < 15) return(NULL)

  clinical_c <- unname(concordance(
    Surv(test$os_time[complete], test$os_event[complete]) ~ clinical_lp[complete],
    reverse = TRUE
  )$concordance)
  gene_c <- unname(concordance(
    Surv(test$os_time[complete], test$os_event[complete]) ~ gene_lp[complete],
    reverse = TRUE
  )$concordance)

  tibble(
    oob_n = sum(complete),
    oob_events = sum(test$os_event[complete]),
    clinical_c_index = clinical_c,
    clinical_gene_c_index = gene_c,
    delta_c_index = gene_c - clinical_c
  )
}

draw_stratified_bootstrap <- function(event) {
  unlist(lapply(sort(unique(event)), function(value) {
    idx <- which(event == value)
    sample(idx, length(idx), replace = TRUE)
  }), use.names = FALSE)
}

set.seed(V2_RESAMPLING$seed + 2L)
repeat_results <- vector("list", V2_RESAMPLING$pipeline_bootstrap_repeats)
repeat_summary <- vector("list", V2_RESAMPLING$pipeline_bootstrap_repeats)
oob_results <- vector("list", V2_RESAMPLING$pipeline_bootstrap_repeats)

detected_cores <- parallel::detectCores(logical = FALSE)
if (!is.finite(detected_cores)) detected_cores <- 2L
bootstrap_workers <- max(1L, min(4L, detected_cores - 1L))
bootstrap_cluster <- NULL
if (bootstrap_workers > 1L) {
  bootstrap_cluster <- parallel::makeCluster(bootstrap_workers)
  local_lib_path <- normalizePath(LOCAL_R_LIB, winslash = "/", mustWork = TRUE)
  parallel::clusterExport(
    bootstrap_cluster,
    "local_lib_path",
    envir = environment()
  )
  parallel::clusterEvalQ(bootstrap_cluster, {
    .libPaths(c(local_lib_path, .libPaths()))
    suppressPackageStartupMessages({
      library(dplyr)
      library(survival)
      library(broom)
    })
  })
  parallel::clusterExport(
    bootstrap_cluster,
    c("vst_mat", "sample_data", "fit_expr"),
    envir = environment()
  )
}

fit_main_screen <- function(gene_ids, indices) {
  worker <- function(gene_id, sampled_indices) {
    out <- fit_expr(gene_id, sampled_indices, c("age", "sex", "stage", "grade"))
    if (is.null(out)) return(NULL)
    mutate(out, tcga_gene_id = gene_id)
  }
  if (is.null(bootstrap_cluster)) {
    return(bind_rows(lapply(gene_ids, worker, sampled_indices = indices)))
  }
  bind_rows(parallel::parLapply(
    bootstrap_cluster,
    gene_ids,
    worker,
    sampled_indices = indices
  ))
}

for (repeat_id in seq_len(V2_RESAMPLING$pipeline_bootstrap_repeats)) {
  inbag_indices <- draw_stratified_bootstrap(sample_data$os_event)
  oob_indices <- setdiff(seq_len(nrow(sample_data)), unique(inbag_indices))
  message(
    "Pipeline bootstrap repeat ",
    repeat_id,
    "/",
    V2_RESAMPLING$pipeline_bootstrap_repeats
  )

  main <- fit_main_screen(gene_table$tcga_gene_id, inbag_indices) |>
    select(-ph_p_value) |>
    mutate(main_fdr = p.adjust(p_value, method = "BH"))

  eligible <- main |>
    filter(
      main_fdr < THRESHOLDS$high_confidence_survival_fdr,
      abs(log_hr) >= THRESHOLDS$high_confidence_abs_log_hr
    ) |>
    inner_join(gene_table, by = "tcga_gene_id") |>
    filter(
      abs(gse40435_log2fc) >= THRESHOLDS$min_geo_abs_log2fc,
      abs(gse53757_log2fc) >= THRESHOLDS$min_geo_abs_log2fc
    )

  hardened <- tibble(
    tcga_gene_id = character(),
    ph_p_value = double(),
    stage_log_hr = double(),
    stage_p_value = double(),
    grade_log_hr = double(),
    grade_p_value = double()
  )
  if (nrow(eligible) > 0) {
    hardened <- bind_rows(hardened, lapply(seq_len(nrow(eligible)), function(i) {
      row <- eligible[i, ]
      main_ph <- fit_expr(row$tcga_gene_id, inbag_indices, c("age", "sex", "stage", "grade"), ph = TRUE)
      stage_fit <- fit_expr(row$tcga_gene_id, inbag_indices, c("age", "sex", "stage"))
      grade_fit <- fit_expr(row$tcga_gene_id, inbag_indices, c("age", "sex", "grade"))
      if (is.null(main_ph) || is.null(stage_fit) || is.null(grade_fit)) return(NULL)
      tibble(
        tcga_gene_id = row$tcga_gene_id,
        ph_p_value = main_ph$ph_p_value,
        stage_log_hr = stage_fit$log_hr,
        stage_p_value = stage_fit$p_value,
        grade_log_hr = grade_fit$log_hr,
        grade_p_value = grade_fit$p_value
      )
    }))
  }

  one_repeat <- main |>
    left_join(gene_table, by = "tcga_gene_id") |>
    left_join(hardened, by = "tcga_gene_id") |>
    mutate(
      repeat_id = repeat_id,
      main_screen_pass = main_fdr < THRESHOLDS$high_confidence_survival_fdr &
        abs(log_hr) >= THRESHOLDS$high_confidence_abs_log_hr,
      geo_effect_support = abs(gse40435_log2fc) >= THRESHOLDS$min_geo_abs_log2fc &
        abs(gse53757_log2fc) >= THRESHOLDS$min_geo_abs_log2fc,
      ph_pass = !is.na(ph_p_value) & ph_p_value >= THRESHOLDS$ph_min_p,
      sensitivity_pass = !is.na(stage_p_value) & !is.na(grade_p_value) &
        stage_p_value < 0.05 & grade_p_value < 0.05 &
        sign(stage_log_hr) == sign(log_hr) & sign(grade_log_hr) == sign(log_hr),
      selected = main_screen_pass & geo_effect_support & ph_pass & sensitivity_pass
    )

  repeat_results[[repeat_id]] <- one_repeat |>
    select(
      repeat_id,
      symbol,
      tcga_gene_id,
      frozen_high_confidence,
      log_hr,
      p_value,
      main_fdr,
      main_screen_pass,
      geo_effect_support,
      ph_pass,
      sensitivity_pass,
      selected
    )

  selected_frozen <- one_repeat |>
    filter(frozen_high_confidence, selected)
  if (nrow(selected_frozen) > 0 && length(oob_indices) > 0) {
    oob_results[[repeat_id]] <- bind_rows(lapply(seq_len(nrow(selected_frozen)), function(i) {
      row <- selected_frozen[i, ]
      out <- fit_oob_increment(row$tcga_gene_id, inbag_indices, oob_indices)
      if (is.null(out)) return(NULL)
      mutate(out, repeat_id = repeat_id, symbol = row$symbol, .before = 1)
    }))
  }

  repeat_summary[[repeat_id]] <- tibble(
    repeat_id,
    inbag_draws = length(inbag_indices),
    unique_inbag_patients = length(unique(inbag_indices)),
    oob_patients = length(oob_indices),
    oob_events = sum(sample_data$os_event[oob_indices]),
    modeled_genes = nrow(main),
    main_screen_genes = sum(one_repeat$main_screen_pass, na.rm = TRUE),
    selected_genes = sum(one_repeat$selected, na.rm = TRUE),
    selected_frozen_candidates = sum(
      one_repeat$selected & one_repeat$frozen_high_confidence,
      na.rm = TRUE
    )
  )
}

if (!is.null(bootstrap_cluster)) parallel::stopCluster(bootstrap_cluster)

all_repeats <- bind_rows(repeat_results)
all_oob <- bind_rows(oob_results)

bootstrap_stability <- all_repeats |>
  group_by(symbol, tcga_gene_id, frozen_high_confidence) |>
  summarise(
    repeats_modeled = n(),
    main_screen_frequency = mean(main_screen_pass, na.rm = TRUE),
    selection_frequency = mean(selected, na.rm = TRUE),
    median_log_hr = median(log_hr, na.rm = TRUE),
    log_hr_ci_low = quantile(log_hr, 0.025, na.rm = TRUE),
    log_hr_ci_high = quantile(log_hr, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(gene_table |> select(symbol, full_main_log_hr), by = "symbol") |>
  mutate(full_fit_direction_agreement = sign(median_log_hr) == sign(full_main_log_hr)) |>
  arrange(desc(selection_frequency), desc(main_screen_frequency))

oob_summary <- all_repeats |>
  filter(frozen_high_confidence) |>
  group_by(symbol) |>
  summarise(
    bootstrap_repeats = n(),
    selected_repeats = sum(selected, na.rm = TRUE),
    selection_frequency = mean(selected, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(
    all_oob |>
      group_by(symbol) |>
      summarise(
        oob_evaluable_repeats = n(),
        mean_oob_clinical_c_index = mean(clinical_c_index),
        mean_oob_clinical_gene_c_index = mean(clinical_gene_c_index),
        mean_oob_delta_c_index = mean(delta_c_index),
        median_oob_delta_c_index = median(delta_c_index),
        oob_delta_ci_low = quantile(delta_c_index, 0.025),
        oob_delta_ci_high = quantile(delta_c_index, 0.975),
        positive_oob_delta_frequency = mean(delta_c_index > 0),
        .groups = "drop"
      ),
    by = "symbol"
  ) |>
  arrange(desc(selection_frequency), desc(mean_oob_delta_c_index))

write_csv_atomic(
  bootstrap_stability,
  file.path(DIRS$tables, "pipeline_bootstrap_all_genes.csv")
)
write_csv_atomic(
  bootstrap_stability |> filter(frozen_high_confidence),
  file.path(DIRS$tables, "pipeline_bootstrap_frozen_candidates.csv")
)
write_csv_atomic(
  bind_rows(repeat_summary),
  file.path(DIRS$tables, "pipeline_bootstrap_repeats.csv")
)
write_csv_atomic(
  all_oob,
  file.path(DIRS$tables, "pipeline_bootstrap_oob_repeats.csv")
)
write_csv_atomic(
  oob_summary,
  file.path(DIRS$tables, "pipeline_bootstrap_oob_summary.csv")
)

message(
  "Pipeline bootstrap complete. Repeats: ",
  V2_RESAMPLING$pipeline_bootstrap_repeats,
  "; frozen candidates assessed: ",
  sum(bootstrap_stability$frozen_high_confidence)
)
