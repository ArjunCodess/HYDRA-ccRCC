# The same clinical-only null as stage 22, with iteration-level checkpoints.
# Load its exact setup/functions without running the outer evaluation or writers.
for (expression in parse("analysis/22_nested_cv.R")) {
  if (identical(expression[[1]], as.name("<-")) && identical(expression[[2]], as.name("predictions"))) break
  eval(expression, envir = .GlobalEnv)
}
n_null <- 200L
rm(se, counts, tumor_counts, geo, geo_a, geo_b, symbols, meta, all_normal, clinical)
gc()
paths <- sprintf("data/processed/nested_cv_fold_checkpoints/repeat_01_fold_%d.rds", 1:5)
first_fold_cache <- lapply(paths, function(path) {
  saved <- readRDS(path)
  stopifnot(!is.null(saved$cache), !is.null(saved$normalization_correction))
  saved$cache
})
signature <- list(n = n_null, r_version = R.version.string,
  files = unname(tools::md5sum(c("analysis/50_clinical_null_refits.R", "analysis/22_nested_cv.R",
    "analysis/functions/frozen_normalization.R", "analysis/00_config.R", "environment/package_versions.csv",
    FILES$tcga_se, FILES$tcga_clinical, paths))))
checkpoint <- file.path(DIRS$processed, "frozen_clinical_null_iterations.rds")
saved <- if (file.exists(checkpoint)) readRDS(checkpoint) else NULL
valid <- !is.null(saved) && identical(saved$signature, signature)
rows <- if (valid) saved$rows else list()
baseline <- coxph(Surv(os_time, os_event) ~ age + sex + stage + grade, data = patients)
null_lp <- as.numeric(predict(baseline, newdata = patients, type = "lp", reference = "zero"))
base_hazard <- tail(basehaz(baseline, centered = FALSE)$hazard, 1)
censor_times <- patients$os_time[patients$os_event == 0L]
set.seed(RESAMPLING$seed + 24L)
if (valid) assign(".Random.seed", saved$rng_after, envir = .GlobalEnv)
if (length(rows)) message("Resumed ", length(rows), " completed clinical-null iterations.")
if (length(rows) < n_null) for (sim in seq.int(length(rows) + 1L, n_null)) {
  simulated_event_time <- -log(runif(nrow(patients))) / (base_hazard * exp(null_lp)) * max(patients$os_time)
  simulated_censor <- sample(censor_times, nrow(patients), replace = TRUE)
  sim_time <- pmin(simulated_event_time, simulated_censor)
  sim_event <- as.integer(simulated_event_time <= simulated_censor)
  fold <- ((sim - 1L) %% 5L) + 1L
  f <- first_fold_cache[[fold]]
  selected <- select_gene(f$fd, f$train, sim_time, sim_event)
  sim_pred <- predict_fold(f$fd, f$train, f$test, selected, sim_time, sim_event)
  rows[[sim]] <- score_repeat(sim_pred) |> mutate(simulation = sim, fold = fold, selected_gene = selected)
  write_rds_atomic(list(signature = signature, rows = rows, rng_after = .Random.seed), checkpoint)
  if (sim %% 10L == 0L) message("Clinical-only null simulation ", sim, "/", n_null)
}
result <- bind_rows(rows)
stopifnot(nrow(result) == 200L, identical(result$simulation, 1:200))
original <- read_csv("results/archive/normalization_20261004/nested_cv_clinical_null.csv", show_col_types = FALSE)
stopifnot(identical(result$selected_gene, original$selected_gene), identical(result$fold, as.integer(original$fold)))
for (column in c("clinical_c", "clinical_brier3", "clinical_calibration_slope")) {
  stopifnot(max(abs(result[[column]] - original[[column]])) < 1e-10)
}
write_csv_atomic(result, file.path(DIRS$tables, "nested_cv_clinical_null.csv"))
message("Completed clinical-only null with unchanged DGP, selector, predictions and seed; checkpointed iterations.")
