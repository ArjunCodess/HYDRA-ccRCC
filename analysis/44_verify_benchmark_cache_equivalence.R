# Recover only proven source-equivalent benchmark checkpoints. No model refit
# is claimed; original files/specifications and all model contents are retained.
source("analysis/00_config.R")
source("analysis/functions/io.R")
suppressPackageStartupMessages({ library(readr); library(dplyr) })

cache_dir <- file.path(DIRS$processed, "nested_benchmark_checkpoints")
paths <- list.files(cache_dir, pattern = "^repeat_[0-9]+_fold_[0-9]+[.]rds$", full.names = TRUE)
if (!length(paths)) {
  write_csv_atomic(tibble(checkpoint = NA_character_, status = "no_cached_benchmark_models; full_refit_required"),
    file.path(DIRS$tables, "review_benchmark_cache_audit.csv"))
  message("No benchmark checkpoints: the normal benchmark stage will refit.")
  quit(status = 0L)
}

git_source <- function(ref, path) {
  x <- system2("git", c("show", paste0(ref, ":", path)), stdout = TRUE)
  if (!is.null(attr(x, "status"))) {
    message("Historical source unavailable; equivalence cannot be established: ", path)
    return("")
  }
  paste(x, collapse = "\n")
}
ast <- function(text) parse(text = text, keep.source = FALSE)
current_source <- function(path) paste(readLines(path, warn = FALSE), collapse = "\n")
old_script_ref <- "fdc07c079638460197c25461dd41a4d4ba0fd6ae"
old_benchmark_ref <- "4f656e84829565e3deb52df702e1a68b5133e973"
old_patient_ref <- "8d67f86cc3be059127bf1e2806368cadde3db399"
script_path <- "analysis/31_nested_selection_benchmark.R"
benchmark_path <- "analysis/functions/nested_benchmark.R"
patient_path <- "analysis/functions/patient_samples.R"

# Normalize precisely the documented output name and the fold-policy replay.
# Reuse additionally requires matching recorded replacement settings below.
normalized_script <- current_source(script_path)
normalized_script <- gsub(", force_replacement_fallback = FALSE", "", normalized_script, fixed = TRUE)
normalized_script <- gsub("replacement_fallback <- force_replacement_fallback",
                          "replacement_fallback <- FALSE", normalized_script, fixed = TRUE)
normalized_script <- gsub("if (!force_replacement_fallback)", "if (attempt == 1L)",
                          normalized_script, fixed = TRUE)
normalized_script <- gsub('    say("DESeq2 trimmed-mean failure in fold, attempt ", attempt,',
  '    replacement_fallback <- TRUE\n    say("DESeq2 trimmed-mean failure in fold, attempt ", attempt,',
  normalized_script, fixed = TRUE)
normalized_script <- gsub("; retrying the recorded primary-fold replacement setting.",
  "; retrying without outlier replacement.", normalized_script, fixed = TRUE)
policy_block <- paste(c(
  "  # Every comparator must use the replacement setting actually used by the",
  "  # primary fit. A transient runtime error must not change that setting.",
  "  expected_fallback <- primary_folds$outlier_replacement_fallback[",
  "    primary_folds$repeat_id == repeat_id & primary_folds$fold == fold",
  "  ]",
  "  if (length(expected_fallback) != 1L || is.na(expected_fallback)) {",
  '    stop("Primary nested-CV replacement setting is missing.")',
  "  }",
  "  fd <- make_fold_data(train_idx, test_idx, force_replacement_fallback = expected_fallback)"
), collapse = "\n")
normalized_script <- gsub(policy_block, "  fd <- make_fold_data(train_idx, test_idx)",
                          normalized_script, fixed = TRUE)
script_equivalent <- identical(
  ast(gsub("nested_selection_benchmark_five_arm.png", "nested_selection_benchmark.png",
           normalized_script, fixed = TRUE)),
  ast(git_source(old_script_ref, script_path)))

# Comments do not occur in the AST. Remove precisely the added validation
# guard, which cannot alter an already successful glmnet inner-CV fit.
guard <- ast('if (length(unique(foldid)) < 2L || any(foldid < 1L)) stop("Ridge inner folds must be a training-only partition.")')[[1]]
benchmark_ast <- ast(current_source(benchmark_path))
removed_guard <- 0L
for (i in seq_along(benchmark_ast)) {
  e <- benchmark_ast[[i]]
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      identical(e[[2]], as.name("fit_ridge_cox"))) {
    fn <- e[[3]]
    body <- as.list(fn[[3]])
    remove <- vapply(body, identical, logical(1), y = guard)
    removed_guard <- sum(remove)
    fn[[3]] <- as.call(body[!remove])
    e[[3]] <- fn
    benchmark_ast[[i]] <- e
  }
}
benchmark_equivalent <- removed_guard == 1L && identical(
  benchmark_ast, ast(git_source(old_benchmark_ref, benchmark_path)))

# The added aliquot helper is unused by this benchmark; all previous patient
# selection functions must have identical parsed bodies.
patient_ast <- ast(current_source(patient_path))
added <- vapply(as.list(patient_ast), function(e) is.call(e) &&
  identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name("apply_aliquot_rule")), logical(1))
patient_equivalent <- sum(added) == 1L &&
  !"apply_aliquot_rule" %in% all.names(ast(current_source(script_path))) &&
  identical(as.expression(as.list(patient_ast)[!added]), ast(git_source(old_patient_ref, patient_path)))
source_equivalent <- script_equivalent && benchmark_equivalent && patient_equivalent

inputs <- c(script_path, benchmark_path, "analysis/00_config.R", patient_path,
            FILES$tcga_se, FILES$tcga_clinical,
            file.path(DIRS$tables, "gse40435_limma_tumor_vs_normal.csv"),
            file.path(DIRS$tables, "gse53757_limma_tumor_vs_normal.csv"))
current_spec <- list(files = unname(tools::md5sum(inputs)),
  screen_n = as.integer(Sys.getenv("HYDRA_SURVIVAL_SCREEN", "2000")),
  seed = RESAMPLING$seed, ridge = "glmnet-5.0-efron-alpha0-lambda.min")
legacy_files <- current_spec$files
legacy_files[c(1, 2, 4)] <- c("a694218eec9cb01e3b82ae97718b7325",
  "878bc1b93c77f168445c39d7f623fa8b", "53e38286aa6f2023a855152f46b03b83")
legacy_spec <- current_spec
legacy_spec$files <- legacy_files
primary_inputs <- c("analysis/22_nested_cv.R", "analysis/00_config.R", patient_path,
  FILES$tcga_se, FILES$tcga_clinical, inputs[7:8])
primary_signature <- unname(tools::md5sum(primary_inputs))
rows <- vector("list", length(paths))
eligible <- logical(length(paths))

for (i in seq_along(paths)) {
  saved <- readRDS(paths[i])
  primary_path <- file.path(DIRS$processed, "nested_cv_fold_checkpoints", basename(paths[i]))
  primary_current <- FALSE; prediction_agreement <- FALSE; replacement_agreement <- FALSE
  ridge <- saved$fold_rows[saved$fold_rows$strategy == "ridge_eligible", ]
  successful_ridge <- nrow(ridge) == 1L && !ridge$fallback_to_clinical &&
    is.finite(ridge$lambda) && ridge$lambda > 0 && ridge$n_model_genes > 0
  if (file.exists(primary_path)) {
    primary <- readRDS(primary_path)
    primary_current <- identical(primary$signature, primary_signature)
    if (primary_current) {
      replacement_agreement <- all(saved$fold_rows$outlier_replacement_fallback ==
        primary$summary$outlier_replacement_fallback)
      a <- primary$prediction
      comparisons <- c(clinical = "clinical", hydra = "gene")
      prediction_agreement <- all(vapply(names(comparisons), function(strategy) {
        b <- saved$predictions[saved$predictions$strategy == strategy, ]
        if (nrow(b) != nrow(a) || anyDuplicated(b$patient_barcode) ||
            !setequal(a$patient_barcode, b$patient_barcode)) return(FALSE)
        b <- b[match(a$patient_barcode, b$patient_barcode), ]
        prefix <- comparisons[[strategy]]
        isTRUE(all.equal(a$os_time, b$os_time, tolerance = 0)) &&
          identical(a$os_event, b$os_event) &&
          isTRUE(all.equal(a[[paste0(prefix, "_lp")]], b$lp, tolerance = 1e-12)) &&
          isTRUE(all.equal(a[[paste0(prefix, "_risk3")]], b$risk3, tolerance = 1e-12)) &&
          (strategy != "hydra" || identical(a$selected_gene, b$selected_gene))
      }, logical(1)))
    }
  }
  previous_audit <- saved$source_equivalence_audit
  legacy_saved <- saved
  if (!is.null(previous_audit)) {
    original <- readRDS(previous_audit$backup)
    model_names <- setdiff(names(original), "spec")
    stopifnot(identical(original[model_names], saved[model_names]))
    legacy_saved <- original
  }
  legacy_match <- identical(legacy_saved$spec, legacy_spec)
  already_current <- identical(saved$spec, current_spec)
  if (already_current && (!primary_current || !prediction_agreement || !replacement_agreement)) {
    stop("Current-signature benchmark checkpoint fails fresh primary-fold checks; refit required: ", paths[i])
  }
  eligible[i] <- source_equivalent && legacy_match && primary_current &&
    prediction_agreement && successful_ridge && replacement_agreement && !already_current
  rows[[i]] <- tibble(checkpoint = paths[i],
    audit_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    audit_script_md5 = unname(tools::md5sum("analysis/44_verify_benchmark_cache_equivalence.R")),
    original_md5 = if (is.null(previous_audit)) unname(tools::md5sum(paths[i])) else previous_audit$original_md5,
    current_md5 = unname(tools::md5sum(paths[i])),
    script_ast_equivalent = script_equivalent, benchmark_ast_equivalent = benchmark_equivalent,
    patient_ast_equivalent = patient_equivalent, legacy_spec_matches = legacy_match,
    current_spec_matches = already_current, primary_signature_current = primary_current,
    primary_predictions_agree = prediction_agreement, successful_ridge = successful_ridge,
    recorded_replacement_setting_agrees = replacement_agreement,
    status = if (already_current && !is.null(previous_audit)) "current_signature; verified_equivalent_cache; not_refitted" else if (already_current) "current_signature" else if (eligible[i]) "eligible_for_verified_migration" else "requires_refit",
    original_source_refs = paste(old_script_ref, old_benchmark_ref, old_patient_ref, sep = ";"))
  if (!is.null(previous_audit)) {
    rows[[i]]$backup <- previous_audit$backup
    rows[[i]]$backup_md5 <- unname(tools::md5sum(previous_audit$backup))
    rows[[i]]$migrated_md5 <- unname(tools::md5sum(paths[i]))
    stopifnot(identical(rows[[i]]$backup_md5, previous_audit$original_md5))
  } else {
    historical_backup <- file.path(DIRS$processed, "benchmark_cache_originals", basename(paths[i]))
    if (file.exists(historical_backup)) {
      rows[[i]]$backup <- historical_backup
      rows[[i]]$backup_md5 <- unname(tools::md5sum(historical_backup))
    }
  }
}

audit_only <- "--audit-only" %in% commandArgs(trailingOnly = TRUE)
if (!audit_only) {
  backup_dir <- file.path(DIRS$processed, "benchmark_cache_originals")
  dir.create(backup_dir, showWarnings = FALSE, recursive = TRUE)
  for (i in which(eligible)) {
    saved <- readRDS(paths[i])
    backup <- file.path(backup_dir, basename(paths[i]))
    if (file.exists(backup)) {
      stopifnot(identical(unname(tools::md5sum(backup)), rows[[i]]$original_md5))
    } else if (!file.copy(paths[i], backup)) stop("Cannot preserve original checkpoint.")
    contents <- saved[!names(saved) %in% c("spec", "source_equivalence_audit")]
    original_spec <- if (is.null(saved$source_equivalence_audit)) saved$spec else saved$source_equivalence_audit$original_spec
    saved$source_equivalence_audit <- list(original_spec = original_spec,
      original_md5 = rows[[i]]$original_md5, backup = backup,
      source_refs = rows[[i]]$original_source_refs,
      audit_script_md5 = unname(tools::md5sum("analysis/44_verify_benchmark_cache_equivalence.R")),
      verified_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      caveat = "Verified cache reuse; not a newly refitted benchmark model.")
    saved$spec <- current_spec
    write_rds_atomic(saved, paths[i])
    check <- readRDS(paths[i])
    stopifnot(identical(check[names(contents)], contents))
    rows[[i]]$status <- "migrated_verified_equivalent; original_preserved; model_contents_unchanged"
    rows[[i]]$backup <- backup
    rows[[i]]$backup_md5 <- unname(tools::md5sum(backup))
    rows[[i]]$migrated_md5 <- unname(tools::md5sum(paths[i]))
  }
}
audit <- bind_rows(rows)
write_csv_atomic(audit, file.path(DIRS$tables, "review_benchmark_cache_audit.csv"))
message("Benchmark cache audit: source equivalence = ", source_equivalent,
        "; eligible = ", sum(eligible), "/", length(paths),
        "; audit-only = ", audit_only,
        ". Unverified caches retain their old signature and require refitting.")
