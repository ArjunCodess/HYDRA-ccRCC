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
candidates <- read_csv(file.path(DIRS$tables, "candidate_gene_evidence_table.csv"), show_col_types = FALSE)

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
  filter(!is.na(os_time), os_time > 0, !is.na(os_event), sample_barcode %in% colnames(vst_mat)) |>
  distinct(patient_barcode, .keep_all = TRUE)

gene_table <- candidates |>
  filter(reproducible_deg, tcga_gene_id %in% rownames(vst_mat)) |>
  select(symbol, tcga_gene_id, gse40435_log2fc, gse53757_log2fc, high_confidence_candidate)

fit_expr <- function(gene_id, sampled_barcodes, covars, ph = FALSE) {
  dat <- sample_data |>
    filter(sample_barcode %in% sampled_barcodes) |>
    mutate(expr = as.numeric(scale(vst_mat[gene_id, sample_barcode]))) |>
    select(os_time, os_event, expr, all_of(covars)) |>
    filter(if_all(everything(), ~ !is.na(.x))) |>
    mutate(across(any_of(c("sex", "stage", "grade")), droplevels))

  factors <- intersect(c("sex", "stage", "grade"), covars)
  if (nrow(dat) < 100 || sum(dat$os_event) < 25 ||
      any(vapply(dat[factors], nlevels, integer(1)) < 2)) return(NULL)

  form <- as.formula(paste("Surv(os_time, os_event) ~ expr +", paste(covars, collapse = " + ")))
  fit <- tryCatch(coxph(form, data = dat), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  term <- tidy(fit) |> filter(term == "expr")
  if (nrow(term) != 1) return(NULL)

  ph_p <- NA_real_
  if (ph) {
    zph <- tryCatch(cox.zph(fit), error = function(e) NULL)
    if (!is.null(zph) && "expr" %in% rownames(zph$table)) ph_p <- zph$table["expr", "p"]
  }

  tibble(log_hr = term$estimate, p_value = term$p.value, ph_p_value = ph_p, n = fit$n, events = fit$nevent)
}

stratified_sample <- function(data, fraction) {
  data |>
    group_by(os_event) |>
    slice_sample(prop = fraction, replace = FALSE) |>
    ungroup() |>
    pull(sample_barcode)
}

set.seed(RESAMPLING$seed)
repeat_results <- vector("list", RESAMPLING$stability_repeats)
repeat_summary <- vector("list", RESAMPLING$stability_repeats)

for (repeat_id in seq_len(RESAMPLING$stability_repeats)) {
  sampled <- stratified_sample(sample_data, RESAMPLING$stability_fraction)
  message("Selection stability repeat ", repeat_id, "/", RESAMPLING$stability_repeats)

  main <- bind_rows(lapply(gene_table$tcga_gene_id, function(gene_id) {
    out <- fit_expr(gene_id, sampled, c("age", "sex", "stage", "grade"))
    if (is.null(out)) return(NULL)
    mutate(out, tcga_gene_id = gene_id)
  })) |>
    select(-ph_p_value) |>
    mutate(main_fdr = p.adjust(p_value, method = "BH"))

  eligible <- main |>
    filter(main_fdr < THRESHOLDS$high_confidence_survival_fdr,
           abs(log_hr) >= THRESHOLDS$high_confidence_abs_log_hr) |>
    inner_join(gene_table, by = "tcga_gene_id") |>
    filter(abs(gse40435_log2fc) >= THRESHOLDS$min_geo_abs_log2fc,
           abs(gse53757_log2fc) >= THRESHOLDS$min_geo_abs_log2fc)

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
      main_ph <- fit_expr(row$tcga_gene_id, sampled, c("age", "sex", "stage", "grade"), ph = TRUE)
      stage_fit <- fit_expr(row$tcga_gene_id, sampled, c("age", "sex", "stage"))
      grade_fit <- fit_expr(row$tcga_gene_id, sampled, c("age", "sex", "grade"))
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
    select(repeat_id, symbol, tcga_gene_id, log_hr, p_value, main_fdr,
           main_screen_pass, geo_effect_support, ph_pass, sensitivity_pass, selected)
  repeat_summary[[repeat_id]] <- tibble(
    repeat_id,
    sampled_patients = length(sampled),
    modeled_genes = nrow(main),
    main_screen_genes = sum(one_repeat$main_screen_pass, na.rm = TRUE),
    selected_genes = sum(one_repeat$selected, na.rm = TRUE)
  )
}

all_repeats <- bind_rows(repeat_results)
stability <- all_repeats |>
  group_by(symbol, tcga_gene_id) |>
  summarise(
    repeats_modeled = n(),
    main_screen_frequency = mean(main_screen_pass, na.rm = TRUE),
    selection_frequency = mean(selected, na.rm = TRUE),
    median_log_hr = median(log_hr, na.rm = TRUE),
    log_hr_iqr = IQR(log_hr, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(candidates |> select(symbol, full_main_log_hr = main_log_hr, frozen_high_confidence = high_confidence_candidate),
            by = "symbol") |>
  mutate(full_fit_direction_agreement = sign(median_log_hr) == sign(full_main_log_hr)) |>
  arrange(desc(selection_frequency), desc(main_screen_frequency))

write_csv_atomic(stability, file.path(DIRS$tables, "selection_stability_all_genes.csv"))
write_csv_atomic(
  stability |> filter(frozen_high_confidence),
  file.path(DIRS$tables, "selection_stability_frozen_candidates.csv")
)
write_csv_atomic(bind_rows(repeat_summary), file.path(DIRS$tables, "selection_stability_repeats.csv"))

message("Selection stability complete. Frozen candidates assessed: ",
        sum(stability$frozen_high_confidence, na.rm = TRUE))
