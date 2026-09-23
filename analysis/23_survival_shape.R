source("analysis/00_config.R")
source("analysis/functions/io.R")
source("analysis/functions/tcga_metadata.R")
source("analysis/functions/patient_samples.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(survival)
  library(splines)
})

vst <- read_required_rds(FILES$tcga_vst)
meta <- read_selected_tcga_coldata(FILES$tcga_coldata, FILES$tcga_counts)
clinical <- read_csv(FILES$tcga_clinical, show_col_types = FALSE) |>
  transmute(patient_barcode = submitter_id, os_time, os_event,
            age = suppressWarnings(as.numeric(age_at_diagnosis)) / 365.25,
            sex = factor(gender), stage = factor(normalize_stage(ajcc_pathologic_stage)),
            grade = factor(case_when(normalize_grade(tumor_grade) %in% c("G1", "G2") ~ "Low grade",
                                     normalize_grade(tumor_grade) %in% c("G3", "G4") ~ "High grade",
                                     TRUE ~ NA_character_)))
dat <- meta |>
  filter(sample_type == "Primary Tumor") |>
  mutate(patient_barcode = substr(sample_barcode, 1, 12)) |>
  inner_join(clinical, by = "patient_barcode") |>
  filter(is.finite(os_time), os_time > 0, os_event %in% c(0L, 1L), is.finite(age),
         !is.na(sex), !is.na(stage), !is.na(grade))
stopifnot(!anyDuplicated(dat$patient_barcode))
candidates <- read_csv(file.path(DIRS$tables, "high_confidence_candidate_genes.csv"),
                       show_col_types = FALSE)

fit_one <- function(i) {
  gene_id <- candidates$tcga_gene_id[i]
  x <- dat |> mutate(expr = as.numeric(scale(vst[gene_id, sample_barcode])))
  linear <- coxph(Surv(os_time, os_event) ~ expr + age + sex + stage + grade, data = x)
  nonlinear <- coxph(Surv(os_time, os_event) ~ ns(expr, df = 3) + age + sex + stage + grade,
                     data = x)
  varying <- tryCatch(coxph(Surv(os_time, os_event) ~ expr + age + sex + stage + grade + tt(expr),
                            data = x, tt = function(z, t, ...) z * log(pmax(t, 1))),
                      error = function(e) NULL)
  nonlinear_lrt <- anova(linear, nonlinear, test = "LRT")
  varying_lrt <- if (is.null(varying)) NA_real_ else
    pchisq(2 * (as.numeric(logLik(varying)) - as.numeric(logLik(linear))),
           df = 1, lower.tail = FALSE)
  tibble(symbol = candidates$symbol[i], n = nrow(x), events = sum(x$os_event),
         linear_log_hr = unname(coef(linear)["expr"]),
         nonlinear_lrt_p = nonlinear_lrt[2, "Pr(>|Chi|)"],
         time_varying_lrt_p = varying_lrt,
         time_varying_coefficient = if (is.null(varying)) NA_real_ else unname(coef(varying)["tt(expr)"]))
}

out <- bind_rows(lapply(seq_len(nrow(candidates)), fit_one)) |>
  mutate(nonlinear_fdr = p.adjust(nonlinear_lrt_p, "BH"),
         time_varying_fdr = p.adjust(time_varying_lrt_p, "BH"))
write_csv_atomic(out, file.path(DIRS$tables, "candidate_survival_shape_sensitivity.csv"))
message("Nonlinear and time-varying survival checks complete: ", nrow(out), " genes.")
