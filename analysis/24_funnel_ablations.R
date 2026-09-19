source("analysis/00_config.R")
source("analysis/functions/io.R")
source("analysis/functions/tcga_metadata.R")
source("analysis/functions/patient_samples.R")

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(Biobase)
  library(readxl)
  library(dplyr)
  library(readr)
  library(stringr)
  library(survival)
})

all_gene <- read_csv(FILES$tcga_apeglm_survival, show_col_types = FALSE) |>
  filter(!is.na(symbol), symbol != "", model_status == "ok") |>
  arrange(global_fdr, gene_id) |>
  distinct(symbol, .keep_all = TRUE)
full <- read_csv(file.path(DIRS$tables, "high_confidence_candidate_genes.csv"),
                 show_col_types = FALSE)
k <- nrow(full)
stopifnot(k > 0L, !anyDuplicated(full$symbol))

list_rows <- list(
  complete_rule = full |> transmute(symbol, gene_id = tcga_gene_id,
                                    tcga_log_hr = main_log_hr),
  survival_only = all_gene |> slice_head(n = k) |>
    transmute(symbol, gene_id, tcga_log_hr = log_hr),
  de_only = all_gene |> arrange(tcga_de_fdr, gene_id) |> slice_head(n = k) |>
    transmute(symbol, gene_id, tcga_log_hr = log_hr)
)

# Refit the omitted-covariate checks for genes admitted by each single-GEO rule.
vst <- read_required_rds(FILES$tcga_vst)
meta <- read_selected_tcga_coldata(FILES$tcga_coldata, FILES$tcga_counts)
clinical <- read_csv(FILES$tcga_clinical, show_col_types = FALSE) |>
  transmute(patient_barcode = submitter_id, os_time, os_event,
            age = suppressWarnings(as.numeric(age_at_diagnosis)) / 365.25,
            sex = factor(gender), stage = factor(normalize_stage(ajcc_pathologic_stage)),
            grade = factor(case_when(normalize_grade(tumor_grade) %in% c("G1", "G2") ~ "Low grade",
                                     normalize_grade(tumor_grade) %in% c("G3", "G4") ~ "High grade",
                                     TRUE ~ NA_character_)))
surv_dat <- meta |> filter(sample_type == "Primary Tumor") |>
  mutate(patient_barcode = substr(sample_barcode, 1, 12)) |>
  inner_join(clinical, by = "patient_barcode") |>
  filter(is.finite(os_time), os_time > 0, os_event %in% c(0L, 1L),
         is.finite(age), !is.na(sex), !is.na(stage), !is.na(grade))
stopifnot(!anyDuplicated(surv_dat$patient_barcode))
passes_sensitivity <- function(gene_id, direction) {
  if (!gene_id %in% rownames(vst)) return(FALSE)
  dat <- surv_dat |> mutate(expr = as.numeric(scale(vst[gene_id, sample_barcode])))
  fits <- lapply(list(c("age", "sex", "stage"), c("age", "sex", "grade")),
    function(vars) tryCatch(coxph(as.formula(paste("Surv(os_time, os_event) ~ expr +",
                                                paste(vars, collapse = " + "))), data = dat),
                            error = function(e) NULL))
  if (any(vapply(fits, is.null, logical(1)))) return(FALSE)
  all(vapply(fits, function(f) {
    z <- summary(f)$coefficients["expr", ]
    is.finite(z["Pr(>|z|)"]) && z["Pr(>|z|)"] < 0.05 &&
      sign(z["coef"]) == sign(direction)
  }, logical(1)))
}
for (kept_geo in c("gse40435", "gse53757")) {
  effect <- paste0(kept_geo, "_log2fc")
  p <- paste0(kept_geo, "_p_value")
  eligible <- all_gene |>
    filter(primary_tcga_deg_gate, is.finite(.data[[effect]]),
           is.finite(.data[[p]]), .data[[p]] < 0.05,
           sign(tcga_mle_log2fc) == sign(.data[[effect]]),
           abs(.data[[effect]]) >= THRESHOLDS$min_geo_abs_log2fc) |>
    mutate(relevant_fdr = p.adjust(p_value, method = "BH")) |>
    filter(relevant_fdr < THRESHOLDS$high_confidence_survival_fdr,
           abs(log_hr) >= THRESHOLDS$high_confidence_abs_log_hr) |>
    arrange(relevant_fdr, desc(abs(log_hr)), gene_id)
  keep <- logical(nrow(eligible))
  for (i in seq_len(nrow(eligible))) {
    keep[i] <- passes_sensitivity(eligible$gene_id[i], eligible$log_hr[i])
    if (sum(keep) >= k) break
  }
  name <- paste0("leave_out_", if (kept_geo == "gse40435") "gse53757" else "gse40435")
  list_rows[[name]] <- eligible[keep, ] |> slice_head(n = k) |>
    transmute(symbol, gene_id, tcga_log_hr = log_hr)
}

# Controls match mean TCGA abundance but fail the genome-wide survival screen.
pool <- all_gene |>
  filter(global_fdr > 0.5, is.finite(base_mean), base_mean > 0,
         !symbol %in% unlist(lapply(list_rows, function(x) x$symbol)))
targets <- full |> left_join(all_gene |> select(gene_id, base_mean),
                            by = c("tcga_gene_id" = "gene_id"))
controls <- vector("list", k)
for (i in seq_len(k)) {
  idx <- which.min(abs(log1p(pool$base_mean) - log1p(targets$base_mean[i])))
  controls[[i]] <- pool[idx, ] |> transmute(symbol, gene_id, tcga_log_hr = log_hr)
  pool <- pool[-idx, ]
}
list_rows$expression_matched_control <- bind_rows(controls)
k_matched <- min(vapply(list_rows, nrow, integer(1)))
if (k_matched < k) message("Matched-list comparison uses ", k_matched,
                           " genes because one leave-one-GEO rule admitted fewer than ", k, ".")
list_rows <- lapply(list_rows, function(x) slice_head(x, n = k_matched))
lists <- bind_rows(lapply(names(list_rows), function(name) {
  list_rows[[name]] |> mutate(rule = name, rank = row_number(), .before = 1)
}))
write_csv_atomic(lists, file.path(DIRS$tables, "funnel_matched_lists.csv"))

fit_external <- function(expr, clinical, cohort) {
  bind_rows(lapply(seq_len(nrow(lists)), function(i) {
    row <- lists[i, ]
    if (!row$symbol %in% rownames(expr)) {
      return(tibble(rule = row$rule, rank = row$rank, symbol = row$symbol,
                    cohort, present = FALSE, log_hr = NA_real_, p_value = NA_real_,
                    same_direction = NA))
    }
    dat <- clinical |> mutate(expression = as.numeric(scale(expr[row$symbol, sample])))
    fit <- tryCatch(coxph(Surv(time, event) ~ expression, data = dat),
                    error = function(e) NULL)
    beta <- if (is.null(fit)) NA_real_ else unname(coef(fit)["expression"])
    p <- if (is.null(fit)) NA_real_ else summary(fit)$coefficients["expression", "Pr(>|z|)"]
    tibble(rule = row$rule, rank = row$rank, symbol = row$symbol,
           cohort, present = TRUE, log_hr = beta, p_value = p,
           same_direction = is.finite(beta) && sign(beta) == sign(row$tcga_log_hr))
  }))
}

gse <- read_required_rds(file.path(DIRS$processed, "gse29609_series_matrix.rds"))
gse_feature <- fData(gse)
symbol_col <- names(gse_feature)[str_detect(tolower(names(gse_feature)), "gene symbol")][1]
gse_symbol <- str_trim(str_split_fixed(as.character(gse_feature[[symbol_col]]),
                                       " /// | // |;|,", 2)[, 1])
gse_expr <- as.data.frame(exprs(gse)) |> mutate(symbol = gse_symbol) |>
  filter(!is.na(symbol), symbol != "") |>
  group_by(symbol) |> summarise(across(where(is.numeric), mean), .groups = "drop") |>
  tibble::column_to_rownames("symbol") |> as.matrix()
gse_pheno <- pData(gse)
gse_clinical <- tibble(sample = rownames(gse_pheno),
  time = as.numeric(gse_pheno[["survival time:ch1"]]),
  event = as.numeric(gse_pheno[["death (1=yes, 0=no):ch1"]]))

em_expr_raw <- read_tsv(file.path(DIRS$raw, "arrayexpress", "E-MTAB-1980",
                                  "ccRCC_exp_log_quantile_normalized.txt"),
                        show_col_types = FALSE, progress = FALSE, name_repair = "minimal")
em_clinical_raw <- readxl::read_excel(file.path(DIRS$raw, "arrayexpress", "E-MTAB-1980",
                                                "supplementary_table_1.xlsx"), sheet = 1, skip = 1)
em_clinical <- em_clinical_raw |>
  transmute(sample = as.character(.data[["sample ID"]]),
            time = as.numeric(.data[["observation period (month)"]]),
            event = as.integer(tolower(as.character(.data[["outcome"]])) == "dead"),
            subtype = as.character(.data[["gene expression profile"]])) |>
  filter(subtype %in% c("ccA", "ccB"), is.finite(time), time > 0)
refseq <- sub("\\.\\d+$", "", as.character(em_expr_raw$SystematicName))
ref_map <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(refseq[!is.na(refseq) & refseq != ""]),
                                  keytype = "REFSEQ", columns = "SYMBOL") |>
  filter(!is.na(SYMBOL)) |> distinct(REFSEQ, SYMBOL)
em_expr <- em_expr_raw |> mutate(refseq = refseq) |>
  inner_join(ref_map, by = c("refseq" = "REFSEQ")) |>
  select(SYMBOL, all_of(em_clinical$sample)) |>
  mutate(across(all_of(em_clinical$sample), as.numeric)) |>
  group_by(SYMBOL) |>
  summarise(across(all_of(em_clinical$sample), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
  tibble::column_to_rownames("SYMBOL") |> as.matrix()

external <- bind_rows(fit_external(gse_expr, gse_clinical, "GSE29609"),
                      fit_external(em_expr, em_clinical, "E-MTAB-1980")) |>
  group_by(rule, cohort) |> mutate(fdr = p.adjust(p_value, method = "BH")) |> ungroup()
write_csv_atomic(external, file.path(DIRS$tables, "funnel_external_gene_results.csv"))
summary <- external |> group_by(rule, cohort) |>
  summarise(list_size = n(), present = sum(present),
            directional_rate = mean(same_direction, na.rm = TRUE),
            same_direction_fdr = sum(same_direction & fdr < 0.05, na.rm = TRUE),
            opposite_direction_fdr = sum(!same_direction & fdr < 0.05, na.rm = TRUE),
            same_direction = sum(same_direction, na.rm = TRUE),
            .groups = "drop")
write_csv_atomic(summary, file.path(DIRS$tables, "funnel_external_summary.csv"))

# Resample the union of genes, preserving pairing between the complete rule
# and survival-only list in each exploratory cohort.
set.seed(RESAMPLING$seed + 25L)
paired <- bind_rows(lapply(c("GSE29609", "E-MTAB-1980"), function(cohort) {
  d <- external |> filter(.data$cohort == .env$cohort, rule %in% c("complete_rule", "survival_only"))
  wide <- d |> select(rule, symbol, same_direction) |>
    tidyr::pivot_wider(names_from = rule, values_from = same_direction)
  genes <- wide$symbol
  draws <- replicate(2000L, {
    x <- wide[match(sample(genes, length(genes), replace = TRUE), genes), ]
    mean(x$complete_rule, na.rm = TRUE) - mean(x$survival_only, na.rm = TRUE)
  })
  tibble(cohort, difference = mean(wide$complete_rule, na.rm = TRUE) -
           mean(wide$survival_only, na.rm = TRUE),
         ci_low = quantile(draws, 0.025, na.rm = TRUE),
         ci_high = quantile(draws, 0.975, na.rm = TRUE),
         acceptance = ci_low > 0)
}))
write_csv_atomic(paired, file.path(DIRS$tables, "funnel_external_paired_bootstrap.csv"))
message("Matched-list funnel ablations complete: ", k_matched, " genes per rule.")
