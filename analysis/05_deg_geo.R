source("analysis/00_config.R")
source("analysis/functions/io.R")

suppressPackageStartupMessages({
  library(Biobase)
  library(limma)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

collapse_to_gene <- function(expr, feature_data, accession) {
  symbol_col <- names(feature_data)[str_detect(tolower(names(feature_data)), "symbol|gene.symbol|gene symbol")]
  if (length(symbol_col) == 0) {
    stop("No gene symbol-like feature column found for ", accession, call. = FALSE)
  }

  symbols <- as.character(feature_data[[symbol_col[1]]])
  symbols <- str_split_fixed(symbols, " /// | // |;|,", 2)[, 1]
  symbols <- str_trim(symbols)
  symbols[symbols == "" | is.na(symbols)] <- NA_character_

  df <- as.data.frame(expr) |>
    mutate(symbol = symbols) |>
    filter(!is.na(symbol))

  df |>
    group_by(symbol) |>
    summarise(across(where(is.numeric), mean), .groups = "drop") |>
    column_to_rownames("symbol") |>
    as.matrix()
}

run_geo_limma <- function(accession) {
  series <- readRDS(file.path(DIRS$processed, paste0(tolower(accession), "_series_matrix.rds")))
  eset <- series[[1]]
  expr <- Biobase::exprs(eset)
  pdata <- Biobase::pData(eset)
  fdata <- Biobase::fData(eset)

  if (max(expr, na.rm = TRUE) > 100) {
    expr <- log2(expr + 1)
  }

  if (accession == "GSE40435") {
    condition <- case_when(
      str_detect(tolower(pdata$source_name_ch1), "adjacent|normal|non[- ]?tumou?r") ~ "normal",
      str_detect(tolower(pdata$source_name_ch1), "tumou?r|carcinoma|ccrcc") ~ "tumor",
      TRUE ~ NA_character_
    )
    patient <- str_match(pdata$title, "patient\\s+([0-9]+)")[, 2]
  } else if (accession == "GSE53757") {
    condition <- case_when(
      str_detect(tolower(pdata$characteristics_ch1), "clear cell renal cell carcinoma") ~ "tumor",
      str_detect(tolower(pdata$characteristics_ch1), "normal kidney") ~ "normal",
      TRUE ~ NA_character_
    )
    patient <- as.character(rep(seq_len(nrow(pdata) / 2), each = 2, length.out = nrow(pdata)))
  } else {
    stop("Unsupported accession: ", accession, call. = FALSE)
  }

  keep_samples <- !is.na(condition) & !is.na(patient)
  expr <- expr[, keep_samples]
  condition <- factor(condition[keep_samples], levels = c("normal", "tumor"))
  patient <- factor(patient[keep_samples])

  gene_expr <- collapse_to_gene(expr, fdata, accession)

  if (nlevels(condition) != 2) {
    stop(accession, " condition parsing failed. Parsed levels: ", paste(levels(droplevels(condition)), collapse = ", "), call. = FALSE)
  }

  design <- model.matrix(~ patient + condition)
  fit <- lmFit(gene_expr, design)
  fit <- eBayes(fit)

  result <- topTable(fit, coef = "conditiontumor", number = Inf, sort.by = "P") |>
    rownames_to_column("symbol") |>
    rename(
      log2FoldChange = logFC,
      pvalue = P.Value,
      padj = adj.P.Val
    ) |>
    mutate(
      accession = accession,
      significant = padj < THRESHOLDS$deg_fdr & abs(log2FoldChange) >= THRESHOLDS$deg_abs_log2fc
    ) |>
    select(accession, symbol, log2FoldChange, AveExpr, t, pvalue, padj, B, significant)

  out_path <- file.path(DIRS$tables, paste0(tolower(accession), "_limma_tumor_vs_normal.csv"))
  write_csv_atomic(result, out_path)

  sample_summary <- tibble(accession = accession, condition = condition, patient = patient) |>
    count(accession, condition, name = "n_samples")
  write_csv_atomic(sample_summary, file.path(DIRS$tables, paste0(tolower(accession), "_sample_summary.csv")))

  message(accession, " limma complete. Significant DEG count: ", sum(result$significant, na.rm = TRUE))
}

args <- commandArgs(trailingOnly = TRUE)
accessions <- if (length(args) > 0) toupper(args) else c("GSE40435", "GSE53757")
for (accession in accessions) {
  run_geo_limma(accession)
  gc()
}
