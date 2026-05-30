source("analysis/00_config.R")
source("analysis/functions/io.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

repro <- read_csv(file.path(DIRS$tables, "reproducible_deg_tcga_gse40435_gse53757.csv"), show_col_types = FALSE)
surv <- read_csv(FILES$tcga_survival, show_col_types = FALSE)

candidates <- repro |>
  left_join(surv, by = c("tcga_gene_id" = "gene_id")) |>
  mutate(
    prognostic = !is.na(fdr) & fdr < THRESHOLDS$survival_fdr,
    candidate = reproducible_deg & prognostic,
    abs_tcga_log2fc = abs(tcga_log2fc),
    abs_log_hr = abs(log_hr)
  ) |>
  arrange(desc(candidate), fdr, desc(abs_tcga_log2fc))

write_csv_atomic(candidates, file.path(DIRS$tables, "candidate_gene_evidence_table.csv"))

summary <- tibble(
  metric = c("reproducible_deg", "prognostic_reproducible_deg"),
  value = c(
    sum(candidates$reproducible_deg, na.rm = TRUE),
    sum(candidates$candidate, na.rm = TRUE)
  )
)
write_csv_atomic(summary, file.path(DIRS$tables, "candidate_summary.csv"))

message("Candidate reproducible + prognostic genes: ", sum(candidates$candidate, na.rm = TRUE))

