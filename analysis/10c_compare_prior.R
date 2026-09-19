source("analysis/00_config.R")
source("analysis/functions/io.R")
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

prior_commit <- "d5a59464006c3a174c0de5d85d82a77eb549189d"
prior_path <- "results/tables/high_confidence_candidate_genes.csv"
prior_text <- system2("git", c("show", paste0(prior_commit, ":", prior_path)),
                      stdout = TRUE, stderr = TRUE)
if (!is.null(attr(prior_text, "status")) || length(prior_text) < 2L) {
  stop("Prior committed candidate table is unavailable in Git history.")
}
prior <- read_csv(I(paste(prior_text, collapse = "\n")), show_col_types = FALSE) |>
  select(symbol, prior_gene_id = tcga_gene_id)
current <- read_csv(file.path(DIRS$tables, "high_confidence_candidate_genes.csv"),
                    show_col_types = FALSE) |>
  select(symbol, current_gene_id = tcga_gene_id)
delta <- full_join(prior, current, by = "symbol") |>
  mutate(status = case_when(is.na(prior_gene_id) ~ "added",
                            is.na(current_gene_id) ~ "removed",
                            TRUE ~ "retained"),
         prior_commit = prior_commit) |>
  arrange(status, symbol)
stopifnot(sum(delta$status == "added") + nrow(prior) -
            sum(delta$status == "removed") == nrow(current))
write_csv_atomic(delta, file.path(DIRS$tables, "prior_candidate_delta.csv"))
message("Relative to the original committed run: ", sum(delta$status == "added"),
        " added, ", sum(delta$status == "removed"), " removed, ",
        sum(delta$status == "retained"), " retained.")
