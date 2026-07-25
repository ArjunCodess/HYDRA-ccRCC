source("analysis/00_config.R")
source("analysis/functions/io.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

raw_dir <- file.path(DIRS$raw, "hpa")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
zip_path <- file.path(raw_dir, "rna_single_cell_type.tsv.zip")
force_download <- identical(Sys.getenv("HYDRA_FORCE_DOWNLOAD"), "1")

if (force_download || !file.exists(zip_path)) {
  message("Downloading Human Protein Atlas single-cell type data.")
  download.file(SOURCE_URLS$hpa_single_cell, zip_path, mode = "wb", quiet = FALSE)
}

members <- unzip(zip_path, list = TRUE)
tsv_member <- members$Name[grepl("\\.tsv$", members$Name, ignore.case = TRUE)][1]
if (is.na(tsv_member)) stop("HPA archive contains no TSV file.")
unzip(zip_path, files = tsv_member, exdir = raw_dir, overwrite = TRUE)
hpa <- read_tsv(file.path(raw_dir, tsv_member), show_col_types = FALSE)

names(hpa) <- make.names(names(hpa))
required <- c("Gene.name", "Cell.type")
if (!all(required %in% names(hpa))) {
  stop("Unexpected HPA single-cell schema: ", paste(names(hpa), collapse = ", "))
}
value_col <- intersect(c("nCPM", "NX"), names(hpa))[1]
if (is.na(value_col)) stop("HPA single-cell table has no nCPM or NX expression column.")

candidates <- read_csv(file.path(DIRS$tables, "candidate_gene_evidence_table.csv"), show_col_types = FALSE) |>
  filter(high_confidence_candidate) |>
  select(symbol, main_log_hr)

compartment <- function(cell_type) {
  case_when(
    grepl("proximal|tubular|collecting duct|urothelial", cell_type, ignore.case = TRUE) ~ "renal_epithelial",
    grepl("endothelial|pericyte|smooth muscle", cell_type, ignore.case = TRUE) ~ "vascular",
    grepl("t-cell|b-cell|macrophage|monocyte|dendritic|neutrophil|mast", cell_type, ignore.case = TRUE) ~ "immune",
    grepl("fibroblast|stromal|mesenchymal", cell_type, ignore.case = TRUE) ~ "stromal",
    TRUE ~ "other"
  )
}

candidate_hpa <- hpa |>
  transmute(
    symbol = Gene.name,
    cell_type = Cell.type,
    expression = suppressWarnings(as.numeric(.data[[value_col]]))
  ) |>
  filter(symbol %in% candidates$symbol, is.finite(expression)) |>
  mutate(compartment = compartment(cell_type))

top_types <- candidate_hpa |>
  group_by(symbol) |>
  arrange(desc(expression), .by_group = TRUE) |>
  mutate(cell_type_rank = row_number()) |>
  filter(cell_type_rank <= 10) |>
  ungroup() |>
  left_join(candidates, by = "symbol")

compartment_summary <- candidate_hpa |>
  group_by(symbol, compartment) |>
  summarise(
    max_expression = max(expression, na.rm = TRUE),
    mean_expression = mean(expression, na.rm = TRUE),
    top_cell_type = cell_type[which.max(expression)],
    .groups = "drop"
  ) |>
  group_by(symbol) |>
  mutate(
    compartment_rank = min_rank(desc(max_expression)),
    dominant_compartment = compartment[which.max(max_expression)]
  ) |>
  ungroup() |>
  left_join(candidates, by = "symbol") |>
  arrange(symbol, compartment_rank)

write_csv_atomic(top_types, file.path(DIRS$tables, "hpa_candidate_top_cell_types.csv"))
write_csv_atomic(compartment_summary, file.path(DIRS$tables, "hpa_candidate_cell_source_summary.csv"))

message("HPA cell-source triangulation complete. Candidates mapped: ",
        n_distinct(candidate_hpa$symbol), "/", nrow(candidates))
