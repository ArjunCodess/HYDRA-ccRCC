source("analysis/00_config.R")
source("analysis/functions/io.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

source_rows <- tribble(
  ~source_id, ~role, ~url,
  "TCGA-KIRC", "discovery expression and survival", SOURCE_URLS$tcga_gdc,
  "GSE40435", "tumor-normal expression validation", SOURCE_URLS$geo_gse40435,
  "GSE53757", "tumor-normal expression validation", SOURCE_URLS$geo_gse53757,
  "GSE29609", "small external survival direction check", SOURCE_URLS$geo_gse29609,
  "E-MTAB-1980", "independent external survival validation", SOURCE_URLS$emtab1980,
  "HPA-v25.1", "independent single-cell type expression", SOURCE_URLS$hpa_single_cell,
  "Aran-2015-CPE", "published consensus TCGA tumor-purity sensitivity", SOURCE_URLS$aran2015_purity_study
) |>
  mutate(
    access_date = format(Sys.Date(), "%Y-%m-%d"),
    retrieval = "scripted public download",
    candidate_definition_role = if_else(
      source_id %in% c("E-MTAB-1980", "HPA-v25.1", "Aran-2015-CPE"),
      "downstream evaluation only; not used to define the reviewer-driven revised candidate set",
      "part of discovery or upstream replication"
    )
  )

write_csv_atomic(source_rows, file.path(DIRS$tables, "source_provenance.csv"))

roots <- c("analysis", "results/tables", "results/figures")
files <- unlist(lapply(roots, function(root) {
  list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
}))
files <- c(files, "README.md", "protocol.md", "run_pipeline.ps1", "environment/sessionInfo.txt")
files <- files[file.exists(files) & !dir.exists(files)]
files <- files[!grepl("run_manifest\\.csv$", files)]
info <- file.info(files)

manifest <- tibble(
  path = gsub("\\\\", "/", files),
  bytes = as.numeric(info$size),
  modified_utc = format(as.POSIXct(info$mtime, tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  md5 = unname(tools::md5sum(files))
) |>
  arrange(path)

write_csv_atomic(manifest, file.path(DIRS$tables, "run_manifest.csv"))
writeLines(
  sub("[[:space:]]+$", "", capture.output(sessionInfo())),
  "environment/sessionInfo.txt"
)

message("Run manifest complete. Files checksummed: ", nrow(manifest))
