source("analysis/00_config.R")
source("analysis/functions/io.R")

suppressPackageStartupMessages({
  library(GEOquery)
  library(readr)
  library(dplyr)
})

manifest_path <- file.path(DIRS$metadata, "geo_manifest.csv")
if (!file.exists(manifest_path)) {
  source("analysis/02_download_geo_manifest.R")
}

geo_manifest <- read_csv(manifest_path, show_col_types = FALSE)

geo_dir <- file.path(DIRS$raw, "geo")
dir.create(geo_dir, recursive = TRUE, showWarnings = FALSE)

download_one_geo <- function(accession) {
  series_path <- file.path(DIRS$processed, paste0(tolower(accession), "_series_matrix.rds"))
  supp_dir <- file.path(geo_dir, accession)
  force_download <- identical(Sys.getenv("HYDRA_FORCE_DOWNLOAD"), "1")

  if (!force_download && file.exists(series_path) && dir.exists(supp_dir) && length(list.files(supp_dir)) > 0) {
    message(accession, " already downloaded. Set HYDRA_FORCE_DOWNLOAD=1 to redownload.")
    return(tibble(
      accession = accession,
      series_rds = series_path,
      supplementary_dir = supp_dir,
      downloaded_at = as.character(Sys.time()),
      status = "cached"
    ))
  }

  message("Downloading GEO Series Matrix for ", accession)
  series <- getGEO(accession, GSEMatrix = TRUE, AnnotGPL = FALSE, destdir = geo_dir)
  write_rds_atomic(series, series_path)

  message("Downloading GEO supplementary files for ", accession)
  dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
  getGEOSuppFiles(accession, baseDir = supp_dir, makeDirectory = FALSE)

  tibble(
    accession = accession,
    series_rds = series_path,
    supplementary_dir = supp_dir,
    downloaded_at = as.character(Sys.time()),
    status = "downloaded"
  )
}

to_download <- geo_manifest |>
  filter(download_now) |>
  arrange(priority) |>
  pull(accession)

download_log <- bind_rows(lapply(to_download, download_one_geo))
write_csv(download_log, file.path(DIRS$metadata, "geo_download_log.csv"))

message("GEO download complete for: ", paste(to_download, collapse = ", "))
