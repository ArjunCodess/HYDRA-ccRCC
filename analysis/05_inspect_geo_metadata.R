source("analysis/00_config.R")

suppressPackageStartupMessages({
  library(Biobase)
  library(readr)
  library(dplyr)
})

summarise_geo <- function(accession) {
  series <- readRDS(file.path(DIRS$processed, paste0(tolower(accession), "_series_matrix.rds")))
  eset <- series[[1]]
  pdata <- Biobase::pData(eset)

  metadata_path <- file.path(DIRS$metadata, paste0(tolower(accession), "_pheno_columns.csv"))
  columns <- tibble(
    column = names(pdata),
    example_1 = vapply(pdata, function(x) as.character(x[1]), character(1)),
    example_2 = vapply(pdata, function(x) as.character(x[min(2, length(x))]), character(1)),
    unique_n = vapply(pdata, function(x) length(unique(as.character(x))), integer(1))
  )
  write_csv(columns, metadata_path)

  message(accession, ": ", nrow(pdata), " samples, ", ncol(Biobase::exprs(eset)), " expression columns")
  message("Wrote metadata column audit: ", metadata_path)
}

summarise_geo("GSE40435")
summarise_geo("GSE53757")

