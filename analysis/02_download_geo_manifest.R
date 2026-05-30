source("analysis/00_config.R")

suppressPackageStartupMessages({
  library(tibble)
  library(readr)
})

geo_manifest <- tribble(
  ~accession, ~role, ~priority, ~download_now, ~notes,
  "GSE40435", "main_validation", 1L, TRUE, "101 paired ccRCC tumor/adjacent non-tumor samples; download Series Matrix and non-normalized supplementary table.",
  "GSE53757", "main_validation", 2L, TRUE, "144 matched tumor/normal samples on Affymetrix GPL570; download Series Matrix and RAW CEL archive."
)

write_csv(geo_manifest, file.path(DIRS$metadata, "geo_manifest.csv"))
message("Saved GEO manifest to data/metadata/geo_manifest.csv")
