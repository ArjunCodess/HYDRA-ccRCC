source("analysis/00_config.R")

suppressPackageStartupMessages({
  library(tibble)
  library(readr)
})

geo_manifest <- tribble(
  ~accession, ~role, ~priority, ~download_now, ~notes,
  "GSE40435", "main_validation", 1L, TRUE, "101 paired ccRCC tumor/adjacent non-tumor samples; download Series Matrix and non-normalized supplementary table.",
  "GSE53757", "main_validation", 2L, TRUE, "144 matched tumor/normal samples on Affymetrix GPL570; download Series Matrix and RAW CEL archive.",
  "GSE46699", "secondary_validation", 3L, FALSE, "Paired ccRCC tumor/normal; exclude documented technical duplicate CEL files after normalization.",
  "GSE36895", "secondary_validation", 4L, FALSE, "Use human primary tumors and normal kidney cortex; exclude tumorgraft/xenograft samples from main validation.",
  "GSE68417", "grade_biology", 5L, FALSE, "Use for grade-biology support, not main tumor-normal validation."
)

write_csv(geo_manifest, file.path(DIRS$metadata, "geo_manifest.csv"))
message("Saved GEO manifest to data/metadata/geo_manifest.csv")
