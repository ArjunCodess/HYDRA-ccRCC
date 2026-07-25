source("analysis/00_config.R")
source("analysis/functions/tcga_metadata.R")
source("analysis/functions/io.R")

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
  library(readr)
})

force_download <- identical(Sys.getenv("HYDRA_FORCE_DOWNLOAD"), "1")
if (!force_download && all(file.exists(c(FILES$tcga_se, FILES$tcga_counts, FILES$tcga_coldata, FILES$tcga_clinical)))) {
  clinical <- read_csv(FILES$tcga_clinical, show_col_types = FALSE)
  if (!"gender" %in% names(clinical) && "sex_at_birth" %in% names(clinical)) {
    clinical <- clinical |> mutate(gender = sex_at_birth)
    write_csv(clinical, FILES$tcga_clinical)
    message("Normalized current GDC sex_at_birth field to legacy gender contract.")
  }
  message("TCGA-KIRC processed files already exist. Set HYDRA_FORCE_DOWNLOAD=1 to redownload.")
  quit(save = "no", status = 0)
}

query <- GDCquery(
  project = PROJECT_ID,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = TCGA_SAMPLE_TYPES
)

GDCdownload(query, directory = DIRS$raw, method = "api")
se <- GDCprepare(query, directory = DIRS$raw)

write_rds_atomic(se, FILES$tcga_se)

counts <- SummarizedExperiment::assay(se, "unstranded")
write_rds_atomic(counts, FILES$tcga_counts)

coldata <- as.data.frame(SummarizedExperiment::colData(se)) |>
  tibble::rownames_to_column("sample_barcode") |>
  mutate(patient_barcode = tcga_patient_barcode(sample_barcode))

write_csv(coldata, FILES$tcga_coldata)

clinical <- GDCquery_clinic(project = PROJECT_ID, type = "clinical") |>
  derive_os() |>
  add_clean_stage_grade() |>
  mutate(submitter_id = as.character(.data$submitter_id))
if (!"gender" %in% names(clinical) && "sex_at_birth" %in% names(clinical)) {
  clinical <- clinical |> mutate(gender = sex_at_birth)
}

write_csv(clinical, FILES$tcga_clinical)
writeLines(capture.output(sessionInfo()), FILES$tcga_session)

message("Saved TCGA SummarizedExperiment: ", FILES$tcga_se)
message("Saved TCGA counts: ", FILES$tcga_counts)
message("Saved TCGA sample metadata: ", FILES$tcga_coldata)
message("Saved TCGA clinical metadata: ", FILES$tcga_clinical)
