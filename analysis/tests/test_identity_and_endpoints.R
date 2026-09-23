source("analysis/00_config.R")
source("analysis/functions/patient_samples.R")
source("analysis/functions/tcga_metadata.R")
suppressPackageStartupMessages(library(dplyr))

meta <- tibble(
  sample_barcode = c("TCGA-AA-0001-01A", "TCGA-AA-0001-01B",
                     "TCGA-AA-0002-01B", "TCGA-AA-0002-01A",
                     "TCGA-AA-0001-11A", "TCGA-AA-0001-11B"),
  sample_type = c(rep("Primary Tumor", 4), rep("Solid Tissue Normal", 2))
)
counts <- matrix(c(2, 2, 4, 4, 1, 9,
                   2, 3, 4, 4, 1, 1), nrow = 2, byrow = TRUE,
                 dimnames = list(c("gene1", "gene2"), meta$sample_barcode))
selected <- select_tcga_patient_samples(meta, counts)
stopifnot(setequal(selected$sample_barcode,
                   c("TCGA-AA-0001-01B", "TCGA-AA-0002-01A", "TCGA-AA-0001-11B")),
          !anyDuplicated(paste(substr(selected$sample_barcode, 1, 12), selected$sample_type)))

clinical <- tibble(vital_status = c("Dead", "Alive", NA_character_),
                   days_to_death = c(10, NA, NA),
                   days_to_last_follow_up = c(20, 30, 40))
surv <- derive_os(clinical)
stopifnot(identical(surv$os_event, c(1L, 0L, NA_integer_)),
          identical(surv$os_time, c(10, 30, NA_real_)))
bad <- clinical
bad$vital_status[1] <- "unknown"
stopifnot(inherits(try(derive_os(bad), silent = TRUE), "try-error"))
message("Synthetic patient identity and endpoint tests passed.")
