select_tcga_patient_samples <- function(coldata, counts) {
  required <- c("sample_barcode", "sample_type")
  if (!all(required %in% names(coldata))) stop("TCGA sample metadata is incomplete.")
  if (anyDuplicated(coldata$sample_barcode)) stop("Duplicate TCGA sample barcodes.")
  if (!all(coldata$sample_barcode %in% colnames(counts))) stop("TCGA counts are missing samples.")

  coldata$patient_barcode <- substr(coldata$sample_barcode, 1L, 12L)
  coldata$library_depth <- as.numeric(colSums(counts[, coldata$sample_barcode, drop = FALSE]))
  selected <- coldata |>
    dplyr::arrange(.data$patient_barcode, .data$sample_type,
                   dplyr::desc(.data$library_depth), .data$sample_barcode) |>
    dplyr::group_by(.data$patient_barcode, .data$sample_type) |>
    dplyr::slice_head(n = 1L) |>
    dplyr::ungroup()
  if (anyDuplicated(paste(selected$patient_barcode, selected$sample_type))) {
    stop("TCGA patient/type selection is not unique.")
  }
  selected
}

read_selected_tcga_coldata <- function(coldata_path, counts_path) {
  coldata <- readr::read_csv(coldata_path, show_col_types = FALSE)
  counts <- readRDS(counts_path)
  select_tcga_patient_samples(coldata, counts)
}
