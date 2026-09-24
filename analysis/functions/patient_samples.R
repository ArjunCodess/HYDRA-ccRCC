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

apply_aliquot_rule <- function(coldata, counts, rule = c("deepest", "first_barcode", "drop_multi", "sum_counts")) {
  rule <- match.arg(rule)
  prepared <- coldata
  prepared$patient_barcode <- substr(prepared$sample_barcode, 1L, 12L)
  prepared$library_depth <- as.numeric(colSums(counts[, prepared$sample_barcode, drop = FALSE]))
  multi_patients <- prepared |>
    dplyr::filter(.data$sample_type == "Primary Tumor") |>
    dplyr::count(.data$patient_barcode) |>
    dplyr::filter(.data$n > 1L) |>
    dplyr::pull(.data$patient_barcode)

  if (rule == "drop_multi") {
    prepared <- prepared |> dplyr::filter(!.data$patient_barcode %in% multi_patients)
    selected <- select_tcga_patient_samples(prepared, counts)
    return(list(coldata = selected, counts = counts[, selected$sample_barcode, drop = FALSE]))
  }

  if (rule == "first_barcode") {
    selected <- prepared |>
      dplyr::arrange(.data$patient_barcode, .data$sample_type, .data$sample_barcode) |>
      dplyr::group_by(.data$patient_barcode, .data$sample_type) |>
      dplyr::slice_head(n = 1L) |>
      dplyr::ungroup()
    return(list(coldata = selected, counts = counts[, selected$sample_barcode, drop = FALSE]))
  }

  if (rule == "sum_counts") {
    normals <- select_tcga_patient_samples(prepared |> dplyr::filter(.data$sample_type != "Primary Tumor"), counts)
    tumors <- prepared |> dplyr::filter(.data$sample_type == "Primary Tumor")
    pieces <- lapply(split(tumors, tumors$patient_barcode), function(rows) {
      chosen <- rows$sample_barcode[which.max(rows$library_depth)]
      summed <- rowSums(counts[, rows$sample_barcode, drop = FALSE])
      list(barcode = chosen, counts = summed, depth = sum(summed), n_aliquots = nrow(rows))
    })
    tumor_counts <- do.call(cbind, lapply(pieces, `[[`, "counts"))
    colnames(tumor_counts) <- vapply(pieces, `[[`, character(1), "barcode")
    tumor_meta <- tumors |>
      dplyr::distinct(.data$patient_barcode, .keep_all = TRUE) |>
      dplyr::mutate(sample_barcode = colnames(tumor_counts)[match(.data$patient_barcode, names(pieces))],
                    library_depth = vapply(pieces, function(x) x$depth, numeric(1))[match(.data$patient_barcode, names(pieces))])
    selected <- dplyr::bind_rows(tumor_meta, normals)
    combined <- cbind(tumor_counts, counts[, normals$sample_barcode, drop = FALSE])
    combined <- combined[, selected$sample_barcode, drop = FALSE]
    return(list(coldata = selected, counts = combined))
  }

  selected <- select_tcga_patient_samples(prepared, counts)
  list(coldata = selected, counts = counts[, selected$sample_barcode, drop = FALSE])
}

read_selected_tcga_coldata <- function(coldata_path, counts_path) {
  coldata <- readr::read_csv(coldata_path, show_col_types = FALSE)
  counts <- readRDS(counts_path)
  select_tcga_patient_samples(coldata, counts)
}
