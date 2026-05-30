read_required_rds <- function(path) {
  if (!file.exists(path)) {
    stop("Required file does not exist: ", path, call. = FALSE)
  }
  readRDS(path)
}

write_rds_atomic <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
  on.exit(unlink(tmp), add = TRUE)

  saveRDS(object, tmp)

  if (file.exists(path)) {
    removed <- unlink(path)
    if (!identical(removed, 0L) || file.exists(path)) {
      stop(
        "Could not replace existing RDS file: ", path,
        "\nClose any R/RStudio session using this project, then rerun the pipeline.",
        call. = FALSE
      )
    }
  }

  if (!file.rename(tmp, path)) {
    stop("Could not move temporary RDS file into place: ", path, call. = FALSE)
  }

  invisible(path)
}

write_csv_atomic <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp")
  readr::write_csv(data, tmp)
  file.rename(tmp, path)
  invisible(path)
}
