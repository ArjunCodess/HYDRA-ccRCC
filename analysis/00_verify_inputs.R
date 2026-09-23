source("analysis/00_config.R")
suppressPackageStartupMessages(library(readr))

path <- file.path(DIRS$tables, "input_manifest.csv")
if (!file.exists(path) || identical(Sys.getenv("HYDRA_FORCE_DOWNLOAD"), "1")) {
  message("Input checksums unavailable or force-download requested; manifest will be regenerated.")
} else {
  locked <- read_csv(path, show_col_types = FALSE)
  stopifnot(all(c("path", "bytes", "md5") %in% names(locked)),
            !anyDuplicated(locked$path), all(file.exists(locked$path)))
  sizes <- file.info(locked$path)$size
  if (any(sizes != locked$bytes)) stop("Cached input size differs from committed manifest.")
  current <- unname(tools::md5sum(locked$path))
  if (any(tolower(current) != tolower(locked$md5))) {
    stop("Cached input checksum differs from committed manifest: ",
         paste(locked$path[tolower(current) != tolower(locked$md5)], collapse = ", "))
  }
  message("Cached input checksums match: ", nrow(locked), " files.")
}
