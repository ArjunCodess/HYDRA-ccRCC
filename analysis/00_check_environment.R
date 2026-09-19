source("analysis/00_config.R")

locked <- readr::read_csv("environment/package_versions.csv", show_col_types = FALSE)
installed <- as.data.frame(installed.packages()[, c("Package", "Version")])
installed <- installed[!duplicated(installed$Package), ]
observed <- installed$Version[match(locked$Package, installed$Package)]
bad <- is.na(observed) | observed != locked$Version
if (any(bad)) {
  stop("R package lock differs: ", paste(locked$Package[bad], collapse = ", "))
}
message("R package versions match the committed lock: ", nrow(locked), " packages.")
