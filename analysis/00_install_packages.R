source("analysis/00_config.R")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

install_if_missing <- function(pkgs, installer) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    installer(missing)
  }
}

install_if_missing(CRAN_PACKAGES, function(pkgs) {
  install.packages(pkgs, repos = "https://cloud.r-project.org")
})

install_if_missing(BIOC_PACKAGES, function(pkgs) {
  BiocManager::install(pkgs, ask = FALSE, update = FALSE)
})

writeLines(capture.output(sessionInfo()), FILES$tcga_session)

