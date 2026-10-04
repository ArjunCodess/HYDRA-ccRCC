# DESeq2 median ratios against a fixed training reference. Do not re-center
# factors using the held-out batch: individual predictions must be batch-invariant.
frozen_size_factors <- function(counts, geo_means) {
  if (nrow(counts) != length(geo_means)) stop("Reference and counts differ.")
  if (any(!is.finite(counts) | counts < 0)) stop("Invalid counts.")
  log_reference <- log(geo_means)
  factors <- apply(counts, 2L, function(column) {
    keep <- is.finite(log_reference) & column > 0
    if (!any(keep)) stop("Sample has no usable reference genes.")
    exp(stats::median(log(column[keep]) - log_reference[keep]))
  })
  if (any(!is.finite(factors) | factors <= 0)) stop("Invalid size factors.")
  factors
}
