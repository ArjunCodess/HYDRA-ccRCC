source("analysis/00_config.R")
source("analysis/functions/frozen_normalization.R")
set.seed(RESAMPLING$seed + 47L)
train <- matrix(rpois(600, 50), 100, 6)
train[1, 1] <- 0
reference <- exp(rowMeans(log(train)))
stopifnot(isTRUE(all.equal(unname(frozen_size_factors(train, reference)),
                          unname(DESeq2::estimateSizeFactorsForMatrix(train)), tolerance = 1e-12)))
test <- matrix(rpois(300, 75), 100, 3)
batch <- frozen_size_factors(test, reference)
alone <- frozen_size_factors(test[, 1, drop = FALSE], reference)
stopifnot(abs(batch[1] - alone[1]) < 1e-12,
          max(abs(frozen_size_factors(3*test, reference) - 3*batch)) < 1e-12)
# Demonstrate the old batch-centering operation has a different deployment rule.
old_batch <- DESeq2::estimateSizeFactorsForMatrix(test, geoMeans = reference)
old_alone <- DESeq2::estimateSizeFactorsForMatrix(test[, 1, drop = FALSE], geoMeans = reference)
stopifnot(old_alone[1] == 1, abs(old_batch[1] - old_alone[1]) > 1e-6)
message("PASS: training equivalence, individual batch invariance, library-depth scaling; old behavior reproduced.")
