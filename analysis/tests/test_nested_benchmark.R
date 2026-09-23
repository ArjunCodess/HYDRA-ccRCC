source("analysis/00_config.R")
source("analysis/functions/nested_benchmark.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(glmnet)
})

set.seed(31)
event <- rep(c(0L, 1L), c(23, 17))
folds <- event_stratified_folds(event, 5L)
event_counts <- tapply(event, folds, sum)
stopifnot(length(unique(folds)) == 5L, max(event_counts) - min(event_counts) <= 1L)

original_folds <- function(event, k) {
  out <- integer(length(event))
  for (e in c(0L, 1L)) {
    idx <- which(event == e)
    out[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  out
}
set.seed(99)
a <- original_folds(event, 5L)
set.seed(99)
b <- event_stratified_folds(event, 5L)
stopifnot(identical(a, b))

bh <- conservative_bh(c(0.001, 0.01), 100)
plain <- p.adjust(c(0.001, 0.01), method = "BH")
stopifnot(all(bh >= plain), all(bh <= 1),
          isTRUE(all.equal(conservative_bh(c(0.01, 0.02), 2),
                           p.adjust(c(0.01, 0.02), method = "BH"))))

de <- tibble(
  gene_id = c("B", "A", "C"),
  padj = c(0.01, 0.01, 0.2),
  log2FoldChange = c(1.2, 2.4, 3)
)
picked <- select_de_only_gene(de, 0.05, 1)
stopifnot(picked$gene_id == "A", picked$n_pass == 2L)
tied <- tibble(gene_id = c("B", "A"), padj = c(0.01, 0.01), log2FoldChange = c(2, 2))
stopifnot(select_de_only_gene(tied, 0.05, 1)$gene_id == "A")

mean_expr <- c(G1 = 1, G2 = 5, G3 = 5.1, G4 = 9)
abs_cor <- c(G1 = 0.01, G2 = 0.02, G3 = 0.2, G4 = 0.9)
stopifnot(select_matched_control(mean_expr, abs_cor, target = 5, excluded = "G4") == "G2")
stopifnot(select_matched_control(mean_expr, abs_cor, target = 5, excluded = c("G1", "G2")) == "G3")

# Cox three-year risk uses the same reference-zero linear predictor as the primary nested CV.
set.seed(7)
n <- 80
x <- rnorm(n)
time <- rexp(n, exp(0.4 * x))
event_i <- rbinom(n, 1, 0.75)
dat <- data.frame(os_time = time, os_event = event_i, x = x)
fit <- coxph(Surv(os_time, os_event) ~ x, data = dat)
direct <- cox_risk(fit, dat, horizon = 0.5)
lp <- as.numeric(predict(fit, type = "lp", reference = "zero"))
stopifnot(max(abs(direct$lp - lp)) < 1e-8, all(direct$risk3 > 0), all(direct$risk3 < 1))

scored <- score_one(dat$os_time, dat$os_event, direct$lp, direct$risk3, horizon = 0.5)
direct_c <- unname(concordance(Surv(dat$os_time, dat$os_event) ~ direct$lp, reverse = TRUE)$concordance)
stopifnot(isTRUE(all.equal(scored$c, direct_c)), is.finite(scored$brier3))

# HYDRA keeps the gene with the survival signal when both pass the expression gates.
set.seed(11)
n <- 160
latent <- rnorm(n)
os_time <- rexp(n, 0.03 * exp(1.3 * latent))
censor <- rexp(n, 0.015)
keep <- os_time < censor
os_event <- as.integer(keep)
os_time <- pmin(os_time, censor)
train <- data.frame(
  os_time = os_time, os_event = os_event,
  age = rnorm(n, 60, 8),
  sex = factor(rep(c("female", "male"), length.out = n)),
  stage = factor(rep(c("Stage I", "Stage II", "Stage III", "Stage IV"), length.out = n)),
  grade = factor(rep(c("Low grade", "High grade"), length.out = n))
)
expr <- rbind(STRONG = latent + rnorm(n, sd = 0.05), NOISE = rnorm(n))
de_h <- tibble(
  gene_id = c("STRONG", "NOISE"),
  log2FoldChange = c(2, 2),
  gse40435_log2fc = c(1, 1),
  gse53757_log2fc = c(1, 1)
)
hydra <- select_hydra_gene(de_h, expr, train, THRESHOLDS)
stopifnot(hydra$gene_id == "STRONG", hydra$n_pass >= 1L)

null_expr <- rbind(expr, NOISE2 = rnorm(n), NOISE3 = rnorm(n))
resid <- residuals(coxph(Surv(os_time, os_event) ~ age + sex + stage + grade, data = train),
                   type = "martingale")
abs_cor <- abs(residual_correlation(null_expr, resid))
survival_pick <- select_survival_only_gene(null_expr, train, abs_cor, THRESHOLDS, n_screen = 4L)
stopifnot(survival_pick$gene_id == "STRONG", survival_pick$screen_rank == 1L)

# Ridge returns finite held-out risks and does not use the test outcomes.
set.seed(3)
n <- 70
clin <- cbind(age = rnorm(n), sex = rep(c(-0.5, 0.5), length.out = n))
genes <- matrix(rnorm(n * 4), n)
signal <- genes[, 1]
time <- rexp(n, exp(0.8 * signal))
event_r <- rbinom(n, 1, 0.7)
x <- cbind(clin, genes)
foldid <- event_stratified_folds(event_r[1:56], 5L)
ridge <- fit_ridge_cox(x[1:56, ], x[57:70, ], time[1:56], event_r[1:56], foldid,
                        n_unpenalized = 2L, horizon = 0.5)
stopifnot(ridge$n_genes == 4L, is.finite(ridge$lambda),
          length(ridge$lp_test) == 14L, all(is.finite(ridge$risk_test)),
          all(ridge$risk_test > 0), all(ridge$risk_test < 1))

message("Nested benchmark helper tests passed.")
