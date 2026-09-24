# Selection and scoring helpers for the nested benchmark.
# Callers attach survival before fitting Cox models.

event_stratified_folds <- function(event, k) {
  event <- as.integer(event)
  if (!all(event %in% c(0L, 1L))) stop("Event indicators must be 0 or 1.")
  out <- integer(length(event))
  for (e in c(0L, 1L)) {
    idx <- which(event == e)
    if (length(idx) == 0L) next
    out[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  out
}

with_random_seed <- function(seed, expr) {
  old <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(expr)
}

conservative_bh <- function(p, n_universe) {
  if (any(!is.finite(p))) stop("Screened p-values must be finite before FDR adjustment.")
  if (n_universe < length(p)) stop("FDR universe is smaller than the screened set.")
  padded <- c(p, rep(1, n_universe - length(p)))
  stats::p.adjust(padded, method = "BH")[seq_along(p)]
}

gene_cox <- function(dat, expr, covariates) {
  dat$expr <- expr
  tryCatch(suppressWarnings(survival::coxph(
    stats::as.formula(paste(
      "Surv(os_time, os_event) ~",
      paste(c("expr", covariates), collapse = "+")
    )),
    data = dat
  )), error = function(e) NULL)
}

cox_coefficient <- function(fit) {
  if (is.null(fit)) return(c(beta = NA_real_, p = NA_real_))
  c(beta = unname(stats::coef(fit)["expr"]),
    p = unname(summary(fit)$coefficients["expr", "Pr(>|z|)"]))
}

select_hydra_gene <- function(de, expr_train, dat, thresholds) {
  if (nrow(de) == 0L) return(list(gene_id = NA_character_, n_pass = 0L))
  if (ncol(expr_train) != nrow(dat)) stop("HYDRA expression and training rows differ.")
  main_model <- t(vapply(de$gene_id, function(id) {
    x <- as.numeric(expr_train[id, ])
    if (!is.finite(stats::sd(x)) || stats::sd(x) == 0) return(c(NA_real_, NA_real_))
    cox_coefficient(gene_cox(dat, as.numeric(scale(x)), c("age", "sex", "stage", "grade")))
  }, numeric(2)))
  colnames(main_model) <- c("main_beta", "main_p")
  de <- dplyr::bind_cols(de, as.data.frame(main_model))
  de$main_fdr <- stats::p.adjust(de$main_p, method = "BH")
  de <- dplyr::filter(
    de,
    is.finite(.data$main_fdr),
    .data$main_fdr < thresholds$high_confidence_survival_fdr,
    abs(.data$main_beta) >= thresholds$high_confidence_abs_log_hr,
    abs(.data$gse40435_log2fc) >= thresholds$min_geo_abs_log2fc,
    abs(.data$gse53757_log2fc) >= thresholds$min_geo_abs_log2fc
  )
  if (nrow(de) == 0L) return(list(gene_id = NA_character_, n_pass = 0L))
  sensitivity <- t(vapply(seq_len(nrow(de)), function(i) {
    x <- as.numeric(scale(expr_train[de$gene_id[i], ]))
    fits <- lapply(list(c("age", "sex", "stage"), c("age", "sex", "grade")),
                   function(vars) gene_cox(dat, x, vars))
    if (any(vapply(fits, is.null, logical(1)))) return(rep(NA_real_, 4L))
    c(vapply(fits, function(f) unname(stats::coef(f)["expr"]), numeric(1)),
      vapply(fits, function(f) summary(f)$coefficients["expr", "Pr(>|z|)"], numeric(1)))
  }, numeric(4)))
  colnames(sensitivity) <- c("stage_beta", "grade_beta", "stage_p", "grade_p")
  de <- dplyr::bind_cols(de, as.data.frame(sensitivity))
  de <- dplyr::filter(
    de,
    sign(.data$stage_beta) == sign(.data$main_beta),
    sign(.data$grade_beta) == sign(.data$main_beta),
    .data$stage_p < 0.05,
    .data$grade_p < 0.05
  )
  de <- dplyr::mutate(
    de,
    evidence_score = -log10(pmax(.data$main_fdr, .Machine$double.xmin)) +
      abs(.data$main_beta) + pmin(abs(.data$log2FoldChange), 5) / 5 +
      pmin(abs(.data$gse40435_log2fc), 3) / 3 +
      pmin(abs(.data$gse53757_log2fc), 3) / 3
  )
  de <- dplyr::arrange(de, dplyr::desc(.data$evidence_score), .data$main_fdr, .data$gene_id)
  list(gene_id = if (nrow(de) == 0L) NA_character_ else de$gene_id[1], n_pass = nrow(de))
}

select_de_only_gene <- function(de_tcga, fdr_cutoff, min_abs_lfc) {
  if (nrow(de_tcga) == 0L) return(list(gene_id = NA_character_, n_pass = 0L))
  passed <- dplyr::filter(
    de_tcga,
    is.finite(.data$padj),
    is.finite(.data$log2FoldChange),
    .data$padj < fdr_cutoff,
    abs(.data$log2FoldChange) >= min_abs_lfc
  )
  passed <- dplyr::arrange(passed, .data$padj, dplyr::desc(abs(.data$log2FoldChange)), .data$gene_id)
  list(gene_id = if (nrow(passed) == 0L) NA_character_ else passed$gene_id[1],
       n_pass = nrow(passed))
}

residual_correlation <- function(expr_train, residuals) {
  if (is.null(rownames(expr_train))) stop("Expression rows need gene ids.")
  if (ncol(expr_train) != length(residuals)) stop("Residual length does not match expression.")
  r <- residuals - mean(residuals)
  r_ss <- sum(r^2)
  if (!is.finite(r_ss) || r_ss <= 0) stop("Clinical martingale residuals have no variation.")
  centered <- expr_train - rowMeans(expr_train)
  ss <- rowSums(centered^2)
  cor <- as.numeric((centered %*% r) / sqrt(ss * r_ss))
  cor[!is.finite(cor) | ss <= 0] <- NA_real_
  stats::setNames(cor, rownames(expr_train))
}

select_survival_only_gene <- function(expr_train, dat, abs_cor, thresholds, n_screen) {
  if (ncol(expr_train) != nrow(dat)) stop("Survival expression and training rows differ.")
  if (!identical(names(abs_cor), rownames(expr_train))) stop("Residual correlations are not aligned to genes.")
  finite <- names(abs_cor)[is.finite(abs_cor)]
  if (length(finite) == 0L) return(list(gene_id = NA_character_, n_pass = 0L, n_main_pass = 0L,
                                        n_universe = 0L, n_screened = 0L, screen_rank = NA_integer_))
  ordered_ids <- finite[order(-abs_cor[finite], finite)]
  screened <- ordered_ids[seq_len(min(n_screen, length(ordered_ids)))]
  main <- t(vapply(screened, function(id) {
    x <- as.numeric(expr_train[id, ])
    cox_coefficient(gene_cox(dat, as.numeric(scale(x)), c("age", "sex", "stage", "grade")))
  }, numeric(2)))
  colnames(main) <- c("main_beta", "main_p")
  stats_df <- data.frame(gene_id = screened, main, stringsAsFactors = FALSE)
  stats_df$main_p[!is.finite(stats_df$main_p)] <- 1
  stats_df$main_fdr <- conservative_bh(stats_df$main_p, length(finite))
  eligible <- stats_df[
    is.finite(stats_df$main_beta) &
      stats_df$main_fdr < thresholds$high_confidence_survival_fdr &
      abs(stats_df$main_beta) >= thresholds$high_confidence_abs_log_hr,
    , drop = FALSE
  ]
  n_main <- nrow(eligible)
  if (n_main == 0L) {
    return(list(gene_id = NA_character_, n_pass = 0L, n_main_pass = 0L,
                n_universe = length(finite), n_screened = length(screened),
                screen_rank = NA_integer_))
  }
  sensitivity <- t(vapply(eligible$gene_id, function(id) {
    x <- as.numeric(scale(expr_train[id, ]))
    fits <- lapply(list(c("age", "sex", "stage"), c("age", "sex", "grade")),
                   function(vars) gene_cox(dat, x, vars))
    if (any(vapply(fits, is.null, logical(1)))) return(rep(NA_real_, 4L))
    c(vapply(fits, function(f) unname(stats::coef(f)["expr"]), numeric(1)),
      vapply(fits, function(f) summary(f)$coefficients["expr", "Pr(>|z|)"], numeric(1)))
  }, numeric(4)))
  colnames(sensitivity) <- c("stage_beta", "grade_beta", "stage_p", "grade_p")
  eligible <- cbind(eligible, sensitivity)
  eligible <- eligible[
    is.finite(eligible$stage_p) & is.finite(eligible$grade_p) &
      sign(eligible$stage_beta) == sign(eligible$main_beta) &
      sign(eligible$grade_beta) == sign(eligible$main_beta) &
      eligible$stage_p < 0.05 & eligible$grade_p < 0.05,
    , drop = FALSE
  ]
  if (nrow(eligible) == 0L) {
    return(list(gene_id = NA_character_, n_pass = 0L, n_main_pass = n_main,
                n_universe = length(finite), n_screened = length(screened),
                screen_rank = NA_integer_))
  }
  eligible$evidence <- -log10(pmax(eligible$main_fdr, .Machine$double.xmin)) + abs(eligible$main_beta)
  eligible <- eligible[order(-eligible$evidence, eligible$main_fdr, eligible$gene_id), , drop = FALSE]
  chosen <- eligible$gene_id[1]
  list(gene_id = chosen, n_pass = nrow(eligible), n_main_pass = n_main,
       n_universe = length(finite), n_screened = length(screened),
       screen_rank = match(chosen, ordered_ids))
}

select_matched_control <- function(mean_expr, abs_cor, target, excluded) {
  if (!identical(names(mean_expr), names(abs_cor))) stop("Control inputs are not aligned.")
  ok <- is.finite(mean_expr) & is.finite(abs_cor)
  if (!any(ok)) return(NA_character_)
  if (!is.finite(target)) target <- stats::median(mean_expr[ok])
  cutoff <- stats::median(abs_cor[ok])
  ids <- names(mean_expr)
  pool <- ids[ok & abs_cor <= cutoff & !ids %in% excluded]
  if (length(pool) == 0L) pool <- ids[ok & !ids %in% excluded]
  if (length(pool) == 0L) return(NA_character_)
  pool <- sort(pool)
  pool[which.min(abs(mean_expr[pool] - target))]
}

# Training-only ridge Cox. The caller must pass training rows only.
# Clinical columns come first and are unpenalized. Gene columns are ridge
# (alpha = 0), so this arm does not zero coefficients or pick a subset.
# lambda.min minimizes inner-fold partial-likelihood deviance, not the C-index
# and not the 1-se rule. glmnet standardizes every column from x_train and
# applies those means and sds to x_test. Test outcomes are not arguments.
fit_ridge_cox <- function(x_train, x_test, time, event, foldid, n_unpenalized, horizon = 1095) {
  if (!requireNamespace("glmnet", quietly = TRUE)) stop("Package glmnet is required for the ridge arm.")
  if (nrow(x_train) != length(time) || nrow(x_train) != length(event) || nrow(x_train) != length(foldid)) {
    stop("Ridge training design is misaligned.")
  }
  if (ncol(x_train) != ncol(x_test)) stop("Ridge design matrices differ.")
  if (n_unpenalized < 1L || n_unpenalized >= ncol(x_train)) stop("Ridge penalty split is invalid.")
  if (length(unique(foldid)) < 2L || any(foldid < 1L)) stop("Ridge inner folds must be a training-only partition.")
  penalty <- c(rep(0, n_unpenalized), rep(1, ncol(x_train) - n_unpenalized))
  fit <- glmnet::cv.glmnet(
    x = x_train,
    y = survival::Surv(time, event),
    family = "cox",
    alpha = 0,
    penalty.factor = penalty,
    foldid = foldid,
    standardize = TRUE,
    type.measure = "deviance",
    cox.ties = "efron"
  )
  lp_test <- as.numeric(stats::predict(fit, newx = x_test, s = "lambda.min", type = "link"))
  curve <- survival::survfit(fit, s = "lambda.min", newx = x_test, x = x_train,
                             y = survival::Surv(time, event))
  surv <- curve$surv
  if (is.null(dim(surv))) surv <- matrix(surv, ncol = 1L)
  if (ncol(surv) != nrow(x_test)) stop("Ridge survival curves do not match the test rows.")
  eligible <- which(curve$time <= horizon)
  risk_test <- if (length(eligible) == 0L) {
    rep(0, ncol(surv))
  } else {
    1 - surv[max(eligible), ]
  }
  list(
    lp_test = lp_test,
    risk_test = as.numeric(risk_test),
    lambda = fit$lambda.min,
    n_genes = ncol(x_train) - n_unpenalized
  )
}

cox_risk <- function(fit, newdata, horizon = 1095) {
  lp <- as.numeric(stats::predict(fit, newdata = newdata, type = "lp", reference = "zero"))
  bh <- survival::basehaz(fit, centered = FALSE)
  h3 <- if (any(bh$time <= horizon)) utils::tail(bh$hazard[bh$time <= horizon], 1) else 0
  list(lp = lp, risk3 = 1 - exp(-h3 * exp(lp)))
}

score_one <- function(os_time, os_event, lp, risk3, horizon = 1095) {
  if (length(os_time) == 0L || any(!is.finite(lp)) || any(!is.finite(risk3))) {
    stop("Cannot score predictions with missing risk.")
  }
  pred <- data.frame(os_time = os_time, os_event = os_event, lp = lp, risk3 = risk3)
  pred$outcome <- survival::Surv(pred$os_time, pred$os_event)
  ci <- unname(survival::concordance(outcome ~ lp, data = pred, reverse = TRUE)$concordance)
  pred$censor_outcome <- survival::Surv(pred$os_time, 1L - pred$os_event)
  cens <- survival::survfit(censor_outcome ~ 1, data = pred)
  censor_surv <- function(t) {
    i <- findInterval(t, cens$time)
    ifelse(i == 0L, 1, cens$surv[pmax(i, 1L)])
  }
  target <- as.numeric(pred$os_time <= horizon & pred$os_event == 1L)
  weight <- ifelse(target == 1, 1 / pmax(censor_surv(pmax(pred$os_time - 1, 0)), 0.01),
                   ifelse(pred$os_time > horizon, 1 / pmax(censor_surv(horizon), 0.01), 0))
  slope <- tryCatch(
    unname(stats::coef(survival::coxph(outcome ~ lp, data = pred))[1]),
    error = function(e) NA_real_
  )
  list(c = ci, brier3 = mean(weight * (target - pred$risk3)^2), calibration_slope = slope)
}
