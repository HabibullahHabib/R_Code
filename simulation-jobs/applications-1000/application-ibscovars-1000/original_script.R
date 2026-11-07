# ============================================================================
# APPLICATION 1: IBScovars
# Shrinkage vs non-shrinkage under controlled missingness
#
# PRIMARY IDEA
# ------------
# IBScovars is complete. We keep ALL 369 observations and use the complete data
# as a benchmark. We then impose MCAR missingness WITHIN dose arm at
# 10%, 20%, 30%, and 40%.
#
# Why this is preferable to taking a smaller subset:
#   * no real information is discarded;
#   * the true complete-data sample contrasts are known;
#   * we can directly measure recovery after missingness/imputation;
#   * 30% is the primary reference-like application scenario;
#   * 10-40% gives a sensitivity analysis analogous to the paper simulations.
#
# MI METHODS
# ----------
# REG: proper Bayesian regression within dose arm, response ~ gender
# PMM: Type-1 PMM within dose arm, K=5, response ~ gender
# PROP: propensity subclasses + ABB, donors restricted within dose arm
#
# Because gender is only a weak outcome predictor in this dataset, the
# missingness experiment should be viewed mainly as a controlled benchmark
# application, not as evidence that IBScovars has a rich MAR structure.
# ============================================================================


# ============================================================================
# COMMON CORE: MANY-TO-ONE MI WITH SHRINKAGE VS NON-SHRINKAGE
# ============================================================================
# This implementation follows the logic used in the thesis/paper simulation:
#   incomplete data -> multiple imputation -> Q^(m), U^(m)
#   -> Rubin pooling -> multivariate-t simultaneous intervals.
#
# IMPORTANT APPLICATION EXTENSION:
# The simulation code used equal n per arm. These real datasets are unbalanced,
# so the James-Stein-type shrinkage intensity is generalized by using the
# arm-specific n_g in Var(s_g^2) ~= 2*s_g^4/(n_g-1).
# ============================================================================

check_packages <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    stop(
      "Please install the following packages first: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
}

fmt_pct <- function(x) paste0(round(100*x), "%")

build_model_matrix <- function(dat, covariates) {
  if (!length(covariates)) {
    return(matrix(1, nrow(dat), 1L,
                  dimnames = list(NULL, "(Intercept)")))
  }
  f <- stats::reformulate(covariates)
  stats::model.matrix(f, data = dat)
}

draw_bayesian_lm <- function(X_obs, y_obs) {
  X_obs <- as.matrix(X_obs)
  y_obs <- as.numeric(y_obs)

  n <- nrow(X_obs)
  p <- ncol(X_obs)

  if (n <= p) {
    stop(sprintf("Too few observed outcomes for Bayesian regression: n=%d, p=%d.", n, p))
  }

  XtX <- crossprod(X_obs)
  if (qr(XtX)$rank < p) {
    stop("Regression design matrix is rank deficient.")
  }

  XtX_inv <- solve(XtX)
  beta_hat <- drop(XtX_inv %*% crossprod(X_obs, y_obs))
  resid <- y_obs - drop(X_obs %*% beta_hat)
  df <- n - p
  sse <- sum(resid^2)

  if (!is.finite(sse) || sse <= 0 || df <= 0) {
    stop("Invalid regression fit for proper MI.")
  }

  sigma2_star <- sse / stats::rchisq(1L, df = df)
  V_beta <- sigma2_star * XtX_inv
  L <- t(chol(V_beta))
  beta_star <- beta_hat + drop(L %*% stats::rnorm(p))

  list(
    beta_hat = beta_hat,
    beta_star = beta_star,
    sigma2_star = sigma2_star
  )
}

impute_regression_once <- function(dat, outcome_covariates, seed) {
  set.seed(seed)
  out <- dat
  X_all <- build_model_matrix(out, outcome_covariates)

  for (g in levels(out$group)) {
    idx <- which(out$group == g)
    obs <- idx[!is.na(out$y[idx])]
    mis <- idx[ is.na(out$y[idx])]

    if (!length(mis)) next

    dr <- draw_bayesian_lm(X_all[obs, , drop = FALSE], out$y[obs])
    mu_mis <- drop(X_all[mis, , drop = FALSE] %*% dr$beta_star)

    out$y[mis] <- stats::rnorm(
      length(mis),
      mean = mu_mis,
      sd = sqrt(dr$sigma2_star)
    )
  }

  if (anyNA(out$y)) stop("Regression imputation left missing outcomes.")
  out
}

impute_pmm_once <- function(dat, outcome_covariates, seed, donors = 5L) {
  set.seed(seed)
  out <- dat
  if (!anyNA(out$y)) return(out)

  X_all <- build_model_matrix(out, outcome_covariates)

  for (g in levels(out$group)) {
    idx <- which(out$group == g)
    obs <- idx[!is.na(out$y[idx])]
    mis <- idx[ is.na(out$y[idx])]

    if (!length(mis)) next
    if (!length(obs)) stop("No observed PMM donors in group ", g, ".")

    dr <- draw_bayesian_lm(X_all[obs, , drop = FALSE], out$y[obs])

    # Type-1 PMM:
    # observed predictions use beta_hat,
    # missing predictions use beta_star.
    pred_obs <- drop(X_all[obs, , drop = FALSE] %*% dr$beta_hat)
    pred_mis <- drop(X_all[mis, , drop = FALSE] %*% dr$beta_star)

    k_use <- min(as.integer(donors), length(obs))

    for (j in seq_along(mis)) {
      dist <- abs(pred_obs - pred_mis[j])
      # Random secondary key prevents row order from deciding exact ties.
      ord <- order(dist, stats::runif(length(dist)))
      pool <- obs[ord[seq_len(k_use)]]
      donor <- sample(pool, 1L)
      out$y[mis[j]] <- out$y[donor]
    }
  }

  if (anyNA(out$y)) stop("PMM left missing outcomes.")
  out
}

approx_bayesian_bootstrap <- function(y_obs, n_missing) {
  if (n_missing == 0L) return(numeric(0))
  if (!length(y_obs)) stop("No ABB donor outcomes available.")

  y_star <- sample(y_obs, size = length(y_obs), replace = TRUE)
  sample(y_star, size = n_missing, replace = TRUE)
}

make_propensity_strata <- function(ps, n_groups = 5L) {
  if (any(!is.finite(ps))) stop("Non-finite propensity score.")

  unique_ps <- sort(unique(ps))
  if (length(unique_ps) == 1L) return(rep.int(1L, length(ps)))

  G <- min(as.integer(n_groups), length(unique_ps))
  score_rank <- rank(unique_ps, ties.method = "average")
  score_group <- ceiling(G * score_rank / length(unique_ps))
  score_group <- pmin(G, pmax(1L, as.integer(score_group)))

  map <- setNames(
    score_group,
    format(unique_ps, digits = 17, scientific = FALSE)
  )

  as.integer(map[format(ps, digits = 17, scientific = FALSE)])
}

impute_propensity_once <- function(dat, propensity_covariates, seed,
                                   n_groups = 5L) {
  set.seed(seed)
  out <- dat
  mis <- is.na(out$y)

  if (!any(mis)) return(out)

  psdat <- out[, c("group", propensity_covariates), drop = FALSE]
  psdat$missing <- as.integer(mis)

  rhs <- c("group", propensity_covariates)
  f <- stats::reformulate(rhs, response = "missing")

  fit <- suppressWarnings(
    stats::glm(f, family = stats::binomial(), data = psdat)
  )

  ps <- as.numeric(stats::predict(fit, newdata = psdat, type = "response"))

  if (any(!is.finite(ps))) {
    stop("Propensity model produced invalid probabilities.")
  }

  for (g in levels(out$group)) {
    idx_g <- which(out$group == g)
    mis_g <- idx_g[is.na(out$y[idx_g])]
    obs_g <- idx_g[!is.na(out$y[idx_g])]

    if (!length(mis_g)) next
    if (!length(obs_g)) stop("No observed propensity donors in group ", g, ".")

    strata <- make_propensity_strata(ps[idx_g], n_groups = n_groups)
    strata_ids <- sort(unique(strata))

    centers <- vapply(
      strata_ids,
      function(s) mean(ps[idx_g[strata == s]]),
      numeric(1)
    )
    names(centers) <- as.character(strata_ids)

    for (s in strata_ids) {
      ids_s <- idx_g[strata == s]
      mis_s <- ids_s[is.na(out$y[ids_s])]
      if (!length(mis_s)) next

      obs_s <- ids_s[!is.na(out$y[ids_s])]

      # If a subclass has no donor, borrow only from the nearest nonempty
      # subclass in the SAME treatment arm.
      if (!length(obs_s)) {
        donor_strata <- strata_ids[
          vapply(
            strata_ids,
            function(ss) {
              ids_ss <- idx_g[strata == ss]
              any(!is.na(out$y[ids_ss]))
            },
            logical(1)
          )
        ]

        if (!length(donor_strata)) {
          stop("No nonempty propensity donor subclass in group ", g, ".")
        }

        nearest <- donor_strata[
          which.min(abs(
            centers[as.character(donor_strata)] -
              centers[as.character(s)]
          ))
        ]

        ids_near <- idx_g[strata == nearest]
        obs_s <- ids_near[!is.na(out$y[ids_near])]
      }

      out$y[mis_s] <- approx_bayesian_bootstrap(
        out$y[obs_s],
        length(mis_s)
      )
    }
  }

  if (anyNA(out$y)) stop("Propensity imputation left missing outcomes.")
  out
}

make_imputations <- function(dat, method, m, base_seed,
                             outcome_covariates,
                             propensity_covariates,
                             pmm_donors = 5L,
                             propensity_ngroups = 5L) {
  out <- vector("list", m)

  for (mm in seq_len(m)) {
    ss <- base_seed + mm - 1L

    out[[mm]] <- switch(
      method,
      regression = impute_regression_once(
        dat, outcome_covariates, seed = ss
      ),
      pmm = impute_pmm_once(
        dat, outcome_covariates, seed = ss, donors = pmm_donors
      ),
      propensity = impute_propensity_once(
        dat, propensity_covariates, seed = ss,
        n_groups = propensity_ngroups
      ),
      stop("Unknown imputation method: ", method)
    )
  }

  out
}

estimate_alpha_hat <- function(s2_vec, n_vec) {
  s2_vec <- as.numeric(s2_vec)
  n_vec <- as.numeric(n_vec)

  G <- length(s2_vec)
  if (length(n_vec) != G) stop("s2_vec and n_vec have different lengths.")
  if (any(n_vec <= 1)) return(0)

  # Unequal-n extension of the thesis simulation formula.
  num <- ((G - 1) / G) * sum(2 * s2_vec^2 / (n_vec - 1))
  den <- sum((s2_vec - mean(s2_vec))^2)

  if (!is.finite(den) || den <= 1e-12) return(0)

  a <- num / den
  min(max(a, 0), 1)
}

shrink_variances <- function(s2_vec, n_vec) {
  a <- estimate_alpha_hat(s2_vec, n_vec)
  target <- mean(s2_vec)
  list(
    s2 = a * target + (1 - a) * s2_vec,
    alpha = a
  )
}

build_U <- function(s2_vec, n_vec, group_levels, control_level) {
  treat_levels <- setdiff(group_levels, control_level)
  labels <- paste0(treat_levels, " vs ", control_level)
  k <- length(treat_levels)

  U <- matrix(
    s2_vec[control_level] / n_vec[control_level],
    nrow = k, ncol = k
  )
  diag(U) <- diag(U) +
    s2_vec[treat_levels] / n_vec[treat_levels]

  rownames(U) <- colnames(U) <- labels
  U
}

make_corr_safe <- function(V) {
  V <- as.matrix(V)
  V <- (V + t(V)) / 2

  ev <- tryCatch(
    eigen(V, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NA_real_
  )

  if (any(!is.finite(ev)) || min(ev, na.rm = TRUE) <= 1e-10) {
    min_ev <- suppressWarnings(min(ev, na.rm = TRUE))
    ridge <- if (is.finite(min_ev)) abs(min_ev) + 1e-8 else 1e-6
    diag(V) <- diag(V) + ridge
  }

  C <- stats::cov2cor(V)
  C[!is.finite(C)] <- 0
  diag(C) <- 1
  C
}

qmvt_critical <- function(corr, df, alpha = 0.05,
                          maxpts = 50000L, abseps = 1e-6,
                          seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  as.numeric(
    mvtnorm::qmvt(
      1 - alpha,
      corr = corr,
      df = df,
      tail = "both.tails",
      algorithm = mvtnorm::GenzBretz(
        maxpts = maxpts,
        abseps = abseps
      )
    )$quantile
  )
}

analyze_complete_many_to_one <- function(dat, use_shrinkage,
                                         control_level,
                                         alpha = 0.05,
                                         qseed = NULL) {
  group_levels <- levels(dat$group)
  treat_levels <- setdiff(group_levels, control_level)
  labels <- paste0(treat_levels, " vs ", control_level)

  if (anyNA(dat$y)) stop("analyze_complete_many_to_one requires complete y.")

  split_y <- split(dat$y, dat$group)
  means <- vapply(split_y, mean, numeric(1))[group_levels]
  s2 <- vapply(split_y, stats::var, numeric(1))[group_levels]
  n_vec <- vapply(split_y, length, integer(1))[group_levels]

  qhat <- means[treat_levels] - means[control_level]
  names(qhat) <- labels

  alpha_hat <- NA_real_
  s2_used <- s2

  if (use_shrinkage) {
    shr <- shrink_variances(s2, n_vec)
    s2_used <- setNames(shr$s2, names(s2))
    alpha_hat <- shr$alpha
  }

  U <- build_U(
    s2_used, n_vec,
    group_levels = group_levels,
    control_level = control_level
  )

  # Same complete-data df convention as the thesis simulation code.
  df <- sum(n_vec) - length(n_vec)
  crit <- qmvt_critical(make_corr_safe(U), df, alpha, seed = qseed)
  se <- sqrt(diag(U))

  tibble::tibble(
    contrast = labels,
    estimate = as.numeric(qhat),
    se = as.numeric(se),
    lower = as.numeric(qhat - crit * se),
    upper = as.numeric(qhat + crit * se),
    width = as.numeric(2 * crit * se),
    significant = (qhat - crit * se > 0) | (qhat + crit * se < 0),
    shrink_alpha = alpha_hat
  )
}

summarize_completions <- function(comps) {
  group_levels <- levels(comps[[1]]$group)
  m <- length(comps)

  means <- matrix(
    NA_real_, m, length(group_levels),
    dimnames = list(NULL, group_levels)
  )
  vars <- means

  for (i in seq_len(m)) {
    means[i, ] <- tapply(comps[[i]]$y, comps[[i]]$group, mean)[group_levels]
    vars[i, ] <- tapply(comps[[i]]$y, comps[[i]]$group, stats::var)[group_levels]
  }

  n_vec <- table(comps[[1]]$group)
  n_vec <- as.numeric(n_vec[group_levels])
  names(n_vec) <- group_levels

  list(means = means, vars = vars, n = n_vec, m = m)
}

pool_mi_many_to_one_both <- function(comps, control_level,
                                     alpha = 0.05,
                                     seed_shrink = NULL,
                                     seed_nonshrink = NULL) {
  summ <- summarize_completions(comps)
  group_levels <- colnames(summ$means)
  treat_levels <- setdiff(group_levels, control_level)
  labels <- paste0(treat_levels, " vs ", control_level)
  k <- length(treat_levels)
  m <- summ$m

  Qmat <- matrix(
    NA_real_, nrow = m, ncol = k,
    dimnames = list(NULL, labels)
  )
  U_ns <- array(NA_real_, dim = c(k, k, m))
  U_sh <- array(NA_real_, dim = c(k, k, m))
  alpha_vec <- numeric(m)

  for (i in seq_len(m)) {
    mu <- summ$means[i, ]
    s2 <- summ$vars[i, ]
    names(s2) <- group_levels

    Qmat[i, ] <- mu[treat_levels] - mu[control_level]

    U_ns[, , i] <- build_U(
      s2, summ$n, group_levels, control_level
    )

    shr <- shrink_variances(s2, summ$n)
    s2_sh <- setNames(shr$s2, group_levels)
    alpha_vec[i] <- shr$alpha

    U_sh[, , i] <- build_U(
      s2_sh, summ$n, group_levels, control_level
    )
  }

  Qbar <- colMeans(Qmat)
  B <- stats::cov(Qmat)

  pool_one <- function(U_array, qseed, variance_method, alpha_value) {
    Ubar <- apply(U_array, c(1, 2), mean)
    Tm <- Ubar + (1 + 1/m) * B
    Tm <- (Tm + t(Tm)) / 2

    rm <- tryCatch(
      (1 + 1/m) *
        (sum(diag(B %*% solve(Ubar))) / k),
      error = function(e) NA_real_
    )

    if (!is.finite(rm) || rm <= 1e-12) {
      nu <- 1e6
    } else {
      nu <- max(1, floor((m - 1) * (1 + 1/rm)^2))
    }

    crit <- qmvt_critical(
      make_corr_safe(Tm), nu, alpha, seed = qseed
    )

    se <- sqrt(diag(Tm))

    tibble::tibble(
      contrast = labels,
      estimate = as.numeric(Qbar),
      se = as.numeric(se),
      lower = as.numeric(Qbar - crit * se),
      upper = as.numeric(Qbar + crit * se),
      width = as.numeric(2 * crit * se),
      significant = (Qbar - crit * se > 0) | (Qbar + crit * se < 0),
      variance_method = variance_method,
      shrink_alpha = alpha_value
    )
  }

  dplyr::bind_rows(
    pool_one(
      U_sh, seed_shrink, "Shrinkage",
      mean(alpha_vec, na.rm = TRUE)
    ),
    pool_one(
      U_ns, seed_nonshrink, "Non-shrinkage",
      NA_real_
    )
  )
}

make_complete_results_both <- function(dat, control_level,
                                       method_label, scenario_label,
                                       missing_rate,
                                       seed_base = 1000L) {
  sh <- analyze_complete_many_to_one(
    dat, TRUE, control_level,
    qseed = seed_base + 1L
  ) |>
    dplyr::mutate(variance_method = "Shrinkage")

  ns <- analyze_complete_many_to_one(
    dat, FALSE, control_level,
    qseed = seed_base + 2L
  ) |>
    dplyr::mutate(
      variance_method = "Non-shrinkage",
      shrink_alpha = NA_real_
    )

  dplyr::bind_rows(sh, ns) |>
    dplyr::mutate(
      scenario = scenario_label,
      missing_rate = missing_rate,
      method = method_label,
      .before = 1
    )
}

make_mi_results <- function(dat_miss, method, method_label,
                            m, control_level,
                            outcome_covariates,
                            propensity_covariates,
                            pmm_donors,
                            propensity_ngroups,
                            scenario_label,
                            missing_rate,
                            seed_base) {
  comps <- make_imputations(
    dat_miss,
    method = method,
    m = m,
    base_seed = seed_base,
    outcome_covariates = outcome_covariates,
    propensity_covariates = propensity_covariates,
    pmm_donors = pmm_donors,
    propensity_ngroups = propensity_ngroups
  )

  pool_mi_many_to_one_both(
    comps,
    control_level = control_level,
    seed_shrink = seed_base + 100000L,
    seed_nonshrink = seed_base + 200000L
  ) |>
    dplyr::mutate(
      scenario = scenario_label,
      missing_rate = missing_rate,
      method = method_label,
      .before = 1
    )
}

impose_mcar_within_arm <- function(dat, rate, seed) {
  if (rate <= 0) return(dat)
  if (rate >= 1) stop("rate must be < 1.")

  set.seed(seed)
  out <- dat

  for (g in levels(out$group)) {
    idx <- which(out$group == g)
    n_remove <- round(rate * length(idx))
    if (n_remove > 0) {
      drop_idx <- sample(idx, n_remove, replace = FALSE)
      out$y[drop_idx] <- NA_real_
    }
  }

  out
}

benchmark_contrasts <- function(dat, control_level) {
  group_levels <- levels(dat$group)
  treat_levels <- setdiff(group_levels, control_level)
  means <- tapply(dat$y, dat$group, mean)
  q <- means[treat_levels] - means[control_level]
  names(q) <- paste0(treat_levels, " vs ", control_level)
  q
}

add_benchmark_metrics <- function(res, truth) {
  res |>
    dplyr::mutate(
      benchmark = unname(truth[contrast]),
      error = estimate - benchmark,
      abs_error = abs(error),
      benchmark_covered = lower <= benchmark & upper >= benchmark
    )
}

summarize_validation <- function(raw) {
  raw |>
    dplyr::group_by(
      missing_rate, method, variance_method
    ) |>
    dplyr::summarise(
      n_rows = dplyr::n(),
      rmse = sqrt(mean(error^2, na.rm = TRUE)),
      mean_abs_error = mean(abs_error, na.rm = TRUE),
      benchmark_recovery_pct = 100 * mean(benchmark_covered, na.rm = TRUE),
      avg_ci_width = mean(width, na.rm = TRUE),
      mean_shrink_alpha = mean(shrink_alpha, na.rm = TRUE),
      .groups = "drop"
    )
}

save_forest_plot <- function(res, filename, title, subtitle = NULL) {
  p <- ggplot2::ggplot(
    res,
    ggplot2::aes(
      y = contrast,
      x = estimate,
      xmin = lower,
      xmax = upper,
      shape = variance_method
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = 2,
      linewidth = 0.5
    ) +
    ggplot2::geom_pointrange(
      position = ggplot2::position_dodge(width = 0.55)
    ) +
    ggplot2::facet_wrap(~ method, scales = "free_x") +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Treatment - control mean difference (simultaneous 95% CI)",
      y = NULL,
      shape = "Variance method"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )

  ggplot2::ggsave(filename, p, width = 12, height = 7, dpi = 300)
}

write_session_info <- function(filename) {
  capture.output(utils::sessionInfo(), file = filename)
}


check_packages(c("DoseFinding", "mvtnorm", "dplyr", "tidyr",
                 "tibble", "ggplot2"))

suppressPackageStartupMessages({
  library(DoseFinding)
  library(mvtnorm)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

# ----------------------------- SETTINGS --------------------------------------
master_seed <- 20261107L
m_imputations <- 40L
pmm_donors <- 5L

# Application arms are much larger than the n=20 simulation arms, so five
# propensity subclasses are stable here.
propensity_ngroups <- 5L

missing_rates <- c(0.10, 0.20, 0.30, 0.40)

# Repeated deletion validation:
# FALSE = run only the main application tables/plots.
# TRUE  = also repeat the random deletion experiment many times.
run_repeated_validation <- TRUE
validation_reps <- 200L

# Change to TRUE for a fast code check before the publication run.
quick_test <- FALSE
if (quick_test) {
  m_imputations <- 10L
  validation_reps <- 10L
}

output_dir <- "IBScovars_application_output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------- LOAD DATA -------------------------------------
data(IBScovars, package = "DoseFinding")

ibs <- IBScovars |>
  transmute(
    group = factor(dose, levels = sort(unique(dose))),
    gender = factor(gender),
    y = as.numeric(resp)
  )

control_level <- levels(ibs$group)[1]   # dose 0
outcome_covariates <- c("gender")
propensity_covariates <- c("gender")

if (anyNA(ibs$y)) {
  stop("IBScovars unexpectedly contains missing response values.")
}

# ----------------------------- DATA SUMMARY ----------------------------------
data_summary <- ibs |>
  group_by(group) |>
  summarise(
    n = n(),
    mean = mean(y),
    sd = sd(y),
    variance = var(y),
    .groups = "drop"
  )

write.csv(
  data_summary,
  file.path(output_dir, "IBScovars_complete_data_summary.csv"),
  row.names = FALSE
)

variance_ratio <- max(data_summary$variance) / min(data_summary$variance)

dataset_info <- tibble(
  dataset = "IBScovars",
  n_total = nrow(ibs),
  n_groups = nlevels(ibs$group),
  natural_missing_pct = 0,
  complete_data_variance_ratio = variance_ratio
)

write.csv(
  dataset_info,
  file.path(output_dir, "IBScovars_dataset_info.csv"),
  row.names = FALSE
)

# Complete-data benchmark contrasts.
truth <- benchmark_contrasts(ibs, control_level)

write.csv(
  tibble(contrast = names(truth), complete_data_contrast = as.numeric(truth)),
  file.path(output_dir, "IBScovars_complete_benchmark_contrasts.csv"),
  row.names = FALSE
)

# Complete-data shrinkage and non-shrinkage simultaneous intervals.
complete_results <- make_complete_results_both(
  ibs,
  control_level = control_level,
  method_label = "Complete-data benchmark",
  scenario_label = "Complete data",
  missing_rate = 0,
  seed_base = master_seed + 1000L
)

# ----------------------- MAIN MISSINGNESS APPLICATION ------------------------
method_map <- c(
  "Regression MI" = "regression",
  "PMM MI" = "pmm",
  "Propensity-score MI" = "propensity"
)

main_results <- list(complete_results)
counter <- 2L

for (r_idx in seq_along(missing_rates)) {
  mr <- missing_rates[r_idx]

  dat_miss <- impose_mcar_within_arm(
    ibs,
    rate = mr,
    seed = master_seed + 10000L * r_idx
  )

  achieved <- mean(is.na(dat_miss$y))
  scenario_label <- paste0(
    "MCAR ", fmt_pct(mr),
    " (achieved ", fmt_pct(achieved), ")"
  )

  # Complete-case analysis is included as a benchmark comparator.
  cc <- dat_miss |>
    filter(!is.na(y)) |>
    droplevels()

  cc_res <- make_complete_results_both(
    cc,
    control_level = control_level,
    method_label = "Complete-case",
    scenario_label = scenario_label,
    missing_rate = mr,
    seed_base = master_seed + 200000L + r_idx * 100L
  )

  main_results[[counter]] <- add_benchmark_metrics(cc_res, truth)
  counter <- counter + 1L

  for (m_idx in seq_along(method_map)) {
    method_label <- names(method_map)[m_idx]
    method_code <- unname(method_map[m_idx])

    mi_res <- make_mi_results(
      dat_miss,
      method = method_code,
      method_label = method_label,
      m = m_imputations,
      control_level = control_level,
      outcome_covariates = outcome_covariates,
      propensity_covariates = propensity_covariates,
      pmm_donors = pmm_donors,
      propensity_ngroups = propensity_ngroups,
      scenario_label = scenario_label,
      missing_rate = mr,
      seed_base = master_seed +
        1000000L +
        r_idx * 100000L +
        m_idx * 10000L
    )

    main_results[[counter]] <- add_benchmark_metrics(mi_res, truth)
    counter <- counter + 1L
  }
}

# Add benchmark columns to complete-data results too.
complete_results2 <- add_benchmark_metrics(complete_results, truth)
main_results[[1]] <- complete_results2

application_results <- bind_rows(main_results) |>
  arrange(missing_rate, method, contrast, variance_method)

write.csv(
  application_results,
  file.path(output_dir, "IBScovars_application_results.csv"),
  row.names = FALSE
)

# ----------------------------- PRIMARY 30% TABLE -----------------------------
primary30 <- application_results |>
  filter(abs(missing_rate - 0.30) < 1e-12) |>
  select(
    scenario, method, variance_method, contrast,
    estimate, lower, upper, width, significant,
    benchmark, error, benchmark_covered, shrink_alpha
  )

write.csv(
  primary30,
  file.path(output_dir, "IBScovars_primary_30pct_table.csv"),
  row.names = FALSE
)

save_forest_plot(
  primary30,
  file.path(output_dir, "IBScovars_primary_30pct_forest.png"),
  title = "IBScovars: 30% random missingness",
  subtitle = "Many-to-one simultaneous 95% intervals; same imputations used for shrinkage and non-shrinkage"
)

# ------------------------ WIDTH TREND GRAPH ----------------------------------
width_trend <- application_results |>
  filter(method %in% names(method_map)) |>
  group_by(missing_rate, method, variance_method) |>
  summarise(
    avg_width = mean(width),
    mean_shrink_alpha = mean(shrink_alpha, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  width_trend,
  file.path(output_dir, "IBScovars_ci_width_trend.csv"),
  row.names = FALSE
)

p_width <- ggplot(
  width_trend,
  aes(
    x = 100 * missing_rate,
    y = avg_width,
    linetype = variance_method,
    shape = variance_method
  )
) +
  geom_line() +
  geom_point(size = 2) +
  facet_wrap(~ method) +
  labs(
    title = "IBScovars: simultaneous CI width as missingness increases",
    x = "Missing rate (%)",
    y = "Average simultaneous CI width",
    linetype = "Variance method",
    shape = "Variance method"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(
  file.path(output_dir, "IBScovars_ci_width_vs_missing_rate.png"),
  p_width, width = 10, height = 6, dpi = 300
)

# ------------------------ SHRINKAGE GAIN TABLE -------------------------------
shrink_gain <- application_results |>
  filter(method %in% names(method_map)) |>
  select(
    missing_rate, method, contrast,
    variance_method, width
  ) |>
  tidyr::pivot_wider(
    names_from = variance_method,
    values_from = width
  ) |>
  mutate(
    width_reduction_pct =
      100 * (`Non-shrinkage` - Shrinkage) / `Non-shrinkage`
  )

write.csv(
  shrink_gain,
  file.path(output_dir, "IBScovars_shrinkage_width_gain.csv"),
  row.names = FALSE
)

# ---------------------- REPEATED VALIDATION ---------------------------------
# This is NOT population coverage. It is a recovery experiment relative to the
# actual complete IBScovars sample contrasts, which are known before deletion.
if (run_repeated_validation) {
  raw_validation <- vector("list",
                           length(missing_rates) *
                             validation_reps *
                             length(method_map))
  cc <- 1L

  for (r_idx in seq_along(missing_rates)) {
    mr <- missing_rates[r_idx]

    for (b in seq_len(validation_reps)) {
      dat_miss <- impose_mcar_within_arm(
        ibs,
        rate = mr,
        seed = master_seed +
          50000000L +
          r_idx * 100000L +
          b
      )

      for (m_idx in seq_along(method_map)) {
        method_label <- names(method_map)[m_idx]
        method_code <- unname(method_map[m_idx])

        res <- make_mi_results(
          dat_miss,
          method = method_code,
          method_label = method_label,
          m = m_imputations,
          control_level = control_level,
          outcome_covariates = outcome_covariates,
          propensity_covariates = propensity_covariates,
          pmm_donors = pmm_donors,
          propensity_ngroups = propensity_ngroups,
          scenario_label = paste0("Repeated MCAR ", fmt_pct(mr)),
          missing_rate = mr,
          seed_base = master_seed +
            100000000L +
            r_idx * 1000000L +
            b * 1000L +
            m_idx * 100L
        ) |>
          add_benchmark_metrics(truth) |>
          mutate(rep = b)

        raw_validation[[cc]] <- res
        cc <- cc + 1L
      }
    }
  }

  validation_raw <- bind_rows(raw_validation)

  validation_summary <- summarize_validation(validation_raw)

  write.csv(
    validation_summary,
    file.path(output_dir, "IBScovars_validation_summary.csv"),
    row.names = FALSE
  )

  # Saving raw validation data can be useful for a paper appendix.
  write.csv(
    validation_raw,
    file.path(output_dir, "IBScovars_validation_raw.csv"),
    row.names = FALSE
  )

  p_recovery <- ggplot(
    validation_summary,
    aes(
      x = 100 * missing_rate,
      y = benchmark_recovery_pct,
      linetype = variance_method,
      shape = variance_method
    )
  ) +
    geom_hline(yintercept = 95, linetype = 3) +
    geom_line() +
    geom_point(size = 2) +
    facet_wrap(~ method) +
    labs(
      title = "IBScovars: recovery of complete-data sample contrasts",
      subtitle = "Repeated random deletion; horizontal line marks 95%",
      x = "Missing rate (%)",
      y = "Benchmark recovery rate (%)",
      linetype = "Variance method",
      shape = "Variance method"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")

  ggsave(
    file.path(output_dir, "IBScovars_validation_recovery.png"),
    p_recovery, width = 10, height = 6, dpi = 300
  )
}

write_session_info(
  file.path(output_dir, "IBScovars_sessionInfo.txt")
)

message("IBScovars application finished.")
message("Main outputs are in: ", normalizePath(output_dir))
