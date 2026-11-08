# ============================================================================
# FINAL APPLICATION CODE 1: IBScovars
#
# PRIMARY RESULT:
#   Complete-data benchmark + ONE primary 30% covariate-driven MAR scenario.
#
# FINAL PAPER CONCEPT:
#   logit P(M=1 | X) = alpha_r + log(2)*X
#   SAME gamma in every dose arm.
#   No MCAR control-arm special case.
#   No -2/+2 differential-gamma scenario.
#
# X for this real dataset:
#   standardized gender indicator (the only fully observed auxiliary covariate).
#
# Appendix only:
#   10%, 20%, 40% MAR sensitivity analyses + group/missingness summaries.
# ============================================================================

rm(list = ls())


# ============================================================================
# COMMON APPLICATION FUNCTIONS
# Final paper concept retained:
#   * m = 40
#   * REG: proper arm-specific Bayesian regression Y ~ X
#   * PMM: arm-specific Type-1 PMM, same-arm donors, K = 5
#   * PROPENSITY: missing ~ group + X, 3 subclasses WITHIN arm,
#                 ABB donors restricted to same arm
#   * SHRINKAGE vs NON-SHRINKAGE use the SAME imputed datasets
#   * Rubin pooling + many-to-one multivariate-t simultaneous 95% CIs
# ============================================================================

check_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Install these packages first: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

safe_solve <- function(M) {
  tryCatch(
    solve(M),
    error = function(e) solve(M + diag(1e-8, nrow(M)))
  )
}

make_corr_safe <- function(V) {
  V <- as.matrix(V)
  V <- (V + t(V)) / 2

  ev <- tryCatch(
    eigen(V, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NA_real_
  )

  if (any(!is.finite(ev)) || min(ev, na.rm = TRUE) <= 1e-8) {
    min_ev <- suppressWarnings(min(ev, na.rm = TRUE))
    ridge <- if (is.finite(min_ev)) abs(min_ev) + 1e-8 else 1e-6
    diag(V) <- diag(V) + ridge
  }

  C <- cov2cor(V)
  C[!is.finite(C)] <- 0
  diag(C) <- 1
  C
}

safe_qmvt_crit <- function(corr, df, alpha = 0.05, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  df_use <- if (!is.finite(df) || df <= 0) 1e6 else max(1, floor(df))

  out <- mvtnorm::qmvt(
    p = 1 - alpha,
    tail = "both.tails",
    df = df_use,
    corr = corr,
    algorithm = mvtnorm::GenzBretz(
      maxpts = 50000L,
      abseps = 1e-6
    )
  )

  crit <- as.numeric(out$quantile)
  if (!is.finite(crit) || crit <= 0) stop("Invalid qmvt critical value.")
  crit
}

# ----------------------------------------------------------------------------
# Unequal-n shrinkage formula retained from the finalized Simulation II code.
# ----------------------------------------------------------------------------
estimate_alpha_hat <- function(s2_vec, n_vec) {
  s2_vec <- as.numeric(s2_vec)
  n_vec <- as.numeric(n_vec)
  G <- length(s2_vec)

  if (length(n_vec) != G) stop("s2_vec and n_vec lengths do not match.")
  if (any(!is.finite(s2_vec)) || any(!is.finite(n_vec)) || any(n_vec <= 1)) {
    stop("Invalid variance/sample-size input for shrinkage.")
  }

  target <- mean(s2_vec)
  num <- ((G - 1) / G) * sum(2 * s2_vec^2 / (n_vec - 1))
  den <- sum((s2_vec - target)^2)

  if (!is.finite(den) || den <= 1e-12) return(1)

  a_hat <- num / den
  min(max(a_hat, 0), 1)
}

shrink_variances <- function(s2_vec, n_vec) {
  alpha_hat <- estimate_alpha_hat(s2_vec, n_vec)
  target <- mean(s2_vec)
  s2_tilde <- alpha_hat * target + (1 - alpha_hat) * s2_vec

  list(
    s2 = setNames(as.numeric(s2_tilde), names(s2_vec)),
    alpha_hat = alpha_hat
  )
}

build_U_unequal <- function(s2_vec, n_vec, group_levels, control_level) {
  treat_levels <- setdiff(group_levels, control_level)
  treat_labels <- paste0(treat_levels, " vs ", control_level)
  k <- length(treat_levels)

  s20 <- s2_vec[control_level]
  n0 <- n_vec[control_level]

  U <- matrix(
    s20 / n0,
    nrow = k,
    ncol = k,
    dimnames = list(treat_labels, treat_labels)
  )

  diag(U) <- vapply(
    treat_levels,
    function(g) s20 / n0 + s2_vec[g] / n_vec[g],
    numeric(1)
  )

  U
}

# ----------------------------------------------------------------------------
# Proper Bayesian arm-specific regression Y ~ X
# ----------------------------------------------------------------------------
draw_arm_regression <- function(x_obs, y_obs) {
  n <- length(y_obs)
  X <- cbind(`(Intercept)` = 1, x = x_obs)
  p <- ncol(X)

  if (n <= p) {
    stop(
      "Too few observed outcomes for Bayesian regression: n_obs=",
      n, ", p=", p
    )
  }

  XtX <- crossprod(X)
  if (qr(XtX)$rank < p) stop("Arm-specific regression matrix is rank deficient.")

  XtX_inv <- solve(XtX)
  beta_hat <- drop(XtX_inv %*% crossprod(X, y_obs))
  resid <- y_obs - drop(X %*% beta_hat)
  df <- n - p
  sse <- sum(resid^2)

  if (!is.finite(sse) || sse <= 0 || df <= 0) {
    stop("Invalid arm-specific regression fit.")
  }

  sigma2_star <- sse / rchisq(1, df = df)
  V_beta <- sigma2_star * XtX_inv
  L <- t(chol(V_beta))
  beta_star <- beta_hat + drop(L %*% rnorm(p))

  list(
    beta_hat = beta_hat,
    beta_star = beta_star,
    sigma2_star = sigma2_star
  )
}

impute_regression_once <- function(dat, seed) {
  set.seed(seed)
  out <- dat

  for (g in levels(out$group)) {
    idx <- which(out$group == g)
    obs <- idx[!is.na(out$y[idx])]
    mis <- idx[ is.na(out$y[idx])]

    if (!length(mis)) next

    dr <- draw_arm_regression(
      x_obs = out$x[obs],
      y_obs = out$y[obs]
    )

    X_mis <- cbind(`(Intercept)` = 1, x = out$x[mis])
    mu_mis <- drop(X_mis %*% dr$beta_star)

    out$y[mis] <- rnorm(
      length(mis),
      mean = mu_mis,
      sd = sqrt(dr$sigma2_star)
    )
  }

  if (anyNA(out$y)) stop("Regression imputation left missing values.")
  out
}

# ----------------------------------------------------------------------------
# Type-1 PMM, same-arm donors, K = 5
# ----------------------------------------------------------------------------
impute_pmm_once <- function(dat, seed, donors = 5L) {
  set.seed(seed)
  out <- dat

  if (!anyNA(out$y)) return(out)

  for (g in levels(out$group)) {
    idx <- which(out$group == g)
    obs <- idx[!is.na(out$y[idx])]
    mis <- idx[ is.na(out$y[idx])]

    if (!length(mis)) next
    if (!length(obs)) stop("No PMM donors in arm ", g, ".")

    dr <- draw_arm_regression(
      x_obs = out$x[obs],
      y_obs = out$y[obs]
    )

    X_obs <- cbind(`(Intercept)` = 1, x = out$x[obs])
    X_mis <- cbind(`(Intercept)` = 1, x = out$x[mis])

    # Type-1 PMM:
    # observed predictions use beta_hat;
    # missing predictions use beta_star.
    pred_obs <- drop(X_obs %*% dr$beta_hat)
    pred_mis <- drop(X_mis %*% dr$beta_star)

    k_use <- min(as.integer(donors), length(obs))

    for (j in seq_along(mis)) {
      d <- abs(pred_obs - pred_mis[j])
      ord <- order(d, runif(length(d)))  # random tie breaking
      donor_pool <- obs[ord[seq_len(k_use)]]
      donor_row <- sample(donor_pool, 1L)
      out$y[mis[j]] <- out$y[donor_row]
    }
  }

  if (anyNA(out$y)) stop("PMM left missing values.")
  out
}

# ----------------------------------------------------------------------------
# Propensity-score subclassification + Approximate Bayesian Bootstrap
# Final concept: missing ~ group + x; 3 subclasses within arm.
# ----------------------------------------------------------------------------
approx_bayesian_bootstrap <- function(y_obs, n_missing) {
  if (n_missing == 0L) return(numeric(0))
  if (!length(y_obs)) stop("No ABB donors.")

  y_star <- sample(y_obs, size = length(y_obs), replace = TRUE)
  sample(y_star, size = n_missing, replace = TRUE)
}

make_propensity_strata <- function(ps, n_groups = 3L) {
  if (any(!is.finite(ps))) stop("Non-finite propensity score.")

  unique_ps <- sort(unique(ps))
  if (length(unique_ps) == 1L) return(rep.int(1L, length(ps)))

  G <- min(as.integer(n_groups), length(unique_ps))
  score_rank <- rank(unique_ps, ties.method = "average")
  score_group <- ceiling(G * score_rank / length(unique_ps))
  score_group <- pmin(G, pmax(1L, as.integer(score_group)))

  mp <- setNames(
    score_group,
    format(unique_ps, digits = 17, scientific = FALSE)
  )

  as.integer(mp[format(ps, digits = 17, scientific = FALSE)])
}

impute_propensity_once <- function(dat, seed, n_groups = 3L) {
  set.seed(seed)
  out <- dat
  mis <- is.na(out$y)

  if (!any(mis)) return(out)

  psdat <- data.frame(
    missing = as.integer(mis),
    group = out$group,
    x = out$x
  )

  fit <- suppressWarnings(
    glm(
      missing ~ group + x,
      family = binomial(),
      data = psdat
    )
  )

  ps <- as.numeric(predict(fit, newdata = psdat, type = "response"))
  if (any(!is.finite(ps))) stop("Invalid propensity scores.")

  for (g in levels(out$group)) {
    idx_g <- which(out$group == g)
    mis_g <- idx_g[is.na(out$y[idx_g])]
    obs_g <- idx_g[!is.na(out$y[idx_g])]

    if (!length(mis_g)) next
    if (!length(obs_g)) stop("No observed propensity donors in arm ", g, ".")

    strata_g <- make_propensity_strata(ps[idx_g], n_groups = n_groups)
    strata_ids <- sort(unique(strata_g))

    centers <- vapply(
      strata_ids,
      function(s) mean(ps[idx_g[strata_g == s]]),
      numeric(1)
    )
    names(centers) <- as.character(strata_ids)

    for (s in strata_ids) {
      ids_s <- idx_g[strata_g == s]
      mis_s <- ids_s[is.na(out$y[ids_s])]

      if (!length(mis_s)) next

      obs_s <- ids_s[!is.na(out$y[ids_s])]

      # Only if a subclass is empty, borrow from nearest nonempty subclass
      # within the SAME arm.
      if (!length(obs_s)) {
        donor_strata <- strata_ids[
          vapply(
            strata_ids,
            function(ss) {
              ids_ss <- idx_g[strata_g == ss]
              any(!is.na(out$y[ids_ss]))
            },
            logical(1)
          )
        ]

        if (!length(donor_strata)) {
          stop("No propensity donor subclass in arm ", g, ".")
        }

        nearest <- donor_strata[
          which.min(
            abs(
              centers[as.character(donor_strata)] -
                centers[as.character(s)]
            )
          )
        ]

        donor_ids <- idx_g[strata_g == nearest]
        obs_s <- donor_ids[!is.na(out$y[donor_ids])]
      }

      out$y[mis_s] <- approx_bayesian_bootstrap(
        out$y[obs_s],
        length(mis_s)
      )
    }
  }

  if (anyNA(out$y)) stop("Propensity imputation left missing values.")
  out
}

make_imputations <- function(dat, method, m = 40L, base_seed,
                             pmm_donors = 5L,
                             propensity_groups = 3L) {
  comps <- vector("list", m)

  for (mm in seq_len(m)) {
    ss <- base_seed + mm - 1L

    comps[[mm]] <- switch(
      method,
      regression = impute_regression_once(dat, seed = ss),
      pmm = impute_pmm_once(dat, seed = ss, donors = pmm_donors),
      propensity = impute_propensity_once(
        dat, seed = ss, n_groups = propensity_groups
      ),
      stop("Unknown imputation method.")
    )
  }

  comps
}

get_completed_stats <- function(comp_list) {
  group_levels <- levels(comp_list[[1]]$group)
  m <- length(comp_list)

  means <- matrix(
    NA_real_, nrow = m, ncol = length(group_levels),
    dimnames = list(NULL, group_levels)
  )
  vars <- means

  for (mm in seq_len(m)) {
    means[mm, ] <- tapply(
      comp_list[[mm]]$y,
      comp_list[[mm]]$group,
      mean
    )[group_levels]

    vars[mm, ] <- tapply(
      comp_list[[mm]]$y,
      comp_list[[mm]]$group,
      var
    )[group_levels]
  }

  n_vec <- table(comp_list[[1]]$group)
  n_vec <- as.numeric(n_vec[group_levels])
  names(n_vec) <- group_levels

  list(means = means, vars = vars, n = n_vec, m = m)
}

# ----------------------------------------------------------------------------
# Rubin + many-to-one pooling, BOTH variance versions from SAME imputations
# ----------------------------------------------------------------------------
pool_many_to_one_both <- function(comp_list, control_level,
                                  alpha = 0.05,
                                  qseed_shrink = NULL,
                                  qseed_nonshrink = NULL) {
  st <- get_completed_stats(comp_list)

  group_levels <- colnames(st$means)
  treat_levels <- setdiff(group_levels, control_level)
  treat_labels <- paste0(treat_levels, " vs ", control_level)
  k <- length(treat_levels)
  m <- st$m

  Qmat <- matrix(
    NA_real_, nrow = m, ncol = k,
    dimnames = list(NULL, treat_labels)
  )

  U_ns <- array(NA_real_, dim = c(k, k, m))
  U_sh <- array(NA_real_, dim = c(k, k, m))
  alpha_hat <- numeric(m)

  for (mm in seq_len(m)) {
    mu <- st$means[mm, ]
    s2 <- st$vars[mm, ]
    names(s2) <- group_levels

    Qmat[mm, ] <- mu[treat_levels] - mu[control_level]

    U_ns[, , mm] <- build_U_unequal(
      s2, st$n, group_levels, control_level
    )

    shr <- shrink_variances(s2, st$n)
    alpha_hat[mm] <- shr$alpha_hat

    U_sh[, , mm] <- build_U_unequal(
      shr$s2, st$n, group_levels, control_level
    )
  }

  Qbar <- colMeans(Qmat)
  B <- cov(Qmat)

  pool_one <- function(U_array, qseed, variance_method, mean_alpha) {
    Ubar <- apply(U_array, c(1, 2), mean)
    Tm <- Ubar + (1 + 1/m) * B
    Tm <- (Tm + t(Tm)) / 2

    if (max(abs(B)) <= 1e-14) {
      nu <- 1e6
    } else {
      rm <- (1 + 1/m) *
        sum(diag(B %*% safe_solve(Ubar))) / k

      if (!is.finite(rm) || rm < 0) {
        stop("Invalid relative increase in variance.")
      }

      nu <- if (rm <= 1e-14) 1e6 else (m - 1) * (1 + 1/rm)^2
    }

    crit <- safe_qmvt_crit(
      make_corr_safe(Tm),
      df = nu,
      alpha = alpha,
      seed = qseed
    )

    se <- sqrt(diag(Tm))
    lower <- Qbar - crit * se
    upper <- Qbar + crit * se

    tibble::tibble(
      contrast = treat_labels,
      estimate = as.numeric(Qbar),
      lower = as.numeric(lower),
      upper = as.numeric(upper),
      width = as.numeric(upper - lower),
      significant = lower > 0 | upper < 0,
      variance_method = variance_method,
      shrink_alpha = mean_alpha
    )
  }

  dplyr::bind_rows(
    pool_one(
      U_sh,
      qseed_shrink,
      "Shrinkage",
      mean(alpha_hat)
    ),
    pool_one(
      U_ns,
      qseed_nonshrink,
      "Non-shrinkage",
      NA_real_
    )
  )
}

# ----------------------------------------------------------------------------
# Complete-data / complete-case many-to-one analysis
# ----------------------------------------------------------------------------
analyze_observed_both <- function(dat, control_level,
                                  method_label,
                                  alpha = 0.05,
                                  seed_base = 1000L) {
  dat2 <- dat[!is.na(dat$y), , drop = FALSE]
  dat2$group <- droplevels(dat2$group)

  group_levels <- levels(dat2$group)
  treat_levels <- setdiff(group_levels, control_level)
  treat_labels <- paste0(treat_levels, " vs ", control_level)

  split_y <- split(dat2$y, dat2$group)
  means <- vapply(split_y, mean, numeric(1))[group_levels]
  s2 <- vapply(split_y, var, numeric(1))[group_levels]
  n_vec <- vapply(split_y, length, integer(1))[group_levels]

  Q <- means[treat_levels] - means[control_level]
  names(Q) <- treat_labels

  make_one <- function(use_shrinkage, variance_method, qseed) {
    s2_use <- s2
    a <- NA_real_

    if (use_shrinkage) {
      shr <- shrink_variances(s2, n_vec)
      s2_use <- shr$s2
      a <- shr$alpha_hat
    }

    U <- build_U_unequal(
      s2_use, n_vec, group_levels, control_level
    )

    df <- sum(n_vec) - length(group_levels)
    crit <- safe_qmvt_crit(
      make_corr_safe(U),
      df = df,
      alpha = alpha,
      seed = qseed
    )

    se <- sqrt(diag(U))
    lower <- Q - crit * se
    upper <- Q + crit * se

    tibble::tibble(
      method = method_label,
      contrast = treat_labels,
      estimate = as.numeric(Q),
      lower = as.numeric(lower),
      upper = as.numeric(upper),
      width = as.numeric(upper - lower),
      significant = lower > 0 | upper < 0,
      variance_method = variance_method,
      shrink_alpha = a
    )
  }

  dplyr::bind_rows(
    make_one(TRUE, "Shrinkage", seed_base + 1L),
    make_one(FALSE, "Non-shrinkage", seed_base + 2L)
  )
}

make_primary_table <- function(long_results) {
  sh <- long_results |>
    dplyr::filter(variance_method == "Shrinkage") |>
    dplyr::select(
      method, contrast,
      estimate,
      lower_sh = lower,
      upper_sh = upper,
      width_sh = width,
      sig_sh = significant
    )

  ns <- long_results |>
    dplyr::filter(variance_method == "Non-shrinkage") |>
    dplyr::select(
      method, contrast,
      lower_ns = lower,
      upper_ns = upper,
      width_ns = width,
      sig_ns = significant
    )

  dplyr::left_join(sh, ns, by = c("method", "contrast")) |>
    dplyr::mutate(
      `Shrinkage 95% CI` = sprintf("[%.3f, %.3f]", lower_sh, upper_sh),
      `Non-shrinkage 95% CI` = sprintf("[%.3f, %.3f]", lower_ns, upper_ns),
      `Conclusion (Sh / NS)` = paste0(
        ifelse(sig_sh, "Significant", "Not significant"),
        " / ",
        ifelse(sig_ns, "Significant", "Not significant")
      )
    ) |>
    dplyr::select(
      Method = method,
      Contrast = contrast,
      Estimate = estimate,
      `Shrinkage 95% CI`,
      `Non-shrinkage 95% CI`,
      `Width Sh` = width_sh,
      `Width NS` = width_ns,
      `Conclusion (Sh / NS)`
    )
}

save_primary_outputs <- function(long_results, primary_table,
                                 output_dir, prefix, title, subtitle) {
  write.csv(
    primary_table,
    file.path(output_dir, paste0(prefix, "_PRIMARY_RESULT_TABLE.csv")),
    row.names = FALSE
  )

  gt_tab <- primary_table |>
    gt::gt() |>
    gt::tab_header(
      title = title,
      subtitle = subtitle
    ) |>
    gt::fmt_number(
      columns = c("Estimate", "Width Sh", "Width NS"),
      decimals = 3
    ) |>
    gt::cols_align(
      align = "center",
      columns = gt::everything()
    )

  gt::gtsave(
    gt_tab,
    file.path(output_dir, paste0(prefix, "_PRIMARY_RESULT_TABLE.html"))
  )

  # Forest plot: this is the ONE main graph.
  p <- ggplot2::ggplot(
    long_results,
    ggplot2::aes(
      x = contrast,
      y = estimate,
      ymin = lower,
      ymax = upper,
      shape = variance_method
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = 2,
      linewidth = 0.5
    ) +
    ggplot2::geom_pointrange(
      position = ggplot2::position_dodge(width = 0.55)
    ) +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~ method, scales = "free_y") +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      y = "Treatment - control mean difference (simultaneous 95% CI)",
      x = NULL,
      shape = "Variance method"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )

  ggplot2::ggsave(
    file.path(output_dir, paste0(prefix, "_PRIMARY_FOREST_PLOT.png")),
    p,
    width = 11,
    height = 7,
    dpi = 300
  )

  print(gt_tab)
  print(p)
}



check_packages(
  c("DoseFinding", "mvtnorm", "dplyr", "tibble", "tidyr", "ggplot2", "gt")
)

suppressPackageStartupMessages({
  library(DoseFinding)
  library(mvtnorm)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(gt)
})

# ============================================================================
# 0. SETTINGS
# ============================================================================
master_seed <- 20261107L
set.seed(master_seed)

alpha <- 0.05
m_imputations <- 40L
pmm_donors <- 5L
propensity_groups <- 3L

# FINAL MAR slope retained from paper simulation.
mar_gamma <- log(2)

# Main application scenario shown in the primary table/figure.
primary_missing_rate <- 0.30

# Appendix only.
appendix_missing_rates <- c(0.10, 0.20, 0.40)

output_dir <- "IBScovars_FINAL_application"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# 1. LOAD AND PREPARE DATA
# ============================================================================
data(IBScovars, package = "DoseFinding")

gender_num <- if (is.factor(IBScovars$gender)) {
  as.numeric(IBScovars$gender)
} else {
  as.numeric(IBScovars$gender)
}

dat_full <- tibble(
  group = factor(
    IBScovars$dose,
    levels = sort(unique(IBScovars$dose))
  ),
  x = as.numeric(scale(gender_num)),
  y = as.numeric(IBScovars$resp)
)

control_level <- levels(dat_full$group)[1]  # dose 0

if (anyNA(dat_full$x)) stop("X contains missing values.")
if (anyNA(dat_full$y)) stop("IBScovars response should be complete.")

# ============================================================================
# 2. FINAL MAR MECHANISM FOR THIS APPLICATION
# ============================================================================
# In the simulation X~N(0,1), so alpha was calibrated under the standard normal.
# Here X is an observed standardized binary covariate. We retain gamma=log(2)
# exactly, but calibrate alpha on the empirical X values so the EXPECTED
# marginal missing rate is the requested rate.

calibrate_empirical_intercept <- function(x, target_rate, gamma = mar_gamma) {
  uniroot(
    function(a) mean(plogis(a + gamma*x)) - target_rate,
    interval = c(-30, 30),
    tol = 1e-12
  )$root
}

impose_final_mar <- function(dat, target_rate, seed) {
  set.seed(seed)
  out <- dat

  alpha_mar <- calibrate_empirical_intercept(
    out$x,
    target_rate = target_rate,
    gamma = mar_gamma
  )

  p_missing <- plogis(alpha_mar + mar_gamma*out$x)
  M <- rbinom(nrow(out), size = 1L, prob = p_missing)

  out$y[M == 1L] <- NA_real_

  attr(out, "alpha_mar") <- alpha_mar
  attr(out, "expected_missing_rate") <- mean(p_missing)
  out
}

# ============================================================================
# 3. COMPLETE-DATA BENCHMARK
# ============================================================================
complete_results <- analyze_observed_both(
  dat_full,
  control_level = control_level,
  method_label = "Complete data",
  alpha = alpha,
  seed_base = master_seed + 1000L
)

# ============================================================================
# 4. PRIMARY 30% MAR APPLICATION
# ============================================================================
dat_primary <- impose_final_mar(
  dat_full,
  target_rate = primary_missing_rate,
  seed = master_seed + 2000L
)

achieved_primary <- mean(is.na(dat_primary$y))

method_map <- c(
  "Regression MI" = "regression",
  "PMM MI" = "pmm",
  "Propensity MI" = "propensity"
)

primary_mi_results <- vector("list", length(method_map))

for (j in seq_along(method_map)) {
  method_label <- names(method_map)[j]
  method_code <- unname(method_map[j])

  comps <- make_imputations(
    dat_primary,
    method = method_code,
    m = m_imputations,
    base_seed = master_seed + 100000L*j,
    pmm_donors = pmm_donors,
    propensity_groups = propensity_groups
  )

  primary_mi_results[[j]] <- pool_many_to_one_both(
    comps,
    control_level = control_level,
    alpha = alpha,
    qseed_shrink = master_seed + 500000L + 1000L*j,
    qseed_nonshrink = master_seed + 600000L + 1000L*j
  ) |>
    mutate(method = method_label, .before = 1)
}

primary_long <- bind_rows(
  complete_results,
  bind_rows(primary_mi_results)
)

primary_table <- make_primary_table(primary_long)

# ============================================================================
# 5. PRIMARY RESULT FIRST
# ============================================================================
cat("\n\n============================================================\n")
cat("IBScovars - PRIMARY RESULT\n")
cat("Target MAR = 30%; achieved missing rate = ",
    sprintf("%.1f%%", 100*achieved_primary), "\n", sep = "")
cat("MAR slope gamma = log(2), SAME in all dose arms.\n")
cat("============================================================\n\n")

save_primary_outputs(
  long_results = primary_long,
  primary_table = primary_table,
  output_dir = output_dir,
  prefix = "IBScovars",
  title = "IBScovars: shrinkage vs non-shrinkage",
  subtitle = paste0(
    "Complete-data benchmark + 30% covariate-driven MAR (achieved ",
    sprintf("%.1f%%", 100*achieved_primary),
    "); m=40; REG, PMM, and propensity MI"
  )
)

# ============================================================================
# 6. APPENDIX ONLY - DO NOT CLUTTER PRIMARY RESULT
# ============================================================================

# A. Group summary for the original complete data.
appendix_group_summary <- dat_full |>
  group_by(group) |>
  summarise(
    n = n(),
    mean = mean(y),
    sd = sd(y),
    variance = var(y),
    .groups = "drop"
  )

write.csv(
  appendix_group_summary,
  file.path(output_dir, "APPENDIX_group_summary.csv"),
  row.names = FALSE
)

# B. Missingness by arm for primary 30% MAR realization.
appendix_primary_missing <- dat_primary |>
  group_by(group) |>
  summarise(
    n = n(),
    n_missing = sum(is.na(y)),
    missing_pct = 100*mean(is.na(y)),
    .groups = "drop"
  )

write.csv(
  appendix_primary_missing,
  file.path(output_dir, "APPENDIX_primary_missingness_by_arm.csv"),
  row.names = FALSE
)

# C. Sensitivity: 10%, 20%, and 40% MAR.
# Same finalized MAR slope gamma=log(2) in every arm.
appendix_sensitivity <- list()
aa <- 1L

for (r_idx in seq_along(appendix_missing_rates)) {
  mr <- appendix_missing_rates[r_idx]

  dat_miss <- impose_final_mar(
    dat_full,
    target_rate = mr,
    seed = master_seed + 1000000L + r_idx*10000L
  )

  for (j in seq_along(method_map)) {
    method_label <- names(method_map)[j]
    method_code <- unname(method_map[j])

    comps <- make_imputations(
      dat_miss,
      method = method_code,
      m = m_imputations,
      base_seed = master_seed +
        2000000L +
        r_idx*100000L +
        j*10000L,
      pmm_donors = pmm_donors,
      propensity_groups = propensity_groups
    )

    res <- pool_many_to_one_both(
      comps,
      control_level = control_level,
      alpha = alpha,
      qseed_shrink = master_seed +
        3000000L +
        r_idx*100000L +
        j*10000L,
      qseed_nonshrink = master_seed +
        4000000L +
        r_idx*100000L +
        j*10000L
    ) |>
      mutate(
        target_missing_rate = mr,
        achieved_missing_rate = mean(is.na(dat_miss$y)),
        method = method_label,
        .before = 1
      )

    appendix_sensitivity[[aa]] <- res
    aa <- aa + 1L
  }
}

appendix_sensitivity_df <- bind_rows(appendix_sensitivity)

write.csv(
  appendix_sensitivity_df,
  file.path(output_dir, "APPENDIX_10_20_40pct_MAR_sensitivity.csv"),
  row.names = FALSE
)

# D. Shrinkage alpha values, primary analysis.
appendix_alpha <- primary_long |>
  filter(variance_method == "Shrinkage") |>
  select(method, contrast, shrink_alpha)

write.csv(
  appendix_alpha,
  file.path(output_dir, "APPENDIX_shrinkage_alpha.csv"),
  row.names = FALSE
)

cat("\nMain files:\n")
cat("  IBScovars_PRIMARY_RESULT_TABLE.html\n")
cat("  IBScovars_PRIMARY_RESULT_TABLE.csv\n")
cat("  IBScovars_PRIMARY_FOREST_PLOT.png\n")
cat("\nAppendix files are saved separately in the same folder.\n")
cat("Output folder: ", normalizePath(output_dir), "\n", sep = "")
