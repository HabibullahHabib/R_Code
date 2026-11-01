# ============================================================================
# THESIS_HABIB - PAPER SIMULATION: 5 TREATMENTS + 1 CONTROL
# UNEQUAL VARIANCES + COVARIATE-DRIVEN MAR + SHRINKAGE VS NON-SHRINKAGE
#
# This is the primary paper-setting extension of the Jiang, Zhou, and Liang
# Simulation-I reconstruction.
#
# PAPER DESIGN RETAINED:
#   * T0 control + T1,...,T5 treatments
#   * n = 20 per arm
#   * variance heterogeneity is now studied at max/min ratios R = 2, 4, 9
#   * Scenario I: all means = 0
#   * Scenario II: means = (0, 1.5, 2.0, 2.5, 3.0, 3.5)
#   * nominal missing rates = 0%, 10%, 20%, 30%, 40%
#   * m = 40 imputations
#   * shrinkage and non-shrinkage many-to-one simultaneous CIs
#
# MAR COVARIATE EXTENSION:
#   * X ~ N(0,1), fully observed
#   * Corr(X,Y | treatment arm) = rho_xy = 0.50
#   * For each arm g:
#       Y = mu_g + beta_x[g] * X + epsilon_g
#       beta_x[g] = rho_xy * sqrt(sig2[g])
#       Var(epsilon_g) = sig2[g] * (1-rho_xy^2)
#     Thus the ORIGINAL marginal arm variance sig2[g] is preserved exactly.
#   * missingness now uses arm-specific gamma values; arm-specific intercepts
#     are calibrated so each arm retains the requested expected missing rate.
#
# IMPUTATION:
#   * REG: proper arm-specific Bayesian linear regression Y ~ X.
#          This respects both unequal arm variances and the MAR covariate.
#   * PMM: arm-specific Type-1 PMM using Y ~ X, same-arm donors, K=5.
#          Observed predictions use beta_hat; missing predictions use beta_star.
#          Exact ties are randomly broken (continuous X makes ties rare).
#   * PROPENSITY: logistic missing ~ group * X; propensity subclasses are made
#          WITHIN treatment arm; ABB donors are restricted to that arm.
#          Primary choice = 3 subclasses because n=20 per arm is too sparse for
#          5 subclasses at 30-40% missingness.
#
# IMPORTANT:
#   The goal is a statistically defensible primary analysis, not parameter
#   tuning to force coverage to a desired number. Diagnostic outputs are saved
#   separately so bias and SE calibration can be checked without changing the
#   primary result tables.
# ============================================================================

suppressPackageStartupMessages({
  library(mvtnorm)
  library(dplyr)
  library(tibble)
  library(gt)
  library(parallel)
})

# =========================================================
# SETTINGS
# =========================================================

master_seed <- 20261002
set.seed(master_seed)

# Final publication run. Set quick_test=TRUE first to verify the code.
n_sim <- 10000
quick_test <- FALSE
if (quick_test) n_sim <- 50

n_treatments <- 5L
n_each <- 20L
m_imputations <- 40L
alpha <- 0.05
missing_rates <- c(0, 0.10, 0.20, 0.30, 0.40)

# Parallel settings.
use_parallel <- TRUE
n_cores <- max(1, parallel::detectCores(logical = TRUE) - 1)

# Numerical settings for multivariate-t critical values.
qmvt_maxpts <- 50000
qmvt_abseps <- 1e-6

# Treatment/control labels.
group_levels <- paste0("T", 0:n_treatments)
treat_groups <- group_levels[-1]
treat_labels <- paste0(treat_groups, "vs.T0")

# Equal sample size per arm.
N_vec <- setNames(rep(n_each, length(group_levels)), group_levels)

# Mean scenarios.
scenario1_label <- sprintf(
  "Scenario I: all means are 0 for T0,T1,...,T%d",
  n_treatments
)

scenario2_last_mean <- 1.5 + 0.5 * (n_treatments - 1L)
scenario2_label <- sprintf(
  "Scenario II: mu0=0, mu1=1.5, mu2=2.0, ..., mu%d=%.1f",
  n_treatments,
  scenario2_last_mean
)

scenarios <- setNames(
  list(
    setNames(rep(0, length(group_levels)), group_levels),
    setNames(
      c(0, seq(1.5, by = 0.5, length.out = n_treatments)),
      group_levels
    )
  ),
  c(scenario1_label, scenario2_label)
)

rho_xy <- 0.50

# ---------------- Professor-feedback simulation settings -------------------
# Stronger heterogeneity: max/min variance ratios R = 2, 4, 9.
# The minimum variance remains 9. T0 has the maximum variance; treatment
# variances increase linearly from 9 toward the control variance.
variance_min <- 9.0
variance_ratios <- c(2, 4, 9)

make_variance_vector <- function(R) {
  if (!is.finite(R) || R < 1) {
    stop("Variance ratio R must be finite and >= 1.")
  }

  max_var <- variance_min * R
  treatment_var <- variance_min +
    (0:(n_treatments - 1L)) / n_treatments *
      (max_var - variance_min)

  setNames(
    c(max_var, treatment_var),
    group_levels
  )
}

# Arm-specific missingness mechanisms requested by the professor.
mar_mechanisms <- c(
  "control_MCAR_treatment_MAR",
  "opposite_sign_MAR"
)

make_mar_gamma_vector <- function(mechanism) {
  out <- setNames(rep(NA_real_, length(group_levels)), group_levels)

  if (mechanism == "control_MCAR_treatment_MAR") {
    out["T0"] <- 0
    out[treat_groups] <- 2
  } else if (mechanism == "opposite_sign_MAR") {
    out["T0"] <- -2
    out[treat_groups] <- 2
  } else {
    stop("Unknown MAR mechanism.")
  }

  out
}

# Active values; refreshed for each condition in the main loop.
sig2 <- make_variance_vector(variance_ratios[1])
mar_gamma_by_group <- make_mar_gamma_vector(mar_mechanisms[1])

if (!is.finite(rho_xy) || abs(rho_xy) >= 1) {
  stop("rho_xy must be strictly between -1 and 1.")
}

beta_x_by_group <- rho_xy * sqrt(sig2)
sigma2_error_by_group <- sig2 * (1 - rho_xy^2)

if (any(!is.finite(sigma2_error_by_group)) ||
    any(sigma2_error_by_group <= 0)) {
  stop("Invalid conditional residual variances.")
}

# ----------------------- Imputation settings --------------------------------
# PMM: K=5 is retained as the primary setting. With n=20, K=10 would be too
# large at 40% missingness (about 12 observed subjects remain per arm).
pmm_donors <- 5L

# Propensity: use 3 subclasses rather than 5 in the PRIMARY n=20 design.
# At 40% missingness, about 12 observed subjects remain per arm, giving about
# 4 observed donors per subclass on average instead of about 2.4 with 5 groups.
propensity_ngroups <- 3L

# Imputation methods in the primary tables.
imp_map <- c(
  "Imputation using linear regression" = "regression",
  "Imputation using predictive mean matching method" = "pmm",
  "Imputation using propensity-score method" = "propensity"
)

fmt_mr <- function(mr) paste0(round(100 * mr), "%")

# ----------------------- MAR calibration ------------------------------------
mean_mar_probability <- function(alpha0, gamma) {
  integrate(
    function(z) plogis(alpha0 + gamma * z) * dnorm(z),
    lower = -Inf,
    upper = Inf,
    rel.tol = 1e-10,
    subdivisions = 1000L
  )$value
}

calibrate_mar_intercept <- function(target_rate, gamma) {
  if (target_rate <= 0 || target_rate >= 1) {
    stop("target_rate must lie strictly between 0 and 1.")
  }

  if (abs(gamma) <= 1e-14) {
    return(qlogis(target_rate))
  }

  uniroot(
    function(a) mean_mar_probability(a, gamma) - target_rate,
    interval = c(-30, 30),
    tol = 1e-12
  )$root
}

build_mar_intercepts <- function(gamma_vec) {
  out <- lapply(
    missing_rates[missing_rates > 0],
    function(r) {
      setNames(
        vapply(
          group_levels,
          function(g) calibrate_mar_intercept(r, gamma_vec[g]),
          numeric(1)
        ),
        group_levels
      )
    }
  )

  names(out) <- fmt_mr(missing_rates[missing_rates > 0])
  out
}

build_mar_calibration <- function(gamma_vec, intercept_list) {
  bind_rows(
    lapply(
      missing_rates[missing_rates > 0],
      function(r) {
        a_vec <- intercept_list[[fmt_mr(r)]]

        tibble(
          missing_rate = r,
          target_percent = 100 * r,
          group = group_levels,
          alpha_mar = as.numeric(a_vec[group_levels]),
          gamma_mar = as.numeric(gamma_vec[group_levels]),
          OR_per_1SD_X = exp(as.numeric(gamma_vec[group_levels])),
          expected_rate = vapply(
            group_levels,
            function(g) mean_mar_probability(a_vec[g], gamma_vec[g]),
            numeric(1)
          )
        )
      }
    )
  )
}

mar_intercepts <- build_mar_intercepts(mar_gamma_by_group)
mar_calibration <- build_mar_calibration(
  mar_gamma_by_group,
  mar_intercepts
)

# =========================================================
# BASIC HELPERS
# =========================================================

build_Uy <- function(s2_vec, n_vec) {
  # Covariance matrix for Q = (mean(T1)-mean(T0), ..., mean(T5)-mean(T0)).
  # Off-diagonal covariance is Var(mean(T0)).
  k <- length(s2_vec) - 1
  U <- matrix(s2_vec[1] / n_vec[1], nrow = k, ncol = k)
  diag(U) <- diag(U) + s2_vec[-1] / n_vec[-1]
  rownames(U) <- colnames(U) <- treat_labels
  U
}

estimate_alpha_hat <- function(s2_vec, n_each) {
  k <- length(s2_vec)
  num <- ((k - 1) / k) * sum(2 * s2_vec^2 / (n_each - 1))
  den <- sum((s2_vec - mean(s2_vec))^2)
  
  if (!is.finite(den) || den <= 1e-12) return(0)
  
  a_hat <- num / den
  min(max(a_hat, 0), 1)
}

shrink_hat <- function(s2_vec, n_each) {
  alpha_hat <- estimate_alpha_hat(s2_vec, n_each = n_each)
  target <- mean(s2_vec)
  s2_tilde <- alpha_hat * target + (1 - alpha_hat) * s2_vec
  
  list(s2 = s2_tilde, alpha_hat = alpha_hat)
}

pooled_hat <- function(s2_vec, n_vec) {
  df_vec <- pmax(as.numeric(n_vec) - 1, 0)
  sp2 <- sum(df_vec * as.numeric(s2_vec)) / sum(df_vec)

  if (!is.finite(sp2) || sp2 <= 0) {
    stop("Invalid pooled variance estimate.")
  }

  setNames(
    rep(sp2, length(s2_vec)),
    names(s2_vec)
  )
}

make_complete_data <- function(mu_vec, data_seed) {
  set.seed(data_seed)

  group <- factor(
    rep(group_levels, times = N_vec),
    levels = group_levels
  )

  # Fully observed baseline covariate.
  x <- rnorm(length(group), mean = 0, sd = 1)

  g_char <- as.character(group)
  mu_row <- unname(mu_vec[g_char])
  beta_row <- unname(beta_x_by_group[g_char])
  resid_sd_row <- sqrt(unname(sigma2_error_by_group[g_char]))

  # This construction preserves every requested MARGINAL arm variance.
  y <- mu_row +
    beta_row * x +
    rnorm(length(group), mean = 0, sd = resid_sd_row)

  tibble(
    group = group,
    x = x,
    y = y
  )
}

# Covariate-driven Bernoulli MAR.
# Missingness depends on fully observed X, not directly on the unobserved Y.
# Because the intercept is calibrated in expectation, an individual replicate
# need not have exactly 10/20/30/40% missing, but averages converge to target.
impose_missingness_mar <- function(dat, mr, miss_seed) {
  if (mr == 0) return(dat)

  set.seed(miss_seed)
  out <- dat

  alpha_vec <- mar_intercepts[[fmt_mr(mr)]]
  g_char <- as.character(out$group)

  p_missing <- plogis(
    unname(alpha_vec[g_char]) +
      unname(mar_gamma_by_group[g_char]) * out$x
  )

  if (any(!is.finite(p_missing)) ||
      any(p_missing < 0) ||
      any(p_missing > 1)) {
    stop("Invalid arm-specific MAR probabilities.")
  }

  missing <- rbinom(
    n = nrow(out),
    size = 1L,
    prob = p_missing
  ) == 1L

  out$y[missing] <- NA_real_

  out
}

# =========================================================
# IMPUTATION ENGINES
# =========================================================

# ----- 1) Proper arm-specific Bayesian regression Y ~ X ---------------------
# The paper design is heteroscedastic. Therefore each arm is allowed its own
# intercept, X slope, and residual variance. This also matches the data-generating
# structure, where the X slope is arm-specific to preserve a common rho_xy while
# retaining the requested marginal variance in each arm.

draw_arm_regression <- function(x_obs, y_obs) {
  n <- length(y_obs)
  X <- cbind(`(Intercept)` = 1, x = x_obs)
  p <- ncol(X)

  if (n <= p + 1L) {
    stop("Too few observed outcomes in an arm for Bayesian regression imputation.")
  }

  XtX <- crossprod(X)
  if (qr(XtX)$rank < p) {
    stop("Arm-specific regression design matrix is rank deficient.")
  }

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
    sigma2_star = sigma2_star,
    df = df
  )
}

impute_regression_once <- function(dat, seed) {
  set.seed(seed)
  out <- dat

  for (g in group_levels) {
    idx <- which(out$group == g)
    obs <- idx[!is.na(out$y[idx])]
    mis <- idx[ is.na(out$y[idx])]

    if (!length(mis)) next

    dr <- draw_arm_regression(
      x_obs = out$x[obs],
      y_obs = out$y[obs]
    )

    X_mis <- cbind(
      `(Intercept)` = 1,
      x = out$x[mis]
    )

    mu_mis <- drop(X_mis %*% dr$beta_star)

    out$y[mis] <- rnorm(
      length(mis),
      mean = mu_mis,
      sd = sqrt(dr$sigma2_star)
    )
  }

  if (anyNA(out$y)) {
    stop("Regression imputation left missing values.")
  }

  out
}

# ----- 2) PMM: arm-specific Type-1 predictive mean matching -----------------
# For each arm separately:
#   observed predictive means = X_obs * beta_hat
#   missing predictive means  = X_mis * beta_star
#   donor search is restricted to the same arm
#   K = pmm_donors
# Continuous X makes exact ties uncommon. If a tie occurs, random ordering is
# used rather than row number so dataset ordering cannot determine donors.

impute_pmm_once <- function(dat, seed, donors = pmm_donors) {
  set.seed(seed)
  out <- dat

  if (!anyNA(out$y)) return(out)

  for (g in group_levels) {
    idx <- which(out$group == g)
    obs <- idx[!is.na(out$y[idx])]
    mis <- idx[ is.na(out$y[idx])]

    if (!length(mis)) next
    if (!length(obs)) {
      stop(paste0("No observed PMM donors in arm ", g, "."))
    }

    dr <- draw_arm_regression(
      x_obs = out$x[obs],
      y_obs = out$y[obs]
    )

    X_obs <- cbind(
      `(Intercept)` = 1,
      x = out$x[obs]
    )

    X_mis <- cbind(
      `(Intercept)` = 1,
      x = out$x[mis]
    )

    pred_obs <- drop(X_obs %*% dr$beta_hat)
    pred_mis <- drop(X_mis %*% dr$beta_star)

    k_use <- min(as.integer(donors), length(obs))

    for (jj in seq_along(mis)) {
      d <- abs(pred_obs - pred_mis[jj])

      # Random secondary key affects only exact/near machine-level ties.
      ord <- order(d, runif(length(d)))

      donor_pool <- obs[ord[seq_len(k_use)]]
      donor_row <- sample(donor_pool, size = 1L)

      out$y[mis[jj]] <- out$y[donor_row]
    }
  }

  if (anyNA(out$y)) {
    stop("PMM imputation left missing values.")
  }

  out
}

# ----- 3) Propensity-score subclassification + ABB ---------------------------

approx_bayesian_bootstrap <- function(y_obs, n_missing) {
  if (n_missing == 0L) return(numeric(0))
  if (!length(y_obs)) {
    stop("No observed donor values available for ABB.")
  }

  # Approximate Bayesian bootstrap:
  # first bootstrap the empirical donor distribution, then draw imputations.
  y_star <- sample(
    y_obs,
    size = length(y_obs),
    replace = TRUE
  )

  sample(
    y_star,
    size = n_missing,
    replace = TRUE
  )
}

make_propensity_strata <- function(ps, n_groups = propensity_ngroups) {
  if (any(!is.finite(ps))) {
    stop("Non-finite propensity score.")
  }

  unique_ps <- sort(unique(ps))

  if (length(unique_ps) == 1L) {
    return(rep.int(1L, length(ps)))
  }

  G <- min(
    as.integer(n_groups),
    length(unique_ps)
  )

  # With continuous X the scores are generally unique. Ranking gives
  # approximately equal-size subclasses and preserves any exact ties.
  score_rank <- rank(unique_ps, ties.method = "average")

  score_group <- ceiling(
    G * score_rank / length(unique_ps)
  )

  score_group <- pmin(
    G,
    pmax(1L, as.integer(score_group))
  )

  mp <- setNames(
    score_group,
    format(unique_ps, digits = 17, scientific = FALSE)
  )

  as.integer(
    mp[
      format(ps, digits = 17, scientific = FALSE)
    ]
  )
}

impute_propensity_once <- function(dat, seed) {
  set.seed(seed)

  out <- dat
  mis <- is.na(out$y)
  propensity_merge_count <- 0L

  if (!any(mis)) {
    attr(out, "propensity_merge_count") <- propensity_merge_count
    return(out)
  }

  # Treatment is included in the propensity model, while treatment arm is also
  # retained as an exact donor stratum. The latter prevents cross-arm borrowing.
  psdat <- data.frame(
    missing = as.integer(mis),
    group = out$group,
    x = out$x
  )

  fit <- suppressWarnings(
    glm(
      missing ~ group * x,
      family = binomial(),
      data = psdat
    )
  )

  ps <- as.numeric(
    predict(
      fit,
      newdata = psdat,
      type = "response"
    )
  )

  if (any(!is.finite(ps))) {
    stop("Propensity model produced invalid fitted probabilities.")
  }

  for (g in group_levels) {
    idx_g <- which(out$group == g)
    mis_g <- idx_g[is.na(out$y[idx_g])]
    obs_g <- idx_g[!is.na(out$y[idx_g])]

    if (!length(mis_g)) next
    if (!length(obs_g)) {
      stop(paste0("Treatment arm ", g, " has no observed ABB donors."))
    }

    strata_g <- make_propensity_strata(
      ps[idx_g],
      n_groups = propensity_ngroups
    )

    strata_ids <- sort(unique(strata_g))

    stratum_center <- vapply(
      strata_ids,
      function(s) {
        mean(ps[idx_g[strata_g == s]])
      },
      numeric(1)
    )
    names(stratum_center) <- as.character(strata_ids)

    for (s in strata_ids) {
      ids_s <- idx_g[strata_g == s]
      mis_s <- ids_s[is.na(out$y[ids_s])]

      if (!length(mis_s)) next

      obs_s <- ids_s[!is.na(out$y[ids_s])]

      # Finite-sample safeguard. Only an EMPTY donor stratum is merged.
      # We do not merge merely because coverage would improve.
      if (!length(obs_s)) {
        propensity_merge_count <- propensity_merge_count + 1L

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
          stop(paste0("No observed propensity donor subclass in arm ", g, "."))
        }

        nearest_s <- donor_strata[
          which.min(
            abs(
              stratum_center[as.character(donor_strata)] -
                stratum_center[as.character(s)]
            )
          )
        ]

        donor_ids <- idx_g[strata_g == nearest_s]
        obs_s <- donor_ids[!is.na(out$y[donor_ids])]
      }

      out$y[mis_s] <- approx_bayesian_bootstrap(
        out$y[obs_s],
        length(mis_s)
      )
    }
  }

  if (anyNA(out$y)) {
    stop("Propensity-score imputation left missing values.")
  }

  attr(out, "propensity_merge_count") <- propensity_merge_count
  out
}

make_imputations <- function(dat, method, base_seed) {
  comps <- vector("list", m_imputations)

  for (mm in seq_len(m_imputations)) {
    ss <- base_seed + mm - 1L

    comps[[mm]] <- switch(
      method,
      regression = impute_regression_once(
        dat,
        seed = ss
      ),
      pmm = impute_pmm_once(
        dat,
        seed = ss,
        donors = pmm_donors
      ),
      propensity = impute_propensity_once(
        dat,
        seed = ss
      ),
      stop("Unknown imputation method.")
    )
  }

  comps
}

summarize_completions <- function(comps) {
  m <- length(comps)
  means_mat <- matrix(NA_real_, m, length(group_levels),
                      dimnames = list(NULL, group_levels))
  vars_mat <- matrix(NA_real_, m, length(group_levels),
                     dimnames = list(NULL, group_levels))
  for (ii in seq_len(m)) {
    means_mat[ii, ] <- tapply(comps[[ii]]$y, comps[[ii]]$group, mean)[group_levels]
    vars_mat[ii, ] <- tapply(comps[[ii]]$y, comps[[ii]]$group, var)[group_levels]
  }
  list(means = means_mat, vars = vars_mat, n = N_vec, m = m)
}

make_corr_safe <- function(V) {
  V <- as.matrix(V)
  V <- (V + t(V)) / 2
  
  ev <- tryCatch(eigen(V, symmetric = TRUE, only.values = TRUE)$values,
                 error = function(e) NA_real_)
  if (any(!is.finite(ev)) || min(ev, na.rm = TRUE) <= 1e-8) {
    ridge <- ifelse(any(is.finite(ev)), abs(min(ev, na.rm = TRUE)) + 1e-8, 1e-6)
    diag(V) <- diag(V) + ridge
  }
  
  corr <- cov2cor(V)
  corr[!is.finite(corr)] <- 0
  diag(corr) <- 1
  corr
}

qmvt_crit_fast <- function(corr, df, qseed = NULL) {
  if (!is.null(qseed)) set.seed(qseed)
  
  as.numeric(
    qmvt(
      1 - alpha,
      corr = corr,
      df = df,
      tail = "both.tails",
      algorithm = GenzBretz(maxpts = qmvt_maxpts, abseps = qmvt_abseps)
    )$quantile
  )
}

# =========================================================
# COMPLETE-DATA CIS
# =========================================================

complete_data_CIs_fast <- function(
  dat,
  alpha = 0.05,
  variance_mode = c("shrinkage", "nonshrinkage", "pooled"),
  qseed = NULL
) {
  variance_mode <- match.arg(variance_mode)

  ys <- lapply(
    group_levels,
    function(g) dat$y[dat$group == g]
  )
  names(ys) <- group_levels

  means <- sapply(ys, mean)
  s2_vec <- sapply(ys, var)
  n_vec <- sapply(ys, length)

  Q <- setNames(
    as.numeric(means[treat_groups] - means["T0"]),
    treat_labels
  )

  # None adjustment is intentionally unchanged across the three variance
  # benchmarks so the existing table format remains exactly the same.
  ci_none_list <- vector("list", length(treat_groups))
  names(ci_none_list) <- treat_labels

  none_se <- setNames(
    rep(NA_real_, length(treat_groups)),
    treat_labels
  )

  for (j in seq_along(treat_groups)) {
    g <- treat_groups[j]
    q <- Q[j]

    se2 <- s2_vec[g] / n_vec[g] +
      s2_vec["T0"] / n_vec["T0"]

    df <- se2^2 / (
      (s2_vec[g] / n_vec[g])^2 / (n_vec[g] - 1) +
        (s2_vec["T0"] / n_vec["T0"])^2 / (n_vec["T0"] - 1)
    )

    se_j <- sqrt(se2)
    crit <- qt(1 - alpha / 2, df = df)

    none_se[j] <- se_j
    ci_none_list[[j]] <- c(
      q - crit * se_j,
      q + crit * se_j
    )
  }

  shrink_alpha <- NA_real_

  if (variance_mode == "shrinkage") {
    shr <- shrink_hat(
      s2_vec,
      n_each = n_each
    )
    s2_for_U <- shr$s2
    shrink_alpha <- shr$alpha_hat
  } else if (variance_mode == "pooled") {
    s2_for_U <- pooled_hat(
      s2_vec,
      n_vec
    )
  } else {
    s2_for_U <- s2_vec
  }

  U <- build_Uy(
    s2_for_U,
    n_vec
  )

  df0 <- sum(n_vec) - length(n_vec)
  corr <- make_corr_safe(U)
  crit <- qmvt_crit_fast(
    corr = corr,
    df = df0,
    qseed = qseed
  )

  many_se <- sqrt(diag(U))

  ci_many_list <- vector("list", length(treat_groups))
  names(ci_many_list) <- treat_labels

  for (j in seq_along(treat_labels)) {
    ci_many_list[[j]] <- c(
      Q[j] - crit * many_se[j],
      Q[j] + crit * many_se[j]
    )
  }

  list(
    ci_none = ci_none_list,
    ci_many = ci_many_list,
    qhat = Q,
    none_se = none_se,
    many_se = many_se,
    shrink_alpha = shrink_alpha
  )
}

# =========================================================
# MI SUMMARY HELPERS
# =========================================================

# Completed-data summaries are produced by summarize_completions() above.

none_MI_CIs_fast <- function(summ, alpha = 0.05) {
  # HETEROSCEDASTIC pairwise MI intervals.
  #
  # Within-imputation variance for Tj-T0:
  #   U_m = s0_m^2 / n0 + sj_m^2 / nj
  #
  # This is consistent with the unequal-variance complete-data Welch analysis.
  m <- summ$m
  means_mat <- summ$means
  vars_mat <- summ$vars
  n_vec <- summ$n

  ci_none <- vector("list", length(treat_groups))
  names(ci_none) <- treat_labels

  qhat <- setNames(
    rep(NA_real_, length(treat_groups)),
    treat_labels
  )
  se <- qhat
  ubar_vec <- qhat
  b_vec <- qhat

  for (j in seq_along(treat_groups)) {
    g <- treat_groups[j]
    lbl <- treat_labels[j]

    Q <- as.numeric(
      means_mat[, g] -
        means_mat[, "T0"]
    )

    U <- as.numeric(
      vars_mat[, "T0"] / n_vec["T0"] +
        vars_mat[, g] / n_vec[g]
    )

    qbar <- mean(Q)
    ubar <- mean(U)
    b <- var(Q)
    tvar <- ubar + (1 + 1 / m) * b

    if (!is.finite(tvar) || tvar <= 0) {
      stop("Invalid scalar Rubin variance.")
    }

    if (b <= 1e-14) {
      nu <- 1e6
    } else {
      r <- (1 + 1 / m) * b / ubar
      nu <- if (
        !is.finite(r) || r <= 1e-14
      ) {
        1e6
      } else {
        (m - 1) * (1 + 1 / r)^2
      }
    }

    se_j <- sqrt(tvar)
    crit <- qt(1 - alpha / 2, df = nu)

    ci_none[[lbl]] <- c(
      qbar - crit * se_j,
      qbar + crit * se_j
    )

    qhat[lbl] <- qbar
    se[lbl] <- se_j
    ubar_vec[lbl] <- ubar
    b_vec[lbl] <- b
  }

  list(
    ci = ci_none,
    qhat = qhat,
    se = se,
    ubar = ubar_vec,
    b = b_vec
  )
}

many_to_one_MI_both_fast <- function(
  summ,
  alpha = 0.05,
  qseed_shrink = NULL,
  qseed_nonshrink = NULL,
  qseed_pooled = NULL
) {
  m <- summ$m
  means_mat <- summ$means
  vars_mat <- summ$vars
  n_vec <- summ$n

  k <- length(treat_labels)

  Qmat <- matrix(
    NA_real_,
    nrow = m,
    ncol = k,
    dimnames = list(NULL, treat_labels)
  )

  U_array_ns <- array(
    NA_real_,
    dim = c(k, k, m),
    dimnames = list(treat_labels, treat_labels, NULL)
  )

  U_array_sh <- array(
    NA_real_,
    dim = c(k, k, m),
    dimnames = list(treat_labels, treat_labels, NULL)
  )

  U_array_pool <- array(
    NA_real_,
    dim = c(k, k, m),
    dimnames = list(treat_labels, treat_labels, NULL)
  )

  shrink_alpha <- rep(NA_real_, m)

  for (ii in seq_len(m)) {
    means <- means_mat[ii, ]
    s2_vec <- vars_mat[ii, ]

    Qmat[ii, ] <- setNames(
      as.numeric(
        means[treat_groups] -
          means["T0"]
      ),
      treat_labels
    )

    U_array_ns[, , ii] <- build_Uy(
      s2_vec,
      n_vec
    )

    shr <- shrink_hat(
      s2_vec,
      n_each = n_each
    )

    shrink_alpha[ii] <- shr$alpha_hat

    U_array_sh[, , ii] <- build_Uy(
      shr$s2,
      n_vec
    )

    U_array_pool[, , ii] <- build_Uy(
      pooled_hat(s2_vec, n_vec),
      n_vec
    )
  }

  B <- cov(Qmat)
  Qm <- colMeans(Qmat)

  pool_many <- function(Ubar, qseed) {
    Tm <- Ubar +
      (1 + 1 / m) * B

    Tm <- (Tm + t(Tm)) / 2

    rm <- tryCatch(
      (1 + 1 / m) *
        (
          sum(
            diag(
              B %*% solve(Ubar)
            )
          ) / k
        ),
      error = function(e) NA_real_
    )

    if (!is.finite(rm) || rm <= 1e-12) {
      nu <- 1e6
    } else {
      nu <- max(
        1,
        floor(
          (m - 1) *
            (1 + 1 / rm)^2
        )
      )
    }

    corr <- make_corr_safe(Tm)

    crit <- qmvt_crit_fast(
      corr = corr,
      df = nu,
      qseed = qseed
    )

    se <- sqrt(diag(Tm))

    ci <- vector("list", k)
    names(ci) <- treat_labels

    for (j in seq_along(treat_labels)) {
      ci[[j]] <- c(
        Qm[j] - crit * se[j],
        Qm[j] + crit * se[j]
      )
    }

    list(
      ci = ci,
      se = se,
      Ubar = Ubar,
      B = B,
      Tm = Tm,
      df = nu
    )
  }

  Ubar_ns <- apply(
    U_array_ns,
    c(1, 2),
    mean
  )

  Ubar_sh <- apply(
    U_array_sh,
    c(1, 2),
    mean
  )

  Ubar_pool <- apply(
    U_array_pool,
    c(1, 2),
    mean
  )

  list(
    qhat = Qm,
    shrink = pool_many(
      Ubar_sh,
      qseed_shrink
    ),
    nonshrink = pool_many(
      Ubar_ns,
      qseed_nonshrink
    ),
    pooled = pool_many(
      Ubar_pool,
      qseed_pooled
    ),
    mean_shrink_alpha = mean(shrink_alpha)
  )
}

make_result_df <- function(
  scen_name,
  method_label,
  mr,
  sim,
  truth,
  ci_none,
  ci_many,
  qhat,
  none_se,
  many_se,
  achieved_missing_rate,
  achieved_missing_by_arm,
  x_mean_missing = NA_real_,
  x_mean_observed = NA_real_,
  propensity_merges = NA_real_,
  shrink_alpha = NA_real_
) {
  none_cov <- vapply(
    treat_labels,
    function(lbl) {
      ci <- ci_none[[lbl]]
      isTRUE(
        ci[1] <= truth[lbl] &&
          ci[2] >= truth[lbl]
      )
    },
    logical(1)
  )

  many_cov <- vapply(
    treat_labels,
    function(lbl) {
      ci <- ci_many[[lbl]]
      isTRUE(
        ci[1] <= truth[lbl] &&
          ci[2] >= truth[lbl]
      )
    },
    logical(1)
  )

  none_len <- vapply(
    treat_labels,
    function(lbl) diff(ci_none[[lbl]]),
    numeric(1)
  )

  many_len <- vapply(
    treat_labels,
    function(lbl) diff(ci_many[[lbl]]),
    numeric(1)
  )

  tibble(
    scenario = scen_name,
    method = method_label,
    mr = mr,
    sim = sim,
    comparison = treat_labels,

    truth = as.numeric(truth[treat_labels]),
    qhat = as.numeric(qhat[treat_labels]),

    none_se = as.numeric(none_se[treat_labels]),
    many_se = as.numeric(many_se[treat_labels]),

    achieved_missing_rate = achieved_missing_rate,
    achieved_mr_T0 = unname(achieved_missing_by_arm["T0"]),
    achieved_mr_treatment = as.numeric(
      achieved_missing_by_arm[treat_groups]
    ),

    x_mean_missing = x_mean_missing,
    x_mean_observed = x_mean_observed,

    propensity_merges = propensity_merges,
    shrink_alpha = shrink_alpha,

    none_covered = as.logical(none_cov),
    many_covered = as.logical(many_cov),
    none_length = as.numeric(none_len),
    many_length = as.numeric(many_len)
  )
}

make_table_block <- function(df_scen, scen_label) {
  scen_row <- df_scen[1, , drop = FALSE]
  scen_row[,] <- NA
  scen_row$MR <- scen_label
  bind_rows(scen_row, df_scen)
}


make_naive_bias_df <- function(
  dat,
  scen_name,
  mr,
  sim,
  truth
) {
  means_cc <- setNames(
    vapply(
      group_levels,
      function(g) {
        vals <- dat$y[dat$group == g]
        vals <- vals[!is.na(vals)]
        if (!length(vals)) NA_real_ else mean(vals)
      },
      numeric(1)
    ),
    group_levels
  )

  qhat_cc <- setNames(
    as.numeric(
      means_cc[treat_groups] -
        means_cc["T0"]
    ),
    treat_labels
  )

  tibble(
    scenario = scen_name,
    mr = mr,
    sim = sim,
    comparison = treat_labels,
    truth = as.numeric(truth[treat_labels]),
    qhat = as.numeric(qhat_cc[treat_labels]),
    bias = as.numeric(qhat_cc[treat_labels] - truth[treat_labels])
  )
}

# =========================================================
# ONE SIMULATION JOB
# =========================================================

run_one_job <- function(i) {
  job <- jobs[i, ]

  scen_idx <- job$scen_idx
  mr_idx <- job$mr_idx
  sim <- job$sim

  scen_name <- names(scenarios)[scen_idx]
  mu_vec <- scenarios[[scen_idx]]
  mr <- missing_rates[mr_idx]

  truth <- setNames(
    as.numeric(
      mu_vec[treat_groups] -
        mu_vec["T0"]
    ),
    treat_labels
  )

  data_seed <- master_seed +
    10000000 +
    scen_idx * 1000000 +
    mr_idx * 100000 +
    sim

  miss_seed <- master_seed +
    20000000 +
    scen_idx * 1000000 +
    mr_idx * 100000 +
    sim

  dat_full <- make_complete_data(
    mu_vec,
    data_seed = data_seed
  )

  out_sh <- list()
  out_ns <- list()
  out_pool <- list()

  if (mr == 0) {
    zero_by_arm <- setNames(
      rep(0, length(group_levels)),
      group_levels
    )

    cis_sh <- complete_data_CIs_fast(
      dat_full,
      alpha = alpha,
      variance_mode = "shrinkage",
      qseed = master_seed +
        30000000 +
        scen_idx * 1000000 +
        sim
    )

    cis_ns <- complete_data_CIs_fast(
      dat_full,
      alpha = alpha,
      variance_mode = "nonshrinkage",
      qseed = master_seed +
        40000000 +
        scen_idx * 1000000 +
        sim
    )

    cis_pool <- complete_data_CIs_fast(
      dat_full,
      alpha = alpha,
      variance_mode = "pooled",
      qseed = master_seed +
        45000000 +
        scen_idx * 1000000 +
        sim
    )

    common_args <- list(
      scen_name = scen_name,
      method_label = "Complete data",
      mr = mr,
      sim = sim,
      truth = truth,
      achieved_missing_rate = 0,
      achieved_missing_by_arm = zero_by_arm,
      x_mean_missing = NA_real_,
      x_mean_observed = mean(dat_full$x),
      propensity_merges = NA_real_
    )

    out_sh[[1]] <- do.call(
      make_result_df,
      c(
        common_args,
        list(
          ci_none = cis_sh$ci_none,
          ci_many = cis_sh$ci_many,
          qhat = cis_sh$qhat,
          none_se = cis_sh$none_se,
          many_se = cis_sh$many_se,
          shrink_alpha = cis_sh$shrink_alpha
        )
      )
    )

    out_ns[[1]] <- do.call(
      make_result_df,
      c(
        common_args,
        list(
          ci_none = cis_ns$ci_none,
          ci_many = cis_ns$ci_many,
          qhat = cis_ns$qhat,
          none_se = cis_ns$none_se,
          many_se = cis_ns$many_se,
          shrink_alpha = NA_real_
        )
      )
    )

    out_pool[[1]] <- do.call(
      make_result_df,
      c(
        common_args,
        list(
          ci_none = cis_pool$ci_none,
          ci_many = cis_pool$ci_many,
          qhat = cis_pool$qhat,
          none_se = cis_pool$none_se,
          many_se = cis_pool$many_se,
          shrink_alpha = NA_real_
        )
      )
    )

    naive <- make_naive_bias_df(
      dat_full,
      scen_name = scen_name,
      mr = mr,
      sim = sim,
      truth = truth
    )

    return(
      list(
        shrink = bind_rows(out_sh),
        nonshrink = bind_rows(out_ns),
        pooled = bind_rows(out_pool),
        naive = naive
      )
    )
  }

  dat_miss <- impose_missingness_mar(
    dat_full,
    mr = mr,
    miss_seed = miss_seed
  )

  naive <- make_naive_bias_df(
    dat_miss,
    scen_name = scen_name,
    mr = mr,
    sim = sim,
    truth = truth
  )

  mis_indicator <- is.na(dat_miss$y)
  achieved_missing_rate <- mean(mis_indicator)

  achieved_missing_by_arm <- setNames(
    vapply(
      group_levels,
      function(g) {
        idx_g <- which(dat_miss$group == g)
        mean(is.na(dat_miss$y[idx_g]))
      },
      numeric(1)
    ),
    group_levels
  )

  x_mean_missing <- if (any(mis_indicator)) {
    mean(dat_miss$x[mis_indicator])
  } else {
    NA_real_
  }

  x_mean_observed <- mean(
    dat_miss$x[!mis_indicator]
  )

  csh <- 1L
  cns <- 1L
  cpool <- 1L

  for (m_idx in seq_along(imp_map)) {
    method_label <- names(imp_map)[m_idx]
    imp_method <- unname(imp_map[m_idx])

    imp_seed <- master_seed +
      100000000 +
      scen_idx * 100000000 +
      mr_idx * 10000000 +
      m_idx * 2000000 +
      sim * 100

    comps <- make_imputations(
      dat_miss,
      method = imp_method,
      base_seed = imp_seed
    )

    summ <- summarize_completions(comps)

    none_res <- none_MI_CIs_fast(
      summ,
      alpha = alpha
    )

    many_both <- many_to_one_MI_both_fast(
      summ,
      alpha = alpha,
      qseed_shrink = master_seed +
        500000000 +
        scen_idx * 100000000 +
        mr_idx * 10000000 +
        m_idx * 2000000 +
        sim * 100,
      qseed_nonshrink = master_seed +
        900000000 +
        scen_idx * 100000000 +
        mr_idx * 10000000 +
        m_idx * 2000000 +
        sim * 100,
      qseed_pooled = master_seed +
        950000000 +
        scen_idx * 100000000 +
        mr_idx * 10000000 +
        m_idx * 2000000 +
        sim * 100
    )

    propensity_merges <- if (imp_method == "propensity") {
      mean(
        vapply(
          comps,
          function(z) {
            val <- attr(
              z,
              "propensity_merge_count"
            )
            if (is.null(val)) 0 else as.numeric(val)
          },
          numeric(1)
        )
      )
    } else {
      NA_real_
    }

    common_args <- list(
      scen_name = scen_name,
      method_label = method_label,
      mr = mr,
      sim = sim,
      truth = truth,
      ci_none = none_res$ci,
      qhat = none_res$qhat,
      none_se = none_res$se,
      achieved_missing_rate = achieved_missing_rate,
      achieved_missing_by_arm = achieved_missing_by_arm,
      x_mean_missing = x_mean_missing,
      x_mean_observed = x_mean_observed,
      propensity_merges = propensity_merges
    )

    out_sh[[csh]] <- do.call(
      make_result_df,
      c(
        common_args,
        list(
          ci_many = many_both$shrink$ci,
          many_se = many_both$shrink$se,
          shrink_alpha = many_both$mean_shrink_alpha
        )
      )
    )
    csh <- csh + 1L

    out_ns[[cns]] <- do.call(
      make_result_df,
      c(
        common_args,
        list(
          ci_many = many_both$nonshrink$ci,
          many_se = many_both$nonshrink$se,
          shrink_alpha = NA_real_
        )
      )
    )
    cns <- cns + 1L

    out_pool[[cpool]] <- do.call(
      make_result_df,
      c(
        common_args,
        list(
          ci_many = many_both$pooled$ci,
          many_se = many_both$pooled$se,
          shrink_alpha = NA_real_
        )
      )
    )
    cpool <- cpool + 1L
  }

  list(
    shrink = bind_rows(out_sh),
    nonshrink = bind_rows(out_ns),
    pooled = bind_rows(out_pool),
    naive = naive
  )
}

# =========================================================
# TABLE BUILDERS - SAME PRESENTATION AS BEFORE
# =========================================================

build_table2 <- function(raw, subtitle_text) {
  none_cols <- paste0(treat_labels, "_none")
  many_cols <- paste0(treat_labels, "_many")
  table_cols <- c("MR", none_cols, "Joint_none", many_cols, "Joint_many")
  
  make_blank_row <- function(label) {
    out <- as.list(rep(NA, length(table_cols)))
    names(out) <- table_cols
    out$MR <- label
    as_tibble(out)
  }
  
  cov_sum <- raw %>%
    group_by(scenario, method, mr, comparison) %>%
    summarise(
      none_cov = mean(none_covered, na.rm = TRUE) * 100,
      many_cov = mean(many_covered, na.rm = TRUE) * 100,
      .groups = "drop"
    )
  
  joint_sum <- raw %>%
    group_by(scenario, method, mr, sim) %>%
    summarise(
      joint_none = all(none_covered),
      joint_many = all(many_covered),
      .groups = "drop"
    ) %>%
    group_by(scenario, method, mr) %>%
    summarise(
      Joint_none = mean(joint_none, na.rm = TRUE) * 100,
      Joint_many = mean(joint_many, na.rm = TRUE) * 100,
      .groups = "drop"
    )
  
  make_cov_row <- function(scen_name, method_label, mr_value) {
    ind <- cov_sum %>%
      filter(scenario == scen_name, method == method_label, mr == mr_value)
    joint <- joint_sum %>%
      filter(scenario == scen_name, method == method_label, mr == mr_value)
    
    out <- as.list(rep(NA_real_, length(table_cols)))
    names(out) <- table_cols
    out$MR <- fmt_mr(mr_value)
    
    for (lbl in treat_labels) {
      tmp <- ind %>% filter(comparison == lbl)
      if (nrow(tmp) == 1) {
        out[[paste0(lbl, "_none")]] <- tmp$none_cov
        out[[paste0(lbl, "_many")]] <- tmp$many_cov
      }
    }
    
    if (nrow(joint) == 1) {
      out$Joint_none <- joint$Joint_none
      out$Joint_many <- joint$Joint_many
    }
    
    as_tibble(out)
  }
  
  table2_rows <- list()
  
  for (scen_name in names(scenarios)) {
    block <- make_cov_row(scen_name, "Complete data", 0)
    
    for (method_label in names(imp_map)) {
      block <- bind_rows(block, make_blank_row(method_label))
      
      for (mr_value in missing_rates[missing_rates > 0]) {
        block <- bind_rows(block, make_cov_row(scen_name, method_label, mr_value))
      }
    }
    
    block <- make_table_block(block, scen_name)
    table2_rows[[length(table2_rows) + 1]] <- block
  }
  
  table2_df <- bind_rows(table2_rows)
  
  col_labs <- c(
    setNames(treat_labels, none_cols),
    Joint_none = "Joint",
    setNames(treat_labels, many_cols),
    Joint_many = "Joint"
  )
  
  table2_df %>%
    gt() %>%
    tab_header(
      title = "Table 2. Coverage (%) of 95% CIs varying with missing rate.",
      subtitle = subtitle_text
    ) %>%
    tab_spanner(
      label = "None Adjustment",
      columns = all_of(c(none_cols, "Joint_none"))
    ) %>%
    tab_spanner(
      label = "Many-to-One Adjustment",
      columns = all_of(c(many_cols, "Joint_many"))
    ) %>%
    cols_label(MR = "MR", .list = col_labs) %>%
    fmt_number(columns = -MR, decimals = 2) %>%
    fmt_missing(everything(), missing_text = "") %>%
    cols_align(align = "center", columns = -MR) %>%
    cols_align(align = "left", columns = MR)
}

build_table3 <- function(raw, subtitle_text) {
  none_cols <- paste0("None_", treat_labels)
  many_cols <- paste0("Many_", treat_labels)
  table_cols <- c("MR", none_cols, many_cols)
  
  make_blank_row <- function(label) {
    out <- as.list(rep("", length(table_cols)))
    names(out) <- table_cols
    out$MR <- label
    as_tibble(out)
  }
  
  fmt_len_er <- function(len, er) {
    if (!is.finite(len)) return("")
    if (!is.finite(er)) return(sprintf("%.3f", len))
    sprintf("%.3f (%.2f)", len, er)
  }
  
  len_sum <- raw %>%
    group_by(scenario, method, mr, comparison) %>%
    summarise(
      None = mean(none_length, na.rm = TRUE),
      Many = mean(many_length, na.rm = TRUE),
      .groups = "drop"
    )
  
  baseline <- len_sum %>%
    filter(method == "Complete data", mr == 0) %>%
    select(
      scenario,
      comparison,
      base_None = None,
      base_Many = Many
    )
  
  len_sum2 <- len_sum %>%
    left_join(baseline, by = c("scenario", "comparison")) %>%
    mutate(
      ER_None = ifelse(mr == 0, NA_real_, (None / base_None - 1) * 100),
      ER_Many = ifelse(mr == 0, NA_real_, (Many / base_Many - 1) * 100)
    )
  
  make_len_row <- function(scen_name, method_label, mr_value) {
    ind <- len_sum2 %>%
      filter(scenario == scen_name, method == method_label, mr == mr_value)
    
    out <- as.list(rep("", length(table_cols)))
    names(out) <- table_cols
    out$MR <- fmt_mr(mr_value)
    
    for (lbl in treat_labels) {
      tmp <- ind %>% filter(comparison == lbl)
      if (nrow(tmp) == 1) {
        out[[paste0("None_", lbl)]] <- fmt_len_er(tmp$None, tmp$ER_None)
        out[[paste0("Many_", lbl)]] <- fmt_len_er(tmp$Many, tmp$ER_Many)
      }
    }
    
    as_tibble(out)
  }
  
  table3_rows <- list()
  
  for (scen_name in names(scenarios)) {
    block <- make_len_row(scen_name, "Complete data", 0)
    
    for (method_label in names(imp_map)) {
      block <- bind_rows(block, make_blank_row(method_label))
      
      for (mr_value in missing_rates[missing_rates > 0]) {
        block <- bind_rows(block, make_len_row(scen_name, method_label, mr_value))
      }
    }
    
    block <- make_table_block(block, scen_name)
    table3_rows[[length(table3_rows) + 1]] <- block
  }
  
  table3_df <- bind_rows(table3_rows)
  
  col_labs <- c(
    setNames(paste0(treat_labels, "  ER (%)"), none_cols),
    setNames(paste0(treat_labels, "  ER (%)"), many_cols)
  )
  
  table3_df %>%
    gt() %>%
    tab_header(
      title = "Table 3. Average lengths of 95% CIs varying with missing rate.",
      subtitle = subtitle_text
    ) %>%
    tab_spanner(label = "None", columns = all_of(none_cols)) %>%
    tab_spanner(label = "Many-to-One", columns = all_of(many_cols)) %>%
    cols_label(MR = "MR", .list = col_labs) %>%
    cols_align(align = "left", columns = MR) %>%
    cols_align(align = "center", columns = -MR) %>%
    tab_source_note(source_note = "NOTE: MR, missing rate; ER, expansion rate. Values are shown as average length (ER %).")
}



# =========================================================
# MAIN SIMULATION: PROFESSOR-FEEDBACK CONDITIONS
# =========================================================
# build_table2() and build_table3() above are unchanged. Therefore the result
# format is preserved exactly. Each variance-ratio x MAR condition is saved
# separately.

all_condition_summaries <- list()
all_bias_summaries <- list()
condition_counter <- 1L

sanitize_tag <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}

for (variance_ratio in variance_ratios) {
  for (mar_mechanism in mar_mechanisms) {

    # Activate the current heterogeneity condition.
    sig2 <- make_variance_vector(variance_ratio)
    beta_x_by_group <- rho_xy * sqrt(sig2)
    sigma2_error_by_group <- sig2 * (1 - rho_xy^2)

    if (any(!is.finite(sigma2_error_by_group)) ||
        any(sigma2_error_by_group <= 0)) {
      stop("Invalid conditional residual variances.")
    }

    # Activate the current arm-specific missingness condition.
    mar_gamma_by_group <- make_mar_gamma_vector(mar_mechanism)
    mar_intercepts <- build_mar_intercepts(mar_gamma_by_group)
    mar_calibration <- build_mar_calibration(
      mar_gamma_by_group,
      mar_intercepts
    )

    ratio_tag <- paste0(
      "R",
      gsub("\\.", "p", format(variance_ratio, trim = TRUE))
    )
    mar_tag <- sanitize_tag(mar_mechanism)
    condition_tag <- paste(ratio_tag, mar_tag, sep = "_")

    message("============================================================")
    message(
      sprintf(
        "Running %d-treatment condition: R=%s | MAR=%s",
        n_treatments,
        format(variance_ratio, trim = TRUE),
        mar_mechanism
      )
    )
    message("Marginal variances:")
    print(sig2)
    message("Arm-specific MAR slopes:")
    print(mar_gamma_by_group)
    message("MAR calibration:")
    print(mar_calibration)

    jobs <- expand.grid(
      scen_idx = seq_along(scenarios),
      mr_idx = seq_along(missing_rates),
      sim = seq_len(n_sim),
      KEEP.OUT.ATTRS = FALSE
    )

    message(
      sprintf(
        "PMM donors K = %d | propensity subclasses = %d",
        pmm_donors,
        propensity_ngroups
      )
    )
    message(sprintf("Total jobs in this condition: %d", nrow(jobs)))
    message(
      sprintf(
        "Parallel: %s | cores used: %d",
        use_parallel,
        ifelse(use_parallel, n_cores, 1)
      )
    )

    if (use_parallel && n_cores > 1) {
      cl <- makeCluster(n_cores)
      clusterEvalQ(cl, {
        library(mvtnorm)
        library(dplyr)
        library(tibble)
      })

      objects_to_export <- setdiff(ls(envir = .GlobalEnv), "cl")
      clusterExport(
        cl,
        varlist = objects_to_export,
        envir = .GlobalEnv
      )

      results <- parLapplyLB(
        cl,
        seq_len(nrow(jobs)),
        run_one_job
      )
      stopCluster(cl)
    } else {
      results <- lapply(
        seq_len(nrow(jobs)),
        run_one_job
      )
    }

    raw_shrink <- bind_rows(lapply(results, `[[`, "shrink"))
    raw_nonshrink <- bind_rows(lapply(results, `[[`, "nonshrink"))
    raw_pooled <- bind_rows(lapply(results, `[[`, "pooled"))
    raw_naive <- bind_rows(lapply(results, `[[`, "naive"))

    # None-adjustment results should be identical because shrinkage/pooling
    # only affect the many-to-one covariance construction.
    print(
      all.equal(
        raw_shrink %>%
          arrange(scenario, method, mr, sim, comparison) %>%
          select(
            scenario, method, mr, sim, comparison,
            none_covered, none_length
          ),
        raw_nonshrink %>%
          arrange(scenario, method, mr, sim, comparison) %>%
          select(
            scenario, method, mr, sim, comparison,
            none_covered, none_length
          ),
        check.attributes = FALSE
      )
    )

    print(
      all.equal(
        raw_shrink %>%
          arrange(scenario, method, mr, sim, comparison) %>%
          select(
            scenario, method, mr, sim, comparison,
            none_covered, none_length
          ),
        raw_pooled %>%
          arrange(scenario, method, mr, sim, comparison) %>%
          select(
            scenario, method, mr, sim, comparison,
            none_covered, none_length
          ),
        check.attributes = FALSE
      )
    )

    # ---------------------------------------------------------
    # Bias comparison requested by the professor
    # ---------------------------------------------------------
    naive_bias_summary <- raw_naive %>%
      group_by(
        scenario,
        mr,
        comparison
      ) %>%
      summarise(
        method = "Naive deletion",
        bias = mean(qhat - truth, na.rm = TRUE),
        empirical_SD = sd(qhat, na.rm = TRUE),
        RMSE = sqrt(
          mean(
            (qhat - truth)^2,
            na.rm = TRUE
          )
        ),
        .groups = "drop"
      )

    # MI point estimates do not depend on shrinkage/non-shrinkage/pooling.
    mi_bias_summary <- raw_shrink %>%
      group_by(
        scenario,
        method,
        mr,
        comparison
      ) %>%
      summarise(
        bias = mean(qhat - truth, na.rm = TRUE),
        empirical_SD = sd(qhat, na.rm = TRUE),
        RMSE = sqrt(
          mean(
            (qhat - truth)^2,
            na.rm = TRUE
          )
        ),
        .groups = "drop"
      )

    bias_comparison <- bind_rows(
      naive_bias_summary,
      mi_bias_summary
    ) %>%
      mutate(
        variance_ratio = variance_ratio,
        mar_mechanism = mar_mechanism,
        gamma_control = unname(mar_gamma_by_group["T0"]),
        gamma_treatment = unname(mar_gamma_by_group[treat_groups[1]])
      ) %>%
      select(
        variance_ratio,
        mar_mechanism,
        gamma_control,
        gamma_treatment,
        everything()
      )

    write.csv(
      bias_comparison,
      paste0(
        "bias_naive_vs_MI_",
        n_treatments,
        "_treatments_",
        condition_tag,
        ".csv"
      ),
      row.names = FALSE
    )

    all_bias_summaries[[condition_counter]] <- bias_comparison

    # ---------------------------------------------------------
    # Existing diagnostics, plus pooled benchmark
    # ---------------------------------------------------------
    diagnostic_long <- bind_rows(
      raw_shrink %>%
        mutate(variance_version = "shrinkage"),
      raw_nonshrink %>%
        mutate(variance_version = "non-shrinkage"),
      raw_pooled %>%
        mutate(variance_version = "pooled")
    )

    diagnostic_summary <- diagnostic_long %>%
      group_by(
        variance_version,
        scenario,
        method,
        mr,
        comparison
      ) %>%
      summarise(
        achieved_mr_pct =
          100 * mean(achieved_missing_rate, na.rm = TRUE),

        achieved_T0_pct =
          100 * mean(achieved_mr_T0, na.rm = TRUE),

        achieved_treatment_pct =
          100 * mean(achieved_mr_treatment, na.rm = TRUE),

        mean_x_missing =
          if (all(is.na(x_mean_missing))) {
            NA_real_
          } else {
            mean(x_mean_missing, na.rm = TRUE)
          },

        mean_x_observed =
          mean(x_mean_observed, na.rm = TRUE),

        bias =
          mean(qhat - truth, na.rm = TRUE),

        empirical_SD =
          sd(qhat, na.rm = TRUE),

        RMSE =
          sqrt(
            mean(
              (qhat - truth)^2,
              na.rm = TRUE
            )
          ),

        mean_none_SE =
          mean(none_se, na.rm = TRUE),

        none_SE_ratio =
          mean(none_se, na.rm = TRUE) /
            sd(qhat, na.rm = TRUE),

        mean_many_SE =
          mean(many_se, na.rm = TRUE),

        many_SE_ratio =
          mean(many_se, na.rm = TRUE) /
            sd(qhat, na.rm = TRUE),

        propensity_merge_replicate_pct =
          if (all(is.na(propensity_merges))) {
            NA_real_
          } else {
            100 *
              mean(
                propensity_merges > 0,
                na.rm = TRUE
              )
          },

        mean_propensity_merges =
          if (all(is.na(propensity_merges))) {
            NA_real_
          } else {
            mean(
              propensity_merges,
              na.rm = TRUE
            )
          },

        mean_shrink_alpha =
          if (all(is.na(shrink_alpha))) {
            NA_real_
          } else {
            mean(
              shrink_alpha,
              na.rm = TRUE
            )
          },

        .groups = "drop"
      ) %>%
      mutate(
        variance_ratio = variance_ratio,
        mar_mechanism = mar_mechanism
      )

    write.csv(
      diagnostic_summary,
      paste0(
        "paper_MAR_diagnostic_summary_",
        n_treatments,
        "_treatments_",
        condition_tag,
        ".csv"
      ),
      row.names = FALSE
    )

    # Raw outputs
    write.csv(
      raw_shrink,
      paste0(
        "raw_shrink_",
        n_treatments,
        "_treatments_",
        condition_tag,
        ".csv"
      ),
      row.names = FALSE
    )

    write.csv(
      raw_nonshrink,
      paste0(
        "raw_nonshrink_",
        n_treatments,
        "_treatments_",
        condition_tag,
        ".csv"
      ),
      row.names = FALSE
    )

    write.csv(
      raw_pooled,
      paste0(
        "raw_pooled_",
        n_treatments,
        "_treatments_",
        condition_tag,
        ".csv"
      ),
      row.names = FALSE
    )

    # ---------------------------------------------------------
    # Monte Carlo standard errors
    # ---------------------------------------------------------
    mcse_individual <- bind_rows(
      raw_shrink %>% mutate(variance_version = "shrink"),
      raw_nonshrink %>% mutate(variance_version = "nonshrink"),
      raw_pooled %>% mutate(variance_version = "pooled")
    ) %>%
      group_by(
        variance_version,
        scenario,
        method,
        mr,
        comparison
      ) %>%
      summarise(
        individual_coverage = mean(many_covered),
        individual_mcse = sqrt(
          individual_coverage *
            (1 - individual_coverage) / n()
        ),
        .groups = "drop"
      )

    mcse_joint <- bind_rows(
      raw_shrink %>% mutate(variance_version = "shrink"),
      raw_nonshrink %>% mutate(variance_version = "nonshrink"),
      raw_pooled %>% mutate(variance_version = "pooled")
    ) %>%
      group_by(
        variance_version,
        scenario,
        method,
        mr,
        sim
      ) %>%
      summarise(
        joint_covered = all(many_covered),
        .groups = "drop"
      ) %>%
      group_by(
        variance_version,
        scenario,
        method,
        mr
      ) %>%
      summarise(
        joint_coverage = mean(joint_covered),
        joint_mcse = sqrt(
          joint_coverage *
            (1 - joint_coverage) / n()
        ),
        .groups = "drop"
      )

    write.csv(
      mcse_individual,
      paste0(
        "coverage_mcse_individual_",
        n_treatments,
        "_treatments_",
        condition_tag,
        ".csv"
      ),
      row.names = FALSE
    )

    write.csv(
      mcse_joint,
      paste0(
        "coverage_mcse_joint_",
        n_treatments,
        "_treatments_",
        condition_tag,
        ".csv"
      ),
      row.names = FALSE
    )

    # ---------------------------------------------------------
    # PRIMARY RESULT TABLES -- FORMAT UNCHANGED
    # ---------------------------------------------------------
    subtitle_common <- sprintf(
      "Variance ratio R=%s; %s; %d treatments",
      format(variance_ratio, trim = TRUE),
      mar_mechanism,
      n_treatments
    )

    table2_shrink_gt <- build_table2(
      raw_shrink,
      paste0(
        "Unequal-variance shrinkage; ",
        subtitle_common
      )
    )

    table3_shrink_gt <- build_table3(
      raw_shrink,
      paste0(
        "Unequal-variance shrinkage; ",
        subtitle_common
      )
    )

    table2_nonshrink_gt <- build_table2(
      raw_nonshrink,
      paste0(
        "Unequal-variance non-shrinkage; ",
        subtitle_common
      )
    )

    table3_nonshrink_gt <- build_table3(
      raw_nonshrink,
      paste0(
        "Unequal-variance non-shrinkage; ",
        subtitle_common
      )
    )

    # Additional pooling benchmark, using the exact same table layout.
    table2_pooled_gt <- build_table2(
      raw_pooled,
      paste0(
        "Common-variance pooled benchmark; ",
        subtitle_common
      )
    )

    table3_pooled_gt <- build_table3(
      raw_pooled,
      paste0(
        "Common-variance pooled benchmark; ",
        subtitle_common
      )
    )

    print(table2_shrink_gt)
    print(table3_shrink_gt)
    print(table2_nonshrink_gt)
    print(table3_nonshrink_gt)
    print(table2_pooled_gt)
    print(table3_pooled_gt)

    gtsave(
      table2_shrink_gt,
      paste0(
        "table2_shrink_",
        n_treatments,
        "_treatments_",
        condition_tag,
        ".html"
      )
    )

    gtsave(
      table3_shrink_gt,
      paste0(
        "table3_shrink_",
        n_treatments,
        "_treatments_",
        condition_tag,
        ".html"
      )
    )

    gtsave(
      table2_nonshrink_gt,
      paste0(
        "table2_nonshrink_",
        n_treatments,
        "_treatments_",
        condition_tag,
        ".html"
      )
    )

    gtsave(
      table3_nonshrink_gt,
      paste0(
        "table3_nonshrink_",
        n_treatments,
        "_treatments_",
        condition_tag,
        ".html"
      )
    )

    gtsave(
      table2_pooled_gt,
      paste0(
        "table2_pooled_",
        n_treatments,
        "_treatments_",
        condition_tag,
        ".html"
      )
    )

    gtsave(
      table3_pooled_gt,
      paste0(
        "table3_pooled_",
        n_treatments,
        "_treatments_",
        condition_tag,
        ".html"
      )
    )

    # Compact cross-condition summary for shrinkage vs non-shrinkage vs pooling.
    joint_summary <- bind_rows(
      raw_shrink %>% mutate(variance_version = "shrinkage"),
      raw_nonshrink %>% mutate(variance_version = "non-shrinkage"),
      raw_pooled %>% mutate(variance_version = "pooled")
    ) %>%
      group_by(
        variance_version,
        scenario,
        method,
        mr,
        sim
      ) %>%
      summarise(
        joint_many = all(many_covered),
        .groups = "drop"
      ) %>%
      group_by(
        variance_version,
        scenario,
        method,
        mr
      ) %>%
      summarise(
        joint_coverage_pct =
          100 * mean(joint_many),
        .groups = "drop"
      ) %>%
      mutate(
        variance_ratio = variance_ratio,
        mar_mechanism = mar_mechanism
      )

    all_condition_summaries[[condition_counter]] <- joint_summary
    condition_counter <- condition_counter + 1L

    rm(
      results,
      raw_shrink,
      raw_nonshrink,
      raw_pooled,
      raw_naive
    )
    invisible(gc())
  }
}

combined_joint_summary <- bind_rows(
  all_condition_summaries
)

combined_bias_summary <- bind_rows(
  all_bias_summaries
)

write.csv(
  combined_joint_summary,
  paste0(
    "professor_feedback_joint_coverage_summary_",
    n_treatments,
    "_treatments.csv"
  ),
  row.names = FALSE
)

write.csv(
  combined_bias_summary,
  paste0(
    "professor_feedback_bias_summary_",
    n_treatments,
    "_treatments.csv"
  ),
  row.names = FALSE
)

message(
  "Done. Existing Table 2/Table 3 format was preserved for every condition. ",
  "Separate pooled benchmark and naive-vs-MI bias diagnostics were also saved."
)
