# =============================================================================
# LITTER DATASET APPLICATION FOR THE PAPER
# =============================================================================
# Dataset:
#   multcomp::litter
#
# Goal:
#   Evaluate the finalized heteroscedastic Dunnett-type procedure with and
#   without variance shrinkage on a real dose-response dataset that has
#   naturally heterogeneous group variances.
#
# Main application:
#   - Control dose = 0
#   - Active doses = 5, 50, 500
#   - Outcome = litter weight
#   - Fully observed MAR covariate X = standardized litter size ("number")
#   - Target expected missingness = 30%
#   - MAR: logit P(M=1|X) = alpha_30 + log(2) X
#   - m = 40 imputations
#   - Primary MI method = proper arm-specific Bayesian regression Y ~ X
#   - Same imputed datasets used for shrinkage and non-shrinkage
#
# Additional validation:
#   Repeatedly impose new 30% MAR patterns on the SAME complete dataset,
#   impute them, and compare shrinkage vs non-shrinkage against the fixed
#   complete-data treatment-control contrasts.
#
# IMPORTANT TERMINOLOGY:
#   The repeated-deletion validation uses the complete sample as an empirical
#   benchmark. Therefore the script calls the resulting quantity
#   "benchmark recovery", NOT population coverage.
#
# Outputs:
#   1. litter_OBSERVED_GROUP_SUMMARY.csv
#   2. litter_COMPLETE_DATA_BENCHMARK.csv
#   3. litter_PRIMARY_REG_APPLICATION_TABLE.csv
#   4. litter_ALL_MI_APPLICATION_RESULTS.csv
#   5. litter_PRIMARY_REFERENCE_STYLE_FOREST.pdf / .png
#   6. litter_PRIMARY_MISSINGNESS_BY_ARM.csv
#   7. litter_PRIMARY_SHRINKAGE_ALPHA.csv
#   8. litter_VALIDATION_SUMMARY.csv
#   9. litter_VALIDATION_BY_CONTRAST.csv
#  10. litter_VALIDATION_SHRINKAGE_MINUS_NONSHRINK.csv
#
# Publication run:
#   First run with publication_run <- FALSE.
#   After confirming that everything works, set publication_run <- TRUE.
# =============================================================================

rm(list = ls())

# =============================================================================
# 0. PACKAGES AND SETTINGS
# =============================================================================

needed <- c("multcomp", "mvtnorm", "dplyr", "ggplot2")
missing_pkg <- needed[
  !vapply(needed, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkg)) {
  stop(
    "Install these packages first: ",
    paste(missing_pkg, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(multcomp)
  library(mvtnorm)
  library(dplyr)
  library(ggplot2)
})

master_seed <- 20261112L
set.seed(master_seed)

alpha_level <- 0.05
target_missing_rate <- 0.30
mar_gamma <- log(2)

m_imputations <- 40L
pmm_donors <- 5L
propensity_groups <- 3L

# Primary application uses regression MI.
primary_mi_method <- "REG"

# Repeated-MAR validation.
run_validation <- TRUE

# FALSE = quick methodological check.
# TRUE  = final publication-grade run.
publication_run <- FALSE

n_validation <- if (publication_run) 10000L else 500L

# Numerical settings for mv-t critical values.
qmvt_maxpts <- 50000L
qmvt_abseps <- 1e-6
qmvt_df_rule <- "floor"

out_dir <- "Litter_Application_Final"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. LOAD AND PREPARE THE DATA
# =============================================================================

data("litter", package = "multcomp")

dat_full <- data.frame(
  group = factor(
    as.character(litter$dose),
    levels = c("0", "5", "50", "500")
  ),
  y = as.numeric(litter$weight),
  number = as.numeric(litter$number),
  gesttime = as.numeric(litter$gesttime)
)

# Primary auxiliary covariate:
# standardized litter size.  Standardize using the COMPLETE dataset once,
# then retain those values in every incomplete version.
dat_full$x <- as.numeric(scale(dat_full$number))

control_level <- "0"
group_levels <- levels(dat_full$group)
treat_levels <- setdiff(group_levels, control_level)
contrast_labels <- paste0(treat_levels, " vs. ", control_level)
k_treat <- length(treat_levels)

stopifnot(!anyNA(dat_full$y))
stopifnot(!anyNA(dat_full$x))
stopifnot(all(group_levels == c("0", "5", "50", "500")))

# =============================================================================
# 2. BASIC DATA SUMMARY
# =============================================================================

observed_group_summary <- dat_full %>%
  group_by(group) %>%
  summarise(
    n = n(),
    mean_y = mean(y),
    sd_y = sd(y),
    variance_y = var(y),
    mean_number = mean(number),
    sd_number = sd(number),
    mean_gesttime = mean(gesttime),
    sd_gesttime = sd(gesttime),
    cor_y_x = cor(y, x),
    .groups = "drop"
  )

variance_ratio <- max(observed_group_summary$variance_y) /
  min(observed_group_summary$variance_y)

write.csv(
  observed_group_summary,
  file.path(out_dir, "litter_OBSERVED_GROUP_SUMMARY.csv"),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("LITTER DATASET: OBSERVED COMPLETE-DATA SUMMARY\n")
cat("============================================================\n")
print(observed_group_summary)
cat(
  "\nLargest/smallest observed variance ratio = ",
  sprintf("%.3f", variance_ratio),
  "\n",
  sep = ""
)

# =============================================================================
# 3. MAR CALIBRATION AND MISSINGNESS
# =============================================================================

# Application calibration is based on the empirical X distribution in THIS
# complete dataset.  This preserves the finalized slope log(2) while making
# the expected overall missing rate exactly 30% for the observed covariates.
calibrate_empirical_mar_intercept <- function(
  x,
  target_rate,
  gamma = mar_gamma
) {
  uniroot(
    function(a) {
      mean(plogis(a + gamma * x)) - target_rate
    },
    interval = c(-30, 30),
    tol = 1e-12
  )$root
}

alpha_mar_30 <- calibrate_empirical_mar_intercept(
  dat_full$x,
  target_rate = target_missing_rate,
  gamma = mar_gamma
)

impose_mar <- function(dat, seed) {
  set.seed(seed)

  out <- dat
  p_miss <- plogis(alpha_mar_30 + mar_gamma * out$x)
  missing <- rbinom(
    n = nrow(out),
    size = 1L,
    prob = p_miss
  )

  out$missing <- missing
  out$p_missing <- p_miss
  out$y[missing == 1L] <- NA_real_

  out
}

# =============================================================================
# 4. IMPUTATION HELPERS
# =============================================================================

draw_arm_regression <- function(x_obs, y_obs) {
  n <- length(y_obs)
  X <- cbind(1, x_obs)
  p <- ncol(X)

  if (n <= p) {
    stop("Too few observed outcomes in an arm for Bayesian regression imputation.")
  }

  XtX <- crossprod(X)
  XtX_inv <- tryCatch(
    solve(XtX),
    error = function(e) {
      solve(XtX + diag(1e-8, nrow(XtX)))
    }
  )

  beta_hat <- drop(XtX_inv %*% crossprod(X, y_obs))

  resid <- y_obs - drop(X %*% beta_hat)
  df_resid <- n - p
  sse <- sum(resid^2)

  if (!is.finite(sse) || sse <= 0) {
    sse <- max(var(y_obs) * df_resid, 1e-8)
  }

  sigma2_star <- sse / rchisq(1L, df = df_resid)

  V_beta <- sigma2_star * XtX_inv
  V_beta <- (V_beta + t(V_beta)) / 2

  beta_star <- beta_hat + drop(
    t(chol(V_beta + diag(1e-12, nrow(V_beta)))) %*%
      rnorm(p)
  )

  list(
    beta_hat = beta_hat,
    beta_star = beta_star,
    sigma2_star = sigma2_star
  )
}

# -------------------------------------------------------------------------
# 4A. Proper arm-specific Bayesian regression MI
# -------------------------------------------------------------------------

impute_reg_once <- function(dat, seed) {
  set.seed(seed)
  out <- dat

  for (g in group_levels) {
    idx_g <- which(out$group == g)

    obs_g <- idx_g[!is.na(out$y[idx_g])]
    mis_g <- idx_g[ is.na(out$y[idx_g])]

    if (!length(mis_g)) next

    dr <- draw_arm_regression(
      out$x[obs_g],
      out$y[obs_g]
    )

    Xmis <- cbind(1, out$x[mis_g])
    mu_mis <- drop(Xmis %*% dr$beta_star)

    out$y[mis_g] <- rnorm(
      length(mis_g),
      mean = mu_mis,
      sd = sqrt(dr$sigma2_star)
    )
  }

  if (anyNA(out$y)) stop("REG imputation left missing values.")
  out
}

# -------------------------------------------------------------------------
# 4B. Type-1 arm-specific PMM
# -------------------------------------------------------------------------

impute_pmm_once <- function(dat, seed, donors = pmm_donors) {
  set.seed(seed)
  out <- dat

  for (g in group_levels) {
    idx_g <- which(out$group == g)

    obs_g <- idx_g[!is.na(out$y[idx_g])]
    mis_g <- idx_g[ is.na(out$y[idx_g])]

    if (!length(mis_g)) next

    dr <- draw_arm_regression(
      out$x[obs_g],
      out$y[obs_g]
    )

    Xobs <- cbind(1, out$x[obs_g])
    Xmis <- cbind(1, out$x[mis_g])

    # Type-1 PMM:
    # observed predictions use beta_hat;
    # missing predictions use beta_star.
    pred_obs <- drop(Xobs %*% dr$beta_hat)
    pred_mis <- drop(Xmis %*% dr$beta_star)

    k_use <- min(as.integer(donors), length(obs_g))

    for (jj in seq_along(mis_g)) {
      distance <- abs(pred_obs - pred_mis[jj])

      # Random tie-breaking.
      ord <- order(
        distance,
        runif(length(distance))
      )

      donor_pool <- obs_g[
        ord[seq_len(k_use)]
      ]

      chosen_donor <- sample(
        donor_pool,
        size = 1L
      )

      out$y[mis_g[jj]] <- out$y[chosen_donor]
    }
  }

  if (anyNA(out$y)) stop("PMM imputation left missing values.")
  out
}

# -------------------------------------------------------------------------
# 4C. Propensity-score subclassification + same-arm ABB
# -------------------------------------------------------------------------

approx_bayesian_bootstrap <- function(y_obs, n_missing) {
  if (n_missing == 0L) return(numeric(0))

  if (!length(y_obs)) {
    stop("No observed donors available for ABB.")
  }

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

make_propensity_strata <- function(ps, n_groups = propensity_groups) {
  # Tie-preserving subclasses:
  # equal propensity values always remain in the same subclass.

  unique_ps <- sort(unique(ps))

  if (length(unique_ps) == 1L) {
    return(rep.int(1L, length(ps)))
  }

  G_use <- min(
    as.integer(n_groups),
    length(unique_ps)
  )

  rnk <- rank(
    unique_ps,
    ties.method = "average"
  )

  subclass_unique <- ceiling(
    G_use * rnk / length(unique_ps)
  )

  subclass_unique <- pmin(
    G_use,
    pmax(1L, as.integer(subclass_unique))
  )

  key <- format(
    unique_ps,
    digits = 17,
    scientific = FALSE
  )

  lookup <- setNames(
    subclass_unique,
    key
  )

  as.integer(
    lookup[
      format(ps, digits = 17, scientific = FALSE)
    ]
  )
}

impute_ps_once <- function(
  dat,
  seed,
  n_groups = propensity_groups
) {
  set.seed(seed)
  out <- dat

  missing <- is.na(out$y)

  if (!any(missing)) return(out)

  psdat <- data.frame(
    missing = as.integer(missing),
    group = out$group,
    x = out$x
  )

  ps_fit <- suppressWarnings(
    glm(
      missing ~ group + x,
      family = binomial(),
      data = psdat
    )
  )

  ps <- as.numeric(
    predict(
      ps_fit,
      newdata = psdat,
      type = "response"
    )
  )

  # Numerical safeguard.
  ps[!is.finite(ps)] <- mean(missing)

  for (g in group_levels) {
    idx_g <- which(out$group == g)

    obs_g <- idx_g[!is.na(out$y[idx_g])]
    mis_g <- idx_g[ is.na(out$y[idx_g])]

    if (!length(mis_g)) next
    if (!length(obs_g)) {
      stop("No observed donors in arm ", g, ".")
    }

    strata_g <- make_propensity_strata(
      ps[idx_g],
      n_groups = n_groups
    )

    strata_ids <- sort(unique(strata_g))

    centers <- vapply(
      strata_ids,
      function(s) {
        mean(ps[idx_g[strata_g == s]])
      },
      numeric(1)
    )

    names(centers) <- as.character(strata_ids)

    for (s in strata_ids) {
      ids_s <- idx_g[strata_g == s]

      mis_s <- ids_s[is.na(out$y[ids_s])]

      if (!length(mis_s)) next

      obs_s <- ids_s[!is.na(out$y[ids_s])]

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

        nearest_s <- donor_strata[
          which.min(
            abs(
              centers[as.character(donor_strata)] -
              centers[as.character(s)]
            )
          )
        ]

        donor_ids <- idx_g[
          strata_g == nearest_s
        ]

        obs_s <- donor_ids[
          !is.na(out$y[donor_ids])
        ]
      }

      out$y[mis_s] <- approx_bayesian_bootstrap(
        out$y[obs_s],
        length(mis_s)
      )
    }
  }

  if (anyNA(out$y)) stop("PS imputation left missing values.")
  out
}

make_imputations <- function(
  dat,
  method,
  m = m_imputations,
  base_seed
) {
  out <- vector("list", m)

  for (mm in seq_len(m)) {
    ss <- base_seed + mm - 1L

    out[[mm]] <- switch(
      method,
      REG = impute_reg_once(dat, seed = ss),
      PMM = impute_pmm_once(dat, seed = ss),
      PS = impute_ps_once(dat, seed = ss),
      stop("Unknown MI method: ", method)
    )
  }

  out
}

# =============================================================================
# 5. SHRINKAGE + MANY-TO-ONE HELPERS
# =============================================================================

estimate_alpha_hat <- function(s2_vec, n_vec) {
  s2_vec <- as.numeric(s2_vec)
  n_vec <- as.numeric(n_vec)

  G <- length(s2_vec)

  if (length(n_vec) != G) {
    stop("s2_vec and n_vec have different lengths.")
  }

  if (
    any(!is.finite(s2_vec)) ||
    any(!is.finite(n_vec)) ||
    any(n_vec <= 1)
  ) {
    stop("Invalid variance/sample-size input for shrinkage.")
  }

  target <- mean(s2_vec)

  numerator <-
    ((G - 1) / G) *
    sum(
      2 * s2_vec^2 / (n_vec - 1)
    )

  denominator <-
    sum(
      (s2_vec - target)^2
    )

  if (
    !is.finite(denominator) ||
    denominator <= 1e-12
  ) {
    return(1)
  }

  a_hat <- numerator / denominator

  min(
    max(a_hat, 0),
    1
  )
}

shrink_variances <- function(s2_vec, n_vec) {
  a_hat <- estimate_alpha_hat(
    s2_vec,
    n_vec
  )

  target <- mean(s2_vec)

  s2_tilde <-
    a_hat * target +
    (1 - a_hat) * s2_vec

  list(
    s2 = setNames(
      as.numeric(s2_tilde),
      names(s2_vec)
    ),
    alpha_hat = a_hat
  )
}

build_U <- function(s2_vec, n_vec) {
  k <- length(treat_levels)

  s20 <- s2_vec[control_level]
  n0 <- n_vec[control_level]

  U <- matrix(
    s20 / n0,
    nrow = k,
    ncol = k,
    dimnames = list(
      contrast_labels,
      contrast_labels
    )
  )

  diag(U) <- vapply(
    treat_levels,
    function(g) {
      s20 / n0 +
      s2_vec[g] / n_vec[g]
    },
    numeric(1)
  )

  U
}

safe_solve <- function(M) {
  tryCatch(
    solve(M),
    error = function(e) {
      solve(
        M + diag(1e-8, nrow(M))
      )
    }
  )
}

make_corr_safe <- function(V) {
  V <- as.matrix(V)
  V <- (V + t(V)) / 2

  ev <- tryCatch(
    eigen(
      V,
      symmetric = TRUE,
      only.values = TRUE
    )$values,
    error = function(e) NA_real_
  )

  if (
    any(!is.finite(ev)) ||
    min(ev, na.rm = TRUE) <= 1e-8
  ) {
    ridge <- if (any(is.finite(ev))) {
      abs(min(ev, na.rm = TRUE)) + 1e-8
    } else {
      1e-6
    }

    diag(V) <- diag(V) + ridge
  }

  corr <- cov2cor(V)

  corr[!is.finite(corr)] <- 0
  diag(corr) <- 1

  corr
}

as_qmvt_df <- function(df) {
  if (!is.finite(df) || df <= 0) {
    return(1000000L)
  }

  if (qmvt_df_rule == "floor") {
    return(
      as.integer(
        max(1, floor(df))
      )
    )
  }

  if (qmvt_df_rule == "round") {
    return(
      as.integer(
        max(1, round(df))
      )
    )
  }

  if (qmvt_df_rule == "ceiling") {
    return(
      as.integer(
        max(1, ceiling(df))
      )
    )
  }

  stop("Unknown qmvt_df_rule.")
}

safe_qmvt_crit <- function(
  corr,
  df,
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  corr <- make_corr_safe(corr)

  out <- mvtnorm::qmvt(
    p = 1 - alpha_level,
    tail = "both.tails",
    df = as_qmvt_df(df),
    corr = corr,
    algorithm = mvtnorm::GenzBretz(
      maxpts = qmvt_maxpts,
      abseps = qmvt_abseps
    )
  )

  crit <- as.numeric(out$quantile)

  if (!is.finite(crit) || crit <= 0) {
    stop("Invalid mv-t critical value.")
  }

  crit
}

# =============================================================================
# 6. COMPLETE / CASE-DELETION ANALYSIS
# =============================================================================

analyze_one_dataset <- function(
  dat,
  use_shrinkage,
  seed
) {
  dat2 <- dat[
    !is.na(dat$y),
    ,
    drop = FALSE
  ]

  split_y <- split(
    dat2$y,
    dat2$group
  )

  # Need every group represented.
  if (!all(group_levels %in% names(split_y))) {
    stop("At least one treatment arm contains no observed outcomes.")
  }

  means <- vapply(
    split_y[group_levels],
    mean,
    numeric(1)
  )

  s2 <- vapply(
    split_y[group_levels],
    var,
    numeric(1)
  )

  n_vec <- vapply(
    split_y[group_levels],
    length,
    integer(1)
  )

  Q <- means[treat_levels] -
    means[control_level]

  names(Q) <- contrast_labels

  a_hat <- NA_real_

  if (use_shrinkage) {
    shr <- shrink_variances(
      s2,
      n_vec
    )

    s2_use <- shr$s2
    a_hat <- shr$alpha_hat
  } else {
    s2_use <- s2
  }

  U <- build_U(
    s2_use,
    n_vec
  )

  # Working complete/case-deletion df used consistently
  # with the paper simulation implementation.
  df_working <- sum(n_vec) -
    length(group_levels)

  crit <- safe_qmvt_crit(
    corr = U,
    df = df_working,
    seed = seed
  )

  se <- sqrt(diag(U))
  lower <- Q - crit * se
  upper <- Q + crit * se

  data.frame(
    contrast = contrast_labels,
    estimate = as.numeric(Q),
    se = as.numeric(se),
    lower = as.numeric(lower),
    upper = as.numeric(upper),
    length = as.numeric(upper - lower),
    alpha_hat = a_hat,
    stringsAsFactors = FALSE
  )
}

# Complete-data empirical benchmark.
complete_benchmark <- analyze_one_dataset(
  dat_full,
  use_shrinkage = FALSE,
  seed = master_seed + 10L
)

write.csv(
  complete_benchmark,
  file.path(out_dir, "litter_COMPLETE_DATA_BENCHMARK.csv"),
  row.names = FALSE
)

# Also report the shrinkage intensity on the fully observed dataset.
complete_s2 <- tapply(
  dat_full$y,
  dat_full$group,
  var
)[group_levels]

complete_n <- table(dat_full$group)[group_levels]

complete_alpha_hat <- estimate_alpha_hat(
  complete_s2,
  complete_n
)

cat(
  "\nComplete-data shrinkage alpha_hat = ",
  sprintf("%.4f", complete_alpha_hat),
  "\n",
  sep = ""
)

# =============================================================================
# 7. RUBIN POOLING: SHRINKAGE AND NON-SHRINKAGE
# =============================================================================

pool_mi_both <- function(
  imputed_list,
  seed_shrink,
  seed_nonshrink
) {
  m <- length(imputed_list)
  k <- length(treat_levels)

  Qmat <- matrix(
    NA_real_,
    nrow = m,
    ncol = k,
    dimnames = list(
      NULL,
      contrast_labels
    )
  )

  U_ns <- array(
    NA_real_,
    dim = c(k, k, m)
  )

  U_sh <- array(
    NA_real_,
    dim = c(k, k, m)
  )

  alpha_vec <- numeric(m)

  for (mm in seq_len(m)) {
    d <- imputed_list[[mm]]

    means <- tapply(
      d$y,
      d$group,
      mean
    )[group_levels]

    s2 <- tapply(
      d$y,
      d$group,
      var
    )[group_levels]

    n_vec <- table(d$group)[group_levels]

    Qmat[mm, ] <-
      means[treat_levels] -
      means[control_level]

    U_ns[, , mm] <- build_U(
      s2,
      n_vec
    )

    shr <- shrink_variances(
      s2,
      n_vec
    )

    alpha_vec[mm] <- shr$alpha_hat

    U_sh[, , mm] <- build_U(
      shr$s2,
      n_vec
    )
  }

  Qbar <- colMeans(Qmat)
  B <- stats::cov(Qmat)

  pool_one <- function(
    U_array,
    version,
    seed
  ) {
    Ubar <- apply(
      U_array,
      c(1, 2),
      mean
    )

    Tmat <-
      Ubar +
      (1 + 1 / m) * B

    Tmat <- (Tmat + t(Tmat)) / 2

    if (max(abs(B)) <= 1e-14) {
      df_mi <- 1000000
    } else {
      r_m <-
        (1 + 1 / m) *
        sum(
          diag(
            B %*%
              safe_solve(Ubar)
          )
        ) /
        k

      if (!is.finite(r_m) || r_m < 0) {
        stop("Invalid relative increase in variance.")
      }

      df_mi <- if (r_m <= 1e-14) {
        1000000
      } else {
        (m - 1) *
          (1 + 1 / r_m)^2
      }
    }

    crit <- safe_qmvt_crit(
      corr = Tmat,
      df = df_mi,
      seed = seed
    )

    se <- sqrt(diag(Tmat))
    lower <- Qbar - crit * se
    upper <- Qbar + crit * se

    data.frame(
      version = version,
      contrast = contrast_labels,
      estimate = as.numeric(Qbar),
      se = as.numeric(se),
      lower = as.numeric(lower),
      upper = as.numeric(upper),
      length = as.numeric(upper - lower),
      df_mi = df_mi,
      stringsAsFactors = FALSE
    )
  }

  list(
    results = bind_rows(
      pool_one(
        U_sh,
        version = "Shrinkage",
        seed = seed_shrink
      ),
      pool_one(
        U_ns,
        version = "Non-shrinkage",
        seed = seed_nonshrink
      )
    ),
    Qmat = Qmat,
    alpha_vec = alpha_vec,
    mean_alpha = mean(alpha_vec)
  )
}

# =============================================================================
# 8. PRIMARY APPLICATION: ONE PRESPECIFIED MAR REALIZATION
# =============================================================================

primary_missing_seed <- master_seed + 1000L

dat_primary <- impose_mar(
  dat_full,
  seed = primary_missing_seed
)

primary_missing_summary <- dat_primary %>%
  group_by(group) %>%
  summarise(
    n_total = n(),
    n_observed = sum(!is.na(y)),
    n_missing = sum(is.na(y)),
    missing_pct = 100 * mean(is.na(y)),
    mean_p_missing = 100 * mean(p_missing),
    .groups = "drop"
  )

write.csv(
  primary_missing_summary,
  file.path(out_dir, "litter_PRIMARY_MISSINGNESS_BY_ARM.csv"),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("PRIMARY 30% MAR APPLICATION\n")
cat("============================================================\n")
cat(
  "Realized overall missingness = ",
  sum(is.na(dat_primary$y)),
  "/",
  nrow(dat_primary),
  " = ",
  sprintf("%.2f%%", 100 * mean(is.na(dat_primary$y))),
  "\n",
  sep = ""
)
print(primary_missing_summary)

# -------------------------------------------------------------------------
# Run all three MI methods so we can inspect robustness.
# -------------------------------------------------------------------------

all_methods <- c("REG")

primary_all_results <- list()
primary_pools <- list()

for (method in all_methods) {
  imp <- make_imputations(
    dat_primary,
    method = method,
    m = m_imputations,
    base_seed =
      master_seed +
      10000L +
      match(method, all_methods) * 1000L
  )

  pooled <- pool_mi_both(
    imp,
    seed_shrink =
      master_seed +
      20000L +
      match(method, all_methods) * 1000L,
    seed_nonshrink =
      master_seed +
      30000L +
      match(method, all_methods) * 1000L
  )

  primary_pools[[method]] <- list(
    imputed = imp,
    pooled = pooled
  )

  rr <- pooled$results
  rr$MI <- method
  rr$mean_alpha <- pooled$mean_alpha

  primary_all_results[[method]] <- rr
}

primary_all_results <- bind_rows(
  primary_all_results
) %>%
  select(
    MI,
    version,
    contrast,
    estimate,
    se,
    lower,
    upper,
    length,
    mean_alpha,
    df_mi
  )

write.csv(
  primary_all_results,
  file.path(out_dir, "litter_ALL_MI_APPLICATION_RESULTS.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# Naive case deletion comparator.
# -------------------------------------------------------------------------

case_deletion <- analyze_one_dataset(
  dat_primary,
  use_shrinkage = FALSE,
  seed = master_seed + 40000L
)

# -------------------------------------------------------------------------
# Primary paper table = REG only.
# -------------------------------------------------------------------------

reg_pooled <- primary_pools[["REG"]]$pooled$results

primary_reg_table <- bind_rows(
  reg_pooled %>%
    mutate(
      method = ifelse(
        version == "Shrinkage",
        "Proposed method with shrinkage",
        "Proposed method without shrinkage"
      )
    ),
  case_deletion %>%
    mutate(
      version = NA_character_,
      df_mi = NA_real_,
      method = "Dunnett-type procedure with naive case deletion"
    )
) %>%
  mutate(
    CI = sprintf(
      "(%.3f, %.3f)",
      lower,
      upper
    )
  ) %>%
  select(
    method,
    contrast,
    estimate,
    CI,
    length
  )

write.csv(
  primary_reg_table,
  file.path(out_dir, "litter_PRIMARY_REG_APPLICATION_TABLE.csv"),
  row.names = FALSE
)

cat("\nPRIMARY REGRESSION-MI APPLICATION TABLE:\n")
print(primary_reg_table)

# -------------------------------------------------------------------------
# Shrinkage alpha output.
# -------------------------------------------------------------------------

primary_alpha <- data.frame(
  imputation = seq_len(m_imputations),
  alpha_hat =
    primary_pools[["REG"]]$pooled$alpha_vec
)

write.csv(
  primary_alpha,
  file.path(out_dir, "litter_PRIMARY_SHRINKAGE_ALPHA.csv"),
  row.names = FALSE
)

# =============================================================================
# 9. REFERENCE-PAPER-STYLE FOREST PLOT
# =============================================================================
# Individual completed datasets are displayed using NON-SHRINKAGE
# within-dataset inference.  The pooled shrinkage and non-shrinkage results are
# both added, together with the complete-data benchmark and case deletion.

reg_imputed <- primary_pools[["REG"]]$imputed

forest_individual <- vector(
  "list",
  m_imputations
)

for (mm in seq_len(m_imputations)) {
  tmp <- analyze_one_dataset(
    reg_imputed[[mm]],
    use_shrinkage = FALSE,
    seed = master_seed + 50000L + mm
  )

  tmp$display <- as.character(mm)

  forest_individual[[mm]] <- tmp
}

forest_individual <- bind_rows(
  forest_individual
)

forest_pooled <- reg_pooled %>%
  mutate(
    display = ifelse(
      version == "Shrinkage",
      "Pooled (shrinkage)",
      "Pooled (non-shrinkage)"
    )
  ) %>%
  select(
    contrast,
    estimate,
    lower,
    upper,
    display
  )

forest_case <- case_deletion %>%
  mutate(
    display = "Case deletion"
  ) %>%
  select(
    contrast,
    estimate,
    lower,
    upper,
    display
  )

forest_complete <- complete_benchmark %>%
  mutate(
    display = "Complete-data benchmark"
  ) %>%
  select(
    contrast,
    estimate,
    lower,
    upper,
    display
  )

forest_df <- bind_rows(
  forest_individual %>%
    select(
      contrast,
      estimate,
      lower,
      upper,
      display
    ),
  forest_pooled,
  forest_case,
  forest_complete
)

forest_levels <- c(
  "Complete-data benchmark",
  "Case deletion",
  "Pooled (shrinkage)",
  "Pooled (non-shrinkage)",
  as.character(
    rev(seq_len(m_imputations))
  )
)

forest_df$display <- factor(
  forest_df$display,
  levels = forest_levels
)

forest_df$contrast <- factor(
  forest_df$contrast,
  levels = contrast_labels
)

p_forest <- ggplot(
  forest_df,
  aes(
    x = estimate,
    y = display,
    xmin = lower,
    xmax = upper
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.4
  ) +
  geom_pointrange(
    linewidth = 0.28
  ) +
  facet_wrap(
    ~ contrast,
    scales = "free_x",
    nrow = 1
  ) +
  labs(
    title = "Litter data: individual imputations and pooled many-to-one inference",
    subtitle = paste0(
      "40 regression imputations; target MAR = 30%; realized missingness = ",
      sprintf("%.1f%%", 100 * mean(is.na(dat_primary$y)))
    ),
    x = "Treatment-control mean difference",
    y = "Analysis / imputed dataset"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 6)
  )

ggsave(
  filename = file.path(
    out_dir,
    "litter_PRIMARY_REFERENCE_STYLE_FOREST.pdf"
  ),
  plot = p_forest,
  width = 12,
  height = 9
)

ggsave(
  filename = file.path(
    out_dir,
    "litter_PRIMARY_REFERENCE_STYLE_FOREST.png"
  ),
  plot = p_forest,
  width = 12,
  height = 9,
  dpi = 300
)

# =============================================================================
# 10. REPEATED-MAR BENCHMARK VALIDATION
# =============================================================================
# This section evaluates whether the incomplete-data interval contains the
# COMPLETE-SAMPLE contrast, which is treated as an empirical benchmark.
#
# It is NOT population coverage because the complete sample itself is random
# relative to an unknown population.

benchmark_Q <- setNames(
  complete_benchmark$estimate,
  complete_benchmark$contrast
)

validation_one <- function(rep_id) {
  incomplete <- impose_mar(
    dat_full,
    seed = master_seed + 1000000L + rep_id
  )

  # Require enough observed cases in each arm for Y~X.
  obs_by_arm <- table(
    incomplete$group[
      !is.na(incomplete$y)
    ]
  )

  if (
    !all(group_levels %in% names(obs_by_arm)) ||
    any(obs_by_arm[group_levels] <= 3)
  ) {
    return(NULL)
  }

  out <- list()

  for (method in all_methods) {
    imp <- tryCatch(
      make_imputations(
        incomplete,
        method = method,
        m = m_imputations,
        base_seed =
          master_seed +
          2000000L +
          rep_id * 10000L +
          match(method, all_methods) * 1000L
      ),
      error = function(e) NULL
    )

    if (is.null(imp)) next

    pooled <- tryCatch(
      pool_mi_both(
        imp,
        seed_shrink =
          master_seed +
          3000000L +
          rep_id * 10000L +
          match(method, all_methods) * 1000L,
        seed_nonshrink =
          master_seed +
          4000000L +
          rep_id * 10000L +
          match(method, all_methods) * 1000L
      ),
      error = function(e) NULL
    )

    if (is.null(pooled)) next

    rr <- pooled$results

    for (vv in c("Shrinkage", "Non-shrinkage")) {
      q <- rr[
        rr$version == vv,
        ,
        drop = FALSE
      ]

      true_vec <- benchmark_Q[
        q$contrast
      ]

      recovered <- (
        q$lower <= true_vec &
        q$upper >= true_vec
      )

      q$replicate <- rep_id
      q$MI <- method
      q$benchmark <- as.numeric(true_vec)
      q$recovered <- recovered
      q$joint_recovered <- all(recovered)
      q$mean_interval_length <- mean(q$length)
      q$mean_alpha <- pooled$mean_alpha
      q$realized_missing_pct <-
        100 * mean(is.na(incomplete$y))

      out[[paste(method, vv)]] <- q
    }
  }

  if (!length(out)) return(NULL)

  bind_rows(out)
}

if (run_validation) {
  cat("\n============================================================\n")
  cat("REPEATED-MAR BENCHMARK VALIDATION\n")
  cat("Replications requested: ", n_validation, "\n", sep = "")
  cat("============================================================\n")

  validation_list <- vector(
    "list",
    n_validation
  )

  for (bb in seq_len(n_validation)) {
    validation_list[[bb]] <- validation_one(bb)

    if (
      bb %% max(1L, floor(n_validation / 20L)) == 0L
    ) {
      cat(
        "Completed ",
        bb,
        " / ",
        n_validation,
        " validation replications\n",
        sep = ""
      )
    }
  }

  validation_raw <- bind_rows(
    validation_list
  )

  if (!nrow(validation_raw)) {
    stop("No successful validation replications.")
  }

  # -----------------------------------------------------------------------
  # Joint benchmark recovery + average interval length.
  # One row per method/version.
  # -----------------------------------------------------------------------

  validation_rep_level <- validation_raw %>%
    group_by(
      replicate,
      MI,
      version
    ) %>%
    summarise(
      joint_recovered = first(joint_recovered),
      mean_interval_length = first(mean_interval_length),
      mean_alpha = first(mean_alpha),
      realized_missing_pct = first(realized_missing_pct),
      .groups = "drop"
    )

  validation_summary <- validation_rep_level %>%
    group_by(
      MI,
      version
    ) %>%
    summarise(
      successful_replicates = n(),
      joint_benchmark_recovery_pct =
        100 * mean(joint_recovered),
      mean_interval_length =
        mean(mean_interval_length),
      sd_interval_length =
        sd(mean_interval_length),
      mean_realized_missing_pct =
        mean(realized_missing_pct),
      mean_alpha =
        ifelse(
          version == "Shrinkage",
          mean(mean_alpha),
          NA_real_
        ),
      .groups = "drop"
    )

  write.csv(
    validation_summary,
    file.path(
      out_dir,
      "litter_VALIDATION_SUMMARY.csv"
    ),
    row.names = FALSE
  )

  # -----------------------------------------------------------------------
  # Contrast-level benchmark recovery + SE calibration.
  # Empirical SD is computed from pooled estimates across repeated MAR runs.
  # Because shrinkage does not change Qbar, empirical SD should be nearly
  # identical between shrinkage and non-shrinkage for a given MI method.
  # -----------------------------------------------------------------------

  validation_by_contrast <- validation_raw %>%
    group_by(
      MI,
      version,
      contrast
    ) %>%
    summarise(
      benchmark = first(benchmark),
      benchmark_recovery_pct =
        100 * mean(recovered),
      mean_estimate =
        mean(estimate),
      bias_to_benchmark =
        mean(estimate - benchmark),
      empirical_sd =
        sd(estimate),
      mean_model_se =
        mean(se),
      se_to_empirical_sd =
        mean(se) / sd(estimate),
      mean_interval_length =
        mean(length),
      .groups = "drop"
    )

  write.csv(
    validation_by_contrast,
    file.path(
      out_dir,
      "litter_VALIDATION_BY_CONTRAST.csv"
    ),
    row.names = FALSE
  )

  # -----------------------------------------------------------------------
  # Direct shrinkage minus non-shrinkage differences.
  # Positive recovery difference favors shrinkage.
  # Negative length difference means shrinkage is shorter.
  # -----------------------------------------------------------------------

  sh <- validation_summary %>%
    filter(version == "Shrinkage") %>%
    select(
      MI,
      sh_recovery =
        joint_benchmark_recovery_pct,
      sh_length =
        mean_interval_length,
      sh_alpha =
        mean_alpha
    )

  ns <- validation_summary %>%
    filter(version == "Non-shrinkage") %>%
    select(
      MI,
      ns_recovery =
        joint_benchmark_recovery_pct,
      ns_length =
        mean_interval_length
    )

  validation_difference <- inner_join(
    sh,
    ns,
    by = "MI"
  ) %>%
    mutate(
      recovery_difference_pp =
        sh_recovery - ns_recovery,
      interval_length_difference =
        sh_length - ns_length,
      relative_length_change_pct =
        100 *
        (sh_length - ns_length) /
        ns_length
    )

  write.csv(
    validation_difference,
    file.path(
      out_dir,
      "litter_VALIDATION_SHRINKAGE_MINUS_NONSHRINK.csv"
    ),
    row.names = FALSE
  )

  cat("\nVALIDATION SUMMARY:\n")
  print(validation_summary)

  cat("\nSHRINKAGE MINUS NON-SHRINKAGE:\n")
  print(validation_difference)
}

# =============================================================================
# 11. FINAL CONSOLE SUMMARY
# =============================================================================

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
cat("Output folder:\n")
cat(normalizePath(out_dir), "\n\n")

cat("Main application files:\n")
cat("  litter_OBSERVED_GROUP_SUMMARY.csv\n")
cat("  litter_COMPLETE_DATA_BENCHMARK.csv\n")
cat("  litter_PRIMARY_REG_APPLICATION_TABLE.csv\n")
cat("  litter_ALL_MI_APPLICATION_RESULTS.csv\n")
cat("  litter_PRIMARY_REFERENCE_STYLE_FOREST.pdf\n")
cat("  litter_PRIMARY_REFERENCE_STYLE_FOREST.png\n")
cat("  litter_PRIMARY_MISSINGNESS_BY_ARM.csv\n")
cat("  litter_PRIMARY_SHRINKAGE_ALPHA.csv\n\n")

if (run_validation) {
  cat("Repeated-MAR validation files:\n")
  cat("  litter_VALIDATION_SUMMARY.csv\n")
  cat("  litter_VALIDATION_BY_CONTRAST.csv\n")
  cat("  litter_VALIDATION_SHRINKAGE_MINUS_NONSHRINK.csv\n")
}

cat("\nFor the final publication run, set:\n")
cat("  publication_run <- TRUE\n")
cat("which changes n_validation from 500 to 10000.\n")
