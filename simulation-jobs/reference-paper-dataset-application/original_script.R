# ============================================================================
# REFERENCE-PAPER DATASET APPLICATION
# Jiang, Zhou & Liang (2023), Table 5 original incomplete dose-response data
#
# PURPOSE
# -------
# Use the EXACT 100-patient incomplete dataset reported in the main reference
# paper and compare:
#
#   1) Proposed many-to-one procedure WITHOUT variance shrinkage
#   2) Proposed many-to-one procedure WITH variance shrinkage
#
# FINAL PAPER CONCEPT RETAINED WHERE THE DATA SUPPORT IT:
#   * m = 40 imputations
#   * same imputed datasets are used for shrinkage and non-shrinkage
#   * Rubin pooling
#   * many-to-one multivariate-t simultaneous 95% confidence intervals
#   * James-Stein-type variance shrinkage toward the mean arm variance
#
# IMPORTANT DATA LIMITATION
# -------------------------
# The published application dataset contains only dose and response.
# It does NOT provide an auxiliary continuous X.
# Therefore, it is not statistically honest to force the finalized Y~X
# regression/PMM application models onto this dataset.
#
# The reference paper used treatment-allocation propensity-score imputation.
# With treatment group as the only reported fully observed predictor, the
# resulting donor structure reduces naturally to treatment-restricted ABB:
# missing responses in a dose arm are imputed from observed responses in the
# SAME dose arm. This is also consistent with the imputed values shown in
# Table 5 of the reference paper.
#
# MAIN OUTPUT FIRST:
#   * one primary result table
#   * one forest plot
#
# APPENDIX:
#   * observed-data group summaries / variances
#   * missingness by arm
#   * shrinkage alpha
#   * reference paper's reported Table 6 numbers for comparison
# ============================================================================

rm(list = ls())

# ============================================================================
# 0. PACKAGES AND SETTINGS
# ============================================================================

needed <- c("mvtnorm", "dplyr", "tibble", "ggplot2", "gt")
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
  library(mvtnorm)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(gt)
})

master_seed <- 20261108L
set.seed(master_seed)

alpha <- 0.05
m_imputations <- 40L

output_dir <- "ReferencePaper_Data_FINAL_Application"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# 1. EXACT ORIGINAL DATA FROM REFERENCE PAPER TABLE 5
# ============================================================================

dose0 <- c(
  NA, 0.14, -0.02, NA, NA, 0.36, 0.31, NA, 0.62, NA,
  -0.31, -0.45, -0.20, NA, NA, NA, -0.16, NA, 0.81, NA
)

dose005 <- c(
  NA, NA, -0.07, 0.58, 0.96, NA, -0.00, 0.52, -0.36, NA,
  NA, 0.44, -0.02, 0.37, 1.53, NA, 0.37, 1.01, 0.28, 0.99
)

dose020 <- c(
  1.20, 1.57, -0.16, 0.21, 1.85, 1.00, 2.45, -0.52, 0.05, 0.63,
  0.53, 0.42, 1.23, 1.87, 1.06, 0.35, NA, 0.48, 0.57, 1.01
)

dose060 <- c(
  1.17, NA, 1.78, 0.31, 0.06, 0.90, 0.74, 0.23, 1.39, 0.91,
  NA, NA, NA, NA, 1.60, 1.58, NA, NA, 2.16, 0.69
)

dose100 <- c(
  2.25, NA, 1.25, 1.86, NA, 1.20, 1.97, 0.63, NA, NA,
  0.89, 0.56, 0.73, NA, 0.49, NA, NA, NA, 1.13, 0.95
)

group_levels <- c("0", "0.05", "0.20", "0.60", "1.00")
control_level <- "0"
treat_levels <- group_levels[-1]
treat_labels <- paste0(treat_levels, " vs 0")

dat <- tibble(
  group = factor(
    rep(group_levels, each = 20L),
    levels = group_levels
  ),
  y = c(dose0, dose005, dose020, dose060, dose100)
)

stopifnot(nrow(dat) == 100L)
stopifnot(sum(is.na(dat$y)) == 32L)

# ============================================================================
# 2. OBSERVED-DATA SUMMARY
# ============================================================================

observed_summary <- dat |>
  group_by(group) |>
  summarise(
    n_total = n(),
    n_observed = sum(!is.na(y)),
    n_missing = sum(is.na(y)),
    missing_pct = 100 * mean(is.na(y)),
    observed_mean = mean(y, na.rm = TRUE),
    observed_sd = sd(y, na.rm = TRUE),
    observed_variance = var(y, na.rm = TRUE),
    .groups = "drop"
  )

variance_ratio <- max(observed_summary$observed_variance) /
  min(observed_summary$observed_variance)

cat("\nObserved-data variance summary:\n")
print(observed_summary)
cat(
  "\nLargest/smallest observed variance ratio = ",
  sprintf("%.3f", variance_ratio),
  "\n",
  sep = ""
)

# ============================================================================
# 3. REFERENCE-PAPER PROPENSITY / SAME-ARM ABB IMPUTATION
# ============================================================================

approx_bayesian_bootstrap <- function(y_obs, n_missing) {
  if (n_missing == 0L) return(numeric(0))
  if (!length(y_obs)) stop("No observed donors available.")

  # Approximate Bayesian bootstrap:
  # first resample the empirical donor distribution,
  # then draw imputed values from that bootstrap distribution.
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

impute_reference_once <- function(dat, seed) {
  set.seed(seed)
  out <- dat

  for (g in group_levels) {
    idx_g <- which(out$group == g)
    obs_g <- idx_g[!is.na(out$y[idx_g])]
    mis_g <- idx_g[ is.na(out$y[idx_g])]

    if (!length(mis_g)) next
    if (!length(obs_g)) {
      stop("No observed donors in dose group ", g, ".")
    }

    out$y[mis_g] <- approx_bayesian_bootstrap(
      out$y[obs_g],
      length(mis_g)
    )
  }

  if (anyNA(out$y)) stop("Imputation left missing responses.")
  out
}

make_imputations <- function(dat, m = m_imputations, base_seed = master_seed) {
  lapply(
    seq_len(m),
    function(mm) {
      impute_reference_once(
        dat,
        seed = base_seed + mm - 1L
      )
    }
  )
}

# ============================================================================
# 4. SHRINKAGE HELPERS
# Finalized unequal-n formula; here n_g = 20 in all completed datasets.
# ============================================================================

estimate_alpha_hat <- function(s2_vec, n_vec) {
  s2_vec <- as.numeric(s2_vec)
  n_vec <- as.numeric(n_vec)
  G <- length(s2_vec)

  if (length(n_vec) != G) stop("s2_vec and n_vec lengths do not match.")
  if (any(!is.finite(s2_vec)) || any(!is.finite(n_vec)) || any(n_vec <= 1)) {
    stop("Invalid variance/sample-size input.")
  }

  target <- mean(s2_vec)

  num <- ((G - 1) / G) *
    sum(2 * s2_vec^2 / (n_vec - 1))

  den <- sum((s2_vec - target)^2)

  if (!is.finite(den) || den <= 1e-12) return(1)

  a_hat <- num / den
  min(max(a_hat, 0), 1)
}

shrink_variances <- function(s2_vec, n_vec) {
  a <- estimate_alpha_hat(s2_vec, n_vec)
  target <- mean(s2_vec)

  list(
    s2 = setNames(
      a * target + (1 - a) * s2_vec,
      names(s2_vec)
    ),
    alpha_hat = a
  )
}

# ============================================================================
# 5. MANY-TO-ONE COVARIANCE MATRIX
# ============================================================================

build_U <- function(s2_vec, n_vec) {
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
    function(g) {
      s20 / n0 + s2_vec[g] / n_vec[g]
    },
    numeric(1)
  )

  U
}

safe_solve <- function(M) {
  tryCatch(
    solve(M),
    error = function(e) {
      solve(M + diag(1e-8, nrow(M)))
    }
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
    ridge <- if (is.finite(min_ev)) {
      abs(min_ev) + 1e-8
    } else {
      1e-6
    }

    diag(V) <- diag(V) + ridge
  }

  C <- cov2cor(V)
  C[!is.finite(C)] <- 0
  diag(C) <- 1
  C
}

safe_qmvt_crit <- function(corr, df, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  df_use <- if (!is.finite(df) || df <= 0) {
    1e6
  } else {
    max(1, floor(df))
  }

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

  if (!is.finite(crit) || crit <= 0) {
    stop("Invalid multivariate-t critical value.")
  }

  crit
}

# ============================================================================
# 6. SUMMARIZE THE SAME 40 IMPUTED DATASETS
# ============================================================================

summarize_completions <- function(comps) {
  m <- length(comps)

  means <- matrix(
    NA_real_,
    nrow = m,
    ncol = length(group_levels),
    dimnames = list(NULL, group_levels)
  )

  vars <- means

  for (mm in seq_len(m)) {
    means[mm, ] <- tapply(
      comps[[mm]]$y,
      comps[[mm]]$group,
      mean
    )[group_levels]

    vars[mm, ] <- tapply(
      comps[[mm]]$y,
      comps[[mm]]$group,
      var
    )[group_levels]
  }

  n_vec <- setNames(
    rep(20L, length(group_levels)),
    group_levels
  )

  list(
    means = means,
    vars = vars,
    n = n_vec,
    m = m
  )
}

# ============================================================================
# 7. RUBIN + MANY-TO-ONE POOLING
#    SHRINKAGE AND NON-SHRINKAGE FROM SAME IMPUTATIONS
# ============================================================================

pool_both <- function(comps) {
  st <- summarize_completions(comps)

  m <- st$m
  k <- length(treat_levels)

  Qmat <- matrix(
    NA_real_,
    nrow = m,
    ncol = k,
    dimnames = list(NULL, treat_labels)
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
    mu <- st$means[mm, ]
    s2 <- st$vars[mm, ]
    names(s2) <- group_levels

    Qmat[mm, ] <- mu[treat_levels] - mu[control_level]

    U_ns[, , mm] <- build_U(
      s2_vec = s2,
      n_vec = st$n
    )

    shr <- shrink_variances(
      s2_vec = s2,
      n_vec = st$n
    )

    alpha_vec[mm] <- shr$alpha_hat

    U_sh[, , mm] <- build_U(
      s2_vec = shr$s2,
      n_vec = st$n
    )
  }

  Qbar <- colMeans(Qmat)
  B <- cov(Qmat)

  pool_one <- function(U_array, variance_method, seed) {
    Ubar <- apply(
      U_array,
      c(1, 2),
      mean
    )

    Tm <- Ubar + (1 + 1/m) * B
    Tm <- (Tm + t(Tm)) / 2

    if (max(abs(B)) <= 1e-14) {
      nu <- 1e6
    } else {
      rm <- (1 + 1/m) *
        sum(diag(B %*% safe_solve(Ubar))) /
        k

      if (!is.finite(rm) || rm < 0) {
        stop("Invalid relative increase in variance.")
      }

      nu <- if (rm <= 1e-14) {
        1e6
      } else {
        (m - 1) * (1 + 1/rm)^2
      }
    }

    crit <- safe_qmvt_crit(
      corr = make_corr_safe(Tm),
      df = nu,
      seed = seed
    )

    se <- sqrt(diag(Tm))
    lower <- Qbar - crit * se
    upper <- Qbar + crit * se

    tibble(
      Procedure = variance_method,
      Comparison = treat_labels,
      Estimate = as.numeric(Qbar),
      Lower = as.numeric(lower),
      Upper = as.numeric(upper),
      Length = as.numeric(upper - lower),
      Significant = lower > 0 | upper < 0
    )
  }

  list(
    results = bind_rows(
      pool_one(
        U_sh,
        variance_method = "With shrinkage",
        seed = master_seed + 100000L
      ),
      pool_one(
        U_ns,
        variance_method = "Without shrinkage",
        seed = master_seed + 200000L
      )
    ),
    mean_alpha = mean(alpha_vec),
    alpha_by_imputation = alpha_vec
  )
}

# ============================================================================
# 8. RUN THE PRIMARY ANALYSIS
# ============================================================================

comps <- make_imputations(
  dat,
  m = m_imputations,
  base_seed = master_seed + 1000L
)

pooled <- pool_both(comps)

results_long <- pooled$results

# Order: without shrinkage first, then with shrinkage,
# matching the style used in the user's previous application result.
results_long$Procedure <- factor(
  results_long$Procedure,
  levels = c("Without shrinkage", "With shrinkage")
)

results_long <- results_long |>
  arrange(Procedure, Comparison)

# ============================================================================
# 9. PRIMARY RESULT TABLE - SHOWN FIRST
# ============================================================================

primary_table <- results_long |>
  mutate(
    `95% simultaneous confidence interval` =
      sprintf("(%.3f, %.3f)", Lower, Upper),
    `Interval length` = round(Length, 3),
    Conclusion = ifelse(
      Significant,
      "Significant",
      "Not significant"
    )
  ) |>
  select(
    Procedure,
    `Treatment-control comparison` = Comparison,
    `Point estimate` = Estimate,
    `95% simultaneous confidence interval`,
    `Interval length`,
    Conclusion
  )

cat("\n\n============================================================\n")
cat("REFERENCE-PAPER DATASET - PRIMARY RESULT\n")
cat("Original missing values retained: 32/100 = 32%\n")
cat("m = 40 treatment-restricted ABB imputations\n")
cat("Shrinkage and non-shrinkage use the SAME imputations\n")
cat("Observed variance ratio = ",
    sprintf("%.3f", variance_ratio), "\n", sep = "")
cat("============================================================\n\n")

print(primary_table)

gt_primary <- primary_table |>
  gt() |>
  tab_header(
    title = "Reference-paper dose-response data",
    subtitle = paste0(
      "Proposed many-to-one procedure: shrinkage vs non-shrinkage; ",
      "40 imputations; original 32% missingness"
    )
  ) |>
  fmt_number(
    columns = c("Point estimate"),
    decimals = 3
  ) |>
  cols_align(
    align = "center",
    columns = everything()
  )

print(gt_primary)

gtsave(
  gt_primary,
  file.path(
    output_dir,
    "ReferencePaper_PRIMARY_RESULT_TABLE.html"
  )
)

write.csv(
  primary_table,
  file.path(
    output_dir,
    "ReferencePaper_PRIMARY_RESULT_TABLE.csv"
  ),
  row.names = FALSE
)

# ============================================================================
# 10. ONE PRIMARY FOREST PLOT
# ============================================================================

p <- ggplot(
  results_long,
  aes(
    x = Comparison,
    y = Estimate,
    ymin = Lower,
    ymax = Upper,
    shape = Procedure
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = 2,
    linewidth = 0.5
  ) +
  geom_pointrange(
    position = position_dodge(width = 0.55)
  ) +
  coord_flip() +
  labs(
    title = "Reference-paper dose-response data",
    subtitle = "Shrinkage vs non-shrinkage simultaneous 95% confidence intervals",
    x = NULL,
    y = "Dose - placebo mean difference",
    shape = "Procedure"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(p)

ggsave(
  file.path(
    output_dir,
    "ReferencePaper_PRIMARY_FOREST_PLOT.png"
  ),
  p,
  width = 9,
  height = 6,
  dpi = 300
)

# ============================================================================
# 11. APPENDIX ONLY
# ============================================================================

# A. Observed-data variance / missingness summary.
write.csv(
  observed_summary,
  file.path(
    output_dir,
    "APPENDIX_observed_group_variances_and_missingness.csv"
  ),
  row.names = FALSE
)

# B. Shrinkage intensity over 40 imputations.
appendix_alpha <- tibble(
  Imputation = seq_len(m_imputations),
  alpha_hat = pooled$alpha_by_imputation
)

write.csv(
  appendix_alpha,
  file.path(
    output_dir,
    "APPENDIX_shrinkage_alpha_by_imputation.csv"
  ),
  row.names = FALSE
)

# C. Reference paper's reported Table 6 values, for comparison only.
# These are NOT generated by this script.
reference_reported <- tibble(
  Method = rep(
    c(
      "Reference paper proposed method",
      "Reference paper naive case deletion"
    ),
    each = 4
  ),
  Comparison = rep(
    c(
      "0.05 vs 0",
      "0.20 vs 0",
      "0.60 vs 0",
      "1.00 vs 0"
    ),
    times = 2
  ),
  Estimate = c(
    0.374, 0.712, 0.951, 1.060,
    0.362, 0.722, 0.929, 1.050
  ),
  Lower = c(
    -0.181, 0.222, 0.370, 0.538,
    -0.270, 0.125, 0.287, 0.396
  ),
  Upper = c(
    0.930, 1.203, 1.532, 1.582,
    0.995, 1.319, 1.572, 1.704
  )
)

write.csv(
  reference_reported,
  file.path(
    output_dir,
    "APPENDIX_reference_paper_reported_Table6.csv"
  ),
  row.names = FALSE
)

# D. One-line dataset summary.
appendix_dataset_summary <- tibble(
  n_total = nrow(dat),
  n_missing = sum(is.na(dat$y)),
  missing_pct = 100*mean(is.na(dat$y)),
  observed_variance_ratio = variance_ratio,
  mean_shrinkage_alpha = pooled$mean_alpha
)

write.csv(
  appendix_dataset_summary,
  file.path(
    output_dir,
    "APPENDIX_dataset_summary.csv"
  ),
  row.names = FALSE
)

cat("\nMain output files:\n")
cat("  ReferencePaper_PRIMARY_RESULT_TABLE.html\n")
cat("  ReferencePaper_PRIMARY_RESULT_TABLE.csv\n")
cat("  ReferencePaper_PRIMARY_FOREST_PLOT.png\n")
cat("\nEverything else is saved as APPENDIX_*.csv\n")
cat("Output folder: ", normalizePath(output_dir), "\n", sep = "")
