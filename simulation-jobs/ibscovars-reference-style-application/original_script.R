
# ============================================================================
# IBScovars APPLICATION - REFERENCE-PAPER STYLE OUTPUT
#
# GOAL
# ----
# Use IBScovars as the application dataset, but present the result in the
# same style as the main reference paper:
#
#   Table:
#     Methods | Dose Comparison | Point Estimate of (mu_i - mu_0) |
#     95% Confidence Interval
#
#   Figure:
#     forest plot showing each imputed dataset result, pooled result(s),
#     and naive case-deletion result, panelled by dose comparison.
#
# FINAL CONCEPT RETAINED
# ----------------------
# Missingness mechanism:
#     logit P(M = 1 | X) = alpha_r + log(2) X
# same across groups
#
# X:
#     standardized gender indicator
#
# Imputation methods available:
#     "propensity"  (default, closest to the reference-paper application)
#     "regression"
#     "pmm"
#
# This script creates:
#   1) one primary result table in reference-paper style
#   2) one forest plot similar to the reference paper
#   3) LaTeX files you can directly include in your manuscript
# ============================================================================

rm(list = ls())

# =========================
# 0. PACKAGES AND SETTINGS
# =========================
needed <- c("DoseFinding", "mvtnorm", "dplyr", "ggplot2")
missing_pkg <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkg)) {
  stop("Please install first: ", paste(missing_pkg, collapse = ", "))
}

suppressPackageStartupMessages({
  library(DoseFinding)
  library(mvtnorm)
  library(dplyr)
  library(ggplot2)
})

master_seed <- 20261108L
set.seed(master_seed)

# Main application choices
target_missing_rate <- 0.30
m_imputations <- 40L
alpha <- 0.05

# Final concept
mar_gamma <- log(2)
pmm_donors <- 5L
propensity_groups <- 3L

# Choose ONE method for the main application presentation.
# This keeps the result simple and matches the reference-paper style.
# Options: "propensity", "regression", "pmm"
imputation_method <- "propensity"

# Output folder
out_dir <- "IBScovars_ReferenceStyle_Output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =========================
# 1. DATA
# =========================
data(IBScovars, package = "DoseFinding")

gender_num <- if (is.factor(IBScovars$gender)) {
  as.numeric(IBScovars$gender)
} else {
  as.numeric(IBScovars$gender)
}

dat_full <- data.frame(
  group = factor(IBScovars$dose,
                 levels = sort(unique(IBScovars$dose))),
  x = as.numeric(scale(gender_num)),
  y = as.numeric(IBScovars$resp)
)

control_level <- levels(dat_full$group)[1]
treat_levels <- setdiff(levels(dat_full$group), control_level)
contrast_labels <- paste0(treat_levels, " vs. ", control_level)

# =========================
# 2. MISSINGNESS MECHANISM
# =========================
calibrate_empirical_intercept <- function(x, target_rate, gamma = mar_gamma) {
  uniroot(
    function(a) mean(plogis(a + gamma * x)) - target_rate,
    interval = c(-30, 30),
    tol = 1e-12
  )$root
}

impose_final_mar <- function(dat, target_rate, seed) {
  set.seed(seed)
  out <- dat

  a0 <- calibrate_empirical_intercept(
    out$x,
    target_rate = target_rate,
    gamma = mar_gamma
  )

  p <- plogis(a0 + mar_gamma * out$x)
  M <- rbinom(nrow(out), size = 1L, prob = p)
  out$y[M == 1L] <- NA_real_

  attr(out, "alpha_mar") <- a0
  attr(out, "expected_missing") <- mean(p)
  out
}

# =========================
# 3. IMPUTATION FUNCTIONS
# =========================
draw_arm_regression <- function(x_obs, y_obs) {
  n <- length(y_obs)
  X <- cbind(1, x_obs)
  p <- ncol(X)

  if (n <= p) stop("Too few observed outcomes for regression imputation.")

  XtX <- crossprod(X)
  XtX_inv <- solve(XtX)
  beta_hat <- drop(XtX_inv %*% crossprod(X, y_obs))
  resid <- y_obs - drop(X %*% beta_hat)
  df <- n - p
  sse <- sum(resid^2)

  sigma2_star <- sse / rchisq(1L, df = df)
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

    dr <- draw_arm_regression(out$x[obs], out$y[obs])
    Xmis <- cbind(1, out$x[mis])
    mu_mis <- drop(Xmis %*% dr$beta_star)
    out$y[mis] <- rnorm(length(mis), mean = mu_mis, sd = sqrt(dr$sigma2_star))
  }

  out
}

impute_pmm_once <- function(dat, seed, donors = 5L) {
  set.seed(seed)
  out <- dat

  for (g in levels(out$group)) {
    idx <- which(out$group == g)
    obs <- idx[!is.na(out$y[idx])]
    mis <- idx[ is.na(out$y[idx])]
    if (!length(mis)) next

    dr <- draw_arm_regression(out$x[obs], out$y[obs])
    Xobs <- cbind(1, out$x[obs])
    Xmis <- cbind(1, out$x[mis])

    # Type-1 PMM
    pred_obs <- drop(Xobs %*% dr$beta_hat)
    pred_mis <- drop(Xmis %*% dr$beta_star)

    k_use <- min(as.integer(donors), length(obs))

    for (j in seq_along(mis)) {
      d <- abs(pred_obs - pred_mis[j])
      ord <- order(d, runif(length(d)))
      donor_pool <- obs[ord[seq_len(k_use)]]
      donor <- sample(donor_pool, 1L)
      out$y[mis[j]] <- out$y[donor]
    }
  }

  out
}

approx_bayesian_bootstrap <- function(y_obs, n_missing) {
  if (n_missing == 0L) return(numeric(0))
  y_star <- sample(y_obs, size = length(y_obs), replace = TRUE)
  sample(y_star, size = n_missing, replace = TRUE)
}

make_propensity_strata <- function(ps, n_groups = 3L) {
  unique_ps <- sort(unique(ps))
  if (length(unique_ps) == 1L) return(rep.int(1L, length(ps)))

  G <- min(as.integer(n_groups), length(unique_ps))
  score_rank <- rank(unique_ps, ties.method = "average")
  score_group <- ceiling(G * score_rank / length(unique_ps))
  score_group <- pmin(G, pmax(1L, as.integer(score_group)))

  mp <- setNames(score_group, format(unique_ps, digits = 17, scientific = FALSE))
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

  fit <- suppressWarnings(glm(missing ~ group + x, family = binomial(), data = psdat))
  ps <- as.numeric(predict(fit, newdata = psdat, type = "response"))

  for (g in levels(out$group)) {
    idx_g <- which(out$group == g)
    mis_g <- idx_g[is.na(out$y[idx_g])]
    obs_g <- idx_g[!is.na(out$y[idx_g])]
    if (!length(mis_g)) next

    strata_g <- make_propensity_strata(ps[idx_g], n_groups = n_groups)
    strata_ids <- sort(unique(strata_g))
    centers <- vapply(strata_ids, function(s) mean(ps[idx_g[strata_g == s]]), numeric(1))
    names(centers) <- as.character(strata_ids)

    for (s in strata_ids) {
      ids_s <- idx_g[strata_g == s]
      mis_s <- ids_s[is.na(out$y[ids_s])]
      if (!length(mis_s)) next

      obs_s <- ids_s[!is.na(out$y[ids_s])]
      if (!length(obs_s)) {
        donor_strata <- strata_ids[
          vapply(strata_ids, function(ss) {
            ids_ss <- idx_g[strata_g == ss]
            any(!is.na(out$y[ids_ss]))
          }, logical(1))
        ]
        nearest <- donor_strata[which.min(abs(
          centers[as.character(donor_strata)] - centers[as.character(s)]
        ))]
        donor_ids <- idx_g[strata_g == nearest]
        obs_s <- donor_ids[!is.na(out$y[donor_ids])]
      }

      out$y[mis_s] <- approx_bayesian_bootstrap(out$y[obs_s], length(mis_s))
    }
  }

  out
}

make_imputations <- function(dat, method, m, base_seed) {
  out <- vector("list", m)

  for (mm in seq_len(m)) {
    ss <- base_seed + mm - 1L
    out[[mm]] <- switch(
      method,
      regression = impute_regression_once(dat, seed = ss),
      pmm = impute_pmm_once(dat, seed = ss, donors = pmm_donors),
      propensity = impute_propensity_once(dat, seed = ss, n_groups = propensity_groups),
      stop("Unknown imputation method")
    )
  }

  out
}

# =========================
# 4. MANY-TO-ONE INFERENCE
# =========================
estimate_alpha_hat <- function(s2_vec, n_vec) {
  G <- length(s2_vec)
  target <- mean(s2_vec)
  num <- ((G - 1) / G) * sum(2 * s2_vec^2 / (n_vec - 1))
  den <- sum((s2_vec - target)^2)
  if (!is.finite(den) || den <= 1e-12) return(1)
  a_hat <- num / den
  min(max(a_hat, 0), 1)
}

shrink_variances <- function(s2_vec, n_vec) {
  a_hat <- estimate_alpha_hat(s2_vec, n_vec)
  target <- mean(s2_vec)
  list(
    s2 = setNames(a_hat * target + (1 - a_hat) * s2_vec, names(s2_vec)),
    alpha_hat = a_hat
  )
}

build_U <- function(s2_vec, n_vec, group_levels, control_level) {
  treat_levels <- setdiff(group_levels, control_level)
  labels <- paste0(treat_levels, " vs. ", control_level)
  k <- length(treat_levels)

  s20 <- s2_vec[control_level]
  n0 <- n_vec[control_level]

  U <- matrix(s20 / n0, nrow = k, ncol = k, dimnames = list(labels, labels))
  diag(U) <- vapply(
    treat_levels,
    function(g) s20 / n0 + s2_vec[g] / n_vec[g],
    numeric(1)
  )
  U
}

safe_solve <- function(M) {
  tryCatch(solve(M), error = function(e) solve(M + diag(1e-8, nrow(M))))
}

make_corr_safe <- function(V) {
  V <- as.matrix(V)
  V <- (V + t(V)) / 2
  ev <- tryCatch(eigen(V, symmetric = TRUE, only.values = TRUE)$values,
                 error = function(e) NA_real_)
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
    algorithm = mvtnorm::GenzBretz(maxpts = 50000L, abseps = 1e-6)
  )
  as.numeric(out$quantile)
}

analyze_single_complete <- function(dat, use_shrinkage, control_level,
                                    alpha = 0.05, seed = NULL) {
  dat2 <- dat[!is.na(dat$y), , drop = FALSE]
  group_levels <- levels(dat2$group)
  treat_levels <- setdiff(group_levels, control_level)
  labels <- paste0(treat_levels, " vs. ", control_level)

  split_y <- split(dat2$y, dat2$group)
  means <- vapply(split_y, mean, numeric(1))[group_levels]
  s2 <- vapply(split_y, var, numeric(1))[group_levels]
  n_vec <- vapply(split_y, length, integer(1))[group_levels]

  Q <- means[treat_levels] - means[control_level]
  names(Q) <- labels

  if (use_shrinkage) {
    shr <- shrink_variances(s2, n_vec)
    s2 <- shr$s2
    a_hat <- shr$alpha_hat
  } else {
    a_hat <- NA_real_
  }

  U <- build_U(s2, n_vec, group_levels, control_level)
  df <- sum(n_vec) - length(group_levels)
  crit <- safe_qmvt_crit(make_corr_safe(U), df = df, alpha = alpha, seed = seed)

  se <- sqrt(diag(U))
  lower <- Q - crit * se
  upper <- Q + crit * se

  data.frame(
    contrast = labels,
    estimate = as.numeric(Q),
    lower = as.numeric(lower),
    upper = as.numeric(upper),
    width = as.numeric(upper - lower),
    significant = lower > 0 | upper < 0,
    shrink_alpha = a_hat
  )
}

pool_many_to_one_both <- function(comp_list, control_level,
                                  alpha = 0.05,
                                  seed_sh = NULL,
                                  seed_ns = NULL) {
  group_levels <- levels(comp_list[[1]]$group)
  treat_levels <- setdiff(group_levels, control_level)
  labels <- paste0(treat_levels, " vs. ", control_level)
  k <- length(treat_levels)
  m <- length(comp_list)

  means <- matrix(NA_real_, nrow = m, ncol = length(group_levels),
                  dimnames = list(NULL, group_levels))
  vars <- means

  for (mm in seq_len(m)) {
    means[mm, ] <- tapply(comp_list[[mm]]$y, comp_list[[mm]]$group, mean)[group_levels]
    vars[mm, ] <- tapply(comp_list[[mm]]$y, comp_list[[mm]]$group, var)[group_levels]
  }

  n_vec <- table(comp_list[[1]]$group)
  n_vec <- as.numeric(n_vec[group_levels])
  names(n_vec) <- group_levels

  Qmat <- matrix(NA_real_, nrow = m, ncol = k, dimnames = list(NULL, labels))
  U_ns <- array(NA_real_, dim = c(k, k, m))
  U_sh <- array(NA_real_, dim = c(k, k, m))
  alpha_vec <- numeric(m)

  for (mm in seq_len(m)) {
    mu <- means[mm, ]
    s2 <- vars[mm, ]
    names(s2) <- group_levels

    Qmat[mm, ] <- mu[treat_levels] - mu[control_level]
    U_ns[, , mm] <- build_U(s2, n_vec, group_levels, control_level)

    shr <- shrink_variances(s2, n_vec)
    alpha_vec[mm] <- shr$alpha_hat
    U_sh[, , mm] <- build_U(shr$s2, n_vec, group_levels, control_level)
  }

  Qbar <- colMeans(Qmat)
  B <- cov(Qmat)

  pool_one <- function(U_array, seed, method_name, mean_alpha) {
    Ubar <- apply(U_array, c(1, 2), mean)
    Tm <- Ubar + (1 + 1/m) * B
    Tm <- (Tm + t(Tm)) / 2

    if (max(abs(B)) <= 1e-14) {
      nu <- 1e6
    } else {
      rm <- (1 + 1/m) * sum(diag(B %*% safe_solve(Ubar))) / k
      nu <- if (rm <= 1e-14) 1e6 else (m - 1) * (1 + 1 / rm)^2
    }

    crit <- safe_qmvt_crit(make_corr_safe(Tm), df = nu, alpha = alpha, seed = seed)
    se <- sqrt(diag(Tm))
    lower <- Qbar - crit * se
    upper <- Qbar + crit * se

    data.frame(
      method = method_name,
      contrast = labels,
      estimate = as.numeric(Qbar),
      lower = as.numeric(lower),
      upper = as.numeric(upper),
      width = as.numeric(upper - lower),
      significant = lower > 0 | upper < 0,
      shrink_alpha = mean_alpha
    )
  }

  bind_rows(
    pool_one(U_sh, seed_sh, "With shrinkage", mean(alpha_vec)),
    pool_one(U_ns, seed_ns, "Without shrinkage", NA_real_)
  )
}

# =========================
# 5. RUN MAIN ANALYSIS
# =========================
dat_incomplete <- impose_final_mar(dat_full, target_missing_rate, seed = master_seed + 5000L)
achieved_missing <- mean(is.na(dat_incomplete$y))

# Same completed datasets used for shrinkage and non-shrinkage
comps <- make_imputations(
  dat_incomplete,
  method = imputation_method,
  m = m_imputations,
  base_seed = master_seed + 100000L
)

# Pooled results
pooled_res <- pool_many_to_one_both(
  comps,
  control_level = control_level,
  alpha = alpha,
  seed_sh = master_seed + 200000L,
  seed_ns = master_seed + 300000L
)

# Complete-case result (naive deletion)
cc_res <- analyze_single_complete(
  dat_incomplete,
  use_shrinkage = FALSE,
  control_level = control_level,
  alpha = alpha,
  seed = master_seed + 400000L
) %>%
  mutate(method = "Dunnett's procedure with naive case deletion", .before = 1)

# For the forest plot, compute each completed dataset result (no shrinkage),
# analogous to the reference paper's per-imputed-dataset display.
per_imp_list <- vector("list", m_imputations)
for (mm in seq_len(m_imputations)) {
  tmp <- analyze_single_complete(
    comps[[mm]],
    use_shrinkage = FALSE,
    control_level = control_level,
    alpha = alpha,
    seed = master_seed + 500000L + mm
  )
  tmp$label <- as.character(mm)
  per_imp_list[[mm]] <- tmp
}
per_imp_df <- bind_rows(per_imp_list)

# Add pooled and case-deletion lines to the figure
plot_df <- bind_rows(
  per_imp_df %>% mutate(label = factor(label, levels = as.character(seq_len(m_imputations)))),
  pooled_res %>%
    mutate(label = ifelse(method == "Without shrinkage", "Pooled (NS)", "Pooled (Sh)")) %>%
    select(contrast, estimate, lower, upper, width, significant, shrink_alpha, label),
  cc_res %>%
    mutate(label = "Case deletion") %>%
    select(contrast, estimate, lower, upper, width, significant, shrink_alpha, label)
)

plot_df$label <- factor(
  plot_df$label,
  levels = c(as.character(seq_len(m_imputations)),
             "Pooled (NS)", "Pooled (Sh)", "Case deletion")
)

# =========================
# 6. PRIMARY TABLE (REFERENCE STYLE)
# =========================
# Main methods shown in the table
main_table_df <- bind_rows(
  pooled_res %>%
    mutate(method = ifelse(method == "With shrinkage",
                           "Proposed method with shrinkage",
                           "Proposed method without shrinkage")),
  cc_res
) %>%
  mutate(
    CI = sprintf("(%.3f, %.3f)", lower, upper),
    dose_comparison = contrast
  ) %>%
  select(method, dose_comparison, estimate, CI)

# Order methods to mimic paper style
main_table_df$method <- factor(
  main_table_df$method,
  levels = c(
    "Proposed method with shrinkage",
    "Proposed method without shrinkage",
    "Dunnett's procedure with naive case deletion"
  )
)

main_table_df <- main_table_df %>%
  arrange(method, dose_comparison)

# CSV for direct use
write.csv(
  main_table_df,
  file.path(out_dir, "IBScovars_ReferenceStyle_Table.csv"),
  row.names = FALSE
)

# Pretty printed version
cat("\n============================================================\n")
cat("IBScovars reference-style application result\n")
cat("Primary missingness: 30% MAR (achieved ",
    sprintf("%.1f%%", 100 * achieved_missing), ")\n", sep = "")
cat("Imputation method shown: ", imputation_method, "\n", sep = "")
cat("============================================================\n\n")
print(main_table_df)

# =========================
# 7. FOREST PLOT (REFERENCE STYLE)
# =========================
plot_title_method <- switch(
  imputation_method,
  propensity = "propensity-score imputation",
  regression = "regression imputation",
  pmm = "PMM imputation"
)

p <- ggplot(
  plot_df,
  aes(
    y = label,
    x = estimate,
    xmin = lower,
    xmax = upper
  )
) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4) +
  geom_pointrange(size = 0.25) +
  facet_wrap(~ contrast, scales = "free_x", ncol = 2) +
  labs(
    title = "Forest plot showing the point estimates and 95% confidence intervals",
    subtitle = paste0(
      "IBScovars, 30% MAR; ",
      plot_title_method,
      "; individual imputations, pooled results, and naive case deletion"
    ),
    x = "95% Confidence Interval",
    y = "Label of Imputed Dataset"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(
  file.path(out_dir, "IBScovars_ReferenceStyle_ForestPlot.png"),
  p,
  width = 10,
  height = 8,
  dpi = 300
)

ggsave(
  file.path(out_dir, "IBScovars_ReferenceStyle_ForestPlot.pdf"),
  p,
  width = 10,
  height = 8
)

# =========================
# 8. LATEX TABLE CODE
# =========================
make_latex_table <- function(df, caption, label) {
  methods <- as.character(df$method)
  comps <- as.character(df$dose_comparison)
  ests <- sprintf("%.3f", df$estimate)
  cis <- as.character(df$CI)

  # Build multirow-style text without requiring external formatting effort.
  out <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\caption{" %+% caption %+% "}",
    "\\label{" %+% label %+% "}",
    "\\begin{tabular}{llc}",
    "\\hline",
    "Methods & Dose Comparison & Point Estimate of $(\\mu_i-\\mu_0)$ \\\\",
    " &  & 95\\% Confidence Interval \\\\",
    "\\hline"
  )

  prev_method <- ""
  for (i in seq_len(nrow(df))) {
    mtxt <- if (methods[i] == prev_method) "" else methods[i]
    row1 <- paste0(
      mtxt, " & ",
      comps[i], " & ",
      ests[i], " \\\\"
    )
    row2 <- paste0(" &  & ", cis[i], " \\\\")
    out <- c(out, row1, row2)
    prev_method <- methods[i]
  }

  out <- c(out, "\\hline", "\\end{tabular}", "\\end{table}")
  paste(out, collapse = "\n")
}

`%+%` <- function(a, b) paste0(a, b)

latex_table_code <- make_latex_table(
  main_table_df,
  caption = paste0(
    "Simultaneous statistical inferences of mean differences in IBScovars dose response ",
    "based on ", plot_title_method, "."
  ),
  label = "tab:ibs_application"
)

writeLines(
  latex_table_code,
  con = file.path(out_dir, "IBScovars_ReferenceStyle_Table.tex")
)

# =========================
# 9. LATEX FIGURE CODE
# =========================
latex_figure_code <- paste(
  "\\begin{figure}[htbp]",
  "\\centering",
  "\\includegraphics[width=0.95\\textwidth]{IBScovars_ReferenceStyle_ForestPlot.pdf}",
  "\\caption{Forest plot showing the point estimates and 95\\% confidence intervals ",
  "for each imputed dataset, the pooled results, and the naive case-deletion result ",
  "for the IBScovars application. The primary missing-data setting is 30\\% MAR and the ",
  "imputation method shown is ", plot_title_method, ".}",
  "\\label{fig:ibs_application_forest}",
  "\\end{figure}",
  sep = "\n"
)

writeLines(
  latex_figure_code,
  con = file.path(out_dir, "IBScovars_ReferenceStyle_Figure.tex")
)

# Combined LaTeX include file
combined_tex <- paste(
  "% Add this line in your preamble if needed:",
  "% \\usepackage{graphicx}",
  "",
  "% Table",
  latex_table_code,
  "",
  "% Figure",
  latex_figure_code,
  sep = "\n"
)

writeLines(
  combined_tex,
  con = file.path(out_dir, "IBScovars_ReferenceStyle_Manuscript_Insert.tex")
)

# =========================
# 10. APPENDIX / SUPPORTING FILES
# =========================
# Small appendix-style support, but not too many files.
group_summary <- dat_full %>%
  group_by(group) %>%
  summarise(
    n = n(),
    mean = mean(y),
    sd = sd(y),
    variance = var(y),
    .groups = "drop"
  )

missing_summary <- dat_incomplete %>%
  group_by(group) %>%
  summarise(
    n = n(),
    missing_n = sum(is.na(y)),
    missing_pct = 100 * mean(is.na(y)),
    .groups = "drop"
  )

write.csv(group_summary,
          file.path(out_dir, "IBScovars_group_summary.csv"),
          row.names = FALSE)

write.csv(missing_summary,
          file.path(out_dir, "IBScovars_primary_missing_summary.csv"),
          row.names = FALSE)

write.csv(
  pooled_res %>% filter(method == "With shrinkage") %>% select(contrast, shrink_alpha),
  file.path(out_dir, "IBScovars_shrinkage_alpha.csv"),
  row.names = FALSE
)

cat("\nCreated files in folder: ", normalizePath(out_dir), "\n", sep = "")
cat("\nMain files:\n")
cat("  IBScovars_ReferenceStyle_Table.csv\n")
cat("  IBScovars_ReferenceStyle_Table.tex\n")
cat("  IBScovars_ReferenceStyle_ForestPlot.png\n")
cat("  IBScovars_ReferenceStyle_ForestPlot.pdf\n")
cat("  IBScovars_ReferenceStyle_Figure.tex\n")
cat("  IBScovars_ReferenceStyle_Manuscript_Insert.tex\n")
