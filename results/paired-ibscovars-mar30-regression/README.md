# Confirmed paired shrinkage versus non-shrinkage study

10,000 repetitions per dataset, 40 regression imputations per repetition. Two entirely separate studies; no cross-dataset analysis, equal-variance comparison or imputation-method comparison. Both procedures share identical data, imputations and point estimates within each repetition. Original corrected application statistical functions are loaded unchanged.

ACTG175 retains all original patients, baseline covariate and natural outcome missingness. Each repetition changes only the imputation draws; no outcomes are additionally masked and no complete-case restriction is imposed. Its graph compares interval widths and imputation variability. Recovery coverage cannot be computed for unknown missing outcomes.

IBScovars starts with the original complete data. Each repetition creates a new covariate-driven Bernoulli MAR pattern: gamma=log(2) in all arms, empirically calibrated to 30% expected missingness. The original complete-data contrasts are the known benchmark. Its graph shows simultaneous benchmark recovery and mean interval width. This fixed-sample recovery is not population-parameter coverage.

Each study uses two workers, a pre-specified independent dataset seed offset, atomic 10-repetition checkpoints, CSV results/logs/warnings/errors, a plain-language comparison graph (PNG/PDF), and a PDF report. Short trials are stored separately from full results. Full runs are restartable using ./run.sh full. No extra procedures or automatic repetition increases are configured.

Completed 10000 attempted repetitions: 10000 successful and 0 failed.

[PDF report](simulation_report.pdf) · [Plain-language graph](comparison_graph.png) · [Graph PDF](comparison_graph.pdf) · [Results ZIP](paired-ibscovars-mar30-regression_Results.zip)

[Summary CSV](comparison_summary.csv) · [Paired comparison CSV](paired_comparison.csv) · [Raw paired repetitions](paired_repetitions.csv) · [Raw intervals](raw_intervals.csv). Protocol and source files document the exact settings.
