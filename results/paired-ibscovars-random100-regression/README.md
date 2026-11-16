# IBScovars: one random subset of 100 patients

One simple random sample without replacement from the 369 original patients, selected once with seed 20261115 before any outcome comparison. selected_patients.csv records original row IDs, dose, gender and response. No searching over subsets or changes based on results. The same 100 patients are retained in all 10000 repetitions.

Each repetition creates a new 30%-expected, gender-driven MAR pattern, then uses 40 proper arm-specific regression imputations shared by shrinkage and non-shrinkage. Original statistical functions and simulation seed schedule are unchanged from the full-cohort run. Gender is standardized within this subset and MAR intercept calibrated over the subset. Four workers, atomic 10-repetition checkpoints, remote backups every 500 new repetitions.

Failed repetitions are recorded rather than replaced or hidden. Benchmark recovery is against fixed complete-subset contrasts, not population-parameter coverage. Narrower intervals alone do not establish superior statistical validity. The script does not use a variance-equality pretest.

Completed 10000 attempted repetitions: 9976 successful and 24 failed.

[PDF report](simulation_report.pdf) · [Plain-language graph](comparison_graph.png) · [Graph PDF](comparison_graph.pdf) · [Results ZIP](paired-ibscovars-random100-regression_Results.zip)

[Summary CSV](comparison_summary.csv) · [Paired comparison CSV](paired_comparison.csv) · [Raw paired repetitions](paired_repetitions.csv) · [Raw intervals](raw_intervals.csv). Protocol and source files document the exact settings.
