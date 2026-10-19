# Completed ten-treatment simulation

**All 100,000 jobs succeeded:** 10,000 replicates for each of two scenarios and five missing-rate settings, with 40 imputations. Original seeds and numerical settings were retained. No simulation jobs failed. One nonfatal `gt::fmt_missing()` deprecation warning was captured during table formatting.

## Download

- [Download the results ZIP](Ten_Treatment_Final_Tables_and_Plots.zip): tables, PDFs, plots, console output, original R script, warning/error logs, and runtime information. Open the file, then choose **Download raw file** on GitHub.
- [Full PDF report](simulation_report.pdf): 24 pages of tables and plots.
- [Tables PDF](main_tables.pdf) and [plots PDF](simulation_plots.pdf).
- [Main results in Markdown](main_results_chat.md).

## CSV tables

- Table 2 coverage: [shrinkage](table2_shrink.csv) and [non-shrinkage](table2_nonshrink.csv).
- Table 3 interval lengths and expansion rates: [shrinkage](table3_shrink.csv) and [non-shrinkage](table3_nonshrink.csv).
- Numerical per-comparison summaries: [coverage](coverage_summary.csv) and [interval length](length_summary.csv).
- [Joint coverage](joint_coverage_summary.csv), [mean interval length across comparisons](mean_length_summary.csv), and [completion counts](completion_summary.csv).

The ZIP contains console output and the original R script. Warning/error records are also provided separately. Supplementary plots were generated from the simulation results; the original script defines tables only.

The two raw per-replicate CSVs each contain 2,600,000 result rows and are approximately 700 MB. They remain saved in the cloud environment and are not included in this ZIP or this Git repository. The CSV summaries and four primary tables above are included in full. The five-treatment simulation is still running; no final five-treatment results are included here.

`SHA256SUMS.txt` records the results ZIP checksum.
