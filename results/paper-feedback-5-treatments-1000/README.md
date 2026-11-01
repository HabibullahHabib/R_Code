# paper-feedback-5-treatments-1000

Completed 60,000 simulation jobs: 1,000 replicates for each of 60 settings, 40 imputations. Original statistical code and seeds preserved; the execution wrapper overrides n_sim to 1000 as requested. Three variance ratios (2,4,9), two MAR mechanisms, two mean scenarios, five missing rates. Shrinkage, non-shrinkage and pooled benchmark included, plus naive deletion bias diagnostics.

[Tables, plots, logs and source ZIP](paper-feedback-5-treatments-1000_Tables_Plots_and_Logs.zip) · [Tables PDF](main_tables.pdf) · [Plots PDF](simulation_plots.pdf)

Original primary tables are provided individually as CSV and HTML. Raw replicate CSVs are gzip compressed and provided separately. Files ending .part001, .part002, etc. are consecutive binary segments of a single gzip file; concatenate in numeric order before decompressing. SHA256SUMS.json verifies downloaded files. No reduced-replicate trial results are included. Supplementary plots average diagnostics across comparisons; CSVs preserve individual comparisons.
