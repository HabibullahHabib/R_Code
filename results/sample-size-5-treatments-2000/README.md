# sample-size-5-treatments-2000

Sample-size sensitivity analysis completed. See completion_summary.csv and run_summary.txt for exact replicate and success/failure counts. Forty imputations. Original statistical code and seeds preserved; the execution wrapper records the requested replicate count. Seven sample sizes (20,50,100,150,200,250,300 per arm), two mean scenarios and five missing rates. Shrinkage and non-shrinkage results are included.

[Tables, plots, logs and source ZIP](sample-size-5-treatments-2000_Tables_Plots_and_Logs.zip) · [Tables PDF](main_tables.pdf) · [Plots PDF](simulation_plots.pdf)

Original primary tables are provided individually as CSV and HTML. Raw replicate CSVs are gzip compressed and provided separately. Files ending .part001, .part002, etc. are consecutive binary segments of a single gzip file; concatenate in numeric order before decompressing. SHA256SUMS.json verifies downloaded files. No reduced-replicate trial results are included. Supplementary plots show coverage and mean interval length across sample sizes; CSVs retain the original detailed tables.
