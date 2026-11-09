# Reference-paper dataset application

Runs the newly uploaded FINAL_Application_ReferencePaper_Dataset.R as supplied. The original 100 dose-response values and all missing-value locations exactly match Jiang, Zhou and Liang (2023), Table 5. Five dose arms have 20 patients each; original 32% missingness is retained.

This is **one real-data application with 40 treatment-restricted ABB imputations**, not a repeated-simulation study. The data contain dose and response only, so the uploaded script uses within-arm approximate Bayesian bootstrap rather than regression with an invented covariate. Both proposed procedures use identical imputations and point estimates; only shrinkage versus non-shrinkage is compared. No equal-variance baseline was added.

The original paper used 20 propensity-score imputations. This dataset transcription is exact, but these 40-imputation numerical results are from the supplied script, not a claim of exact reproduction of the original imputation implementation or reported Table 6. The APPENDIX_reference_paper_reported_Table6.csv contains transcribed published numbers, clearly separate from computed results.

The primary table follows the thesis Table 4.2 layout and contains actual simultaneous 95% confidence intervals. The combined PDF contains the table and forest plot. Console output, warnings/errors, all original CSV/HTML/PNG outputs, precise computed intervals, dataset verification and an RDS snapshot of all 40 imputations are retained.

[Full PDF report](simulation_report.pdf) · [Primary table PDF](main_tables.pdf) · [Forest plot](comparison_graph.png) · [Results ZIP](ReferencePaper_Dataset_Application_Results.zip)

[Full-precision computed results](computed_results_full_precision.csv) · [Dataset verification](DATA_VERIFICATION.json).
