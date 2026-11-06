# Sample-size studies: 2,000 total jobs per script

The user explicitly selected **2,000 total jobs per script**, not 2,000 per setting. Each script retains all 70 settings (7 sample sizes × 2 scenarios × 5 missing rates), allocated 28 or 29 replicates each according to allocation_plan.csv. The allocation is balanced to 1,000 jobs per scenario and 400 per missingness rate. Forty imputations, original scientific functions and original replicate seeds are retained.

The queue runs 5 treatments, then 10 treatments, each with four workers and restartable checkpoints. Short trials validate CSV/HTML tables and PDF exports. Exact function comparisons gate full execution. These files are execution setup, not completed results. Reports and final CSVs will appear under results/sample-size-5-treatments-2000 and results/sample-size-10-treatments-2000 when ready. Small per-setting replicate counts make coverage estimates preliminary.
