# application-actg175-1000

Completed **1,000 repeated-deletion experiments total**. ACTG175 uses 1,000 experiments at 30% MCAR on the observed-outcome benchmark subset; its original natural-missingness primary analysis is also retained. IBScovars uses 250 experiments at each of 10%, 20%, 30% and 40% missingness. Every experiment evaluates all three MI methods with 40 imputations and both shrinkage and non-shrinkage intervals.

[Results ZIP](application-actg175-1000_Tables_Plots_and_Logs.zip) · [Full PDF report](simulation_report.pdf) · [Tables PDF](main_tables.pdf) · [Plots PDF](simulation_plots.pdf)

CSV tables and raw validation results are included, together with source, logs, warnings/errors and session information. Benchmark recovery is measured against fixed complete-sample contrasts, not population parameters. Both variance versions use the same imputations. All original statistical functions and seeds are preserved; only validation_reps is overridden. No increase beyond 1,000 experiments is enabled.
