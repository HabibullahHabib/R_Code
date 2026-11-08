# Corrected application-actg175-corrected: 10,000 imputations

Completed the corrected application script with **10,000 imputations per method in each application analysis**, as explicitly confirmed by the user. No repeated-deletion simulation loop was added. ACTG175 retains natural outcome missingness; IBScovars retains its primary 30% MAR analysis and appendix 10%, 20%, 40% MAR analyses.

[Results ZIP](application-actg175-corrected_Tables_Plots_and_Logs.zip) · [Full PDF report](simulation_report.pdf) · [Tables PDF](main_tables.pdf) · [Plots PDF](simulation_plots.pdf)

Original source, seeds, statistical formulas, MAR mechanisms, donors and propensity groups are preserved. Only m_imputations is overridden to 10000; plot subtitles are updated accordingly. Two worker processes per application generate seeded imputations with atomic 100-imputation checkpoints. Every method/analysis verifies exact equality of its first three worker-generated imputations against the original serial function. Both variance versions share all 10,000 imputations.

CSV tables, original HTML tables and PNG/PDF figures, full-precision primary results, completion counts, logs, warnings/errors and session information are included. Imputed-data checkpoints remain in the cloud workspace and local backups rather than the results archive. No automatic further increase is configured.
