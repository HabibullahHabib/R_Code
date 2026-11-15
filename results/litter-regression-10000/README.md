# Litter: 10000 repetitions, 40 regression imputations per repetition

Extends the previous 500-repetition regression-only study to 10000 total repetitions, using identical repetition-specific seeds and statistical functions. The first 500 verified checkpoints are reused; repetitions 501–10000 are newly computed. Primary application remains one 40-imputation analysis. Four arms: control and three active doses.

The uploaded source is preserved. Execution changes: REG only per user preference; publication_run TRUE; scalar summary-export condition and correct SD evaluation order to prevent duplicated rows and a 100-million-row join. All statistical function definitions must match the previous execution exactly. Both proposed procedures share the same imputations.

Local checkpoints are sealed every 10 repetitions. Remote checkpoint snapshots are saved every 500 new repetitions when repository access works. On completion the controller builds PDF/CSV/plots and publishes the results separately. Environment shutdown stops computation; snapshots support resuming.

[Full PDF report](simulation_report.pdf) · [Primary table PDF](main_tables.pdf) · [Forest plot](comparison_graph.png) · [Results ZIP](Litter_Regression_Application_Results.zip)

[Full-precision primary results](computed_pooled_results_full_precision.csv) · [Original output tables](Litter_Application_Final).
