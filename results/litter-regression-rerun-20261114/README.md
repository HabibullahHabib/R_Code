# Litter regression application: fresh rerun

The newly uploaded script is byte-for-byte identical to the previous Litter script. Original source is preserved. Regression only, 40 imputations, 500 repeated-MAR validation repetitions as specified by the uploaded default. Four groups: control and three active doses; this does not meet the separate request for placebo plus four active arms.

Both proposed procedures use the same imputations and retain unequal group variances before shrinkage. The existing summary-export bug is repaired from unchanged raw validation results; the original exports and warnings are retained. Highlighted forest rows show Rubin-combined with-shrinkage and without-shrinkage results.

The first launcher attempt completed the primary analysis but could not create local worker sockets under the initial network sandbox. Its console/error records are retained under launch_attempt_1. The run was restarted after socket access was enabled.

[Full PDF report](simulation_report.pdf) · [Primary table PDF](main_tables.pdf) · [Forest plot](comparison_graph.png) · [Results ZIP](Litter_Regression_Application_Results.zip)

[Full-precision primary results](computed_pooled_results_full_precision.csv) · [Original output tables](Litter_Application_Final).
