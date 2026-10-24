# Simulation II: ten treatments

Completed 40,000 attempted jobs: **40000 successful**, **0 failed**. Each of four N/scenario settings used 10,000 Monte Carlo replicates, 40 imputations and 30% expected covariate-driven MAR missingness. N=40 and N=80 per arm, ten treatments and one control. The original statistical script, settings and seeds were preserved. See completion_summary.csv and warning/error logs for actual denominators and diagnostics.

- [Download results ZIP](Simulation_II_10_Treatments_Final_Results.zip): open the file and choose **Download raw file** on GitHub.
- [Full PDF report](simulation_report.pdf), [primary tables PDF](main_tables.pdf), [plots PDF](simulation_plots.pdf).
- Original primary CSV tables: [shrinkage](paper_sim2_10treat_shrinkage_N40_N80_table.csv), [non-shrinkage](paper_sim2_10treat_nonshrinkage_N40_N80_table.csv).
- [Numeric results](numerical_summary.csv), [replicate rejection indicators](raw_replicates.csv), [completion counts](completion_summary.csv).

Table cells for MI methods use **PS / PMM / RE** (propensity / predictive mean matching / regression). Scenario I reports joint family-wise error rate (at least one rejection); Scenario II reports individual power and joint power (all ten nulls rejected). Supplementary figures were generated from results; the supplied script defines tables only. The ZIP also contains the original script and console output. SHA256SUMS.txt verifies the ZIP.
