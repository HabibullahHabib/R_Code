# application-ibscovars-thesis-aligned: 40-imputation real-data adaptation

The user requested 40 imputations after review of the uploaded thesis and reference article. Both documents use 20 for their original real-data application; this count difference is explicit. Original dataset applications differ from ACTG175 and IBScovars, so these results are adaptations, not exact replications.

[Full PDF](simulation_report.pdf) · [Tables PDF](main_tables.pdf) · [Plots PDF](simulation_plots.pdf) · [Results ZIP](application-ibscovars-thesis-aligned_Tables_Plots_and_Logs.zip)

[Thesis-primary Regression MI comparison](THESIS_PRIMARY_regression_equal_vs_shrinkage.csv): equal-variance versus unequal-variance shrinkage intervals. [Reference-method Propensity MI comparison](REFERENCE_METHOD_propensity_equal_vs_shrinkage.csv) is reported separately. PMM and the original unshrunk unequal-variance comparison remain supplementary. The unshrunk unequal-variance procedure is not the common-variance reference baseline.

See METHOD_REVIEW.md and METHOD_ALIGNMENT.txt for source sections, formulas, dataset differences and user overrides. No repeated-deletion simulations were added. Both variance constructions share the same 40 imputations and point estimates. Existing verified seeded imputation prefixes were reused when available, retaining original data, seeds and imputation functions. Session information, logs, warnings/errors, detailed CSVs and completion counts are included.
