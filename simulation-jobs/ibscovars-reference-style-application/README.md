# IBScovars reference-style application

Runs the uploaded IBScovars_ReferenceStyle_Application.R as supplied: one application analysis with **40 propensity-score imputations**, target 30% empirical covariate-driven MAR (gamma=log(2), standardized gender covariate). The realized missingness is **118/369 = 31.98%**. No repeated-simulation loop or imputation-method comparison was added.

The primary table contains the supplied script's three procedures: proposed with shrinkage, proposed without shrinkage, and naive case deletion. Both proposed procedures share the same imputed datasets and point estimates. Original scientific calculations, seeds and settings were retained.

The combined PDF contains the reference-style table and forest plot. The forest plot shows 40 individual completed-dataset results, pooled shrinkage/non-shrinkage, and case deletion in four dose-comparison panels. Its readable copy has a white background and greater height; values and plot data are unchanged. Original figure outputs are preserved in IBScovars_ReferenceStyle_Output.

The original LaTeX files compiled successfully. Use latex-ready/Compile_Application.tex for a standalone manuscript example; latex-ready contains unchanged original table/figure include code and a copy of the readable figure. Compile_Application.pdf is supplied. No scientific calculations were changed to prepare these files.

Full-precision pooled and case-deletion results, all 160 individual-imputation intervals, 172 forest-plot rows, input data, 40-imputation RDS snapshot, logs, warnings/errors, session information and original source are included. The console's initial unused socket-probe message did not affect the serial R analysis; the current serial wrapper omits that unnecessary probe.
