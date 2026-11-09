# Corrected IBScovars application: regression imputation

**This is the corrected result for the user's regression-only instruction.** The earlier propensity-default application is superseded for this request and retained for audit.

The executed script changes exactly one analysis parameter in the uploaded source: imputation_method is regression. All 40 imputations, patient data, the exact 118/369 missing-value locations, MAR settings and seeds are retained. Same-patient/missingness input equality was independently verified. This is one real-data application analysis, not a new repeated-simulation experiment.

The reference article's Figure 1 (p.13) highlights the Pooled row: the final MI estimate and interval combining imputations with Rubin's rules. It does not identify a selected individual imputation or imply pooled/equal treatment variances. Our figure labels the two highlighted rows directly as **Proposed with shrinkage** and **Proposed without shrinkage**. Orange diamonds and blue squares identify these procedures; grey numbered rows are individual regression imputations; black triangles are the case-deletion comparator included in the source script. Colour does not indicate statistical significance or superiority.

Both proposed procedures share the same 40 imputations and point estimates. No equal-variance comparison or imputation-method comparison was added. Original source figure files are retained in IBScovars_ReferenceStyle_Output; the report's highlighted plot is derived from the same 172 plot rows without changing estimates or intervals. The regression LaTeX files compiled successfully.

Use executed_script.R to reproduce the corrected analysis. original_script.R preserves the uploaded source, including its conflicting propensity default, for traceability only. The wrapper, full-precision CSVs, individual intervals, input snapshots, correction metadata, logs, source and compiled LaTeX are included.
