# Verified document methods and user-requested application adaptation

The PDF named **Variance Shrinkage in Dunnett-Type Multiple Comparisons with Miss.pdf** is Md Habibullah's thesis, *Variance Shrinkage in Dunnett-Type Multiple Comparisons with Missing Data* (2026). The PDF named **Thesis Paper.pdf** is the Jiang, Zhou and Liang article, *An Improved Dunnett’s Procedure for Comparing Multiple Treatments with a Control in the Presence of Missing Observations*, Mathematics 2023, 11, 3233, DOI 10.3390/math11143233. File names do not identify the documents reliably; their title pages were checked.

| Document and section | Verified procedure |
|---|---|
| Thesis §4.1, printed pp.43–44, PDF pp.51–52 | Real-data application: regression imputation, 20 completed datasets; equal-variance versus unequal-variance shrinkage intervals. |
| Reference article §6, printed p.12 | Real-data application: propensity-score imputation, 20 completed datasets. |
| Thesis §3.1 and §3.2, printed pp.28 and39 | Simulation studies: 5,000 independent runs, 20 imputations within each incomplete dataset. |
| Reference article §5.1, pp.7–8 | Simulation I: 10,000 independent runs, 40 imputations per incomplete dataset. |
| Reference article §5.2, p.10 | Simulation II: 15,000 simulations; distinct from the single real-data application. |

The original real-data application is a 100-patient, five-arm dose-response dataset. It is not ACTG175 or IBScovars. The corrected uploaded scripts use different data, covariates and an unequal-arm-size extension; those features are retained and identified as adaptations.

After document review, the user explicitly requested **40 imputations**. This overrides the documents' count of 20 and the previous request for 10,000 imputations. No repeated-deletion experiments are added. Primary thesis-style output uses Regression MI and compares common-variance with unequal-variance shrinkage. Propensity MI is reported separately to show the reference article's imputation choice on the new datasets. PMM and unshrunk unequal-variance intervals are retained as supplementary comparisons. The latter is not the reference article's equal-variance comparator.

Common-variance covariance uses within-arm degrees-of-freedom-weighted pooled variance (article §3), followed by Rubin pooling and the trace-based multivariate-t degrees of freedom. Point estimates and imputations are shared across variance versions. Original corrected-script shrinkage and non-shrinkage calculations are preserved. First 40 original seeded imputations may be reused from verified 10,000-imputation checkpoints after checking source/data/parameters, with exact original-function comparisons.
