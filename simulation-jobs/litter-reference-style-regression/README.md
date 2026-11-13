# Litter application: regression imputation only

Runs the uploaded Litter_Final_Application_Shrinkage_vs_NonShrinkage.R with the user's standing regression-only instruction. The only analysis-parameter change is all_methods <- c("REG") instead of c("REG", "PMM", "PS"). The primary REG seeds, formulas, missingness model, primary application, 40 imputations and publication_run=FALSE are preserved. The supplied default validation therefore requests 500 repetitions, not 10,000.

The dataset is multcomp::litter: doses 0, 5, 50, 500; litter weight is the outcome, standardized litter size (number) is the MAR covariate. Target expected missingness is 30% with gamma=log(2), same across arms. The two proposed procedures share identical imputed datasets and point estimates.

The primary reference-style table shows shrinkage, non-shrinkage and the case-deletion comparator included in the supplied script. The figure retains individual-imputation and complete-data benchmark rows as supplied, while identifying the two proposed procedures directly: blue square = without shrinkage; orange diamond = with shrinkage. Highlighted results combine all 40 regression imputations with Rubin's rules; this is not an equal-variance assumption and colour does not indicate superiority.

Validation uses new MAR patterns on the fixed complete dataset. Its benchmark recovery is not population coverage. The original safe validation function and NULL/omission behaviour are retained. Parallel execution adds four workers and atomic 10-repetition checkpoints, without changing the REG calculations; repetition 1 is checked against the original serial function. Failed/omitted validation IDs and actual successful counts are retained. Only regression is run.

CSV tables, complete validation rows, original figures, highlighted figure, combined PDF, logs, warning/error capture, input data, session information and a primary-imputation RDS snapshot are retained. The data/results for other applications are separate.

## Validation export repair

The supplied summary used a vector-valued ifelse inside summarise, duplicating procedure rows and causing a many-to-many join. It also computed a standard deviation after masking the original column with its mean, giving NA. repair_validation_summary.R regenerates the two-row summary and one-row comparison from unchanged raw results using distinct per-repetition aliases. Original as-generated exports and the two corresponding warnings are preserved. No primary estimates, imputations, intervals or repetition-level results were changed. Both procedures have 500 successful validation repetitions.
