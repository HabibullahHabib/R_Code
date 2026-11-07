# Application studies: 1,000 experiments per script

ACTG175: retain the original natural-missingness primary clinical analysis, then run 1,000 repeated 30% MCAR deletion experiments on the observed-outcome benchmark subset. IBScovars: retain the complete-data benchmark and original main applications, then run 250 repeated experiments at each of 10%, 20%, 30%, 40% missingness (1,000 total). All three MI methods and 40 imputations are retained. Both variance methods share imputations.

Original uploaded sources and scientific functions are preserved. The runner overrides validation_reps only, checkpoints every five experiments and uses four local R workers for repeated validation. Independent checks showed exact original-serial versus worker-output equality. Reports include CSV tables, raw results, logs, warnings/errors, PDF tables and plots. Recovery of complete-sample benchmark contrasts is not population-parameter coverage. These files are setup, not final results; final folders will be results/application-actg175-1000 and results/application-ibscovars-1000. No automatic increase in repetitions is enabled.

Runtime data packages: BART 2.9.10 (official CRAN source mirror) and DoseFinding 1.3-1 (Debian trixie), with ggplot2 3.5.1, mvtnorm, dplyr, tidyr and tibble. Prepared runtime and immutable checkpoints are included in ongoing local backup archives.
