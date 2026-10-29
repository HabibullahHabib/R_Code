# Professor feedback simulation jobs

These are source scripts and resumable runners, **not completed results**. Each original script requests 600,000 jobs: 10,000 replicates for each variance-ratio × MAR-mechanism × scenario × missing-rate setting; 40 imputations.

The equal-variance simulation is paused at the user’s request. The execution queue runs 5 treatments, then 10 treatments, each with four workers, and automatically resumes the equal-variance simulation afterward. Short trial runs and exact comparisons against the original scientific functions gate full execution. Checkpoints and automatic local backups are enabled. Running processes and checkpoints are not stored in this source-only Git commit. An environment shutdown stops its processes.

Original seeds, condition levels and statistical functions are preserved. Instrumentation adds atomic 20-job checkpoints, console/warning/error capture, lossless gzip CSV export, CSV/HTML copies of the original tables, and PDF tables and diagnostic plots. Full-run results will be published separately under results/paper-feedback-5-treatments and results/paper-feedback-10-treatments after completion.

To restore on a compatible machine, place each run folder under /workspace, restore the prepared R runtime and packages and available checkpoint backups, then run ./run.sh full. The runtime requires R, mvtnorm, dplyr, tibble, gt and parallel; PDF merging uses pdfunite. The queue helper resumes the paused equal-variance process when available, or restarts it from saved checkpoints.
