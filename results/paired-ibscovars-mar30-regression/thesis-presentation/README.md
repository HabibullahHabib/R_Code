# IBScovars: thesis-style presentation

Reformats the existing completed results without rerunning or changing any analysis.

- [PDF tables](Thesis_Style_Results.pdf): page 1 gives descriptive averages across all 10,000 repetitions; page 2 gives actual simultaneous 95% intervals from the first prespecified 40-imputation repetition, matching the application interpretation of thesis Table 4.2.
- [10,000-repetition averages CSV](thesis_style_10000_averages.csv). Mean lower/upper bounds are averages of separate intervals; they are **not** a single 95% confidence interval.
- [Actual 40-imputation example CSV](thesis_style_40_imputation_example.csv), repetition 1, not selected for its results.

Only shrinkage and non-shrinkage are shown. Regression imputation is fixed throughout. A new 30% expected MAR pattern is imposed in each repetition on the complete dataset.

Original source: ../raw_intervals.csv; SHA256 6323b7d092c381278d5df2929b2d671784892395ac5a4767ce1aea9fb69c3e99. All original analyses, figures and logs remain unchanged.
