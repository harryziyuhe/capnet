# Reproducing the paper experiments

Run `bash tools/reproduce_paper.sh` from the repository root. It installs the
checkout into an isolated temporary R library, sets requested BLAS/OpenMP thread
limits to one, runs the experiments, regenerates the manuscript tables and
application figures, renders the three experiment reports, and rebuilds
`manuscript/main.pdf`. It requires R and the package's suggested vignette
packages, Pandoc, and a LaTeX installation. Dependencies are not downloaded by
the script.

The source of each experiment is in `vignettes/paper/experiments.R` and
`vignettes/paper/simulation.R`; the reports are `example_simulation.Rmd`,
`example_elan.Rmd`, and `example_wiki.Rmd`. By default the reports read the latest
saved run, rather than performing a second fit just to typeset it. Render with
`params = list(rerun = TRUE)` to recompute an individual report's fits. Use the
shell command above to regenerate all manuscript artifacts together.

For individual computations after installing the checkout, use
`Rscript tools/run_paper_experiments.R TASK`, where TASK is `simulation`,
`finance`, `wiki`, `runtime`, or `export`. Set `R_LIBS_USER` to the isolated
library when doing this manually. Runtime should be run when the other
experiments have finished; the all-task runner does this sequentially.

## Saved outputs

- `simulation.rds`: configuration, per-replication losses and convergence
  diagnostics, mean/Monte Carlo standard error, and paired-difference summaries.
  The two simulation CSVs expose the same data without requiring R.
- `finance.rds`, `wiki.rds`: complete predictions, feature contributions,
  monthly lambda selections, requested/effective penalties, and numerical
  diagnostics. A walk-forward retry may reduce gamma; these reductions are
  recorded rather than hidden.
- `*-scores.csv`: pooled squared-error out-of-sample R-squared for every
  specification. `wiki-2024-scores.csv` restricts the saved predictions to 2024
  and resets the expanding benchmark within that subperiod.
- `*-maximum-contributions.csv`: the exact nonnegative row-wise maxima used in
  the corrected plots, computed independently for each model.
- `*-snapshot.csv`, `*-moments.csv`, `*-plot-summary.csv`: snapshot contributions,
  descriptive statistics, and the chosen illustrative specifications.
- `runtime.rds`, `runtime-replications.csv`, `runtime-summary.csv`: three serial
  full-backtest timings per method, with an active gamma of 0.1 for both joint
  and walk-forward capnet fits.
- `*-session.txt`: R/package versions, platform, requested thread settings,
  run completion time, and input/source MD5 hashes.

`export_results.R` generates the LaTeX tabular fragments in `manuscript/tables`
and the four application figures in `manuscript/figs`. The manuscript includes
those fragments directly, avoiding manual copies of numerical results.

## Interpretation

Tables 1 and 2 use a new, declared 10-replication experiment, not a recovered
version of the earlier undocumented simulation. It has five nonzero signal
coefficients; features 1 and 2 have proxy-noise standard deviation 10 with
spike probability 0.10 in training and 0.20 in the test batch. The mixture is
an explicit test-time stress shift. The homoskedastic control and all penalty
configurations are retained, even where caps worsen error. Grid winners are
not independently validated deployment choices.

The finance protocol preserves the original vignette's 25-year windows,
correlation ceilings, random inner CV, and one-row walk steps. Proprietary
feature release times are not established by the supplied CSV.

The web protocol preserves the original five-row step for speed and treats it
as a week-like **conditional batch**. It conditions on all five daily covariate
rows, retains daily next-day targets and scoring, restarts slices each month,
and uses no evaluation responses during fitting. It is not a five-day-ahead
forecast based solely on the first day's information. Its rolling score uses
180 observations, summed directly within each window to avoid cumulative-sum
cancellation after extreme errors. Poisson is a working mean-model loss; the exported residual
summaries do not assert a calibrated count distribution.
