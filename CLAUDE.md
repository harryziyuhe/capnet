# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`capnet` is an R package implementing contribution-capped elastic net regularization. Unlike standard elastic net which penalizes coefficients, capnet penalizes per-feature contributions (i.e., `beta_i * X_i`), preventing any single feature from dominating predictions. Developed for quantitative finance but general-purpose.

The optimization objective is:
```
min_beta  ||y - Xβ||² + λ(α||β||₁ + (1-α)||β||²) + γ·Σ max(0, |β_i·Z_i| - L_i)²
```

OWL-QN is used when `alpha > 0` (L1 component); L-BFGS otherwise. Both via the `lbfgs` package.

## R package commands

```r
# Install dependencies and the package itself
devtools::install_deps()
devtools::install()

# Regenerate documentation from roxygen2 comments
devtools::document()

# Run R CMD CHECK (build + check)
devtools::check()

# Build source tarball
devtools::build()

# Build vignettes
devtools::build_vignettes()
```

Tests live under `tests/testthat/` (run via `devtools::test()`). `tests/testthat/helper-gradient-check.R` provides `.loss_gradient_check()` for numerical gradient verification against the analytical gradients in `objective_factory.R`; `tests/testthat/test-gradients.R` exercises it across all four supported families and includes a regression test for the `capnet_fit.R` contribution-cap gradient scale-invariance bug.

## Architecture

The public API consists of three estimators and several utilities:

| Public function | Description |
|---|---|
| `capnet()` | Single-fit estimator; returns a `"capnet"` object |
| `cv_capnet()` | Grid search over `alpha` × `lambda` with K-fold CV; returns `"cv_capnet"` |
| `walk_capnet()` | Walk-forward refitting over an evaluation matrix `newx`; returns `"walk_capnet"` |
| `coef_path()` | Fits models along one hyperparameter dimension; returns a `"capnet_path"` data frame |
| `capnet_violations()` | Checks which (row, feature) pairs exceed their cap `L` |
| `plot_redistribution()` | Plots how coefficients/contributions shift as caps tighten |

### Internal pipeline (shared by all three estimators)

Every fit goes through this four-step internal pipeline:

1. **`.capnet_spec()`** (`capnet_prepare.R`) — validates and normalizes all inputs into a canonical `spec` list: coerces `X`, `y`, `L`, `multiplier` to correct shapes; resolves `family` via `normalize_family()` and `validate_family_supported()`; sets default `lower.limits`/`upper.limits`.

2. **`.capnet_standardize_train()`** (`capnet_prepare.R`) — standardizes `X` and `y` on the training rows (if `standardize = TRUE`), stores `X_center` / `X_scale` for later rescaling. Returns a `train` list.

3. **`.capnet_cap_context()`** (`capnet_prepare.R`) — slices `newx` and `multiplier` to the evaluation rows, bundles with `L` into a `cap` list.

4. **`.capnet_fit()`** (`capnet_fit.R`) — builds the composite loss (likelihood + ridge + cap penalty) and gradient, then calls `lbfgs::lbfgs()` with `orthantwise_c = lambda * alpha` for the L1 term. The cap penalty is always computed in the **original (unstandardized) scale** even when the model is fit on standardized data.

5. **`.capnet_output()`** (`capnet_output.R`) — maps standardized coefficients back to the original scale, computes `feature_contributions`, and wraps everything into the S3 class object.

### Loss / gradient factory

`objective_factory.R` exports `make_likelihood(family, X, y)`, which returns a list `{loss, gradient}` closed over `(X, y)`. Supported families:

- `gaussian` — identity link, MSE loss
- `binomial` — logit link, cross-entropy loss (numerically stable via `log1pexp`)
- `poisson` — log link
- `Gamma` — log link only. `normalize_family()` special-cases the string `"gamma"`/`"Gamma"` to build `stats::Gamma(link = "log")` (bypassing `stats::Gamma()`'s default `inverse` link, which does not match this package's Gamma loss/gradient); `validate_family_supported()` accordingly requires `log` link for Gamma.

### S3 methods (`methods.R`)

- `coef.capnet`, `predict.capnet` — standard extraction
- `coef.walk_capnet`, `predict.walk_capnet` — return the stored path/prediction matrices
- `plot.cv_capnet` — heatmap or slice plot of CV errors
- `plot.capnet_path` — coefficient paths along a hyperparameter
- `print.capnet_violations` — sparse matrix display via `Matrix::printSpMatrix`

### Key parameter relationships

- `L` is always in the **original feature-contribution scale**, not the standardized scale. Internally, caps are enforced on `newx * beta_raw` where `beta_raw = beta_standardized / X_scale`.
- `multiplier` scales contributions row-wise during cap enforcement only; it does not affect the likelihood.
- `walk_capnet` retries with `gamma / 10` up to `max_gamma_tries` times if the optimizer fails to converge — useful when `gamma` is too large for a particular evaluation window.
- `cv_capnet` uses `gamma` fixed and searches over `alpha × lambda`; `gamma` and `L` must be set before running CV.
