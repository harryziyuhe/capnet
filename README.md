# Feature Contribution Cap Elastic Net Model (CAPNET)

**capnet** is an R package for fitting regularized linear models with feature contribution caps. It uses the OWL-QN optimizer to minimize a custom penalized objective function and supports user-defined per-feature contribution caps.

## Installation

You can install the development version from GitHub:

```r
remotes::install_github("harryziyuhe/capnet")
```

## Features
- Optimized with OWL-QN for smooth convergence
- Differentiable customized contribution cap penalty
- Built in cross-validation and walk-forward tools
- Allows specified upper and lower limit for each feature
- Can be used as a drop-in placeholder for `glmnet`

## Example Usage
```r
library(capnet)

set.seed(123)
X <- matrix(rnorm(100 * 5), 100, 5)
y <- rnorm(100)
caps <- c(0.5, 1, 0.2, 0.3, 0.6)

fit <- capnet(X, y, lambda = 0.1, alpha = 0.5, mu = 1, L = caps)

summary(fit)
```

## Methodology

The model solves the following penalized optimization problem:

```math
\arg\min_{\beta}\|Y-X\beta\|_2^2+\lambda\left(\alpha\|\beta\|_1+(1-\alpha)\|\beta\|_2^2\right)+\mu\sum_{i=1}^k\left\|\max\left(0,|\beta_iZ_i|-L\right)\right\|_2^2
```

Where:
- $k$ is the number of covariates and $L$ is the contribution cap vector.
- $\mu$ is a user-defined regularization strength
- The penalty is differentiable to support gradient based optimization via OWL-QN 

## Function Reference

| Function                | Description
|-------------------------|------------------------------
| `capnet()`              | Fit a linear model with contribution caps
| `cv.capnet()`           | Cross-validation to select lambda/alpha/mu value
| `capnet.walk()`         | Walk forward to fit contribution caps over multi-row evaluation matrix
| `coef_path()`           | Fit models along one hyperparameter path and return the coefficients
| `plot_redistribution()` | Plot redistribution of feature coefficients or contributions
| `capnet.violations()`   | Check contribution violations

## Dependencies

- R >= 4.1.0
- lbfgs, ggplot2

## License

MIT License © 2025 Harry He
