# CAPNET
*A Contribution-Capped Elastic Net for Robust and Interpretable Prediction

## Why
### Quantitative Finance Use Case
In financial modeling, predictive features often exhibit unstable variance and regime-dependent scaling--whether derived from volatility measures, fundamentals, technical ratios, or alternative data. These shifts often reflect real changes in the market structure rather than mere noise.

Traditional regularization methods with L1 and L2 penalties constrain coefficients, not contributions, and therefore cannot directly control how much any single feature drives model predictions. In fact, the coefficient-shrinking behavior of these penalties can inadvertently magnify the dominance of the most predictive or volatile features. When one variable carries disproportionately large variance--especially when such cases are underrepresented in the training sample--it can dominate the prediction output, leading to undesired idiosyncratic exposure. Conceptually, this mirrors porftolio construction: just as excess reliance on a single asset increases portfolio risk, excess reliance on one predictive feature increases model risk.

Other approaches to this problem include winsorization and outlier filtering, but these methods can unintentionally suppress valuable information. Extreme observations, absent measurement errors, often capture regime transitions or structural signals whose magnitudes carry predictive meaning. Moreover, when feature distributions are unknown or non-stationary, fixed winsorization thresholds offer no theoretical guarantee of mitigating single-feature dominance. In practice, experiments I conducted have shown that manually tuning these thresholds consistently performan worse than a principled contribution-based regularization implemented with **capnet**.

**capnet** introduces a new regularization framework that constrains per-feature contributions rather than raw coefficients. This allows a model to learn from the training data while ensuring that no single feature in the evaluation set disproportionately drives predictions, achieving a more robust balance between information retention and risk control.

Mathematically, capnet augments the elastic net objective with contribution caps that impose a tunable soft ceiling on each variable's influence on fitted values. The resulting estimator adapts naturally to heteroskedastic and correlated feature spaces, making it particularly effective for:
 - Cross-sectional return prediction under shifting feature volatility
 - Portfolio construction pipelines seeking to avoid concentration from single predictive factors

Developed in collaboration with Hull Tactical Asset Allocation, **capnet** has demonstrated measurable improvements in signal stability and generalization when implemented in medium-frequency market-timing models. It has also demonstrated superior performances when compared to other approaches such as winsorization and outlier filtering.

### General Use Cases
Originally developed for financial modeling, **capnet** generalizes to any predictive setting where a small number of features can disproportionately influence outcomes.

By explicitly constraining feature contributions, **capnet** promotes fairer and more balanced models, ensuring that predictive power is not overly concentrated in a handful of variables. This property is valuable in contexts such as:
 - **Social Science**: where hterogeneous covariates can dominate due to skewed distributions
 - **Health and Biomedical Modeling**: where rare conditions or biomarkers have disproportionate effects but must be handled without arbitrary truncation
 - **General Statistical Learning**: where observed features are imperfect or noisy realizations of underlying causal variables with high and unstable variance

Across these settings, capnet offers a principled alternative to ad hoc preprocessing methods by controlling feature influence directly within the optimization framework.

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
