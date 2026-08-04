# Numerically verifies the analytic gradient returned by make_likelihood()
# against a central-difference approximation of its own loss function.
# Returns a structured result rather than warning(), so callers (tests) can
# assert on it directly.
.loss_gradient_check <- function(X, y, family, betas = NULL, eps = 1e-6) {
  family <- normalize_family(family)

  ll <- make_likelihood(family, X, y)
  n <- nrow(X)
  p <- ncol(X)

  if (is.null(betas)) {
    betas <- rnorm(p + 1)
  }
  beta_0 <- betas[1]
  beta <- betas[-1]

  loss <- ll$loss
  grad <- ll$gradient

  b0_loss <- (loss(beta_0 + eps, beta) - loss(beta_0 - eps, beta)) / (2 * eps)
  b0_grad <- grad(beta_0, beta)$g0
  intercept_diff <- abs(b0_loss - b0_grad)

  coef_diffs <- numeric(p)
  for (i in seq_len(p)) {
    beta_up <- beta_down <- beta
    beta_up[i] <- beta_up[i] + eps
    beta_down[i] <- beta_down[i] - eps
    beta_loss <- (loss(beta_0, beta_up) - loss(beta_0, beta_down)) / (2 * eps)
    beta_grad <- grad(beta_0, beta)$gb[i]
    coef_diffs[i] <- abs(beta_loss - beta_grad)
  }

  list(
    intercept_diff = intercept_diff,
    coef_diffs = coef_diffs,
    max_diff = max(intercept_diff, coef_diffs)
  )
}

# Numerically verifies the analytic gradient of the *full* composite objective
# built by .capnet_objective() (likelihood + ridge + contribution-cap penalty,
# including the standardization chain-rule rescaling) against a
# central-difference approximation of its own loss function, at a specific
# beta supplied by the caller. Unlike .loss_gradient_check(), this exercises
# the whole closure .capnet_fit() actually hands to lbfgs(), so it can catch
# bugs in the ridge/cap-penalty/chain-rule wiring that a likelihood-only check
# cannot.
.objective_gradient_check <- function(train, cap, params, betas, eps = 1e-6) {
  obj <- .capnet_objective(train, cap, params)
  p <- length(betas) - 1L

  beta_0 <- betas[1]
  beta <- betas[-1]

  b0_loss <- (obj$loss(c(beta_0 + eps, beta)) - obj$loss(c(beta_0 - eps, beta))) / (2 * eps)
  b0_grad <- obj$gradient(betas)[1]
  intercept_diff <- abs(b0_loss - b0_grad)

  coef_diffs <- numeric(p)
  for (i in seq_len(p)) {
    beta_up <- beta_down <- beta
    beta_up[i] <- beta_up[i] + eps
    beta_down[i] <- beta_down[i] - eps
    beta_loss <- (obj$loss(c(beta_0, beta_up)) - obj$loss(c(beta_0, beta_down))) / (2 * eps)
    beta_grad <- obj$gradient(betas)[i + 1L]
    coef_diffs[i] <- abs(beta_loss - beta_grad)
  }

  list(
    intercept_diff = intercept_diff,
    coef_diffs = coef_diffs,
    max_diff = max(intercept_diff, coef_diffs)
  )
}
