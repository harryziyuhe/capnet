test_that("capnet respects coefficient lower and upper bounds", {
  set.seed(31)
  n <- 120
  p <- 3
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  y <- 0.5 + drop(X %*% c(-2, 2, 1)) + rnorm(n, sd = 0.1)
  lower <- c(0, -Inf, -0.25)
  upper <- c(Inf, 0, 0.25)

  fit <- capnet(
    X, y, L = 1, lambda = 0, gamma = 0,
    lower.limits = lower, upper.limits = upper
  )

  expect_true(all(fit$beta >= lower - 1e-8))
  expect_true(all(fit$beta <= upper + 1e-8))
  expect_equal(fit$beta[1], lower[1], tolerance = 1e-5)
  expect_equal(fit$beta[2], upper[2], tolerance = 1e-5)
  expect_equal(fit$beta[3], upper[3], tolerance = 1e-5)
})

test_that("bounded L1 fits apply soft-thresholding and respect bounds", {
  set.seed(32)
  n <- 80
  p <- 3
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  X <- sweep(X, 2, colMeans(X), "-")
  X <- sqrt(n) * qr.Q(qr(X))
  beta_true <- c(1.2, 0.1, -1)
  y <- 0.7 + drop(X %*% beta_true)

  fit <- capnet(
    X, y, L = 1, lambda = 0.2, alpha = 1, gamma = 0,
    lower.limits = c(0, -Inf, -0.6),
    upper.limits = c(0.85, Inf, Inf),
    standardize = FALSE,
    par = c(0, 0, 0.2, 0.4)
  )

  # X is centered and has X'X / n = I, so the unconstrained lasso solution is
  # soft-threshold(beta_true, lambda). The box then clips coefficients 1 and 3.
  expect_equal(fit$convergence, 0)
  expect_equal(fit$a0, 0.7, tolerance = 1e-5)
  expect_equal(fit$beta, c(0.85, 0, -0.6), tolerance = 1e-5)
  expect_true(all(fit$beta >= c(0, -Inf, -0.6) - 1e-8))
  expect_true(all(fit$beta <= c(0.85, Inf, Inf) + 1e-8))
})
