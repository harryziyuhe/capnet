test_that("uncapped gaussian capnet agrees with ordinary least squares", {
  set.seed(11)
  n <- 120
  p <- 5
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  beta <- c(1.5, -0.8, 0.4, 0, 0.2)
  y <- 0.7 + drop(X %*% beta) + rnorm(n, sd = 0.5)

  fit <- capnet(
    X, y, L = 1, lambda = 0, alpha = 0, gamma = 0,
    standardize = FALSE
  )
  baseline <- stats::lm(y ~ X)

  expect_equal(
    c(fit$a0, fit$beta), unname(stats::coef(baseline)),
    tolerance = 1e-6
  )
  expect_equal(
    as.numeric(predict(fit, newdata = X)),
    as.numeric(stats::predict(baseline)),
    tolerance = 1e-6
  )
})

test_that("uncapped gaussian elastic net agrees with glmnet", {
  skip_if_not_installed("glmnet")
  set.seed(12)
  n <- 160
  p <- 7
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  beta <- c(1.2, -0.9, 0.5, 0, 0, 0.3, 0)
  y <- -0.4 + drop(X %*% beta) + rnorm(n, sd = 0.7)
  alpha <- 0.8
  lambda <- 0.03

  fit <- capnet(
    X, y, L = 1, lambda = lambda, alpha = alpha, gamma = 0,
    standardize = FALSE
  )
  baseline <- glmnet::glmnet(
    X, y, alpha = alpha, lambda = lambda,
    standardize = FALSE, intercept = TRUE, thresh = 1e-12
  )
  baseline_coef <- as.matrix(stats::coef(baseline))[, 1]

  # Allow modest coefficient and prediction differences between OWL-QN and
  # glmnet's coordinate-descent implementation and Gaussian scaling choices.
  expect_equal(
    c(fit$a0, fit$beta), unname(baseline_coef),
    tolerance = 5e-3
  )
  expect_equal(
    as.numeric(predict(fit, newdata = X)),
    as.numeric(predict(baseline, newx = X)),
    tolerance = 2e-2
  )
})
