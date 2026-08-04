test_that("capnet_violations detects known cap violations from a fitted object", {
  set.seed(1)
  n <- 60; p <- 4
  X <- matrix(rnorm(n * p), n, p)
  beta <- c(3, -3, 0, 0)
  y <- as.numeric(X %*% beta + rnorm(n))

  # gamma = 0 (uncapped) plus a very tight L guarantees violations exist.
  fit <- capnet(X, y, L = 0.05, alpha = 0, lambda = 0, gamma = 0)

  v <- capnet_violations(fit)

  expect_s3_class(v, "capnet_violations")
  expect_true(length(v$columns) > 0)
  expect_true(length(v$rows) > 0)
  expect_equal(dim(v$excess), dim(fit$feature_contributions))
  expect_true(all(v$excess >= 0))
  expect_true(all(v$excess[, v$columns] > 0 | v$excess[, v$columns] == 0))
})

test_that("capnet_violations reports no violations when L is generous", {
  set.seed(2)
  n <- 60; p <- 4
  X <- matrix(rnorm(n * p), n, p)
  beta <- c(0.01, -0.01, 0, 0)
  y <- as.numeric(X %*% beta + rnorm(n))

  fit <- capnet(X, y, L = 1000, alpha = 0, lambda = 0, gamma = 0)

  expect_message(result <- capnet_violations(fit), "No cap violations detected")
  expect_null(result)
})

test_that("capnet_violations evaluates a fresh z matrix when supplied", {
  set.seed(3)
  n <- 50; p <- 3; m <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(X %*% c(2, -2, 0) + rnorm(n))
  z <- matrix(rnorm(m * p) * 5, m, p)  # deliberately larger-scale evaluation rows

  fit <- capnet(X, y, L = 0.1, alpha = 0, lambda = 0, gamma = 0)

  v <- capnet_violations(fit, z = z)
  expect_equal(dim(v$excess), c(m, p))
})
