test_that("coef_path traces a lambda path and returns a tidy capnet_path data.frame", {
  set.seed(1)
  n <- 50; p <- 6
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", 1:p)
  beta <- c(1.5, 0.8, 0.2, rep(0, p - 3))
  y <- as.numeric(X %*% beta + rnorm(n))

  lambda <- exp(seq(1, -3, length.out = 10))
  path <- coef_path(X, y, L = 0.5, alpha = 0.5, lambda = lambda, gamma = 1)

  expect_s3_class(path, "capnet_path")
  expect_equal(nrow(path), length(lambda))
  expect_equal(colnames(path), c("lambda (log)", colnames(X)))
  expect_equal(path[["lambda (log)"]], log(lambda))
})

test_that("coef_path traces an alpha path when alpha has length > 1", {
  set.seed(2)
  n <- 40; p <- 4
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(X %*% c(1, -1, 0, 0) + rnorm(n))

  alpha <- seq(0, 1, length.out = 6)
  path <- coef_path(X, y, L = 1, alpha = alpha, lambda = 0.1, gamma = 0.5)

  expect_equal(nrow(path), length(alpha))
  expect_equal(colnames(path)[1], "alpha")
  expect_equal(path[["alpha"]], alpha)
})

test_that("coef_path errors when more than one of lambda/alpha/gamma varies", {
  set.seed(3)
  n <- 30; p <- 3
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(X %*% c(1, 0, -1) + rnorm(n))

  expect_error(
    coef_path(X, y, L = 1, alpha = c(0, 1), lambda = c(0.1, 0.2)),
    "constant value"
  )
})

test_that("coef_path's standardize argument no longer collides with ... when passed explicitly", {
  set.seed(4)
  n <- 30; p <- 3
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(X %*% c(1, 0, -1) + rnorm(n))

  expect_error(
    coef_path(X, y, L = 1, alpha = 0.5, lambda = c(0.1, 0.2), gamma = 1,
              standardize = TRUE),
    NA
  )
})
