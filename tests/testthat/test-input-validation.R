test_that("capnet stops on malformed core inputs", {
  set.seed(21)
  X <- matrix(rnorm(30), nrow = 10, ncol = 3)
  y <- rnorm(10)
  z <- matrix(rnorm(12), nrow = 4, ncol = 3)
  fit_args <- list(lambda = 0, gamma = 0, standardize = FALSE)

  expect_error(
    do.call(capnet, c(list(X = replace(X, 1, NA_real_), y = y, L = 1), fit_args)),
    "contains NA"
  )
  expect_error(
    do.call(capnet, c(list(X = X, y = y[-1], L = 1), fit_args)),
    "y must have length nrow\\(X\\)"
  )
  expect_error(
    do.call(capnet, c(list(X = X, y = y, L = 1, z = z[, 1:2]), fit_args)),
    "same number of columns"
  )
  expect_error(
    do.call(capnet, c(
      list(X = X, y = y, L = 1, z = z, multiplier = c(1, 2)),
      fit_args
    )),
    "multiplier must be length 1 or length nrow\\(z\\)"
  )
  expect_error(
    do.call(capnet, c(list(X = X, y = y, L = NULL), fit_args)),
    "L must be provided"
  )
  expect_error(
    do.call(capnet, c(list(X = X, y = y, L = c(1, 2)), fit_args)),
    "L must be length 1 or length ncol\\(X\\)"
  )
  expect_error(
    do.call(capnet, c(list(X = X, y = y, L = 1, par = numeric(2)), fit_args)),
    "par must be length 1 \\+ ncol\\(X\\)"
  )
})

test_that("capnet rejects unsupported family links", {
  X <- matrix(seq_len(20), nrow = 10, ncol = 2)
  y <- seq_len(10)

  expect_error(
    capnet(X, y, L = 1, family = stats::gaussian(link = "log")),
    "Unsupported gaussian link 'log'"
  )
})

test_that("capnet_violations warns when a varying multiplier is omitted", {
  set.seed(22)
  X <- matrix(rnorm(60), nrow = 20, ncol = 3)
  y <- rnorm(20)
  z <- matrix(rnorm(12), nrow = 4, ncol = 3)
  fit <- capnet(
    X, y, L = 0.1, z = z, multiplier = seq(0.5, 2, length.out = 4),
    lambda = 0, gamma = 0, standardize = FALSE
  )

  expect_warning(
    result <- capnet_violations(fit, z = z),
    "Did not supply a multiplier"
  )
  expect_s3_class(result, "capnet_violations")
})
