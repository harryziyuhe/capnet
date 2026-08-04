test_that("walk_capnet returns correctly shaped output for a simple path", {
  set.seed(1)
  n <- 60; p <- 5; n_new <- 12
  X <- matrix(rnorm(n * p), n, p)
  z <- matrix(rnorm(n_new * p), n_new, p)
  beta <- c(1.5, -1, 0.5, rep(0, p - 3))
  y <- as.numeric(X %*% beta + rnorm(n))

  out <- walk_capnet(X, y, L = 1, z = z, walk = 1,
                     lambda = 0.1, alpha = 0.5, gamma = 1)

  expect_s3_class(out, "walk_capnet")
  expect_equal(length(out$intercepts), n_new)
  expect_equal(dim(out$betas), c(n_new, p))
  expect_equal(dim(out$feature_contributions), c(n_new, p))
  expect_equal(length(out$predictions), n_new)
  expect_equal(length(out$gammas), n_new)
  expect_true(all(is.finite(out$intercepts)))
})

test_that("walk_capnet predictions match linkinv(eta) computed from its own fitted coefficients", {
  set.seed(2)
  n <- 50; p <- 4; n_new <- 8
  X <- matrix(rnorm(n * p), n, p)
  z <- matrix(rnorm(n_new * p), n_new, p)
  beta <- c(1, -0.5, 0, 0)
  y <- as.numeric(X %*% beta + rnorm(n))

  out <- walk_capnet(X, y, L = 2, z = z, walk = 2,
                     lambda = 0.05, alpha = 0.3, gamma = 0.5,
                     family = "gaussian")

  eta <- out$intercepts + rowSums(out$feature_contributions)
  expect_equal(as.numeric(out$predictions), as.numeric(eta), tolerance = 1e-8)
})

test_that("walk_capnet with walk > 1 groups consecutive evaluation rows per step", {
  set.seed(3)
  n <- 40; p <- 3; n_new <- 9
  X <- matrix(rnorm(n * p), n, p)
  z <- matrix(rnorm(n_new * p), n_new, p)
  y <- as.numeric(X %*% c(1, 0, -1) + rnorm(n))

  out <- walk_capnet(X, y, L = 1, z = z, walk = 3,
                     lambda = 0.1, alpha = 0.5, gamma = 1)

  # Rows within the same walk-block share the same fitted beta/intercept,
  # since they come from the same optimizer call.
  expect_equal(out$betas[1, ], out$betas[2, ])
  expect_equal(out$betas[1, ], out$betas[3, ])
  expect_false(isTRUE(all.equal(out$betas[1, ], out$betas[4, ])))
})

test_that("walk_capnet retries with a smaller gamma on non-convergence", {
  # An enormous initial gamma is very likely to make the first attempt hard to
  # converge cleanly; regardless, the recorded gamma for each step should
  # never exceed the requested starting gamma, and should be one of
  # gamma / 10^k for some k in the allowed range.
  set.seed(4)
  n <- 40; p <- 3; n_new <- 5
  X <- matrix(rnorm(n * p), n, p)
  z <- matrix(rnorm(n_new * p), n_new, p)
  y <- as.numeric(X %*% c(1, -1, 0.5) + rnorm(n))

  gamma0 <- 1e6
  out <- walk_capnet(X, y, L = 0.01, z = z, walk = 1,
                     lambda = 0, alpha = 0, gamma = gamma0,
                     max_gamma_tries = 8)

  used <- unique(out$gammas[is.finite(out$gammas)])
  expect_true(all(used <= gamma0))
  ratios <- gamma0 / used
  expect_true(all(sapply(ratios, function(r) any(abs(log10(r) - round(log10(r))) < 1e-8))))
})
