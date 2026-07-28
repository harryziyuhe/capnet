test_that("analytic likelihood gradient matches numerical gradient for gaussian", {
  set.seed(1)
  d <- simulate_capnet_data(n = 150, p = 6, family = "gaussian")
  result <- .loss_gradient_check(d$X, d$y, "gaussian")
  expect_lt(result$max_diff, 1e-6)
})

test_that("analytic likelihood gradient matches numerical gradient for binomial", {
  set.seed(1)
  d <- simulate_capnet_data(n = 150, p = 6, family = "binomial")
  result <- .loss_gradient_check(d$X, d$y, "binomial")
  expect_lt(result$max_diff, 1e-6)
})

test_that("analytic likelihood gradient matches numerical gradient for poisson", {
  set.seed(1)
  d <- simulate_capnet_data(n = 150, p = 6, family = "poisson")
  result <- .loss_gradient_check(d$X, d$y, "poisson")
  expect_lt(result$max_diff, 1e-6)
})

test_that("analytic likelihood gradient matches numerical gradient for gamma", {
  set.seed(1)
  d <- simulate_capnet_data(n = 150, p = 6, family = "gamma")
  result <- .loss_gradient_check(d$X, d$y, "gamma")
  expect_lt(result$max_diff, 1e-6)
})

test_that("likelihood gradient check is repeatable across random beta draws", {
  # .loss_gradient_check() draws a random beta when none is supplied; run it
  # a few times per family so a fluke beta near a kink/asymptote can't mask a
  # real gradient bug.
  for (fam in c("gaussian", "binomial", "poisson", "gamma")) {
    set.seed(2)
    d <- simulate_capnet_data(n = 150, p = 6, family = fam)
    for (i in 1:5) {
      result <- .loss_gradient_check(d$X, d$y, fam)
      expect_lt(result$max_diff, 1e-6, label = sprintf("%s draw %d", fam, i))
    }
  }
})

test_that("contribution-cap gradient (capnet_fit.R) is scale-invariant across mismatched-scale columns", {
  # Regression test for the capnet_fit.R bug where the cap-penalty gradient
  # was multiplied by scaling_factor instead of divided by it, causing the
  # gradient to blow up as X_scale^2 for large-raw-scale covariates instead
  # of staying O(1) (contributions Z %*% beta are scale-invariant by
  # construction, so their gradient should be too). See capnet_fit.R for the
  # fixed chain-rule derivation.
  set.seed(1)
  n <- 200; p <- 3
  X <- cbind(rnorm(n, sd = 1), rnorm(n, sd = 10), rnorm(n, sd = 100))
  beta_true <- c(0.01, 0.001, 0.0001)
  eta <- as.numeric(X %*% beta_true)
  y <- rpois(n, exp(eta))

  for (g in c(0, 0.01, 1, 10, 100)) {
    fit <- capnet(X, y, L = rep(0.05, p), family = "poisson",
                  lambda = 0, alpha = 0, gamma = g)
    expect_equal(fit$convergence, 0, info = sprintf("gamma = %s", g))
  }

  # At a nontrivial gamma, the cap penalty should meaningfully shrink the
  # largest-scale column's contribution rather than leaving it unaffected
  # (which is what happens if the scale-dependent gradient bug regresses).
  fit_uncapped <- capnet(X, y, L = rep(0.05, p), family = "poisson",
                         lambda = 0, alpha = 0, gamma = 0)
  fit_capped <- capnet(X, y, L = rep(0.05, p), family = "poisson",
                       lambda = 0, alpha = 0, gamma = 100)
  expect_lt(
    max(abs(fit_capped$feature_contributions)),
    max(abs(fit_uncapped$feature_contributions))
  )
})
