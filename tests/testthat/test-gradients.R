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

# The tests below check the *full* composite objective built by
# .capnet_objective() -- likelihood + ridge + contribution-cap penalty,
# including the standardize/original-scale chain rule -- rather than just the
# base likelihood covered above. This is the exact closure .capnet_fit() hands
# to lbfgs(), so it exercises the ridge and cap-penalty gradient terms (and
# their interaction with standardization and a non-uniform multiplier) that
# the likelihood-only checks above cannot reach.
.build_objective_fixture <- function(family, standardize, multiplier = 1,
                                      n = 120, p = 5, m = 15, seed = 1) {
  set.seed(seed)
  d <- simulate_capnet_data(n = n, p = p, m = m, family = family)
  L <- d$L

  spec <- .capnet_spec(
    X = d$X, y = d$y, L = L, family = family,
    standardize = standardize, z = d$Z, multiplier = multiplier
  )
  train <- .capnet_standardize_train(spec)
  cap <- .capnet_cap_context(spec)
  list(train = train, cap = cap, p = p)
}

test_that("full objective gradient (capnet_fit.R) matches numerical gradient across families", {
  for (fam in c("gaussian", "binomial", "poisson", "gamma")) {
    for (standardize in c(TRUE, FALSE)) {
      fx <- .build_objective_fixture(fam, standardize)
      params <- list(alpha = 0.3, lambda = 0.7, gamma = 5)

      set.seed(3)
      betas <- rnorm(fx$p + 1, sd = 0.3)
      result <- .objective_gradient_check(fx$train, fx$cap, params, betas)
      expect_lt(
        result$max_diff, 1e-5,
        label = sprintf("family=%s standardize=%s", fam, standardize)
      )
    }
  }
})

test_that("full objective gradient matches numerical gradient with a non-uniform multiplier", {
  fx <- .build_objective_fixture(
    "gaussian", standardize = TRUE,
    multiplier = seq(0.5, 2, length.out = 15)
  )
  params <- list(alpha = 0.2, lambda = 0.4, gamma = 8)

  set.seed(4)
  betas <- rnorm(fx$p + 1, sd = 0.3)
  result <- .objective_gradient_check(fx$train, fx$cap, params, betas)
  expect_lt(result$max_diff, 1e-5)
})

test_that("full objective gradient is well-behaved at the cap-penalty kink (|beta*Z| == L)", {
  # The hinge-squared cap penalty is continuously differentiable everywhere,
  # including exactly at the activation boundary |beta_j * Z_ij| == L_j (the
  # derivative of max(0,u)^2 is 2*max(0,u)*u', which is 0 -- not discontinuous
  # -- at u=0). Construct beta so a contribution lands almost exactly on L_j
  # and confirm the analytic gradient still matches the numerical one there.
  set.seed(5)
  n <- 100; p <- 3; m <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(X %*% c(1, -1, 0.5) + rnorm(n))
  Z <- matrix(rnorm(m * p), m, p)

  # Choose L so that column 1's cap is exactly hit by a mid-sized beta.
  beta_probe <- c(0.4, -0.2, 0.1)
  contributions <- abs(sweep(Z, 2, beta_probe, "*"))
  L <- apply(contributions, 2, function(col) sort(col, decreasing = TRUE)[m %/% 2])

  spec <- .capnet_spec(X = X, y = y, L = L, family = "gaussian",
                       standardize = TRUE, z = Z)
  train <- .capnet_standardize_train(spec)
  cap <- .capnet_cap_context(spec)
  params <- list(alpha = 0, lambda = 0.1, gamma = 20)

  betas <- c(0, beta_probe * train$scaling_factor)
  result <- .objective_gradient_check(train, cap, params, betas)
  expect_lt(result$max_diff, 1e-5)
})
