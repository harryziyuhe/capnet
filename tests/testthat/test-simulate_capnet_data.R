test_that("simulate_capnet_data() is reproducible given the same seed (basic mode, all families)", {
  for (fam in c("gaussian", "binomial", "poisson", "gamma")) {
    set.seed(2026)
    d1 <- simulate_capnet_data(n = 60, p = 8, m = 10, family = fam)
    set.seed(2026)
    d2 <- simulate_capnet_data(n = 60, p = 8, m = 10, family = fam)

    expect_identical(d1$X, d2$X, info = fam)
    expect_identical(d1$y, d2$y, info = fam)
    expect_identical(d1$beta_true, d2$beta_true, info = fam)
    expect_identical(d1$L, d2$L, info = fam)
    expect_identical(d1$Z, d2$Z, info = fam)
  }
})

test_that("simulate_capnet_data() is reproducible given the same seed (proxy mode)", {
  set.seed(2026)
  d1 <- simulate_capnet_data(mode = "proxy", n = 60, p = 8, m = 10)
  set.seed(2026)
  d2 <- simulate_capnet_data(mode = "proxy", n = 60, p = 8, m = 10)

  expect_identical(d1$X, d2$X)
  expect_identical(d1$y, d2$y)
  expect_identical(d1$X_latent, d2$X_latent)
  expect_identical(d1$Z, d2$Z)
})

test_that("simulate_capnet_data() also accepts the seed= argument directly", {
  d1 <- simulate_capnet_data(n = 60, p = 8, m = 10, seed = 7)
  d2 <- simulate_capnet_data(n = 60, p = 8, m = 10, seed = 7)
  expect_identical(d1$X, d2$X)
  expect_identical(d1$y, d2$y)
})

test_that("simulate_capnet_data() output matches a pinned golden snapshot (gaussian, basic mode)", {
  # This is the exact call used in the JSS paper's software walkthrough
  # (Section 4) and vignettes. If this test ever fails, either a real
  # regression was introduced, or the change was intentional and the
  # paper's numbers need to be regenerated together with this snapshot.
  set.seed(123)
  d <- simulate_capnet_data(n = 100, p = 10, m = 20)

  expect_equal(
    unname(d$X[1:3, 1:3]),
    matrix(c(-0.5604756, -0.2301775, 1.5587083,
              0.1176466, -0.9474746, -0.4905574,
             -0.7886220, -0.5021987, 1.4960607),
           nrow = 3),
    tolerance = 1e-6
  )
  expect_equal(
    d$beta_true,
    c(0.6303808, 0.9429246, -0.1396230, 0, 0, 1.4544213, 0.7000023, 0, -1.3541284, -1.9516388),
    tolerance = 1e-6
  )
  expect_equal(d$y[1:5],
    c(-2.9712089, -6.4706397, 0.3439867, 0.1267018, 1.1501935),
    tolerance = 1e-6
  )
  expect_equal(d$L, pmax(abs(d$beta_true), 0.1), tolerance = 1e-6)
})

test_that("simulate_capnet_data() preserves positional-argument compatibility", {
  # Guards against silently reordering the parameter list: family/poisson_mean/
  # gamma_mean/gamma_shape/binomial_prob must stay appended after the original
  # 14 parameters (n, p, m, mode, sparsity, beta_range, sigma, seed, xstar_sd,
  # a_sd, unstable_idx, unstable_frac, sigma_low, sigma_high, high_prob), or
  # any code calling this function positionally silently breaks.
  d_named <- simulate_capnet_data(
    n = 100, p = 10, m = 20, mode = "basic", sparsity = 0.7,
    beta_range = c(-2, 2), sigma = 1, seed = 42, xstar_sd = 1, a_sd = 0.3,
    unstable_idx = NULL, unstable_frac = 0.2, sigma_low = 1.0,
    sigma_high = 8.0, high_prob = 0.1
  )
  d_positional <- simulate_capnet_data(
    100, 10, 20, "basic", 0.7, c(-2, 2), 1, 42, 1, 0.3, NULL, 0.2, 1.0, 8.0, 0.1
  )

  expect_identical(d_named$X, d_positional$X)
  expect_identical(d_named$y, d_positional$y)
  expect_identical(d_named$beta_true, d_positional$beta_true)
})

test_that("simulate_capnet_data() generates family-appropriate responses", {
  set.seed(1)
  d_gauss <- simulate_capnet_data(n = 200, p = 6, family = "gaussian")
  expect_true(is.numeric(d_gauss$y))

  d_binom <- simulate_capnet_data(n = 200, p = 6, family = "binomial")
  expect_true(all(d_binom$y %in% c(0, 1)))

  d_pois <- simulate_capnet_data(n = 200, p = 6, family = "poisson")
  expect_true(all(d_pois$y >= 0))
  expect_true(all(d_pois$y == round(d_pois$y)))

  d_gamma <- simulate_capnet_data(n = 200, p = 6, family = "gamma")
  expect_true(all(d_gamma$y > 0))
})

test_that("family = <non-gaussian> is rejected in proxy mode", {
  expect_error(
    simulate_capnet_data(mode = "proxy", family = "poisson"),
    "only supported in mode"
  )
})

test_that("simulated data fits cleanly through capnet() for every family", {
  for (fam in c("gaussian", "binomial", "poisson", "gamma")) {
    set.seed(3)
    d <- simulate_capnet_data(n = 150, p = 6, m = 15, family = fam)
    fit <- capnet(d$X, d$y, d$L, z = d$Z, family = fam, alpha = 0.5, lambda = 0.05, gamma = 1)
    expect_equal(fit$convergence, 0, info = fam)
  }
})
