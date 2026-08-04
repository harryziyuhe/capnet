test_that("cv_capnet returns correctly shaped output and a valid best pair", {
  set.seed(1)
  d <- simulate_capnet_data(n = 80, p = 6, sparsity = 0.5)

  alpha <- c(0, 0.5, 1)
  lambda <- exp(seq(0, -3, length.out = 5))

  cv <- cv_capnet(d$X, d$y, gamma = 1, L = d$L, alpha = alpha, lambda = lambda, K = 4)

  expect_s3_class(cv, "cv_capnet")
  expect_equal(dim(cv$cv_errors), c(length(alpha), length(lambda), 4))
  expect_equal(dim(cv$mean_errors), c(length(alpha), length(lambda)))
  expect_true(cv$best_alpha %in% alpha)
  expect_true(cv$best_lambda %in% lambda)
  expect_equal(cv$best_error, min(cv$mean_errors, na.rm = TRUE))
  expect_equal(cv$metric, "mse")
})

test_that("cv_capnet honors explicit splits and infers K from them", {
  set.seed(2)
  d <- simulate_capnet_data(n = 60, p = 5, sparsity = 0.5)
  splits <- rep(1:3, length.out = 60)

  cv <- cv_capnet(d$X, d$y, gamma = 0.5, L = d$L,
                  alpha = c(0, 1), lambda = c(0.1, 0.5),
                  splits = splits)

  expect_equal(dim(cv$cv_errors)[3], 3)
  expect_equal(cv$splits, as.integer(splits))
})

test_that("cv_capnet picks a family-appropriate default metric", {
  set.seed(3)
  d_binom <- simulate_capnet_data(n = 80, p = 5, family = "binomial", binomial_prob = 0.4)

  cv <- cv_capnet(d_binom$X, d_binom$y, gamma = 0.5, L = d_binom$L,
                  family = "binomial", alpha = c(0, 1), lambda = c(0.1, 0.5), K = 3)

  expect_equal(cv$metric, "logloss")
})
