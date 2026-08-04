test_that("coef.capnet returns a labeled intercept + coefficient column", {
  set.seed(1)
  n <- 40; p <- 4
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("v", 1:p)
  y <- as.numeric(X %*% c(1, -1, 0, 0) + rnorm(n))
  fit <- capnet(X, y, L = 1, alpha = 0.5, lambda = 0.1, gamma = 1)

  cf <- coef(fit)
  expect_equal(dim(cf), c(p + 1, 1))
  expect_equal(rownames(cf), c("(Intercept)", colnames(X)))
  expect_equal(cf["(Intercept)", 1], fit$a0)
})

test_that("predict.capnet gives link and response scale predictions consistent with the fit", {
  set.seed(2)
  n <- 40; p <- 3
  X <- matrix(rnorm(n * p), n, p)
  y_lin <- as.numeric(X %*% c(1, -1, 0.5))
  y <- rbinom(n, 1, plogis(y_lin))
  fit <- capnet(X, y, L = 1, family = "binomial", alpha = 0.2, lambda = 0.05, gamma = 0.5)

  eta <- predict(fit, newdata = X, type = "link")
  mu <- predict(fit, newdata = X, type = "response")

  expect_equal(as.numeric(eta), as.numeric(fit$a0 + X %*% fit$beta), tolerance = 1e-8)
  expect_equal(as.numeric(mu), plogis(as.numeric(eta)), tolerance = 1e-8)
  expect_true(all(mu >= 0 & mu <= 1))
})

test_that("predict.capnet errors on mismatched newdata dimensions", {
  set.seed(3)
  n <- 30; p <- 3
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(X %*% c(1, 0, -1) + rnorm(n))
  fit <- capnet(X, y, L = 1, alpha = 0.5, lambda = 0.1, gamma = 1)

  expect_error(predict(fit, newdata = X[, 1:2]), "ncol")
})

test_that("coef.walk_capnet and predict.walk_capnet expose the stored path", {
  set.seed(4)
  n <- 40; p <- 3; n_new <- 6
  X <- matrix(rnorm(n * p), n, p)
  z <- matrix(rnorm(n_new * p), n_new, p)
  y <- as.numeric(X %*% c(1, -1, 0) + rnorm(n))
  out <- walk_capnet(X, y, L = 1, z = z, walk = 1, lambda = 0.1, alpha = 0.5, gamma = 1)

  cf <- coef(out)
  expect_equal(dim(cf), c(n_new, p + 1))
  expect_equal(as.numeric(cf[, 1]), as.numeric(out$intercepts))

  cf_subset <- coef(out, index = 1:3)
  expect_equal(nrow(cf_subset), 3)

  expect_equal(predict(out), out$predictions)
})

test_that("plot.cv_capnet produces a ggplot for the heatmap, alpha-slice, and lambda-slice views", {
  skip_if_not_installed("ggplot2")
  set.seed(5)
  d <- simulate_capnet_data(n = 60, p = 5, sparsity = 0.5)
  cv <- cv_capnet(d$X, d$y, gamma = 1, L = d$L,
                  alpha = c(0, 0.5, 1), lambda = c(0.1, 0.5), K = 3)

  p1 <- plot(cv)
  expect_s3_class(p1, "ggplot")

  p2 <- plot(cv, alpha = 0.5)
  expect_s3_class(p2, "ggplot")

  p3 <- plot(cv, lambda = 0.1)
  expect_s3_class(p3, "ggplot")
})

test_that("plot.capnet_path produces a ggplot", {
  skip_if_not_installed("ggplot2")
  set.seed(6)
  n <- 40; p <- 4
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", 1:p)
  y <- as.numeric(X %*% c(1, -1, 0, 0) + rnorm(n))
  path <- coef_path(X, y, L = 0.5, alpha = 0.5,
                    lambda = exp(seq(1, -3, length.out = 8)), gamma = 1)

  p <- plot(path)
  expect_s3_class(p, "ggplot")
})

test_that("print.capnet_violations runs without error on a violations object", {
  set.seed(7)
  n <- 50; p <- 3
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(X %*% c(3, -3, 0) + rnorm(n))
  fit <- capnet(X, y, L = 0.05, alpha = 0, lambda = 0, gamma = 0)
  v <- capnet_violations(fit)

  expect_output(print(v))
})
