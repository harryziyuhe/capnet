# Builds the smooth loss/gradient closures (likelihood + ridge + contribution
# cap penalty) that .capnet_fit() hands to lbfgs(). Factored out so tests can
# finite-difference-check the exact composite objective actually optimized,
# not just the base likelihood.
.capnet_objective <- function(train, cap, params) {
  X <- train$X
  y <- train$y
  newx <- cap$z

  # Loss function combining (only the smooth part)
  # 1. Likelihood loss
  # 2. Ridge penalty (LASSO penalty is non-smooth and therefore excluded)
  # 3. Contribution cap penalty

  lik <- make_likelihood(train$family, X, y)

  loss <- function(beta) {
    beta_0 <- beta[1]
    beta_rest <- beta[-1]

    # Likelihood loss
    ll <- lik$loss(beta_0, beta_rest)

    # Ridge penalty
    ridge_penalty <- ((1 - params$alpha) / 2) * sum(beta_rest^2)

    # Contribution cap penalty
    # Contribution cap needs to be computed in the original scale, so if the
    # data is standardized, transform contributions back to the original scale
    m <- nrow(newx)
    beta_raw <- beta_rest / train$scaling_factor

    feature_contribution <- sweep(
      sweep(newx, 2, beta_raw, "*"), 1,
      cap$multiplier, "*"
    )
    excess_contribution <- sweep(
      pmax(sweep(abs(feature_contribution), 2, cap$L, "-"), 0),
      1, cap$multiplier, "/"
    )
    excess_penalty <- sum(excess_contribution^2) / m

    # Total loss (smooth components)
    total <- ll + params$lambda * ridge_penalty + params$gamma * excess_penalty

    total
  }

  # Gradient function (only the smooth part)
  gradient <- function(beta) {
    beta_0 <- beta[1]
    beta_rest <- beta[-1]

    # Likelihood gradient
    g_ll <- lik$gradient(beta_0, beta_rest)
    gradient_0 <- g_ll$g0
    gradient_rest <- g_ll$gb

    # Ridge penalty gradient
    gradient_rest <- gradient_rest +
      (1 - params$alpha) * params$lambda * beta_rest

    # Contribution cap gradient
    # The optimizer operates on standardized beta_rest, but the cap penalty is
    # evaluated on raw-scale contributions (beta_raw = beta_rest /
    # scaling_factor).
    # By the chain rule, d(beta_raw_j)/d(beta_rest_j) = 1 / scaling_factor_j, so
    # that factor is applied once, at the end, to the raw-scale gradient below.

    m <- nrow(newx)
    beta_raw <- beta_rest / train$scaling_factor

    feature_contribution <- sweep(
      sweep(newx, 2, beta_raw, "*"), 1,
      cap$multiplier, "*"
    )
    excess_contribution <- sweep(
      pmax(sweep(abs(feature_contribution), 2, cap$L, "-"), 0),
      1, cap$multiplier, "/"
    )
    d_excess <- sweep(
      excess_contribution * sign(feature_contribution) * newx,
      2, train$scaling_factor, "/"
    )
    gradient_cap <- (2 * params$gamma / m) * colSums(d_excess)

    gradient_rest <- gradient_rest + gradient_cap

    c(gradient_0, gradient_rest)
  }

  list(loss = loss, gradient = gradient)
}

# Engine for fitting capnet models
.capnet_fit <- function(train, cap, params, verbose = 0) {
  obj <- .capnet_objective(train, cap, params)

  if (any(is.finite(train$lower.limits)) ||
      any(is.finite(train$upper.limits))) {
    return(.capnet_fit_box(train, obj, params))
  }

  fit <- lbfgs(obj$loss, obj$gradient, train$par,
    orthantwise_c = params$lambda * params$alpha,
    orthantwise_start = 1,
    max_iterations = train$maxit,
    invisible = 1 - verbose
  )

  list(
    a0 = fit$par[1],
    beta = fit$par[-1],
    convergence = fit$convergence,
    message = fit$message,
    value = fit$value
  )
}

# L-BFGS/OWL-QN has no box-bound support. Use a projected proximal-gradient
# solve when finite coefficient bounds are present so every iterate remains
# feasible, including steps that would otherwise overshoot a bound.
.capnet_fit_box <- function(train, obj, params) {
  lower <- c(-Inf, train$lower.limits)
  upper <- c(Inf, train$upper.limits)
  x <- train$par

  if (any(x < lower | x > upper)) {
    stop("Initial parameter values must satisfy coefficient bounds.",
      call. = FALSE
    )
  }

  l1 <- params$lambda * params$alpha
  smooth_value <- obj$loss
  total_value <- function(beta) {
    smooth_value(beta) + l1 * sum(abs(beta[-1]))
  }
  proximal <- function(value, step) {
    out <- value
    out[-1] <- sign(value[-1]) * pmax(abs(value[-1]) - l1 * step, 0)
    pmax(lower, pmin(upper, out))
  }

  y <- x
  momentum <- 1
  lipschitz <- 1
  value_x <- total_value(x)
  converged <- FALSE

  for (iteration in seq_len(train$maxit)) {
    restarted <- FALSE
    repeat {
      grad <- obj$gradient(y)
      value_y <- smooth_value(y)
      trial <- proximal(y - grad / lipschitz, 1 / lipschitz)
      delta <- trial - y
      majorizer <- value_y + sum(grad * delta) +
        lipschitz * sum(delta^2) / 2

      if (smooth_value(trial) <= majorizer + 1e-12) break
      lipschitz <- lipschitz * 2
    }

    value_trial <- total_value(trial)
    if (!restarted && value_trial > value_x) {
      y <- x
      momentum <- 1
      restarted <- TRUE
      next
    }

    if (max(abs(trial - x)) <= 1e-8 * (1 + max(abs(x)))) {
      x <- trial
      value_x <- value_trial
      converged <- TRUE
      break
    }

    next_momentum <- (1 + sqrt(1 + 4 * momentum^2)) / 2
    y <- trial + ((momentum - 1) / next_momentum) * (trial - x)
    x <- trial
    momentum <- next_momentum
    value_x <- value_trial
  }

  list(
    a0 = x[1],
    beta = x[-1],
    convergence = if (converged) 0L else -1L,
    message = if (converged) "Projected optimization converged." else
      "Projected optimization reached maxit.",
    value = value_x
  )
}
