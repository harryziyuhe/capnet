# Builds the smooth loss/gradient closures (likelihood + ridge + contribution
# cap penalty) that .capnet_fit() hands to lbfgs(). Factored out so tests can
# finite-difference-check the exact composite objective actually optimized,
# not just the base likelihood.
.capnet_objective <- function(train, cap, params) {

  X <- train$X
  y <- train$y
  n <- nrow(X)
  scaling_factor <- train$scaling_factor

  alpha <- params$alpha
  lambda <- params$lambda
  gamma <- params$gamma

  newx <- cap$z
  multiplier <- cap$multiplier
  L <- cap$L
  m <- nrow(newx)

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
    ridge_penalty <- ((1 - params$alpha) / 2) * sum(beta_rest ^ 2)
    
    # Contribution cap penalty
    # Contribution cap needs to be computed in the original scale, so if the
    # data is standardized, transform contributions back to the original scale
    m <- nrow(newx)
    beta_raw <- beta_rest / train$scaling_factor

    feature_contribution <- sweep(sweep(newx, 2, beta_raw, "*"), 1, cap$multiplier, "*")
    excess_contribution <- sweep(
      pmax(sweep(abs(feature_contribution), 2, cap$L, "-"), 0),
      1, cap$multiplier, "/"
    )
    excess_penalty <- sum(excess_contribution ^ 2) / m
    
    # Total loss (smooth components)
    total <- ll + params$lambda * ridge_penalty + params$gamma * excess_penalty
    
    return(total)
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
    gradient_rest <- gradient_rest + (1 - params$alpha) * params$lambda * beta_rest
    
    # Contribution cap gradient
    # The optimizer operates on standardized beta_rest, but the cap penalty is
    # evaluated on raw-scale contributions (beta_raw = beta_rest / scaling_factor).
    # By the chain rule, d(beta_raw_j)/d(beta_rest_j) = 1 / scaling_factor_j, so
    # that factor is applied once, at the end, to the raw-scale gradient below.

    m <- nrow(newx)
    beta_raw <- beta_rest / train$scaling_factor

    feature_contribution <- sweep(sweep(newx, 2, beta_raw, "*"), 1, cap$multiplier, "*")
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
    
    gradient_rest[beta_rest <= train$lower.limits + train$tol & gradient_rest > 0] <- 0
    gradient_rest[beta_rest >= train$upper.limits - train$tol & gradient_rest < 0] <- 0

    return(c(gradient_0, gradient_rest))
  }

  list(loss = loss, gradient = gradient)
}

# Engine for fitting capnet models
.capnet_fit <- function(train, cap, params, verbose = 0) {

  obj <- .capnet_objective(train, cap, params)

  fit <- lbfgs(obj$loss, obj$gradient, train$par,
               orthantwise_c = params$lambda * params$alpha,
               orthantwise_start = 1,
               max_iterations = train$maxit,
               invisible = 1 - verbose)

  list(
    a0 = fit$par[1],
    beta = fit$par[-1],
    convergence = fit$convergence,
    message = fit$message,
    value = fit$value
  )
}
