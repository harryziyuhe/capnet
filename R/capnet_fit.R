# Engine for fitting capnet models
.capnet_fit <- function(train, cap, params, verbose = 0) {
  
  X <- train$X
  y <- train$y
  n <- nrow(X)
  scaling_factor <- train$scaling_factor
  
  alpha <- params$alpha
  lambda <- params$lambda
  mu <- params$mu
  
  newx <- cap$newx
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
    m <- nrow(cap$newx)
    beta_raw <- beta_rest / train$scaling_factor
    
    feature_contribution <- sweep(sweep(cap$newx, 2, beta_raw, "*"), 1, cap$multiplier, "*")
    excess_contribution <- sweep(
      sweep(
        pmax(sweep(abs(feature_contribution), 2, cap$L, "-"), 0), 
        2, train$scaling_factor, "*"
      ), 
      1, cap$multiplier, "/"
    )
    excess_penalty <- sum(excess_contribution ^ 2) / m
    
    # Total loss (smooth components)
    total <- ll + params$lambda * ridge_penalty + params$mu * excess_penalty
    
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
    # Because contribution cap penalty is calculated in the original scale,
    # so if the data is standardized, scaling factors also show up in the gradient
    
    m <- nrow(cap$newx)
    beta_raw <- beta_rest / train$scaling_factor
    
    feature_contribution <- sweep(sweep(cap$newx, 2, beta_raw, "*"), 1, cap$multiplier, "*")
    excess_contribution <- sweep(
      sweep(
        pmax(sweep(abs(feature_contribution),  2, cap$L, "-"), 0), 
        2, train$scaling_factor, "*"
      ), 
      1, cap$multiplier, "/"
    )
    d_excess <- excess_contribution * sign(feature_contribution) * cap$newx
    gradient_cap <- (2 * params$mu / m) * colSums(d_excess)

    gradient_rest <- gradient_rest + gradient_cap
    
    gradient_rest[beta_rest <= train$lower.limits + train$tol & gradient_rest > 0] <- 0
    gradient_rest[beta_rest >= train$upper.limits - train$tol & gradient_rest < 0] <- 0
    
    return(c(gradient_0, gradient_rest))
  }
  
  fit <- lbfgs(loss, gradient, train$par,
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
