# Finalize capnet outputs
.capnet_output <- function(train, cap, fit, params, call = NULL) {
  a0 <- fit$a0
  beta <- fit$beta
  
  value <- fit$value + params$lambda * params$alpha * sum(abs(beta))
  
  sf <- train$scaling_factor
  beta_raw <- beta / sf
  if (train$intercept) {
    a0_raw <- a0 - sum((train$X_center / train$X_scale) * beta)
  } else {
    a0_raw <- 0
  }

  contributions <- sweep(sweep(cap$newx, 2, beta_raw, "*"), 1, cap$multiplier, "*")
  
  structure(list(
    a0 = a0_raw,
    beta = beta_raw,
    value = value,
    feature_contributions = contributions,
    newx = cap$newx,
    multiplier = cap$multiplier,
    L = cap$L,
    convergence = fit$convergence,
    message = fit$message,
    alpha = params$alpha,
    lambda = params$lambda,
    mu = params$mu,
    family = train$family,
    intercept = train$intercept,
    standardize = train$standardize,
    X_center = train$X_center,
    X_scale = train$X_scale,
    call = call
  ), class = "capnet")
}
