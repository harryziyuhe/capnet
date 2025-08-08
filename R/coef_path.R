#' Calculate coefficients along a parameter path
#' 
#' Calculate the coefficient values from the capnet model along the path
#' of one hyperparameter values.
#' 
#' @param X input matrix
#' @param y response variable
#' @param lambda the strength of the elasticnet regularizer
#' @param alpha the elasticnet mixing parameter
#' @param mu the strength of the contribution cap regularizer
#' @param standardize flag for standardization of the `x` and `y` variables.
#' The coefficients are always returned on the original scale. Default is 
#' `standardize=TRUE`
#' 
#' @return A data.frame of the coefficients for the hyperparameter sets

coef_path <- function(X, y, L,
                      lambda = exp(seq(1, -5, length.out = 50)), 
                      alpha = 0.5, 
                      mu = 1,
                      standardize = TRUE,
                      ...) {
  # Stop if there is any NA value in data
  if (anyNA(x) || anyNA(y)) {
    stop("X or y contains NA values")
  }
  
  betas <- NULL
  
  if (standardize) {
    X_scaled <- scale(X)
    X_center <- attr(X_scaled, "scaled:center")
    X_scale <- attr(X_scaled, "scaled:scale")
    
    y_scaled <- scale(y)
    y_center <- attr(y_scaled, "scaled:center")
    y_scale <- attr(y_scaled, "scaled:scale")
  } else {
    X_scaled <- X
    X_center <- rep(0, p)
    X_scale <- rep(1, p)
    
    y_scaled <- y
    y_center <- 0
    y_scale <- 1
  }
  
  if (length(lambda) > 1) {
    if (length(alpha) != 1) {
      stop("alpha must be a constant value when evaluating along lambda path")
    }
    if (length(mu) != 1) {
      stop("mu must be a constant value when evaluating along lambda path")
    }
    param = log(lambda)
    path = "λ (log)"
    alpha = rep(alpha, length(lambda))
    mu = rep(mu, length(lambda))
  } else if (length(mu) > 1) {
    if (length(alpha) != 1) {
      stop("alpha must be a constant value when evaluating along mu path")
    }
    param = log(mu)
    path = "μ (log)"
    lambda = rep(lambda, length(mu))
    alpha = rep(alpha, length(mu))
  } else if (length(alpha) > 1) {
    param = alpha
    path = "α"
    lambda = rep(lambda, length(alpha))
    mu = rep(mu, length(alpha))
  } else {
    stop("at least one of lambda, alpha, and mu must have length more than 1")
  }
  
  for (i in seq_along(lambda)) {
    l = lambda[i]
    a = alpha[i]
    m = mu[i]
    capnet_results <- capnet(X_scaled,
                             y_scaled,
                             lambda = l,
                             alpha = a,
                             mu = m,
                             L = L,
                             standardize = FALSE,
                             ...)
    betas <- rbind(betas, capnet_results$beta)
  }
  betas <- data.frame(cbind(param, betas))
  colnames(betas) <- c(path, colnames(X))
  
  structure(betas, class = c("capnet_path", "data.frame"))
}
