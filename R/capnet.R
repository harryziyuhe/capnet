#' Fit a linear model with elasticnet and contribution cap regularization
#' 
#' Fit a linear model with via penalized maximum likelihood. The function uses 
#' limited memory BFGS (L-BFGS) when L1 regularization is absent and Orthant-Wise
#' Limited-memory Quasi-Newton (OWL-QN) when L1 regularization is present.
#' 
#' @import lbfgs
#' 
#' @param X input matrix
#' @param y response variable
#' @param lambda the strength of the elasticnet regularizer
#' @param alpha the elasticnet mixing parameter. `alpha=1` collapses to the 
#' LASSO penalty and `alpha=0` the ridge penalty
#' @param mu the strength of the contribution cap regularizer
#' @param L contribution ceiling. If set to a single value the same ceiling is
#' applied to all parameters
#' @param newx data used to evaluate and apply the contribution caps. If unspecified
#' default to the input matrix. 
#' @param par A numeric vector containing the initial values for all variables. 
#' If unspecified default to zero initialization.
#' @param multiplier A numeric value or vector used for scaling feature contribution
#' @param intercept flag for whether the intercept should be fitted. Default to
#' `intercept=TRUE`
#' @param standardize flag for standardization of the `x` and `y` variables.
#' The coefficients are always returned on the original scale. Default is 
#' `standardize=TRUE`
#' @param lower.limits A numeric value or vector containing lower bounds for fitted
#' parameters. The initialized values must be above the `lower.limits`.
#' @param upper.limits A numeric value or vector containing upper bounds for fitted
#' parameters. The initialized values must be below the `upper.limits`.
#' @param tol The tolerance parameter for gradient masking when `lower.limits` or
#' `upper.limits` are specified.
#' @param maxit Maximum number of passes over the data. Default is 10^5
#' @param check.finite Catch and report errors in loss and gradient calculation. 
#' Default to `check.finite=FALSE`
#' @param verbose Defaults to 0. Set to 1 to show console outputs during optimization.
#' 
#' @return A list with the following components:
#' \item{a0}{Best intercept value.}
#' \item{beta}{A numerical array. The best set of parameters found.}
#' \item{convergence}{An integer code. Zero indicates that convergence was
#' reached without issues. Negative values indicate errors in the execution of 
#' the OWL-QN routine.}
#' \item{message}{A character object dealing with execution error. Only returned
#' if the convergence code is different from zero.}
#' \item{value}{The minimized value of the objective function.}
#' \item{feature_contributions}{A numerical matrix recording the contributions of
#' each covariate.}
#' \item{alpha}{The elastic net mixing parameter}
#' \item{lambda}{The elastic net penalty strength specified in the model.}
#' \item{mu}{The contribution cap penalty strength specified in the model.}
#' \item{L}{The contribution cap specified in the model.}
#' \item{multiplier}{The scalar multiplier specified when calculating feature contributions.}
#' \item{call}{The call that produced this object.}
#' 
#' @export


# Build intercept part later
capnet <- function(X, y, lambda, alpha, 
                   mu, L, newx = NULL, par = NULL, multiplier = 1,
                   intercept = TRUE, standardize = TRUE,
                   lower.limits = NULL, upper.limits = NULL,
                   tol = 1e-8, maxit = 10000,
                   check.finite = FALSE, verbose = 0) {
  
  if (anyNA(X) || anyNA(y) || anyNA(newx)) {
    stop("X or y or newx contains NA values")
  }
  
  X_raw <- X
  y_raw <- y
  n <- nrow(X)
  p <- ncol(X)
  
  if (is.null(par)) {
    par <- rep(0, p+1)
  }
  
  if (is.null(newx)) {
    newx <- X_raw
  } else {
    newx <- as.matrix(newx)
  }
  
  
  if (is.data.frame(L) || is.matrix(L)) L <- as.numeric(unlist(L))
  if (is.list(L)) L <- unlist(L)
  if (length(L) == 1) L <- rep(L, p)
  if (length(L) != p) stop("L must be length 1 or match the number of covariates (p = ", p, ")")
  multiplier = as.numeric(multiplier)
  if (length(multiplier) == 1) multiplier <- rep(multiplier, nrow(newx))
  if (length(multiplier) != nrow(newx))stop("multiplier must be length 1 or match the number of observations")
  
  # Standardize data
  if (standardize) {
    y <- scale(y_raw)
    y_center <- attr(y, "scaled:center")
    y_scale <- attr(y, "scaled:scale")
    
    X <- scale(X_raw)
    X_center <- attr(X, "scaled:center")
    X_scale <- attr(X, "scaled:scale")
    
    scaling_factor <- X_scale / y_scale
  } else {
    scaling_factor <- rep(1, p)
  }
  
  
  if (is.null(lower.limits)) {
    lower.limits <- rep(-Inf, p)
  } else if (length(lower.limits) == 1) {
    lower.limits <- rep(lower.limits, p)
  } else if (length(lower.limits) != p) {
    stop("Dimension mismatch.")
  }
  
  if (is.null(upper.limits)) {
    upper.limits <- rep(Inf, p)
  } else if (length(upper.limits) == 1) {
    upper.limits <- rep(upper.limits, p)
  } else if (length(upper.limits) != p) {
    stop("Dimension mismatch.")
  }
  
  # Loss function combining (only the smooth part)
  # 1. Least squares loss
  # 2. Ridge penalty (LASSO penalty is non-smooth and therefore excluded)
  # 3. Concentration cap penalty
  loss <- function(beta) {
    beta_0 <- beta[1]
    beta_rest <- beta[-1]
    
    # Least squares loss
    pred <- beta_0 + X %*% beta_rest
    rss <- sum((y - pred) ^ 2) / (2 * n)
    
    # Ridge penalty
    ridge_penalty <- ((1 - alpha) / 2) * sum(beta_rest ^ 2)
    
    # Concentration cap penalty
    # - Concentration cap needs to be computed in the original scale, so if the
    #   data is standardized, transform contributions back to original scale
    beta_raw <- if (standardize) beta_rest / scaling_factor else beta_rest
    
    feature_contributions <- sweep(sweep(newx, 2, beta_raw, "*"), 1, multiplier, "*")
    
    excess_concentration <- sweep(sweep(sweep(abs(feature_contributions), 2, L, "-"), 2, scaling_factor, "*"), 1, multiplier, "/")

    concentration_penalty <- sum(pmax(0, excess_concentration)^2) / nrow(newx)

    # Total loss (smooth components)
    total <- rss + lambda * ridge_penalty + mu * concentration_penalty
    return(total)
  }
  
  # Gradient function calculation (only the smooth part)
  gradient <- function(beta) {
    beta_0 <- beta[1]
    beta_rest <- beta[-1]
    
    pred <- beta_0 + X %*% beta_rest
    
    # Intercept gradient
    gradient_0 <- -sum(y - pred) / n
    
    # Covariates gradient
    # Least squares loss gradient
    gradient_rss <- (-t(X) %*% (y - pred)) / n
    
    # Ridge penalty gradient
    gradient_ridge <- (1 - alpha) * lambda * beta_rest
    
    # Combining the first two components
    gradient_rest <- gradient_rss + gradient_ridge
    
    # Contribution cap gradient
    # - Because concentration cap penalty is calculated in the original scale,
    #   so if the data is standardized, scaling factors also show up in the gradient
    beta_raw <- if (standardize) beta_rest / scaling_factor else beta_rest
    for (j in seq_along(beta_raw)) {
      Xj <- newx[, j]
      contribution_j <- Xj * beta_raw[j]
      mask <- abs(contribution_j) * multiplier > L[j]
      
      if (any(mask)) {
        excess <- abs(contribution_j[mask]) * multiplier[mask] - L[j]
        signs <- sign(contribution_j[mask])
        penalty_grad_j <- sum(2 * mu * excess * signs * Xj[mask] / multiplier[mask]) / nrow(newx) * scaling_factor[j]
        gradient_rest[j] <- gradient_rest[j] + penalty_grad_j
      }
    }
    gradient_rest[beta_rest <= lower.limits + tol & gradient_rest > 0] <- 0
    gradient_rest[beta_rest >= upper.limits - tol & gradient_rest < 0] <- 0
    return(c(gradient_0, gradient_rest))
  }
  
  if (check.finite) {
    safe_loss <- function(beta) {
      tryCatch({
        val <- loss(beta)
        if (!is.finite(val)) stop("Loss is non-finite")
        val
      }, error = function(e) {
        print(beta)
        stop(e)
      })
    }
    
    safe_grad <- function(beta) {
      tryCatch({
        grad <- gradient(beta)
        if (any(!is.finite(grad))) stop("Gradient has non-finite values")
        grad
      }, error = function(e) {
        print(beta)
        stop(e)
      })
    }
    
    model <- lbfgs(safe_loss, safe_grad, par,
                   orthantwise_c = lambda * alpha,
                   orthantwise_start = 1,
                   max_iterations = maxit,
                   invisible = 1 - verbose)
  } else {
    model <- lbfgs(loss, gradient, par,
                   orthantwise_c = lambda * alpha,
                   orthantwise_start = 1,
                   max_iterations = maxit,
                   invisible = 1)
  }
  
  # Extract coefficients
  a0 <- model$par[1]
  beta_rest <- model$par[-1]
  
  # Calculate minimized value of the total objective function (including the L1 term)
  value <- model$value + lambda * alpha * sum(abs(beta_rest))
  
  if (standardize) {
    beta_rest <- beta_rest / X_scale * y_scale
    a0 <- y_center - sum(beta_rest * X_center) + a0 * y_scale
  }
  
  # Calculate feature contributions for final model
  contributions <- sweep(sweep(newx, 2, beta_rest, "*"), 1, multiplier, "*")
  
  # Return results
  structure(list(
    a0 = a0,
    beta = beta_rest,
    newx = newx,
    convergence = model$convergence,
    message = model$message,
    value = value,
    feature_contributions = contributions,
    alpha = alpha,
    lambda = lambda,
    mu = mu,
    L = L,
    multiplier = multiplier,
    call = match.call()
  ), class = c("capnet"))
}
