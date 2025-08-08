#' Perform Walk Foward Model Fit
#' 
#' Fit a linear elastic net model with optimized contribution cap.
#' Apply walk forward on the evaluation set.
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
#' @param walk number of rows to evaluate and fit models at each time. Default to
#' `walk = 1`
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
#' 
#' @return A list with the following components:
#' \item{intercepts}{Best intercept values.}
#' \item{betas}{The best sets of parameters found.}
#' \item{feature_contributions}{A numerical matrix recording the contributions of
#' each covariate.}
#' \item{predictions}{A numerical matrix recording the model predictions.}
#' \item{mus}{The contribution cap penalty strengths specified in the model.}
#' 
#' @export

capnet.walk <- function(X, y, lambda, alpha, mu, L, newx,
                        walk = 1, par = NULL, multiplier = 1,
                        standardize = TRUE,
                        lower.limits = NULL, upper.limits = NULL,
                        tol = 1e-8, maxit = 10000) {
  
  # Stop if there is any NA values in data
  if (anyNA(X) || anyNA(y) || anyNA(newx)) {
    stop("X or y or newx contains NA values")
  }
  if (length(multiplier) == 1) multiplier <- rep(multiplier, nrow(newx))
  if (length(multiplier) != nrow(newx)) stop("multiplier must be length 1 or match the number of observations")
  
  n <- nrow(newx)
  intercepts <- numeric(n)
  betas <- NULL
  contributions <- NULL
  predictions <- numeric(n)
  mu_values <- numeric(n)
  
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
  
  fitpoints <- c(seq(1, n, walk),n+1)
  for (i in 1:(length(fitpoints)-1)) {
    start_row <- fitpoints[i]
    end_row <- fitpoints[i+1]
    m <- end_row - start_row
    newx.tmp <- matrix(newx[start_row:(end_row-1), ], nrow = m)
    multiplier.tmp <- as.numeric(multiplier)[start_row:(end_row-1)]
    convergence.tmp <- -999
    while (convergence.tmp < 0) {
      capnet_results <- capnet(X_scaled, y_scaled, lambda = lambda, alpha = alpha,
                               mu = mu.tmp, L = L, newx = newx.tmp, 
                               standardize = FALSE,
                               multiplier = multiplier.tmp,
                               lower.limits = lower.limits, 
                               upper.limits = upper.limits, tol = tol,
                               maxit = maxit)
      convergence.tmp <- capnet_results$convergence
      mu.tmp <- mu.tmp / 10
    }
    intercepts[start_row:(end_row-1)] <- capnet_results$a0
    betas <- rbind(betas, matrix(rep(capnet_results$beta, m), nrow = m))
    contributions <- rbind(contributions, capnet_results$contributions)
    predictions[start_row:(end_row-1)] <- rowSums(capnet_results$contributions) + capnet_results$a0
    mu_values[start_row:(end_row-1)] <- capnet_results$mu
  }
  
  if (is.xts(newx)) {
    intercepts <- as.xts(intercepts, order.by = index(newx))
    betas <- as.xts(betas, order.by = index(newx))
    contributions <- as.xts(contributions, order.by = index(newx))
    predictions <- as.xts(predictions, order.by = index(newx))
    mu_values <- as.xts(mu_values, order.by = index(newx))
    colnames(intercepts) <- "intercept"
    colnames(betas) <- colnames(contributions) <- colnames(X)
    colnames(predictions) <- "prediction"
    colnames(mu_vaues) <- "mu"
  } else {
    intercept <- matrix(intercept, ncol = 1)
    predictions <- matrix(predictions, ncol = 1)
    mu_values <- matrix(mu_values, ncol = 1)
    colnames(intercepts) <- "intercept"
    colnames(betas) <- colnames(contributions) <- colnames(X)
    colnames(predictions) <- "prediction"
    colnames(mu_vaues) <- "mu"
  }
  
  structure(list(
    intercepts = intercepts,
    betas = betas,
    feature_contributions = contributions,
    predictions = predictions,
    mus = mu_values,
    call = match.call()
  ), class = c("capnet_walk"))
}