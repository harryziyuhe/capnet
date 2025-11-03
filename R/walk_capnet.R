#' Perform walk-forward model fitting with contribution caps
#' 
#' Fits a linear elastic net model with an additional contribution-cap penalty
#' by repeatedly refitting over expanding data and generating out-of-sample
#' predictions in a walk-forward manner. At each step, the model is trained on
#' all rows up to the current cut point and evaluated on the next \code{walk}
#' rows.
#' 
#' @import xts
#' 
#' @param X Numeric predictor matrix of shape \eqn{n\times p}. Columns are
#'  features and rows are observations.
#' @param y Numeric response vector of length \eqn{n}.
#' @param lambda Nonnegative numeric scalar; overall strength of the elastic net penalty. 
#' @param alpha Numeric scalar in \eqn{[0,1]}; elastic net mixing parameter.
#'  \code{alpha = 1} is Lasso, \code{alpha=0} is Ridge.
#' @param mu Nonnegative numeric scalar; strength of the contribution-cap penalty
#' @param L Nonnegative numeric scalar or length-\eqn{p} vector giving the 
#'  contribution ceiling(s). If scalar, the same ceiling is applied to all
#'  coefficients
#' @param newx Optional numeric matrix with \eqn{p} columns used to evaluate and 
#'  enforce contribution caps. If \code{NULL}, defaults to \code{X}.
#' @param walk Integer; number of consecutive rows predicted at each step before
#'  advancing the window. Default \code{1}.
#' @param par Optional numeric vector of length \eqn{p} with initial coefficient
#'  values. If \code{NULL}, uses zero initialization. 
#' @param multiplier Optional numeric scalar or length-\eqn{n} vector used to
#'  scale feature contributions during the capping step. Defaults to 1.
#' @param intercept Logical; should an intercept be fitted? Default \code {TRUE}.
#' @param standardize Logical; if \code{TRUE}, columns of \code{X} and \code{y}
#'  are standardized for fitting; coefficients are returned on the original scale.
#'  Default \code{TRUE}.
#' @param lower.limits lower.limits Optional numeric scalar or length-\eqn{p} vector of lower
#'  bounds on coefficients. Initial values must satisfy the bounds.
#' @param upper.limits Optional numeric scalar or length-\eqn{p} vector of upper
#'  bounds on coefficients. Initial values must satisfy the bounds. 
#' @param tol Nonnegative numeric tolerance used for gradient masking when 
#'  \code{lower.limits} or \code{upper.limits} are specified. Default \code{1e-8}.
#' @param maxit Integer; maximum number of quasi-Newton iterations. Default 
#'  \code{1e5}.
#' 
#' @return An object of class \code{"capnet_walk"} with components:
#' \itemize{
#'  \item \code{intercepts} Numeric vector of length \eqn{S} with fitted 
#'    intercepts for each step.
#'  \item \code{betas} Numeric matrix of shape \eqn{p\times S} with fitted
#'    coefficients per step
#'  \item \code{feature_contributions} Numeric matrix of shape 
#'    \eqn{S\times p} giving per-row, per-feature contributions 
#'    stacked across all evaluation rows in order of prediction.
#'  \item \code{predictions} Numeric matrix of out-of-sample predictions for 
#'    the evaluation rows; shape \eqn{S\times p}.
#'  \item \code{mus} Numeric vector of length \eqn{n_\mathrm{new}} for the 
#'    \code{mu} used at each step.
#' }
#' 
#' @details
#' Given the evaluation matrix \code{newx} of size \eqn{S\times p}. At step 
#' \eqn{s=1,\dots,S}, the model is trained on \code{X} and \code{y} and the 
#' contribution caps are evaluated on rows 
#' \eqn{(s-1)\times\text{walk}:s\times\text{walk}}. Each fit calls 
#' \code{capnet()} internally with the provided hyperparameters and constraints.
#' 
#' @seealso [capnet()], [predict.capnet()], [coef.capnet()]
#' 
#' @examples
#' set.seed(1)
#' n <- 60; p <- 6; n_new <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' newx <- matrix(rnorm(n_new * p), n_new, p)
#' beta <- c(2.5, 1.5, 0.8, rep(0, p - 3))
#' y <- as.numeric(X %*% beta + rnorm(n))
#' out <- capnet.walk(X, y, lambda = 0.1, alpha = 0.5, mu = 1, L = 0.5, newx = newx, walk = 1)
#' 
#' @export

walk_capnet <- function(X, y, lambda, alpha, mu, L, newx,
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
    p <- ncol(X)
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
    mu.tmp = mu
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
    contributions <- rbind(contributions, capnet_results$feature_contributions)
    predictions[start_row:(end_row-1)] <- rowSums(capnet_results$feature_contributions) + capnet_results$a0
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
    colnames(mu_values) <- "mu"
  } else {
    intercepts <- matrix(intercepts, ncol = 1)
    predictions <- matrix(predictions, ncol = 1)
    mu_values <- matrix(mu_values, ncol = 1)
    colnames(intercepts) <- "intercept"
    colnames(betas) <- colnames(contributions) <- colnames(X)
    colnames(predictions) <- "prediction"
    colnames(mu_values) <- "mu"
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
