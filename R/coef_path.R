#' Calculate coefficients along a penalty path
#' 
#' Computes fitted coefficients from \code{capnet()} along a path of 
#' \code{lambda} values (with \code{alpha} held fixed). This is useful for
#' inspecting coefficient shrinkage and entry/exit as regularization increases.
#' 
#' @param X Numeric predictor matrix of shape \eqn{n\times p}.
#' @param y Numeric response vector of length \eqn{n}.
#' @param L Nonnegative numeric scalar or length-\eqn{p} vector giving the 
#'  contribution ceiling(s). If scalar, the same ceiling is applied to all
#'  coefficients.
#' @param lambda Numeric vector (default \code{exp(seq(1, -5, length.out = 50))})
#'  of nonnegative elastic-net penalties to trace.
#' @param alpha Numeric scalar in \eqn{[0,1]} (default \code{0.5}); the 
#'  elastic-net mixing parameter. \code{alpha = 1} is Lasso, \code{alpha=0} is 
#'  Ridge. 
#' @param mu Nonnegative numeric scalar; strength of the contribution-cap 
#' penalty (default \code{1}).
#' @param standardize Logical; if \code{TRUE}, columns of \code{X} and \code{y}
#'  are standardized for fitting; coefficients are returned on the original scale.
#'  Default \code{TRUE}.
#' @param ... Additional arguments forwarded to \code{capnet()}, e.g.,
#'  \code{newx}, \code{par}, \code{multiplier}, \code{intercept},
#'  \code{lower.limits}, \code{upper.limits}, \code{tol}, \code{maxit},
#'  \code{check.finite}, \code{verbose}.
#' 
#' @return A \code{data.frame} in long/tidy format with one row per
#'  (\code{lambda}, coefficient) pair.
#'  
#' @details
#' For each \code{lambda} in the supplied vector, \code{capnet()} is fit with
#' fixed \code{alpha}, \code{mu}, and \code{L}, and the resulting coefficients
#' are collected. If feature names are available from \code{colnames(X)}, they
#' are used as column names.
#' 
#' @seealso [capnet()], [plot.capnet_path()]
#'  
#' @examples
#' set.seed(1)
#' n <- 50; p <- 8
#' X <- matrix(rnorm(n * p), n, p)
#' colnames(X) <- paste0("x", 1:p)
#' beta <- c(1.5, 0.8, 0.2, rep(0, p - 3))
#' y <- as.numeric(X %*% beta + rnorm(n))
#' path <- coef_path(X, y, L = 0.5, alpha = 0.5,
#'                   lambda = exp(seq(1, -5, length.out = 50)), mu = 1)
#' @export

coef_path <- function(X, y, L,
                      lambda = exp(seq(1, -5, length.out = 50)), 
                      alpha = 0.5, 
                      mu = 1,
                      standardize = TRUE,
                      ...) {
  # Stop if there is any NA value in data
  if (anyNA(X) || anyNA(y)) {
    stop("X or y contains NA values")
  }
  
  betas <- NULL
  p <- ncol(X)
  
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
    path = "lambda (log)"
    alpha = rep(alpha, length(lambda))
    mu = rep(mu, length(lambda))
  } else if (length(mu) > 1) {
    if (length(alpha) != 1) {
      stop("alpha must be a constant value when evaluating along mu path")
    }
    param = log(mu)
    path = "mu (log)"
    lambda = rep(lambda, length(mu))
    alpha = rep(alpha, length(mu))
  } else if (length(alpha) > 1) {
    param = alpha
    path = "alpha"
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
