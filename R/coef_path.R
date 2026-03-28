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
#' @param gamma Nonnegative numeric scalar; strength of the contribution-cap 
#' penalty (default \code{1}).
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
#' fixed \code{alpha}, \code{gamma}, and \code{L}, and the resulting coefficients
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
#'                   lambda = exp(seq(1, -5, length.out = 50)), gamma = 1)
#' @export

coef_path <- function(X, y, L,
                      lambda = exp(seq(1, -5, length.out = 50)), 
                      alpha = 0.5, 
                      gamma = 1,
                      ...) {
  # Stop if there is any NA value in data
  if (anyNA(X) || anyNA(y)) {
    stop("X or y contains NA values")
  }
  
  betas <- NULL
  
  if (length(lambda) > 1) {
    if (length(alpha) != 1) {
      stop("alpha must be a constant value when evaluating along lambda path")
    }
    if (length(gamma) != 1) {
      stop("gamma must be a constant value when evaluating along lambda path")
    }
    param = log(lambda)
    path = "lambda (log)"
    alpha = rep(alpha, length(lambda))
    gamma = rep(gamma, length(lambda))
  } else if (length(gamma) > 1) {
    if (length(alpha) != 1) {
      stop("alpha must be a constant value when evaluating along gamma path")
    }
    param = log(gamma)
    path = "gamma (log)"
    lambda = rep(lambda, length(gamma))
    alpha = rep(alpha, length(gamma))
  } else if (length(alpha) > 1) {
    param = alpha
    path = "alpha"
    lambda = rep(lambda, length(alpha))
    gamma = rep(gamma, length(alpha))
  } else {
    stop("at least one of lambda, alpha, and gamma must have length more than 1")
  }
  
  for (i in seq_along(lambda)) {
    l = lambda[i]
    a = alpha[i]
    g = gamma[i]
    capnet_results <- capnet(X,
                             y,
                             lambda = l,
                             alpha = a,
                             gamma = g,
                             L = L,
                             standardize = FALSE,
                             ...)
    betas <- rbind(betas, capnet_results$beta)
  }
  betas <- data.frame(cbind(param, betas))
  colnames(betas) <- c(path, colnames(X))
  
  structure(betas, class = c("capnet_path", "data.frame"))
}
