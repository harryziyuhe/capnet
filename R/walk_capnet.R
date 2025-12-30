#' Perform walk-forward model fitting with contribution caps
#' 
#' Fits a linear elastic net model with an additional contribution-cap penalty
#' by repeatedly refitting over expanding data and generating out-of-sample
#' predictions in a walk-forward manner. At each step, the model is trained on
#' all rows up to the current cut point and evaluated on the next \code{walk}
#' rows.
#' 
#' @import xts
#' @importFrom zoo index
#' 
#' @param X Numeric predictor matrix of shape \eqn{n\times p}. Columns are
#'  features and rows are observations.
#' @param y Numeric response vector of length \eqn{n}.
#' @param L Nonnegative numeric scalar or length-\eqn{p} vector giving the 
#'  contribution ceiling(s). If scalar, the same ceiling is applied to all
#'  coefficients
#' @param newx Optional numeric matrix with \eqn{p} columns used to evaluate and 
#'  enforce contribution caps. If \code{NULL}, defaults to \code{X}.
#' @param family Optional character scalar (e.g. "binomial"), function (e.g. 
#' \code{stats::binomial}), or family object (e.g. \code{stats::binomial()}).
#' @param intercept Logical; should an intercept be fitted? Default \code{TRUE}.
#' @param standardize Logical; if \code{TRUE}, columns of \code{X} and \code{y}
#'  are standardized for fitting; coefficients are returned on the original scale.
#'  Default \code{TRUE}.
#' @param multiplier Optional numeric scalar or length-\eqn{n} vector used to
#'  scale feature contributions during the capping step. Defaults to 1.
#' @param walk Integer; number of consecutive rows predicted at each step before
#'  advancing the window. Default \code{1}.
#' @param lambda Nonnegative numeric scalar; overall strength of the elastic net penalty. 
#' @param alpha Numeric scalar in \eqn{[0,1]}; elastic net mixing parameter.
#'  \code{alpha = 1} is Lasso, \code{alpha=0} is Ridge.
#' @param mu Nonnegative numeric scalar; strength of the contribution-cap penalty
#' @param ... Additional arguments used in fitting. See [capnet()] for more details.
#' 
#' @return An object of class \code{"walk_capnet"} with components:
#'  \item{\code{intercepts}}{Numeric vector of length \eqn{S} with fitted 
#'    intercepts for each step.}
#'  \item{\code{betas}}{Numeric matrix of shape \eqn{p\times S} with fitted
#'    coefficients per step.}
#'  \item{\code{feature_contributions}}{Numeric matrix of shape 
#'    \eqn{S\times p} giving per-row, per-feature contributions 
#'    stacked across all evaluation rows in order of prediction.}
#'  \item{\code{predictions}}{Numeric matrix of out-of-sample predictions for 
#'    the evaluation rows; shape \eqn{S\times p}.}
#'  \item{\code{mus}}{Numeric vector of length \eqn{n_\mathrm{new}} for the 
#'    \code{mu} used at each step.}
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
#' out <- walk_capnet(X, y, lambda = 0.1, alpha = 0.5, mu = 1, L = 0.5, newx = newx, walk = 1)
#' 
#' @export

walk_capnet <- function(X, y, L, newx,
                        family = "gaussian",
                        intercept = TRUE, 
                        standardize = TRUE,
                        multiplier = 1,
                        walk = 1,
                        lambda = 0,
                        alpha = 0,
                        mu = 0,
                        max_mu_tries = 6,
                        ...) {
  
  # Stop if there is any NA values in data
  if (anyNA(X) || anyNA(y) || anyNA(newx)) {
    stop("X or y or newx contains NA values")
  }
  
  spec <- .capnet_spec(
    X = X, y = y, L = L,
    family = family,
    intercept = intercept,
    standardize = standardize,
    newx = newx,
    multiplier = multiplier,
    ...
  )
  
  train <- .capnet_standardize_train(spec)
  m <- nrow(spec$newx)
  p <- spec$p
  
  intercepts <- numeric(m)
  predictions <- numeric(m)
  mu_values <- numeric(m)
  
  betas <- matrix(NA_real_, nrow = m, ncol = p)
  contributions <- matrix(NA_real_, nrow = m, ncol = p)
  
  fitpoints <- unique(c(seq(1, m, by = walk), m + 1))
  
  for (i in seq_len(length(fitpoints) - 1L)) {
    start_row <- fitpoints[i]
    end_row <- fitpoints[i + 1L]
    idx_cap <- start_row:(end_row - 1L)
    m <- length(idx_cap)
    
    cap <- .capnet_cap_context(spec, idx_cap = idx_cap)
    
    mu_try <- mu
    fit <- NULL
    for (k in seq_len(max_mu_tries)) {
      params <- list(alpha = alpha, lambda = lambda, mu = mu_try)
      
      fit_k <- tryCatch(
        .capnet_fit(train, cap, params),
        error = function(e) NULL
      )
      
      if (!is.null(fit_k) && (fit_k$convergence >= 0)) {
        fit <- fit_k
        break
      }
      
      mu_try <- mu_try / 10
    }
    
    if (is.null(fit)) {
      warning(sprintf(
        "walk_capnet: failed to converge for cap slice [%d, %d]; storing NA outputs.",
        start_row, end_row - 1L
      ))
      intercepts[idx_cap] <- NA_real_
      predictions[idx_cap] <- NA_real_
      mu_values[idx_cap] <- NA_real_
      next
    }
    
    model <- .capnet_output(train, cap, fit, params)
    
    intercepts[idx_cap] <- model$a0
    mu_values[idx_cap] <- mu_try
    betas[idx_cap, ] <- matrix(rep(model$beta, m), nrow = m, byrow = TRUE)
    contributions[idx_cap, ] <- model$feature_contributions
    predictions[idx_cap] <- rowSums(model$feature_contributions) + model$a0
  }
  
  if (xts::is.xts(newx)) {
    ord <- zoo::index(newx)
    
    intercepts <- xts::as.xts(intercepts, order.by = ord)
    betas <- xts::as.xts(betas, order.by = ord)
    contributions <- xts::as.xts(contributions, order.by = ord)
    predictions <- xts::as.xts(predictions, order.by = ord)
    mu_values <- xts::as.xts(mu_values, order.by = ord)
    
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
  ), class = "walk_capnet")
}
