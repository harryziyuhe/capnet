#' Perform walk-forward model fitting with contribution caps
#'
#' Fits a linear elastic net model with an additional contribution-cap penalty
#' once on the full training data \code{X}, \code{y}, then walks the
#' contribution-cap penalty forward across the evaluation matrix \code{z},
#' re-solving the objective at each step with the penalty applied only to the
#' next \code{walk} row(s) of \code{z} (rather than to all of \code{z} at
#' once, as a single call to \code{\link{capnet}} would do). The training data
#' is not expanded or refit incrementally as \code{z} is walked; this mirrors
#' a deployment setting where \code{X}, \code{y} are a data snapshot (e.g. a
#' rolling training window) already prepared upstream, and repeating that
#' preprocessing for every new evaluation row would be prohibitively
#' expensive. To simulate a fully expanding training window instead, call
#' \code{\link{capnet}} directly in a loop with an updated \code{X}, \code{y}
#' at each step.
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
#' @param z Numeric matrix with \eqn{p} columns used to evaluate and
#'  enforce contribution caps (required).
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
#' @param gamma Nonnegative numeric scalar; strength of the contribution-cap penalty.
#' @param max_gamma_tries Optional numeric scalar; maximum number of tries of different
#'  gamma values for convergence. See \link[=walk_capnet]{Details}.
#' @param parallel Logical; if \code{TRUE}, walk-forward steps are dispatched
#'  across a \code{parallel::makeCluster()} PSOCK cluster via
#'  \code{parallel::parLapply()}. Each step fits independently, so this does
#'  not change results. Default \code{FALSE}.
#' @param workers Optional integer; number of parallel workers to use when
#'  \code{parallel = TRUE}. Defaults to \code{parallel::detectCores() - 1}.
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
#'    the evaluation rows, on the response scale (i.e. \code{family$linkinv(eta)});
#'    shape \eqn{S\times p}.}
#'  \item{\code{gammas}}{Numeric vector of length \eqn{n_\mathrm{new}} for the 
#'    \code{gamma} used at each step.}
#' 
#' @details
#' Given the evaluation matrix \code{z} of size \eqn{S\times p}, \code{X} and
#' \code{y} are standardized (if requested) once, up front, and held fixed
#' across all steps. At step \eqn{s=1,\dots,S}, \code{capnet()} is called
#' internally on this fixed training fit, with the contribution-cap penalty
#' restricted to rows \eqn{(s-1)\times\text{walk}:s\times\text{walk}} of
#' \code{z}; because the cap penalty (and its gradient) depend on the
#' evaluation slice, coefficients are genuinely re-optimized at each step, not
#' merely re-evaluated. Penalizing a smaller, local slice of \code{z} instead
#' of the whole evaluation matrix at once (as plain \code{capnet()} would)
#' avoids the over-penalization that comes from forcing a single fit to
#' satisfy contribution caps uniformly across every evaluation row
#' simultaneously.
#'
#' @seealso [capnet()], [predict.capnet()], [coef.capnet()]
#' 
#' @examples
#' set.seed(1)
#' n <- 60; p <- 6; n_new <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' z <- matrix(rnorm(n_new * p), n_new, p)
#' beta <- c(2.5, 1.5, 0.8, rep(0, p - 3))
#' y <- as.numeric(X %*% beta + rnorm(n))
#' out <- walk_capnet(X, y, L = 0.5, z = z, lambda = 0.1, alpha = 0.5, gamma = 1, walk = 1)
#' 
#' @export

walk_capnet <- function(X, y, L, z,
                        family = "gaussian",
                        intercept = TRUE, 
                        standardize = TRUE,
                        multiplier = 1,
                        walk = 1,
                        lambda = 0,
                        alpha = 0,
                        gamma = 0,
                        max_gamma_tries = 6,
                        parallel = FALSE,
                        workers = NULL,
                        ...) {
  
  # Stop if there is any NA values in data
  if (anyNA(X) || anyNA(y) || anyNA(z)) {
    stop("X, y, or z contains NA values")
  }
  
  spec <- .capnet_spec(
    X = X, y = y, L = L,
    family = family,
    intercept = intercept,
    standardize = standardize,
    z = z,
    multiplier = multiplier,
    ...
  )
  
  train <- .capnet_standardize_train(spec)
  m <- nrow(spec$z)
  p <- spec$p
  
  intercepts <- numeric(m)
  predictions <- numeric(m)
  gamma_values <- numeric(m)
  
  betas <- matrix(NA_real_, nrow = m, ncol = p)
  contributions <- matrix(NA_real_, nrow = m, ncol = p)
  
  fitpoints <- unique(c(seq(1, m, by = walk), m + 1))

  # Each walk-forward step fits independently from the shared `train` object
  # (no warm start carried between steps), so steps can run in parallel.
  run_step <- function(i) {
    start_row <- fitpoints[i]
    end_row <- fitpoints[i + 1L]
    idx_cap <- start_row:(end_row - 1L)
    m_i <- length(idx_cap)

    cap <- .capnet_cap_context(spec, idx_cap = idx_cap)

    gamma_try <- gamma
    fit <- NULL
    for (k in seq_len(max_gamma_tries)) {
      params <- list(alpha = alpha, lambda = lambda, gamma = gamma_try)

      fit_k <- tryCatch(
        .capnet_fit(train, cap, params),
        error = function(e) NULL
      )

      # -1001 is libLBFGS's LBFGSERR_ROUNDING_ERROR: the line search could not
      # find further improvement to floating-point precision. This is a benign
      # near-convergence stop (not divergence) and is common on problems with
      # highly collinear features; treat it as an acceptable fit rather than
      # retrying with a smaller gamma.
      if (!is.null(fit_k) && (fit_k$convergence >= 0 || fit_k$convergence == -1001)) {
        fit <- fit_k
        break
      }

      gamma_try <- gamma_try / 10
    }

    if (is.null(fit)) {
      warning(sprintf(
        "walk_capnet: failed to converge for cap slice [%d, %d]; storing NA outputs.",
        start_row, end_row - 1L
      ))
      return(list(idx_cap = idx_cap, ok = FALSE))
    }

    model <- .capnet_output(train, cap, fit, params)

    list(
      idx_cap = idx_cap,
      ok = TRUE,
      a0 = model$a0,
      gamma = gamma_try,
      beta = model$beta,
      feature_contributions = model$feature_contributions,
      m_i = m_i
    )
  }

  step_results <- .capnet_lapply(
    seq_len(length(fitpoints) - 1L),
    run_step,
    parallel = parallel,
    workers = workers
  )

  for (res in step_results) {
    idx_cap <- res$idx_cap

    if (!res$ok) {
      intercepts[idx_cap] <- NA_real_
      predictions[idx_cap] <- NA_real_
      gamma_values[idx_cap] <- NA_real_
      next
    }

    intercepts[idx_cap] <- res$a0
    gamma_values[idx_cap] <- res$gamma
    betas[idx_cap, ] <- matrix(rep(res$beta, res$m_i), nrow = res$m_i, byrow = TRUE)
    contributions[idx_cap, ] <- res$feature_contributions
    eta <- rowSums(res$feature_contributions) + res$a0
    predictions[idx_cap] <- train$family$linkinv(eta)
  }

  if (xts::is.xts(z)) {
    ord <- zoo::index(z)
    
    intercepts <- xts::as.xts(intercepts, order.by = ord)
    betas <- xts::as.xts(betas, order.by = ord)
    contributions <- xts::as.xts(contributions, order.by = ord)
    predictions <- xts::as.xts(predictions, order.by = ord)
    gamma_values <- xts::as.xts(gamma_values, order.by = ord)
    
    colnames(intercepts) <- "intercept"
    colnames(betas) <- colnames(contributions) <- colnames(X)
    colnames(predictions) <- "prediction"
    colnames(gamma_values) <- "gamma"
  } else {
    intercepts <- matrix(intercepts, ncol = 1)
    predictions <- matrix(predictions, ncol = 1)
    gamma_values <- matrix(gamma_values, ncol = 1)
    
    colnames(intercepts) <- "intercept"
    colnames(betas) <- colnames(contributions) <- colnames(X)
    colnames(predictions) <- "prediction"
    colnames(gamma_values) <- "gamma"
  }
  
  structure(list(
    intercepts = intercepts,
    betas = betas,
    feature_contributions = contributions,
    predictions = predictions,
    gammas = gamma_values,
    call = match.call()
  ), class = "walk_capnet")
}
