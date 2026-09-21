#' Fit a linear model with elastic net and contribution cap regularization
#' 
#' Fits a penalized linear model by maximum likelihood with an elastic net penalty
#' and an additional penalty that caps per-feature contributions. The optimizer
#' uses limited memory BFGS (L-BFGS) when there is no L1 component and
#' Orthant-Wise Limited-Memory Quasi-Newton (OWL-QN) when L1 regularization is
#' present. 
#' 
#' @importFrom lbfgs lbfgs
#' @importFrom stats sd
#' 
#' @param X Numeric predictor matrix of shape \eqn{n\times p}. Columns are
#'  features and rows are observations.
#' @param y Numeric response vector of length \eqn{n}.
#' @param L Nonnegative numeric scalar or length-\eqn{p} vector giving the
#'  contribution ceiling(s). \eqn{L_j} caps the absolute value of feature
#'  \eqn{j}'s contribution \eqn{|z_{ij}\beta_j|} for each evaluation row \eqn{i}.
#'  If scalar, the same ceiling is applied to all features. When \code{gamma = 0},
#'  \code{L} has no effect on the fit.
#' @param family Optional character scalar (e.g. \code{"binomial"}), function
#'  (e.g. \code{stats::binomial}), or family object (e.g.
#'  \code{stats::binomial()}). Supported families: \code{"gaussian"},
#'  \code{"binomial"}, \code{"poisson"}, \code{"Gamma"} (log link only).
#' @param intercept Logical; should an intercept be fitted? Default \code{TRUE}.
#' @param standardize Logical; if \code{TRUE}, columns of \code{X} are
#'  standardized (zero mean, unit SD) before fitting and coefficients are
#'  returned on the original scale. The cap penalty is always evaluated on the
#'  original scale regardless of this setting. Default \code{TRUE}.
#' @param z Optional numeric matrix with \eqn{p} columns used to evaluate the
#'  contribution cap penalty. The cap is enforced on \eqn{|z_{ij}\beta_j|}, not
#'  on training rows. If \code{NULL}, defaults to \code{X} (caps are enforced on
#'  the training data). Pass a held-out evaluation set to enforce caps on
#'  out-of-sample rows. Supply \code{z} on the same scale as \code{X}.
#' @param multiplier Optional numeric scalar or length-\eqn{m} vector (where
#'  \eqn{m = \mathrm{nrow}(z)}) used to scale feature contributions during the
#'  capping step. Defaults to 1. Useful when each evaluation row requires a
#'  different scaling of its contribution budget.
#' @param lambda Nonnegative numeric scalar; overall strength of the elastic-net
#'  penalty. When \code{lambda = 0}, no elastic-net penalty is applied.
#' @param alpha Numeric scalar in \eqn{[0,1]}; elastic net mixing parameter.
#'  \code{alpha = 1} is LASSO (pure L1), \code{alpha = 0} is Ridge (pure L2).
#'  When \code{alpha > 0} and \code{lambda > 0}, OWL-QN is used; otherwise
#'  L-BFGS is used.
#' @param gamma Nonnegative numeric scalar; strength of the contribution-cap
#'  penalty. When \code{gamma = 0}, no cap penalty is applied and the model
#'  reduces to standard elastic net.
#' @param lower.limits Optional numeric scalar or length-\eqn{p} vector of lower
#'  bounds on coefficients. Initial parameter values must satisfy the bounds;
#'  the default zero initialization satisfies any symmetric bound.
#' @param upper.limits Optional numeric scalar or length-\eqn{p} vector of upper
#'  bounds on coefficients. Initial parameter values must satisfy the bounds.
#' @param tol Nonnegative numeric tolerance used for gradient masking when
#'  \code{lower.limits} or \code{upper.limits} are specified. Default \code{1e-8}.
#' @param maxit Integer; maximum number of quasi-Newton iterations. Default
#'  \code{1e5}.
#' @param par Optional numeric vector of length \eqn{p+1} (intercept first,
#'  then slopes) with initial parameter values. If \code{NULL}, uses zero
#'  initialization.
#' @param ... Additional arguments used in fitting. Currently unused.
#' 
#' @return An object of class \code{"capnet"} with components:
#'  \item{\code{a0}}{Best intercept (numeric scalar).}
#'  \item{\code{beta}}{Numeric vector (length \eqn{p}); fitted coefficients.}
#'  \item{\code{value}}{Numeric; minimized objective value.}
#'  \item{\code{feature_contributions}}{Numeric matrix of shape
#'    \eqn{\mathrm{nrow}(z)\times p} giving per-feature contributions
#'    \eqn{z_{ij}\hat\beta_j} evaluated on \code{z}.}
#'  \item{\code{z}}{The evaluation matrix used for cap enforcement.}
#'  \item{\code{convergence}}{Integer convergence code from the optimizer:
#'    \code{0} = successful convergence; \code{-1001} = rounding-error stop
#'    (benign near-convergence, common with highly collinear features);
#'    \code{-998} = line search failure (try reducing \code{gamma} or
#'    \code{lambda}); other negative values indicate optimizer errors.
#'    The function does \strong{not} raise an error on convergence failure;
#'    always inspect this field after fitting.}
#'  \item{\code{message}}{Character string from the optimizer (may be
#'    \code{NULL} or uninformative for some failure codes).}
#'  \item{\code{alpha}}{alpha value passed in input.}
#'  \item{\code{lambda}}{lambda value passed in input.}
#'  \item{\code{gamma}}{gamma value passed in input.}
#'  \item{\code{L}}{L value passed in input.}
#'  \item{\code{multiplier}}{multiplier value passed in input.}
#'  \item{\code{family}}{model family passed in input.}
#'  \item{\code{call}}{The matched call.}
#' 
#' @details
#' When \code{alpha > 0} and \code{lambda > 0}, OWL-QN is used to handle the
#' L1 component; otherwise L-BFGS is used. Box constraints are enforced via
#' gradient masking with tolerance \code{tol}. If \code{standardize = TRUE},
#' the model is fit on standardized \code{X} and coefficients are mapped back
#' to the original scale on return; the cap penalty is always computed on the
#' original scale.
#'
#' \strong{Soft constraints}: the contribution cap is a \emph{penalty}, not a
#' hard constraint. The returned \code{feature_contributions} may still exceed
#' \code{L} when \code{gamma} is small or the optimizer does not fully converge.
#' Use \code{capnet_violations()} to check how much each cap is exceeded.
#'
#' \strong{Standardization with multiplier}: standardization is not recommended
#' when a non-unit \code{multiplier} is supplied, as the scaling step may
#' distort the intended contribution budget.
#'
#' @seealso [cv_capnet()], [walk_capnet()], [predict.capnet()],
#'   [coef.capnet()], [capnet_violations()]
#' 
#' @examples
#' set.seed(1)
#' n <- 40; p <- 8
#' X <- matrix(rnorm(n * p), n, p)
#' beta <- c(2, 1.5, rep(0, p - 2))
#' y <- as.numeric(X %*% beta + rnorm(n))
#' L <- rep(1.5, p)
#'
#' # Standard elastic net (gamma = 0, no cap)
#' fit0 <- capnet(X, y, lambda = 0.1, alpha = 0.5, gamma = 0, L = L)
#'
#' # With contribution cap penalty
#' fit1 <- capnet(X, y, lambda = 0.1, alpha = 0.5, gamma = 1, L = L)
#' fit1$convergence  # check optimizer status; 0 = success
#'
#' # Enforce caps on a held-out evaluation set, not training rows
#' X_new <- matrix(rnorm(10 * p), 10, p)
#' fit2 <- capnet(X, y, lambda = 0.1, alpha = 0.5, gamma = 1, L = L, z = X_new)
#'
#' # Box constraints (non-negative coefficients)
#' fit3 <- capnet(X, y, lambda = 0.1, alpha = 0.5, gamma = 1, L = L, lower.limits = 0)
#' 
#' @export


# Build intercept part later
capnet <- function(X, y, L,
                   family = "gaussian",
                   intercept = TRUE,
                   standardize = TRUE,
                   z = NULL,
                   multiplier = 1,
                   lambda = 0,
                   alpha = 0,
                   gamma = 0,
                   lower.limits = NULL,
                   upper.limits = NULL,
                   tol = 1e-8,
                   maxit = 10000L,
                   par = NULL,
                   ...) {
  
  if (anyNA(X) || anyNA(y) || anyNA(z)) {
    stop("X, y, or z contains NA values")
  }
  
  call <- match.call()
  
  spec <- .capnet_spec(
    X = X, y = y, L = L,
    family = family,
    intercept = intercept,
    standardize = standardize,
    z = z,
    multiplier = multiplier,
    lower.limits = lower.limits,
    upper.limits = upper.limits,
    tol = tol,
    maxit = maxit,
    par = par,
    ...
  )
  
  train <- .capnet_standardize_train(spec)
  cap <- .capnet_cap_context(spec)
  
  params <- list(alpha = alpha, lambda = lambda, gamma = gamma)
  
  fit <- .capnet_fit(train, cap, params)
  
  output <- .capnet_output(train, cap, fit, params, call = call)
  
  return(output)
}
