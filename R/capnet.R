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
#'  contribution ceiling(s). If scalar, the same ceiling is applied to all
#'  coefficients
#' @param family Optional character scalar (e.g. "binomial"), function (e.g. 
#' \code{stats::binomial}), or family object (e.g. \code{stats::binomial()}).
#' @param intercept Logical; should an intercept be fitted? Default \code{TRUE}.
#' @param standardize Logical; if \code{TRUE}, columns of \code{X} and \code{y}
#'  are standardized for fitting; coefficients are returned on the original scale.
#'  Default \code{TRUE}.
#' @param newx Optional numeric matrix with \eqn{p} columns used to evaluate and 
#'  enforce contribution caps. If \code{NULL}, defaults to \code{X}.
#' @param multiplier Optional numeric scalar or length-\eqn{n} vector used to
#'  scale feature contributions during the capping step. Defaults to 1.
#' @param lambda Nonnegative numeric scalar; overall strength of the elastic-net 
#'  penalty. 
#' @param alpha Numeric scalar in \eqn{[0,1]}; elastic net mixing parameter.
#'  \code{alpha = 1} is Lasso, \code{alpha=0} is Ridge.
#' @param mu Nonnegative numeric scalar; strength of the contribution-cap penalty
#' @param lower.limits Optional numeric scalar or length-\eqn{p} vector of lower
#'  bounds on coefficients. Initial values must satisfy the bounds.
#' @param upper.limits Optional numeric scalar or length-\eqn{p} vector of upper
#'  bounds on coefficients. Initial values must satisfy the bounds. 
#' @param tol Nonnegative numeric tolerance used for gradient masking when 
#'  \code{lower.limits} or \code{upper.limits} are specified. Default \code{1e-8}.
#' @param maxit Integer; maximum number of quasi-Newton iterations. Default 
#'  \code{1e5}.
#' @param par Optional numeric vector of length \eqn{p} with initial coefficient
#'  values. If \code{NULL}, uses zero initialization. 
#' @param ... Additional arguments used in fitting. Currently unused.
#' 
#' @return An object of class \code{"capnet"} with components:
#'  \item{\code{a0}}{Best intercept (numeric scalar).}
#'  \item{\code{beta}}{Numeric vector (length \eqn{p}); fitted coefficients.}
#'  \item{\code{value}}{Numeric; minimized objective value.}
#'  \item{\code{feature_contributions}}{Numeric matrix of shape
#'    \eqn{n_{\mathrm{new}}\times p} giving per-feature contributions evaluated
#'    on \code{newx} (rows) for each feature (columns).}
#'  \item{\code{newx}}{The evaluation matrix.}
#'  \item{\code{convergence}}{Integer code; \code{0} indicates successful
#'    convergence, negative values indicate OWL-QN/L-BFGS execution errors.}
#'  \item{\code{message}}{Character string describing any optimizer message 
#'    (present if \code{convergence != 0})}
#'  \item{\code{alpha}}{alpha value passed in input.}
#'  \item{\code{lambda}}{lambda value passed in input.}
#'  \item{\code{mu}}{mu value passed in input.}
#'  \item{\code{L}}{L value passed in input.}
#'  \item{\code{multiplier}}{multiplier value passed in input.}
#'  \item{\code{family}}{model family passed in input.}
#'  \item{\code{call}}{The matched call.}
#' 
#' @details
#' When \code{alpha > 0} and \code{lambda > 0}, OWL-QN is used to handle the L1
#' component; otherwise L-BFGS is used. Box constraints are enforced via masked
#' gradients with tolerance \code{tol}. If \code{standardize = TRUE}, the model
#' is fit on standardized variables and coefficients are mapped back to the
#' original scale on return.
#' 
#' Note that standardization is not recommended when a non-unit \code{multiplier}
#' is supplied, as the scaling step may distort the intended effect of the 
#' multiplier on feature contributions.
#' 
#' @seealso [predict.capnet()], [coef.capnet()]
#' 
#' @examples
#' set.seed(1)
#' n <- 40; p <- 8
#' X <- matrix(rnorm(n * p), n, p)
#' beta <- c(2, 1.5, rep(0, p-2))
#' y <- X %*% beta + rnorm(n)
#' fit1 <- capnet(X, as.numeric(y), lambda = 0.1, alpha = 0.5, mu = 1, L = 2)
#' 
#' # Box constraints
#' fit2 <- capnet(X, as.numeric(y), lambda = 0.1, alpha = 0.5, mu = 1, L = 2, lower.limits = 0)
#' 
#' @export


# Build intercept part later
capnet <- function(X, y, L,
                   family = "gaussian",
                   intercept = TRUE,
                   standardize = TRUE,
                   newx = NULL,
                   multiplier = 1,
                   lambda = 0,
                   alpha = 0,
                   mu = 0,
                   lower.limits = NULL,
                   upper.limits = NULL,
                   tol = 1e-8,
                   maxit = 10000L,
                   par = NULL,
                   ...) {
  
  if (anyNA(X) || anyNA(y) || anyNA(newx)) {
    stop("X or y or newx contains NA values")
  }
  
  call <- match.call()
  
  spec <- .capnet_spec(
    X = X, y = y, L = L,
    family = family,
    intercept = intercept,
    standardize = standardize,
    newx = newx,
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
  
  params <- list(alpha = alpha, lambda = lambda, mu = mu)
  
  fit <- .capnet_fit(train, cap, params)
  
  output <- .capnet_output(train, cap, fit, params, call = call)
  
  return(output)
}
