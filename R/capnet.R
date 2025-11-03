#' Fit a linear model with elastic net and contribution cap regularization
#' 
#' Fits a penalized linear model by maximum likelihood with an elastic net penalty
#' and an additional penalty that caps per-feature contributions. The optimizer
#' uses limited memory BFGS (L-BFGS) when there is no L1 component and
#' Orthant-Wise Limited-Memory Quasi-Newton (OWL-QN) when L1 regularization is
#' present. 
#' 
#' @importFrom lbfgs lbfgs
#' 
#' @param X Numeric predictor matrix of shape \eqn{n\times p}. Columns are
#'  features and rows are observations.
#' @param y Numeric response vector of length \eqn{n}.
#' @param lambda Nonnegative numeric scalar; overall strength of the elastic-net 
#'  penalty. 
#' @param alpha Numeric scalar in \eqn{[0,1]}; elastic net mixing parameter.
#'  \code{alpha = 1} is Lasso, \code{alpha=0} is Ridge.
#' @param mu Nonnegative numeric scalar; strength of the contribution-cap penalty
#' @param L Nonnegative numeric scalar or length-\eqn{p} vector giving the 
#'  contribution ceiling(s). If scalar, the same ceiling is applied to all
#'  coefficients
#' @param newx Optional numeric matrix with \eqn{p} columns used to evaluate and 
#'  enforce contribution caps. If \code{NULL}, defaults to \code{X}.
#' @param par Optional numeric vector of length \eqn{p} with initial coefficient
#'  values. If \code{NULL}, uses zero initialization. 
#' @param multiplier Optional numeric scalar or length-\eqn{n} vector used to
#'  scale feature contributions during the capping step. Defaults to 1.
#' @param intercept Logical; should an intercept be fitted? Default \code {TRUE}.
#' @param standardize Logical; if \code{TRUE}, columns of \code{X} and \code{y}
#'  are standardized for fitting; coefficients are returned on the original scale.
#'  Default \code{TRUE}.
#' @param lower.limits Optional numeric scalar or length-\eqn{p} vector of lower
#'  bounds on coefficients. Initial values must satisfy the bounds.
#' @param upper.limits Optional numeric scalar or length-\eqn{p} vector of upper
#'  bounds on coefficients. Initial values must satisfy the bounds. 
#' @param tol Nonnegative numeric tolerance used for gradient masking when 
#'  \code{lower.limits} or \code{upper.limits} are specified. Default \code{1e-8}.
#' @param maxit Integer; maximum number of quasi-Newton iterations. Default 
#'  \code{1e5}.
#' @param check.finite Logical; if \code{TRUE}, detect and report non-finite
#'  values encountered in the loss/gradient. Default \code{FALSE}.
#' @param verbose Integer; \code{0} for silent, \code{1} to print optimization
#'  progress. Default \code{0}.
#' 
#' @return An object of class \code{"capnet"} with components:
#' \itemize{
#'  \item \code{a0} Best intercept (numeric scalar).
#'  \item \code{beta} Numeric vector (length \eqn{p}); fitted coefficients.
#'  \item \code{value} Numeric; minimized objective value.
#'  \item \code{feature_contributions} Numeric matrix of shape
#'    \eqn{n_{\mathrm{new}}\times p} giving per-feature contributions evaluated
#'    on \code{newx} (rows) for each feature (columns).
#'  \item \code{newx} The evaluation matrix.
#'  \item \code{convergence} Integer code; \code{0} indicates successful
#'    convergence, negative values indicate OWL-QN/L-BFGS execution errors.
#'  \item \code{message} Character string describing any optimizer message 
#'    (present if \code{convergence != 0})
#'  \item \code{alpha}, \code{lambda}, \code{mu}, \code{L}, \code{multiplier}
#'    Echoed tuning parameters.
#'  \item \code{call} The matched call.
#' }
#' 
#' @details
#' When \code{alpha > 0} and \code{\lambda > 0}, OWL-QN is used to handle the L1
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
    value = value,
    feature_contributions = contributions,
    newx = newx,
    convergence = model$convergence,
    message = model$message,
    alpha = alpha,
    lambda = lambda,
    mu = mu,
    L = L,
    multiplier = multiplier,
    call = match.call()
  ), class = c("capnet"))
}
