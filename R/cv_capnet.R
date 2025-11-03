#' Perform cross-validation for \code{capnet}
#'
#' Runs a grid search over \code{alpha} and \code{lambda} with K-fold
#' cross-validation for the \code{capnet()} model. Returns the searched grid,
#' fold-wise errors, mean errors, and the best hyperparameters.
#' 
#' @param X Numeric predictor matrix of shape \eqn{n\times p}. Columns are
#'  features and rows are observations.
#' @param y Numeric response vector of length \eqn{n}.
#' @param mu Nonnegative numeric scalar; strength of the contribution-cap penalty
#' @param L Nonnegative numeric scalar or length-\eqn{p} vector giving the 
#'  contribution ceiling(s). If scalar, the same ceiling is applied to all
#'  coefficients
#' @param lambda Numeric vector (default \code{exp(seq(1, -5, length.out = 50))})
#'  of nonnegative elastic-net penalty strengths to search.
#' @param alpha Numeric vector (default \code{seq(0, 1, length.out = 5)}) with
#'  values in \eqn{[0,1]} to search (elastic-net mixing parameters;
#'  \code{alpha = 1} is Lasso, \code{alpha = 0} is Ridge).
#' @param newx Optional numeric matrix with \eqn{p} columns used to evaluate and 
#'  apply contribution caps. If \code{NULL}, defaults to \code{X}.
#' @param multiplier Optional numeric scalar or length-\eqn{n} vector used to
#'  scale feature contributions during the capping step; defaults to 1.
#' @param standardize Logical; if \code{TRUE}, columns of \code{X} and \code{y}
#'  are standardized for fitting; coefficients are returned on the original scale.
#'  Default \code{TRUE}.
#' @param splits Optional integer vector of length \eqn{n} with values in
#'  \code{1:K} giving the fold assignment for each row. If \code{NULL},
#'  K-fold CV is created internally using \code{K}.
#' @param K Integer \eqn{\ge 2}; number of folds used when \code{splits} is
#'  \code{NULL}. Default \code{5}. If \code{splits} is provided, \code{K} is
#'  set to \code{length(unique(splits))}.
#' @param metric Character string; CV scoring metric: \code{"mse"} (lower is 
#'  better) or \code{"rsq"} (higher is better). Default \code{"mse"}.
#' @param verbose Integer; \code{0} for silent, \code{1} to print progress.
#' @param ... Additional arguments forwarded to \code{capnet()}, e.g.,
#'  \code{lower.limits}, \code{upper.limits}, \code{tol}, \code{maxit}.
#' 
#' @return An object of class \code{"cv_capnet"} with components:
#' \itemize{
#'  \item \code{alpha} Numeric vector of alpha values searched.
#'  \item \code{lambda} Numeric vector of lambda values searched.
#'  \item \code{cv_errors} Numeric array of shape
#'    \eqn{A\times L\times K} with fold-wise errors for each
#'    \code{alpha} (\eqn{A}) and \code{lambda}(\eqn{L}) pair across folds.
#'  \item \code{mean_errors} Numeric matrix of shape \eqn{A\times L} with
#'    mean CV error across folds for each parameter pair.
#'  \item \code{best_alpha} Numeric; the selected alpha.
#'  \item \code{best_lambda} Numeric; the selected lambda.
#'  \item \code{best_error} Numeric; CV score at the selected pair (mean over
#'    folds, consistent with \code{metric}).
#' }
#' 
#' @details
#' If \code{split} is \code{NULL}, folds are generated as
#' \code{sample(rep(1:K, length.out = n))}. Set a seed beforehand to reproduce
#' the random fold allocation. If \code{splits} is supplied as an integer vector
#' of fold IDs, \code{K} is inferred as \code{length(unique(splits))}.
#' 
#' For each (\code{alpha}, \code{\lambda}) pair, \code{capnet()} is fit on the 
#' training portion of each fold and evaluated on the held-out rows using 
#' \code{metric}.
#' 
#' Note: The default random fold creation assumes i.i.d. rows. For time-ordered
#' data, supply \code{splits} explicitly to ensure leakage-free evaluation.
#' 
#' @seealso [capnet()], [plot.cv_capnet]
#' 
#' @examples
#' set.seed(1)
#' n <- 80; p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' beta <- c(1.2, 0.7, 0.5, rep(0, p - 3))
#' y <- as.numeric(X %*% beta + rnorm(n))
#' cv <- cv_capnet(X, y, mu = 1, L = 1.5)
#' cv$best_alpha; cv$best_lambda; cv$best_error
#' 
#' @export

cv_capnet <- function(X, y,
                      mu, L,
                      lambda = exp(seq(1, -5, length.out = 50)),
                      alpha = seq(0, 1, length.out = 5),
                      newx = NULL, multiplier = 1,
                      standardize = TRUE, 
                      splits = NULL, K = 5,
                      metric = "mse", verbose = 0,
                      ...) {
  n <- nrow(X)
  p <- ncol(X)
  
  # Default to K-folds split if not specified
  if (is.null(splits)) {
    splits <- sample(rep(1:K, length.out = n))
  } else {
    K <- length(unique(splits))
  }
  
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
  
  cv_errors <- array(NA, dim = c(length(alpha), length(lambda), K),
                    dimnames = list(
                      paste0("alpha_", alpha),
                      paste0("lambda_", format(lambda, digits = 3)),
                      paste0("fold_", 1:K)
                     ))
  
  for (split in 1:K) {
    if (verbose >= 1) message("Fold", split, "/", K)
    
    train_idx <- which(splits != split)
    test_idx <- which(splits == split)
    
    X_train <- X_scaled[train_idx, ,drop = FALSE]
    y_train <- y_scaled[train_idx]
    X_test <- X_scaled[test_idx, ,drop = FALSE]
    y_test <- y_scaled[test_idx]
    
    for (a in seq_along(alpha)) {
      alpha_value <- alpha[a]
      par <- NULL
      
      for (l in seq_along(lambda)) {
        lambda_value <- lambda[l]
        fit <- tryCatch({
          capnet(X_train, y_train,
                 alpha = alpha_value,
                 lambda = lambda_value,
                 mu = mu, L = L,
                 par = par,
                 standardize = FALSE,
                 multiplier = multiplier,
                 newx = newx)
        }, error = function(e) {
          warning(sprintf("alpha = %.3f, lambda = %.4f failed: %s", alpha_value, lambda_value, e$message))
          return(NULL)
        })

        if (!is.null(fit)) {
          par <- c(fit$a0, fit$beta)
          preds <- fit$a0 + as.numeric(X_test %*% fit$beta)
          
          error <- switch(metric,
                          mse = mean((y_test - preds)^ 2),
                          rsq = 1 - sum((y_test - preds)^2) / sum((y_test - mean(y_test)) ^ 2),
                          stop("Unsupported metric"))
          
          cv_errors[a, l, split] <- error
        }
      }
    }
  }
  
  mean_errors <- apply(cv_errors, c(1, 2), mean, na.rm = TRUE)
  best_idx <- which(mean_errors == min(mean_errors, na.rm = TRUE), arr.ind = TRUE)[1,]
  
  structure(list(
    alpha = alpha,
    lambda = lambda,
    cv_errors = cv_errors,
    mean_errors = mean_errors,
    best_alpha = alpha[best_idx[1]],
    best_lambda = lambda[best_idx[2]],
    best_error = mean_errors[best_idx[1], best_idx[2]]
  ), class = "cv_capnet")
}


