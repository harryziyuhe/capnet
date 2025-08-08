#' Perform Cross-Validation for capnet
#'
#' Does grid-search hyperparameter search with cross-validation for capnet. Return
#' the parameter grid and the best parameters pair. Tunes alpha and lambda values.
#' 
#' @param X input matrix
#' @param y response variable
#' @param mu the strength of the contribution cap regularizer
#' @param L contribution ceiling. If set to a single value the same ceiling is
#' applied to all parameters
#' @param lambda the strength of the elasticnet regularizer
#' @param alpha the elasticnet mixing parameter. `alpha=1` collapses to the 
#' LASSO penalty and `alpha=0` the ridge penalty
#' @param newx data used to evaluate and apply the contribution caps. If unspecified
#' default to the input matrix. 
#' @param multiplier A numeric value or vector used for scaling feature contribution
#' @param standardize flag for standardization of the `x` and `y` variables.
#' The coefficients are always returned on the original scale. Default is 
#' `standardize=TRUE`
#' @param splits the user-specified splits index for cross-validation. If
#' unspecified default to k-fold cross-validation
#' @param K the number of folds to split the data into for k-fold cross-validation
#' when no splits index is specified.
#' @param metric the metric used to evaluate cross-validation results. Can be 
#' "mse" or "rsq". Default is `metric="mse"`
#' @param verbose Defaults to 0. Set to 1 to show console outputs cross-validation.
#' 
#' @return A list of the following components:
#' \item {alpha}{Alpha values searched}
#' \item {lambda}{Lambda values searched}
#' \item {cv_errors}{Matrix of cross-validation error for each parameter pair and each fold}
#' \item {mean_errors}{Matrix of average cross-validation error for ech parameter pair}
#' \item {best_alpha}{The best alpha value}
#' \item {best_lambda}{The best lambda value}
#' \item {best_error}{Validation error for the best parameter pair
#' @export

cv.capnet <- function(X, y,
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
  mean_glmnet_errors <- apply(glmnet_errors, c(1, 2), mean, na.rm = TRUE)
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


