#' Extract coefficients from a fitted \code{capnet} model
#'
#' Returns the estimated intercept and coefficients from a fitted
#' \code{capnet} model object.
#'
#' @param object A fitted object of class \code{"capnet"}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A \eqn{(1+p)\times 1} numeric matrix. The first row is the
#'  intercept (named \code{"(Intercept)"}); the remaining \eqn{p} rows are
#'  the fitted slopes, named from \code{colnames(X)} when available.
#'
#' @seealso [capnet()], [predict.capnet()]
#'
#' @examples
#' set.seed(1)
#' n <- 40; p <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' y <- as.numeric(X %*% c(1, -1, 0.5, 0, 0) + rnorm(n))
#' fit <- capnet(X, y, lambda = 0.1, alpha = 0.5, gamma = 1, L = 1)
#' coef(fit)          # (1 + p) x 1 matrix
#' coef(fit)[-1, ]    # slopes only
#'
#' @export
#' @method coef capnet
coef.capnet <- function(object, ...) {
  coefs <- matrix(c(object$a0, object$beta), ncol = 1)
  beta_names <- names(object$beta) %||% paste0("V", seq_along(object$beta))
  rownames(coefs) <- c("(Intercept)", beta_names)
  colnames(coefs) <- "beta"
  coefs
}

#' Extract coefficient paths from a walk-forward \code{capnet} fit
#'
#' @param object A fitted object of class \code{"walk_capnet"} returned by
#'  \code{walk_capnet()}.
#' @param index Optional integer vector of step indices (1 to \code{nrow(z)})
#'  to subset. If \code{NULL}, all steps are returned.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A \eqn{\mathrm{nrow}(z)\times(1+p)} matrix (or \code{xts} if the
#'  input \code{z} was \code{xts}). The first column is the intercept path
#'  (named \code{"intercept"}); the remaining \eqn{p} columns are coefficient
#'  paths named from \code{colnames(X)}. Rows within the same walk step share
#'  identical values. \code{NA} rows indicate steps that failed to converge.
#'  If \code{index} is supplied, only those rows are returned.
#'
#' @seealso [walk_capnet()], [predict.walk_capnet()]
#'
#' @examples
#' set.seed(1)
#' n <- 50; p <- 5; n_new <- 8
#' X <- matrix(rnorm(n * p), n, p)
#' z <- matrix(rnorm(n_new * p), n_new, p)
#' y <- as.numeric(X %*% c(1, -1, 0.5, 0, 0) + rnorm(n))
#' fit <- walk_capnet(X, y, L = 1, z = z, lambda = 0.1, alpha = 0.5, gamma = 1)
#' coef(fit)           # nrow(z) x (1 + p) path
#' coef(fit, index = 1:3)  # first three steps only
#'
#' @export
#' @method coef walk_capnet
coef.walk_capnet <- function(object, index = NULL, ...) {
  coefs <- cbind(object$intercepts, object$betas)
  if (!is.null(index)) {
    coefs <- coefs[index,]
  }
  coefs
}

#' Predict from a fitted \code{capnet} model
#'
#' @param object A fitted object of class \code{"capnet"}.
#' @param newdata Optional numeric matrix with \eqn{p} columns for prediction.
#'  If \code{NULL}, uses \code{object$z} (the evaluation matrix stored at fit
#'  time).
#' @param type Character; \code{"link"} returns the linear predictor
#'  \eqn{\hat\eta = \hat\beta_0 + X\hat\beta}; \code{"response"} applies the
#'  inverse link function and returns the fitted mean \eqn{\hat\mu}. For
#'  Gaussian models the two are identical. Default \code{"link"}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Numeric vector of length \code{nrow(newdata)} (or \code{nrow(z)} if
#'  \code{newdata} is \code{NULL}).
#'
#' @details
#' When \code{newdata = NULL}, predictions are made on \code{object$z}, which
#' is the evaluation matrix used when the model was fit. This is convenient for
#' inspecting fitted contributions without re-specifying the data. When
#' \code{newdata} is supplied it must have exactly \eqn{p} columns (matching
#' \code{length(object$beta)}) but may have any number of rows.
#'
#' @seealso [capnet()], [coef.capnet()]
#'
#' @examples
#' set.seed(1)
#' n <- 40; p <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' y <- as.numeric(X %*% c(1, -1, 0.5, 0, 0) + rnorm(n))
#' fit <- capnet(X, y, lambda = 0.1, alpha = 0.5, gamma = 1, L = 1)
#' predict(fit, type = "response")            # fitted values on training z
#' predict(fit, newdata = X, type = "link")   # linear predictor on new data
#'
#' @export
#' @method predict capnet
predict.capnet <- function(object, newdata  = NULL, type = c("link", "response"), ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    newx <- object$z
  } else {
    newx <- as.matrix(newdata)
  }
  
  if (is.null(object$beta) || is.null(object$a0)) {
    stop("capnet object missing coefficients.")
  }
  if (ncol(newx) != length(object$beta)) {
    stop("newdata must have ncol equal to length(coef slopes).")
  }
  
  eta <- object$a0 + newx %*% object$beta
  
  if (type == "link") return(eta)
  
  family <- object$family
  if (is.character(family)) {
    family <- normalize_family(family)
  }
  if (is.null(family$linkinv)) {
    stop("Family object does not have linkinv().")
  }
  
  mu <- family$linkinv(eta)
  mu
}

#' Predict from a walk-forward \code{capnet} fit
#'
#' @param object A fitted object of class \code{"walk_capnet"}.
#' @param ... Currently unused.
#'
#' @return A \eqn{\mathrm{nrow}(z)\times 1} matrix (or \code{xts} if the
#'  input \code{z} was \code{xts}) of predictions on the response scale,
#'  as computed during the walk-forward evaluation. \code{NA} entries indicate
#'  steps where the optimizer failed to converge. To extract the full
#'  coefficient path, use \code{coef(object)}.
#'
#' @seealso [walk_capnet()], [coef.walk_capnet()]
#'
#' @examples
#' set.seed(1)
#' n <- 50; p <- 5; n_new <- 8
#' X <- matrix(rnorm(n * p), n, p)
#' z <- matrix(rnorm(n_new * p), n_new, p)
#' y <- as.numeric(X %*% c(1, -1, 0.5, 0, 0) + rnorm(n))
#' fit <- walk_capnet(X, y, L = 1, z = z, lambda = 0.1, alpha = 0.5, gamma = 1)
#' predict(fit)
#'
#' @export
#' @method predict walk_capnet
predict.walk_capnet <- function(object, ...) {
  object$predictions
}

#' Plot cross-validation results for \code{cv_capnet}
#'
#' @param x An object of class \code{"cv_capnet"} returned by 
#'  \code{cv_capnet()}.
#' @param alpha Optional numeric; if provided, show CV errors vs \code{lambda} 
#'  at this alpha.
#' @param lambda Optional numeric; if provided, show CV errors vs \code{alpha} 
#'  at this lambda.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return The \code{ggplot} object, returned invisibly. As a side effect the
#'  plot is printed to the current graphics device.
#'
#' @details
#' Three display modes depending on which arguments are supplied:
#' \describe{
#'   \item{Neither \code{alpha} nor \code{lambda}}{Heatmap of ranked mean CV
#'     errors across the full \eqn{\alpha\times\lambda} grid. Lighter fill =
#'     better rank. Useful for identifying the most promising region before
#'     inspecting slices.}
#'   \item{\code{alpha} specified}{Line plot of mean CV error (± 1 SE across
#'     folds) vs \eqn{\log\lambda} at the given \eqn{\alpha}. A dashed
#'     vertical line marks the best \eqn{\lambda}.}
#'   \item{\code{lambda} specified}{Line plot of mean CV error (± 1 SE) vs
#'     \eqn{\alpha} at the given \eqn{\lambda}.}
#' }
#' Supply at most one of \code{alpha} and \code{lambda}; both must be values
#' that appear in the searched grid.
#'
#' @examples
#' \donttest{
#'   set.seed(1)
#'   n <- 80; p <- 10
#'   X <- matrix(rnorm(n * p), n, p)
#'   y <- as.numeric(X %*% c(rep(1, 3), rep(0, p - 3)) + rnorm(n))
#'   cv <- cv_capnet(X, y, gamma = 1, L = 1.5)
#'   plot(cv)                        # full heatmap
#'   plot(cv, alpha = cv$best_alpha) # error vs lambda slice
#'   plot(cv, lambda = cv$best_lambda) # error vs alpha slice
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_tile labs theme_minimal scale_fill_viridis_c
#' @importFrom ggplot2 stat_summary geom_vline mean_se
#' @importFrom rlang .data
#' @export
#' @method plot cv_capnet
plot.cv_capnet <- function(x, alpha = NULL, lambda = NULL, ...) {
  object <- x
  cv_errors  <- object$cv_errors        # A x L x K array
  mean_errors <- object$mean_errors     # A x L matrix
  alphas     <- object$alpha
  lambdas    <- object$lambda
  
  dims <- dim(cv_errors)
  A <- dims[1]; L <- dims[2]; K <- dims[3]
  
  if (is.null(alpha) && is.null(lambda)) {
    # heatmap of mean errors (ranked)
    df <- expand.grid(alpha = alphas, lambda = log10(lambdas))
    df$value <- as.vector(mean_errors)
    r <- rank(df$value, ties.method = "average")
    df$rank_scaled <- (r - min(r)) / (max(r) - min(r))
    p <- ggplot(df, aes(x = .data$alpha, 
                        y = .data$lambda, 
                        fill = .data$rank_scaled)) +
      geom_tile() +
      scale_fill_viridis_c() +
      labs(title = "Ranked Mean CV Errors",
           x = expression(alpha), y = expression(log[10](lambda)), fill = "Rank") +
      theme_minimal()
    print(p)
    return(invisible(p))
  }
  
  if (!is.null(alpha)) {
    if (!(alpha %in% alphas)) stop("alpha value not evaluated during cross validation")
    ia <- which(alphas == alpha)
    errs <- cv_errors[ia, , , drop = TRUE]  # L x K
    tot <- rowSums(errs)
    lambda.min <- lambdas[which.min(tot)]
    bar_width <- (log(max(lambdas)) - log(min(lambdas))) / (length(lambdas) - 1)
    
    df <- data.frame(
      lambda = rep(log(lambdas), times = K),
      fold   = rep(seq_len(K), each = L),
      error  = as.vector(errs)
    )
    
    p <- ggplot(df, aes(x = .data$lambda, y = .data$error)) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = bar_width) +
      stat_summary(fun = mean, geom = "point", size = 2) +
      geom_vline(xintercept = log(lambda.min), linetype = 2) +
      labs(
        title = bquote("Cross-validation errors for " * alpha == .(alpha)),
        x = expression(log(lambda)), y = "CV error"
      ) +
      theme_minimal()
    print(p)
    return(invisible(p))
  }
  
  if (!is.null(lambda)) {
    if (!(lambda %in% lambdas)) stop("lambda value not evaluated during cross validation")
    il <- which(lambdas == lambda)
    errs <- cv_errors[, il, , drop = TRUE]  # A x K
    tot <- rowSums(errs)
    alpha.min <- alphas[which.min(tot)]
    bar_width <- (max(alphas) - min(alphas)) / (length(alphas) - 1)
    
    df <- data.frame(
      alpha = rep(alphas, times = K),
      fold  = rep(seq_len(K), each = A),
      error = as.vector(errs)
    )
    
    p <- ggplot(df, aes(x = .data$alpha, 
                        y = .data$error)) +
      stat_summary(fun.data = mean_se, geom = "errorbar", width = bar_width) +
      stat_summary(fun = mean, geom = "point", size = 2) +
      geom_vline(xintercept = alpha.min, linetype = 2) +
      labs(
        title = bquote("Cross-validation errors for " * lambda == .(lambda)),
        x = expression(alpha), y = "CV error"
      ) +
      theme_minimal()
    print(p)
    return(invisible(p))
  }
  
  invisible(object)
}

#' Plot coefficient paths along a single hyperparameter
#'
#' @param x An object of class \code{"capnet_path"} returned by
#'  \code{coef_path()}, or any \code{data.frame} with the same structure:
#'  first column is the path parameter, remaining columns are coefficients.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return The \code{ggplot} object, returned invisibly. As a side effect the
#'  plot is printed to the current graphics device.
#'
#' @details
#' Produces a line plot with the path parameter on the x-axis (log scale for
#' \code{lambda} and \code{gamma} paths, linear for \code{alpha}) and
#' coefficient values on the y-axis. Each feature is drawn as a separate
#' colored line, making it easy to see which coefficients enter, shrink, or
#' change sign along the path.
#'
#' @examples
#' set.seed(1)
#' n <- 50; p <- 6
#' X <- matrix(rnorm(n * p), n, p)
#' colnames(X) <- paste0("x", seq_len(p))
#' y <- as.numeric(X %*% c(1.5, -1, 0.5, 0, 0, 0) + rnorm(n))
#' path <- coef_path(X, y, L = 0.5, alpha = 0.5,
#'                   lambda = exp(seq(1, -5, length.out = 30)), gamma = 1)
#' plot(path)
#'
#' @importFrom ggplot2 ggplot aes geom_line labs theme_minimal theme_bw theme element_text
#' @importFrom rlang .data
#' @export
#' @method plot capnet_path
plot.capnet_path <- function(x, ...) {
  object <- x
  if (!is.data.frame(object)) object <- as.data.frame(object)
  param_name <- names(object)[1]
  n <- nrow(object)
  feat_names <- names(object)[-1]
  df_long <- data.frame(
    param  = rep(object[[1]], times = length(feat_names)),
    feature = rep(feat_names, each = n),
    coef   = as.vector(as.matrix(object[-1]))
  )
  p <- ggplot(df_long, aes(x = .data$param, 
                           y = .data$coef, 
                           color = .data$feature)) +
    geom_line(linewidth = 0.8) +
    labs(x = param_name, y = "Coefficient") +
    theme_bw() +
    theme(text = element_text(family = "serif", size = 14))
  print(p)
  invisible(p)
}

#' Print the violation matrix from a \code{capnet_violations} object
#'
#' @param x An object of class \code{"capnet_violations"} returned by
#'  \code{capnet_violations()}.
#' @param ... Optional arguments forwarded to \code{Matrix::printSpMatrix()},
#'  e.g. \code{col.names}, \code{digits}, \code{align}.
#'
#' @return \code{x}, invisibly.
#'
#' @details
#' Displays the excess-contribution matrix as a sparse matrix: only non-zero
#' entries (i.e., actual cap violations) are printed. Each non-zero entry
#' gives the amount by which \eqn{|z_{ij}\hat\beta_j|} exceeded \eqn{L_j}
#' for that row-feature combination.
#'
#' @examples
#' set.seed(1)
#' n <- 40; p <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' y <- as.numeric(X %*% c(2, -1.5, 0, 0, 0) + rnorm(n))
#' fit <- capnet(X, y, lambda = 0.05, alpha = 0.5, gamma = 0.1, L = 0.5)
#' v <- capnet_violations(fit)
#' if (!is.null(v)) print(v)         # sparse display of excess contributions
#' if (!is.null(v)) v$excess[1:5, ]  # raw access
#'
#' @seealso [capnet_violations()], [capnet()]
#'
#' @importFrom Matrix Matrix printSpMatrix
#' @export
#' @method print capnet_violations
print.capnet_violations <- function(x, ...) {
  sparse_x <- Matrix(x$excess, sparse = TRUE)
  
  # Default arguments
  defaults <- list(
    col.names = TRUE,
    align     = "right",
    digits    = 4
  )
  
  # Capture user arguments
  user_args <- list(...)
  
  # Let user override defaults
  defaults[names(user_args)] <- user_args
  
  # Call printSpMatrix safely
  do.call(printSpMatrix, c(list(x = sparse_x), defaults))
}




















