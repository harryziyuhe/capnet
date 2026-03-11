#' Extract coefficients from a fitted \code{capnet} model
#' 
#' Returns the estimated intercept and coefficients from a fitted
#' \code{capnet} model object.
#' 
#' @param x A fitted object of class \code{"capnet"}.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return 
#' A named numeric vector (or single-column matrix) containing the intercept
#' followed by the fitted coefficients, labeled with their corresponding
#' variable names.
#' 
#' @seealso [capnet()], [predict()] [predict.capnet()]
#' 
#' @export
#' @method coef capnet
coef.capnet <- function(x, ...) {
  coefs <- matrix(c(x$a0, x$beta), ncol = 1)
  rownames(coefs) <- c("(Intercept)", colnames(x$newx))
  colnames(coefs) <- "beta"
  coefs
}

#' Extract coefficient paths from a walk-forward \code{capnet} fit
#' 
#' @param x A fitted object of class \code{"walk_capnet"} returned by
#'  \code{walk_capnet()}
#' @param index Optional integer vector of rows/steps to return
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return A numeric matrix whose first column is the intercept path and the
#'  remaining columns are coefficient paths.
#' 
#' @seealso [walk_capnet()], [predict()], [predict.walk_capnet()]
#' 
#' @export
#' @method coef walk_capnet
coef.walk_capnet <- function(x, index = NULL, ...) {
  coefs <- cbind(x$intercepts, x$betas)
  if (!is.null(index)) {
    coefs <- coefs[index,]
  }
  coefs
}

#' Predict from a fitted \code{capnet} model
#' 
#' @param object A fitted object of class \code{"capnet"}
#' @param newdata Optional numeric matrix for prediction. If \code{NULL}, uses
#'  \code{object$newx}
#' @param type "link" returns linear predictor eta; "response" returns mean mu.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return Numeric vector of predictions (length = number of rows in 
#'  \code{newx}).
#'  
#' @seealso [capnet()], [predict()]  
#'  
#' @export
#' @method predict capnet
predict.capnet <- function(object, newdata  = NULL, type = c("link", "response"), ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    newx <- object$newx
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
#' @param x A fitted object of class \code{"walk_capnet"}.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @return Numeric vector (or matrix) of predictions stored in x
#' 
#' @seealso [walk_capnet()], [predict()]
#' 
#' @export
#' @method predict walk_capnet
predict.walk_capnet <- function(x, ...) {
  x$predictions
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
#' @return The \code{x} (invisibly). The function prints a ggplot.
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
    return(invisible(object))
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
    return(invisible(object))
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
    return(invisible(object))
  }
  
  invisible(object)
}

#' Plot coefficient paths along a single hyperparameter
#'
#' @param x A data.frame or matrix where the first column is the path
#'   parameter (e.g., \code{lambda}) and the remaining columns are coefficients.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return The \code{ggplot} is printed; the plot x is returned invisibly.
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

#' Print the violation matrix from \code{capnet_violations} object
#' 
#' @param x A fitted object of class \code{"capnet_violations"}.
#' @param ... Further arguments passed to or from other methods.
#' 
#' @importFrom Matrix Matrix printSpMatrix
#' 
#' @return Sparse numeric matrix of contribution cap violations in fitted object x
#' 
#' @seealso [capnet_violations()], [print()]
#' 
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




















