#' Perform cross-validation for \code{capnet}
#'
#' Runs a grid search over \code{alpha} and \code{lambda} with K-fold
#' cross-validation for the \code{capnet()} model. Returns the searched grid,
#' fold-wise errors, mean errors, and the best hyperparameters.
#' 
#' @importFrom stats sd predict
#'
#' @param X Numeric predictor matrix of shape \eqn{n\times p}. Columns are
#'  features and rows are observations.
#' @param y Numeric response vector of length \eqn{n}.
#' @param gamma Nonnegative numeric scalar; strength of the contribution-cap
#'  penalty. Unlike \code{alpha} and \code{lambda}, \code{gamma} is held fixed
#'  during the CV grid search. To compare different \code{gamma} values, call
#'  \code{cv_capnet()} separately for each.
#' @param L Nonnegative numeric scalar or length-\eqn{p} vector giving the
#'  contribution ceiling(s). If scalar, the same ceiling is applied to all
#'  features.
#' @param family Optional character scalar (e.g. "binomial"), function (e.g. 
#' \code{stats::binomial}), or family object (e.g. \code{stats::binomial()}).
#' @param intercept Logical; should an intercept be fitted? Default \code{TRUE}.
#' @param lambda Numeric vector (default \code{exp(seq(1, -5, length.out = 50))})
#'  of nonnegative elastic-net penalty strengths to search.
#' @param alpha Numeric vector (default \code{seq(0, 1, length.out = 5)}) with
#'  values in \eqn{[0,1]} to search (elastic-net mixing parameters;
#'  \code{alpha = 1} is Lasso, \code{alpha = 0} is Ridge).
#' @param z Optional numeric matrix with \eqn{p} columns used to evaluate and
#'  apply contribution caps. If \code{NULL}, defaults to \code{X}.
#' @param multiplier Optional numeric scalar or length-\eqn{n} vector used to
#'  scale feature contributions during the capping step; defaults to 1.
#' @param standardize Logical; if \code{TRUE}, columns of \code{X} and \code{y}
#'  are standardized for fitting; coefficients are returned on the original scale.
#'  Default \code{TRUE}.
#' @param splits Optional integer vector of length \eqn{n} giving the fold
#'  assignment for each row. Each unique integer value defines one fold; rows
#'  with \code{splits == k} are held out when fold \code{k} is the test fold.
#'  Folds need not be contiguous or equal-sized. If \code{NULL}, folds are
#'  assigned randomly using \code{K}. For time-ordered data, supply
#'  \code{splits} explicitly to avoid leakage --- for example:
#'  \code{splits <- rep(1:5, each = n / 5)} assigns the first fifth of rows to
#'  fold 1, the second fifth to fold 2, and so on.
#' @param K Integer \eqn{\ge 2}; number of folds used when \code{splits} is
#'  \code{NULL}. Default \code{5}. If \code{splits} is provided, \code{K} is
#'  set to \code{length(unique(splits))}.
#' @param metric Character string; CV scoring metric. When omitted, the default
#'  is chosen by family: \code{"mse"} for Gaussian (lower is better),
#'  \code{"logloss"} for binomial (lower is better), \code{"deviance"} for
#'  Poisson and Gamma (lower is better). Available options: \code{"mse"},
#'  \code{"rsq"} (higher is better), \code{"logloss"}, \code{"brier"},
#'  \code{"deviance"}.
#' @param verbose Integer; \code{0} for silent, \code{1} to print progress.
#' @param parallel Logical; if \code{TRUE}, the \code{(fold, alpha)} grid is
#'  dispatched across a \code{parallel::makeCluster()} PSOCK cluster via
#'  \code{parallel::parLapply()}. The \code{lambda} path within each
#'  \code{(fold, alpha)} pair is always fit sequentially, warm-starting from
#'  the previous \code{lambda}. Default \code{FALSE}.
#' @param workers Optional integer; number of parallel workers to use when
#'  \code{parallel = TRUE}. Defaults to \code{min(n_tasks, detectCores() - 1,
#'  90)}: capped at the number of (fold x alpha) tasks (idle workers add no
#'  value) and at 90 to avoid exhausting R's 128-connection limit.
#' @param ... Additional arguments forwarded to \code{capnet()}, e.g.,
#'  \code{lower.limits}, \code{upper.limits}, \code{tol}, \code{maxit}.
#' 
#' @return An object of class \code{"cv_capnet"} with components:
#'  \item{\code{alpha}}{Numeric vector of alpha values searched.}
#'  \item{\code{lambda}}{Numeric vector of lambda values searched.}
#'  \item{\code{gamma}}{gamma value passed in input.}
#'  \item{\code{L}}{L value passed in input.}
#'  \item{\code{metric}}{Metric used to evaluate performance.}
#'  \item{\code{splits}}{Cross-validation splits.}
#'  \item{\code{cv_errors}}{Numeric array of shape \eqn{A\times L\times K};
#'    \code{cv_errors[a, l, k]} is the fold-\code{k} error for
#'    \code{alpha[a]}, \code{lambda[l]}. Units match \code{metric}.}
#'  \item{\code{mean_errors}}{Numeric matrix of shape \eqn{A\times L} with
#'    mean CV error across folds for each parameter pair.}
#'  \item{\code{best_alpha}}{Numeric; the selected alpha.}
#'  \item{\code{best_lambda}}{Numeric; the selected lambda.}
#'  \item{\code{best_error}}{Numeric; mean CV score at \code{(best_alpha,
#'    best_lambda)}. For \code{"mse"}, \code{"logloss"}, \code{"deviance"},
#'    and \code{"brier"} this is the minimum; for \code{"rsq"} it is the
#'    maximum.}
#' 
#' @details
#' If \code{splits} is \code{NULL}, folds are generated as
#' \code{sample(rep(1:K, length.out = n))}. Set a seed beforehand to reproduce
#' the random fold allocation. If \code{splits} is supplied as an integer vector
#' of fold IDs, \code{K} is inferred as \code{length(unique(splits))}.
#'
#' For each (\code{alpha}, \code{lambda}) pair, \code{capnet()} is fit on the
#' training portion of each fold and evaluated on the held-out rows using
#' \code{metric}.
#'
#' \strong{Time-series data}: the default random fold creation assumes i.i.d.
#' rows. For time-ordered data supply \code{splits} explicitly so that each
#' test fold lies strictly after its training fold:
#' \preformatted{
#'   # Walk-forward folds: fold k tests the k-th block of rows
#'   splits <- rep(seq_len(K), each = ceiling(n / K))[seq_len(n)]
#' }
#' 
#' @seealso [capnet()], [walk_capnet()], [plot.cv_capnet()]
#'
#' @examples
#' \donttest{
#'   set.seed(1)
#'   n <- 80; p <- 10
#'   X <- matrix(rnorm(n * p), n, p)
#'   beta <- c(1.2, 0.7, 0.5, rep(0, p - 3))
#'   y <- as.numeric(X %*% beta + rnorm(n))
#'
#'   # Basic cross-validation
#'   cv <- cv_capnet(X, y, gamma = 1, L = 1.5)
#'   cv$best_alpha; cv$best_lambda; cv$best_error
#'
#'   # Visualize the CV error surface
#'   plot(cv)
#'   plot(cv, alpha = cv$best_alpha)
#'
#'   # Time-ordered data: use sequential folds to avoid leakage
#'   K <- 5
#'   splits <- rep(seq_len(K), each = ceiling(n / K))[seq_len(n)]
#'   cv_ts <- cv_capnet(X, y, gamma = 1, L = 1.5, splits = splits)
#' }
#' 
#' @export

cv_capnet <- function(X, y,
                      gamma, L,
                      family = "gaussian",
                      lambda = exp(seq(1, -5, length.out = 50)),
                      alpha = seq(0, 1, length.out = 5),
                      z = NULL,
                      multiplier = 1,
                      intercept = TRUE,
                      standardize = TRUE, 
                      splits = NULL,
                      K = 5,
                      metric = c("mse", "rsq", "logloss", "brier", "deviance"),
                      verbose = 0,
                      parallel = FALSE,
                      workers = NULL,
                      ...) {
  metric_missing <- missing(metric) || is.null(metric)
  metric <- match.arg(metric)

  X <- as.matrix(X)
  y <- as.numeric(y)
  
  n <- nrow(X)
  
  # Default to K-folds split if not specified
  if (is.null(splits)) {
    splits <- sample(rep(seq_len(K), length.out = n))
  } else {
    splits <- as.integer(splits)
    K <- length(unique(splits))
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
  
  if (metric_missing) {
    f <- tolower(spec$family$family)
    metric <- switch(f,
      gaussian = "mse",
      binomial = "logloss",
      poisson = "deviance",
      gamma = "deviance",
      "mse"
    )
  }
  
  cv_errors <- array(
    NA_real_,
    dim = c(length(alpha), length(lambda), K),
    dimnames = list(
      paste0("alpha_", format(alpha, digits = 3)),
      paste0("lambda_", format(lambda, digits = 3)),
      paste0("fold_", seq_len(K))
    )
  )
  
  # Precompute per-fold contexts once; shared read-only across that fold's
  # alpha tasks below.
  fold_contexts <- lapply(seq_len(K), function(fold) {
    idx_train <- which(splits != fold)
    idx_test <- which(splits == fold)

    list(
      fold = fold,
      X_test = spec$X[idx_test, , drop = FALSE],
      y_test = spec$y[idx_test],
      train_fold = .capnet_standardize_train(spec, idx_train = idx_train),
      cap_fold = .capnet_cap_context(spec)
    )
  })

  # (fold, alpha) pairs are independent; the lambda path within each pair is
  # warm-started sequentially (lambda to lambda), same as before.
  tasks <- unlist(lapply(seq_len(K), function(fold) {
    lapply(seq_along(alpha), function(a) list(fold = fold, a = a))
  }), recursive = FALSE)

  if (verbose >= 1) {
    message(sprintf("cv_capnet: fitting %d (fold x alpha) tasks...", length(tasks)))
  }

  run_task <- function(task) {
    fc <- fold_contexts[[task$fold]]
    train_fold <- fc$train_fold
    cap_fold <- fc$cap_fold
    a <- task$a

    errs <- rep(NA_real_, length(lambda))
    par0 <- train_fold$par

    for (l in seq_along(lambda)) {
      train_run <- train_fold
      train_run$par <- par0

      params <- list(alpha = alpha[a], lambda = lambda[l], gamma = gamma)

      fit <- tryCatch(
        .capnet_fit(train_run, cap_fold, params),
        error = function(e) NULL
      )

      if (is.null(fit)) {
        warning(sprintf(
          "cv_capnet: failed (fold=%d, alpha=%.3f, lambda=%.4g); storing NA.",
          fc$fold, alpha[a], lambda[l]
        ))
        next
      }

      model <- .capnet_output(train_fold, cap_fold, fit, params)

      preds <- predict(model, newdata = fc$X_test, type = "response")
      errs[l] <- .cv_capnet_error(fc$y_test, preds, train_fold$family, metric)

      par0 <- c(fit$a0, fit$beta)
    }

    list(fold = fc$fold, a = a, errs = errs)
  }

  results <- .capnet_lapply(
    tasks,
    run_task,
    parallel = parallel,
    workers = workers
  )

  for (res in results) {
    cv_errors[res$a, , res$fold] <- res$errs
  }

  mean_errors <- apply(cv_errors, c(1, 2), mean, na.rm = TRUE)
  
  if (metric == "mse") {
    best <- which(mean_errors == min(mean_errors, na.rm = TRUE), arr.ind = TRUE)[1, ]
  } else {
    best <- which(mean_errors == max(mean_errors, na.rm = TRUE), arr.ind = TRUE)[1, ]
  }
  
  structure(list(
    alpha = alpha,
    lambda = lambda,
    gamma = gamma,
    L = L,
    metric = metric,
    splits = splits,
    cv_errors = cv_errors,
    mean_errors = mean_errors,
    best_alpha = alpha[best[1]],
    best_lambda = lambda[best[2]],
    best_error = mean_errors[best[1], best[2]],
    call = match.call()
  ), class = "cv_capnet")
}

.cv_capnet_error <- function(y, pred, family, metric) {
  f <- tolower(family$family)
  
  if (metric == "mse") {
    return(mean((y - pred)^2))
  }
  
  if (metric == "rsq") {
    return(1 - sum((y - pred)^2) / sum((y - mean(y))^2))
  }
  
  if (metric == "logloss") {
    eps <- 1e-15
    p <- pmin(pmax(pred, eps), 1 - eps)
    return(-mean(y * log(p) + (1 - y) * log(1 - p)))
  }
  
  if (metric == "brier") {
    return(mean((y - pred) ^ 2))
  }
  
  if (metric == "deviance") {
    if (!is.null(family$dev.resids)) {
      d <- family$dev.resids(y, pred, wt = rep(1, length(y)))
      return(mean(d))
    }
    stop("Family does not provide dev.resids; cannot compute deviance.")
  }
  
  stop("Unsupported metric: ", metric)
}
