#' Check contribution-cap violations for a fitted \code{capnet} model
#' 
#' Identifies where feature contributions exceed the specified contribution cap.
#' Returns the indices of offending rows and columns along with the excess
#' contribution matrix.
#' 
#' @param object A fitted \code{capnet} model returned by [capnet()].
#' @param z Optional numeric matrix of new data on which to evaluate
#'  contribution caps. If \code{NULL}, uses the contribution matrix
#'  stored in \code{object$feature_contributions}.
#' @param multiplier Optional numeric scalar or length-\eqn{n} vector (where
#'  \eqn{n} is \code{nrow(z)}) to override the multiplier used when the model
#'  was fit. If \code{NULL}, defaults to \code{object$multiplier}.
#' 
#' @return
#' A list with the following components:
#'  \item{\code{columns}}{Integer vector of column indices (features) with any
#'    violations of the contribution cap.}
#'  \item{\code{rows}}{Integer vector of row indices (observations) where at
#'    least one feature exceeded its cap.}
#'  \item{\code{excess}}{Numeric matrix of the same dimension as \code{z}
#'    (or \code{object$feature_contributions}) containing the amount by which
#'    each absolute contribution exceeded its cap. Zero entries indicate no 
#'    violation.}
#' 
#' If no violations are detected (the sum of \code{excess} is zero), the function
#' prints a message and returns \code{NULL} invisibly.
#' 
#' @details
#' Contributions are computed as
#' \eqn{|X_{ij} \beta_j|\times multiplier_i},
#' compared against the cap \eqn{L_j} stored in the fitted object. The returned
#' \code{excess} matrix gives
#' \eqn{\max(|X_{ij} \beta_j| \times multiplier_i - L_j, 0)}.
#' 
#' @seealso [capnet()]
#' 
#' @examples
#' set.seed(1)
#' n <- 40; p <- 8
#' X <- matrix(rnorm(n * p), n, p)
#' beta <- c(2, 1.5, rep(0, p-2))
#' y <- X %*% beta + rnorm(n)
#' fit <- capnet(X, y, lambda = 0.1, alpha = 0.5, gamma = 1, L = 1)
#' capnet_violations(fit)
#' 
#' @export

capnet_violations <- function(object, z = NULL, multiplier = NULL) {
  if (!is.null(z)) {
    m <- nrow(z)
    if (is.null(multiplier)) {
      if (length(unique(object$multiplier)) != 1) {
        warning("Did not supply a multiplier. Original model has varying multiplier values. Defaulting to the first observed multiplier.")
      }
      multiplier <- rep(unique(object$multiplier)[1], m)
    }
    beta <- object$beta
    contribution <- sweep(sweep(z, 2, beta, "*"), 1, multiplier, "*")
  } else {
    contribution <- object$feature_contributions
  }
  
  excess_contribution <- pmax(sweep(abs(contribution), 2, object$L, "-"), 0)
  if (sum(excess_contribution) == 0){
    message("No cap violations detected.")
    return(invisible(NULL))
  }
  
  structure(list(
    columns = which(colSums(excess_contribution) > 0),
    rows = which(rowSums(excess_contribution) > 0),
    excess = excess_contribution
  ), class = "capnet_violations")
}