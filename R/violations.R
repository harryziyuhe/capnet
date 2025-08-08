#' Check contribution cap violations for a capnet model.
#' 
#' Check contribution cap violations for a fitted capnet model. Returns the column
#' and row indices as well as the excess contribution matrix.
#' 
#' @param object A fitted capnet model
#' @param newx New data. If unspecified use `newx` originally passed to the fitted
#' capnet
#' @param multiplier override multiplier. If unspecified use `multiplier` 
#' originally passed to the fitted capnet
#' 
#' @return A list with the following components:
#' \item {columns}{The column index containing violations of contribution cap.}
#' \item {rows}{The row index containing violations of contribution cap.}
#' \item {excess}{The matrix containing violations of contribution cap.}

check_cap_violations <- function(object, newx = NULL, multiplier = NULL) {
  if (!is.null(newx)) {
    if (is.null(multiplier)) {
      multiplier <- object$multiplier
    }
    beta <- object$beta
    contribution <- sweep(sweep(newx, 2, beta, "*"), 2, multiplier, "*")
  } else {
    contribution <- object$feature_contributions
  }
  
  excess_contribution <- pmax(sweep(abs(contribution), 2, L, "-"), 0)
  if (sum(excess_contribution) == 0){
    message("No cap violations detected.")
    return(invisible(NULL))
  }
  
  list(
    columns = which(colSums(excess_contribution) > 0),
    rows = which(rowSums(excess_contribution) > 0),
    excess = excess_contribution
  )
}