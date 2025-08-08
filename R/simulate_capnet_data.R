#' Generate synthetic data for testing capnet
#'
#' Creates a synthetic regression dataset with some true coefficients,
#' additive noise, and suggested contribution caps.
#'
#' @param n Number of observations (default = 100)
#' @param p Number of predictors (default = 20)
#' @param sparsity Proportion of non-zero coefficients (default = 0.7)
#' @param sigma Noise level (default = 1)
#' @param seed Random seed for reproducibility (default = NULL)
#' @return A list containing:
#'   \item{X}{Design matrix (n x p)}
#'   \item{y}{Response vector}
#'   \item{beta_true}{True beta coefficients}
#'   \item{L}{Suggested contribution cap vector}
#' @export

simulate_capnet_data <- function(n = 100, p = 20,
                                 sparsity = 0.7,
                                 sigma = 1,
                                 seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("X", 1:p)
  
  beta_true <- numeric(p)
  nonzero_idx <- sample(1:p, size = floor(sparsity * p))
  beta_true[nonzero_idx] <- runif(length(nonzero_idx), -2, 2)
  
  y <- X %*% beta_true + rnorm(n, sd = sigma)
  y <- drop(y)
  
  L <- pmax(abs(1 * beta_true), 0.1)
  
  list(X = X, y = y, beta_true = beta_true, L = L)
}
