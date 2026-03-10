#' Generate synthetic data for testing \code{capnet}
#'
#' Creates synthetic regression data for two scenarios with a unified return
#' format. In both cases the function returns \code{X}, \code{y}, the ground-truth
#' coefficients \code{beta_true}, suggested caps \code{L}, and a matrix \code{Z}
#' of latent drivers. In \code{mode = "basic"}, there are no latent drivers, so
#' \code{Z} is set to \code{X}.
#' 
#' @importFrom stats rnorm runif rbinom
#'
#' @param n Integer; number of observations. Default \code{100}.
#' @param p Integer; number of predictors. Default \code{20}.
#' @param mode Character; one of \code{"basic"} or \code{"proxy"}. Default \code{"basic"}.
#' @param sparsity Numeric in \eqn{[0,1]}; fraction of nonzero coefficients
#'   in \code{beta_true}. Default \code{0.7}
#' @param beta_range Numeric length-2; uniform range for nonzero \code{beta_true}
#'   entries. Default \code{c(-2, 2)}.
#' @param sigma Positive numeric; noise sd for \code{y}. Default \code{1}.
#' @param seed Optional integer for reproducibility. Default \code{NULL}.
#'
#' @section Proxy/unstable-feature controls (used when \code{mode = "proxy"}):
#' These arguments control the generation of unstable proxy features when
#' simulating data under \code{mode = "proxy"}.
#' 
#' @param z_sd Positive numeric; sd of latent drivers \code{Z}. Default \code{1}.
#' @param a_sd Nonnegative numeric; sd of log-scales \code{a_j = exp(N(0, a_sd^2))}.
#'   Default \code{0.3}.
#' @param unstable_idx Optional integer vector of columns of \code{X} with
#'   unstable variance. If \code{NULL}, \code{round(unstable_frac * p)} columns
#'   are sampled at random.
#' @param unstable_frac Numeric in \eqn{[0,1]}; fraction of features to label
#'   unstable when \code{unstable_idx} is \code{NULL}. Default \code{0.25}.
#' @param sigma_low Positive numeric; base noise sd for \code{X}. Default \code{1}.
#' @param sigma_high Positive numeric; intermittent high noise sd for unstable
#'   features. Default \code{8}.
#' @param high_prob Numeric in \eqn{[0,1]}; per-row probability an unstable
#'   feature uses \code{sigma_high}. Default \code{0.10}.
#'
#' @return A list with the same core components for both modes:
#'   \item{\code{X}}{\eqn{(n \times p)} observed design matrix (colnames 
#'    auto-filled as \code{"X1"},… if missing).}
#'   \item{\code{y}}{Numeric response vector of length \eqn{n}.}
#'   \item{\code{beta_true}}{Numeric length-\eqn{p} vector of ground-truth 
#'    coefficients.}
#'   \item{\code{L}}{Length-\eqn{p} vector of suggested contribution caps 
#'    (nonnegative).}
#'   \item{\code{Z}}{\eqn{(n \times p)} latent drivers. In 
#'    \code{mode = "basic"}, \code{Z = X}.}
#' For \code{mode = "proxy"}, additional convenience fields are included:
#' \code{a} (feature scales), \code{unstable_idx}, \code{stable_idx}.
#'
#' @details
#' In \strong{basic} mode, nonzeros in \code{beta_true} are placed uniformly at
#' random and drawn from \code{runif(beta_range[1], beta_range[2])}. The suggested
#' caps are \code{L_j = pmax(|beta_true_j|, 0.1)} as a simple starting point.
#'
#' In \strong{proxy} mode, \code{y = Z \%*\% beta_true + N(0, sigma^2)} and
#' \code{X = diag(a) Z + varepsilon}. For columns in \code{unstable_idx}, each
#' row independently flips to \code{sigma_high} with probability \code{high_prob}.
#'
#' @examples
#' set.seed(42)
#' d_basic <- simulate_capnet_data(n = 100, p = 10, sparsity = 0.5)
#' names(d_basic)
#'
#' d_proxy <- simulate_capnet_data(mode = "proxy", n = 300, p = 8,
#'                                 sparsity = 0.6, unstable_frac = 0.25,
#'                                 sigma_low = 1, sigma_high = 6, high_prob = 0.2)
#' names(d_proxy)
#'
#' @export
simulate_capnet_data <- function(
    n = 100, p = 20,
    mode = c("basic", "proxy"),
    sparsity = 0.7,
    beta_range = c(-2, 2),
    sigma = 1,
    seed = NULL,
    z_sd = 1,
    a_sd = 0.3,
    unstable_idx = NULL,
    unstable_frac = 0.2,
    sigma_low = 1.0,
    sigma_high = 8.0,
    high_prob = 0.1
) {
  mode <- match.arg(mode)
  if (!is.null(seed)) set.seed(seed)
  
  if (mode == "basic") {
    X <- matrix(rnorm(n * p), n, p)
    if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(p))
    
    beta_true <- numeric(p)
    k <- floor(sparsity * p)
    if (k > 0) {
      nz <- sort(sample(seq_len(p), size = k, replace = FALSE))
      beta_true[nz] <- runif(k, min = beta_range[1], max = beta_range[2])
    }
    
    y <- drop(X %*% beta_true + rnorm(n, sd = sigma))
    L <- pmax(abs(beta_true), 0.1)
    
    return(list(X = X, y = y, beta_true = beta_true,
                L = L, Z = X))
  }
  
  Z <- matrix(rnorm(n * p, sd = z_sd), n, p)
  beta_true <- rnorm(p)
  k0 <- floor((1 - sparsity) * p)
  if (k0 > 0) beta_true[sample(seq_len(p), k0)] <- 0
  y <- drop(Z %*% beta_true + rnorm(n, sd = sigma))
  
  a <- exp(rnorm(p, mean = 0, sd = a_sd))
  mean_part <- sweep(Z, 2, a, "*")
  
  if (is.null(unstable_idx)) {
    k_unstable <- max(1, round(unstable_frac * p))
    unstable_idx <- sort(sample(seq_len(p), k_unstable, replace = FALSE))
  }
  stable_idx <- setdiff(seq_len(p), unstable_idx)
  
  Eta <- matrix(rnorm(n * p, sd = sigma_low), n, p)
  if (length(unstable_idx) > 0) {
    flips <- matrix(rbinom(n * length(unstable_idx), 1, high_prob), n,
                    length(unstable_idx))
    sd_mat <- sigma_low + flips * (sigma_high - sigma_low)
    Eta[, unstable_idx] <- matrix(rnorm(n * length(unstable_idx)), n, 
                                  length(unstable_idx)) * sd_mat
  }
  
  X <- mean_part + Eta
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(p))
  
  L <- pmax(abs(cor(X,y) * 5), 0.1)
  
  list(
    X = X, y = y,beta_true = beta_true,
    L = L, Z = Z
  )
}
