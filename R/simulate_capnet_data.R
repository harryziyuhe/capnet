#' Generate synthetic data for testing \code{capnet}
#'
#' Creates synthetic regression data for two scenarios with a unified return
#' format. In both cases the function returns \code{X}, \code{y}, the ground-truth
#' coefficients \code{beta_true}, and suggested caps \code{L}. When \code{m > 0},
#' an evaluation matrix \code{Z} of \code{m} rows is also returned, generated from
#' the same data-generating process and ready to pass as the \code{z} argument of
#' \code{capnet()} or \code{walk_capnet()}.
#'
#' @importFrom stats rnorm runif rbinom cor rpois rgamma qlogis plogis sd
#'
#' @param n Integer; number of training observations. Default \code{100}.
#' @param p Integer; number of predictors. Default \code{20}.
#' @param m Integer; number of evaluation rows to generate. When \code{m > 0}
#'   the function returns a \code{Z} matrix of \code{m} rows drawn from the same
#'   DGP as \code{X}. Default \code{0} (no evaluation set returned).
#' @param mode Character; one of \code{"basic"} or \code{"proxy"}. Default \code{"basic"}.
#' @param sparsity Numeric in \eqn{[0,1]}; fraction of nonzero coefficients
#'   in \code{beta_true}. Default \code{0.7}
#' @param beta_range Numeric length-2; uniform range for nonzero \code{beta_true}
#'   entries. Default \code{c(-2, 2)}.
#' @param sigma Positive numeric; noise sd for \code{y}. Only used when
#'   \code{family = "gaussian"}. Default \code{1}.
#' @param seed Optional integer for reproducibility. Default \code{NULL}.
#'
#' @section Proxy/unstable-feature controls (used when \code{mode = "proxy"}):
#' These arguments control the generation of unstable proxy features when
#' simulating data under \code{mode = "proxy"}.
#'
#' @param xstar_sd Positive numeric; sd of latent true signal \eqn{X^*}. Default \code{1}.
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
#' @param family Character; response family, one of \code{"gaussian"},
#'   \code{"binomial"}, \code{"poisson"}, or \code{"gamma"}. Default
#'   \code{"gaussian"}. Only supported in \code{mode = "basic"};
#'   \code{mode = "proxy"} always generates a Gaussian response, since that
#'   mode's measurement-error structure is orthogonal to response family. For
#'   every non-Gaussian family, the slope part of the linear predictor,
#'   \code{X \%*\% beta_true}, is rescaled to unit standard deviation and
#'   \code{beta_true} is updated to match, then a data-generating intercept
#'   (not returned, analogous to the \code{a0} that \code{capnet()} fits
#'   internally) shifts the linear predictor so the response's mean matches
#'   \code{poisson_mean}/\code{gamma_mean}/\code{binomial_prob} as applicable.
#'   This keeps the response in a numerically sane range regardless of
#'   \code{p} or \code{beta_range}.
#' @param poisson_mean Positive numeric; target average response when
#'   \code{family = "poisson"}. Default \code{10}.
#' @param gamma_mean Positive numeric; target average response when
#'   \code{family = "gamma"}. Default \code{10}.
#' @param gamma_shape Positive numeric; fixed Gamma shape parameter used when
#'   \code{family = "gamma"} (rate is set to \code{gamma_shape / mu} so the
#'   mean is \code{mu}). Default \code{2}.
#' @param binomial_prob Numeric in \eqn{(0,1)}; target average success
#'   probability when \code{family = "binomial"}. Default \code{0.5}.
#'
#' @return A list with the following components:
#'   \item{\code{X}}{\eqn{(n \times p)} training design matrix.}
#'   \item{\code{y}}{Numeric response vector of length \eqn{n}.}
#'   \item{\code{beta_true}}{Numeric length-\eqn{p} vector of ground-truth
#'    coefficients.}
#'   \item{\code{L}}{Length-\eqn{p} vector of suggested contribution caps
#'    (nonnegative).}
#'   \item{\code{Z}}{(\eqn{m \times p}) evaluation matrix, returned only when
#'    \code{m > 0}. Generated from the same DGP as \code{X} and ready to pass
#'    as \code{z} in \code{capnet()}.}
#'  For \code{mode = "proxy"}, additional fields are included:
#'   \item{\code{X_latent}}{\eqn{(n \times p)} latent true-signal matrix
#'    \eqn{X^*} for the training rows.}
#'   \item{\code{a}}{Feature scales.}
#'   \item{\code{unstable_idx}}{Integer vector of unstable feature indices.}
#'   \item{\code{stable_idx}}{Integer vector of stable feature indices.}
#'
#' @details
#' In \strong{basic} mode, nonzeros in \code{beta_true} are placed uniformly at
#' random and drawn from \code{runif(beta_range[1], beta_range[2])}. When
#' \code{family = "gaussian"}, \code{y = X \%*\% beta_true + N(0, sigma^2)}.
#' For the other families (see the rescaling description under the
#' \code{family} argument above), \code{y ~ Poisson(exp(eta))} for
#' \code{"poisson"}, \code{y ~ Gamma(shape = gamma_shape, rate = gamma_shape /
#' exp(eta))} for \code{"gamma"}, and \code{y ~ Bernoulli(plogis(eta))} for
#' \code{"binomial"}. In all cases, the suggested caps are
#' \code{L_j = pmax(|beta_true_j|, 0.1)} as a simple starting point. When
#' \code{m > 0}, \code{Z} is an \eqn{m \times p} matrix of i.i.d. standard
#' normal draws from the same pool as \code{X}.
#'
#' In \strong{proxy} mode, \code{y = X^* \%*\% beta_true + N(0, sigma^2)} and
#' \code{X = diag(a) X^* + varepsilon}. For columns in \code{unstable_idx}, each
#' row independently uses \code{sigma_high} with probability \code{high_prob}.
#' When \code{m > 0}, \code{Z} is generated by applying the same noisy proxy
#' process to \code{m} fresh latent rows, giving the observed evaluation matrix
#' a practitioner would actually have available at evaluation time.
#'
#' @examples
#' set.seed(42)
#' d_basic <- simulate_capnet_data(n = 100, p = 10, m = 20, sparsity = 0.5)
#' names(d_basic)  # includes Z
#'
#' d_pois <- simulate_capnet_data(n = 100, p = 10, m = 20, sparsity = 0.5,
#'                                 family = "poisson", poisson_mean = 15)
#' range(d_pois$y)
#'
#' d_binom <- simulate_capnet_data(n = 100, p = 10, m = 20, sparsity = 0.5,
#'                                  family = "binomial", binomial_prob = 0.3)
#' table(d_binom$y)
#'
#' d_gamma <- simulate_capnet_data(n = 100, p = 10, m = 20, sparsity = 0.5,
#'                                  family = "gamma", gamma_mean = 5)
#' range(d_gamma$y)
#'
#' d_proxy <- simulate_capnet_data(mode = "proxy", n = 300, p = 8, m = 50,
#'                                 sparsity = 0.6, unstable_frac = 0.25,
#'                                 sigma_low = 1, sigma_high = 6, high_prob = 0.2)
#' names(d_proxy)  # includes Z and X_latent
#'
#' @export
simulate_capnet_data <- function(
    n = 100, p = 20, m = 0,
    mode = c("basic", "proxy"),
    sparsity = 0.7,
    beta_range = c(-2, 2),
    sigma = 1,
    seed = NULL,
    xstar_sd = 1,
    a_sd = 0.3,
    unstable_idx = NULL,
    unstable_frac = 0.2,
    sigma_low = 1.0,
    sigma_high = 8.0,
    high_prob = 0.1,
    family = c("gaussian", "binomial", "poisson", "gamma"),
    poisson_mean = 10,
    gamma_mean = 10,
    gamma_shape = 2,
    binomial_prob = 0.5
) {
  mode <- match.arg(mode)
  family <- match.arg(family)
  if (family != "gaussian" && mode == "proxy") {
    stop(sprintf("family = \"%s\" is only supported in mode = \"basic\".", family))
  }
  if (!is.null(seed)) set.seed(seed)

  if (mode == "basic") {
    total <- n + m
    X_all <- matrix(rnorm(total * p), total, p)
    colnames(X_all) <- paste0("X", seq_len(p))

    X <- X_all[seq_len(n), , drop = FALSE]

    beta_true <- numeric(p)
    k <- floor(sparsity * p)
    if (k > 0) {
      nz <- sort(sample(seq_len(p), size = k, replace = FALSE))
      beta_true[nz] <- runif(k, min = beta_range[1], max = beta_range[2])
    }

    if (family == "gaussian") {
      y <- drop(X %*% beta_true + rnorm(n, sd = sigma))
    } else {
      # Rescale so the slope part of the linear predictor has unit sd, then
      # shift (via a data-generating intercept, not returned as part of
      # beta_true) so the response's mean matches the family-specific target.
      target_eta0 <- switch(family,
        poisson  = log(poisson_mean),
        gamma    = log(gamma_mean),
        binomial = qlogis(binomial_prob)
      )
      eta_slope <- drop(X %*% beta_true)
      s <- sd(eta_slope)
      if (s > 0) {
        beta_true <- beta_true / s
        eta_slope <- eta_slope / s
      }
      eta <- eta_slope + (target_eta0 - mean(eta_slope))

      y <- switch(family,
        poisson  = rpois(n, exp(eta)),
        gamma    = rgamma(n, shape = gamma_shape, rate = gamma_shape / exp(eta)),
        binomial = rbinom(n, size = 1, prob = plogis(eta))
      )
    }
    L <- pmax(abs(beta_true), 0.1)

    out <- list(X = X, y = y, beta_true = beta_true, L = L)
    if (m > 0) {
      out$Z <- X_all[(n + 1):(n + m), , drop = FALSE]
    }
    return(out)
  }

  # proxy mode: generate n + m latent rows from the same X* distribution
  total <- n + m
  X_star_all <- matrix(rnorm(total * p, sd = xstar_sd), total, p)

  beta_true <- rnorm(p)
  k0 <- floor((1 - sparsity) * p)
  if (k0 > 0) beta_true[sample(seq_len(p), k0)] <- 0

  X_latent <- X_star_all[seq_len(n), , drop = FALSE]
  y <- drop(X_latent %*% beta_true + rnorm(n, sd = sigma))

  a <- exp(rnorm(p, mean = 0, sd = a_sd))

  if (is.null(unstable_idx)) {
    k_unstable <- max(1, round(unstable_frac * p))
    unstable_idx <- sort(sample(seq_len(p), k_unstable, replace = FALSE))
  }
  stable_idx <- setdiff(seq_len(p), unstable_idx)

  .make_proxy <- function(X_star_rows) {
    nr <- nrow(X_star_rows)
    mean_part <- sweep(X_star_rows, 2, a, "*")
    Eta <- matrix(rnorm(nr * p, sd = sigma_low), nr, p)
    if (length(unstable_idx) > 0) {
      flips <- matrix(rbinom(nr * length(unstable_idx), 1, high_prob),
                      nr, length(unstable_idx))
      sd_mat <- sigma_low + flips * (sigma_high - sigma_low)
      Eta[, unstable_idx] <- matrix(rnorm(nr * length(unstable_idx)),
                                    nr, length(unstable_idx)) * sd_mat
    }
    mean_part + Eta
  }

  X <- .make_proxy(X_latent)
  colnames(X)        <- paste0("X", seq_len(p))
  colnames(X_latent) <- paste0("X", seq_len(p))

  L <- pmax(abs(cor(X, y) * 5), 0.1)

  out <- list(
    X = X, y = y, beta_true = beta_true,
    L = L, X_latent = X_latent,
    a = a, unstable_idx = unstable_idx, stable_idx = stable_idx
  )

  if (m > 0) {
    X_star_eval <- X_star_all[(n + 1):(n + m), , drop = FALSE]
    Z <- .make_proxy(X_star_eval)
    colnames(Z) <- paste0("X", seq_len(p))
    out$Z <- Z
  }

  out
}
