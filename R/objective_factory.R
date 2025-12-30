make_likelihood <- function(family, X, y) {
  
  n <- length(y)
  
  if (family$family == "gaussian") {
    loss <- function(beta0, beta) {
      eta <- as.vector(beta0 + X %*% beta)
      sum((y - eta)^2) / (2 * n)
    }
    gradient <- function(beta0, beta) {
      eta <- as.vector(beta0 + X %*% beta)
      r <- (eta - y)                       # note sign for convenience
      g0 <- sum(r) / n
      gb <- as.vector(crossprod(X, r) / n)
      list(g0 = g0, gb = gb)
    }
    return(list(loss = loss, gradient = gradient, family = family))
  }
  
  if (family$family == "binomial") {
    loss <- function(beta0, beta) {
      eta <- as.vector(beta0 + X %*% beta)
      mean(log1pexp(eta) - y * eta)
    }
    gradient <- function(beta0, beta) {
      eta <- as.vector(beta0 + X %*% beta)
      mu <- sigmoid(eta)
      r <- (mu - y)
      g0 <- mean(r)
      gb <- as.vector(crossprod(X, r) / n)
      list(g0 = g0, gb = gb)
    }
    return(list(loss = loss, gradient = gradient, family = family))
  }
  
  if (family$family == "poisson") {
    loss <- function(beta0, beta) {
      eta <- as.vector(beta0 + X %*% beta)
      mu <- safe_exp(eta)
      mean(mu - y * eta)
    }
    gradient <- function(beta0, beta) {
      eta <- as.vector(beta0 + X %*% beta)
      mu <- safe_exp(eta)
      r <- (mu - y)
      g0 <- mean(r)
      gb <- as.vector(crossprod(X, r) / n)
      list(g0 = g0, gb = gb)
    }
    return(list(loss = loss, gradient = gradient, family = family))
  }
  
  if (family$family == "Gamma") {
    # log link only, constants dropped
    loss <- function(beta0, beta) {
      eta <- as.vector(beta0 + X %*% beta)
      mu <- safe_exp(eta)
      mean(y / mu + log(mu))              # = mean(y*exp(-eta) + eta)
    }
    gradient <- function(beta0, beta) {
      eta <- as.vector(beta0 + X %*% beta)
      mu <- safe_exp(eta)
      r <- (1 - y / mu)                   # d/deta
      g0 <- mean(r)
      gb <- as.vector(crossprod(X, r) / n)
      list(g0 = g0, gb = gb)
    }
    return(list(loss = loss, gradient = gradient, family = family))
  }
  
  stop("Unsupported family.", call. = FALSE)
}
