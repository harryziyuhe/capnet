.loss_gradient_check <- function(X, y, family, betas = NULL, eps = 1e-6) {
  family <- normalize_family(family)
  
  ll <- make_likelihood(family, X, y)
  n <- nrow(X)
  p <- ncol(X)
  
  if (is.null(betas)) {
    betas <- rnorm(p + 1)
  }
  beta_0 <- betas[1]
  beta <- betas[-1]
  
  loss <- ll$loss
  grad <- ll$grad
  
  b0_loss <- (loss(beta_0 + eps, beta) - loss(beta_0 - eps, beta)) / (2 * eps)
  b0_grad <- grad(beta_0, beta)$g0
  b0_diff <- abs(b0_loss - b0_grad)
  if (b0_diff > eps) {
    warning("Error detected in calculating intercept gradient.")
  }
  
  for (i in seq_len(p)) {
    beta_up <- beta_down <- beta
    beta_up[i] <- beta_up[i] + eps
    beta_down[i] <- beta_down[i] - eps
    beta_loss <- (loss(beta_0, beta_up) - loss(beta_0, beta_down)) / (2 * eps)
    beta_grad <- grad(beta_0, beta)$gb[i]
    beta_diff <- abs(beta_loss - beta_grad)
    if (beta_diff > eps) {
      warning("Error detected in calculating coefficient %d", i)
    }
  }
}
