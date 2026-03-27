.capnet_spec <- function(
  X, y, L,
  family = "gaussian",
  intercept = TRUE,
  standardize = TRUE,
  newx = NULL,
  multiplier = 1,
  lower.limits = NULL,
  upper.limits = NULL,
  tol = 1e-8,
  maxit = 10000L,
  par = NULL,
    ...
) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  
  n <- nrow(X)
  p <- ncol(X)
  if (length(y) != n) stop("y must have length nrow(X).")
  
  family <- normalize_family(family)
  validate_family_supported(family)
  
  newx <- newx %||% X
  newx <- as.matrix(newx)
  if (ncol(newx) != p) stop("newx must have same number of columns as X.")
  
  m <- nrow(newx)
  
  multiplier <- as.numeric(multiplier)
  if (length(multiplier) == 1) {
    multiplier <- rep(multiplier, m)
  } else if (length(multiplier) != m) {
    stop("multiplier must be length 1 or length nrow(newx).")
  }
  
  if (is.null(L)) stop("L must be provided (scalar or length p).")
  L <- as.numeric(L)
  if (length(L) == 1) {
    L <- rep(L, p)
  } else if (length(L) != p) {
    stop("L must be length 1 or length ncol(X).")
  }
  
  # Set default lower and upper bounds
  if (is.null(lower.limits)) {
    lower.limits <- rep(-Inf, p)
  } else if (length(lower.limits) == 1) {
    lower.limits <- rep(lower.limits, p)
  } else if (length(lower.limits) != p) {
    stop("Lower limits dimension mismatch.")
  }
  
  if (is.null(upper.limits)) {
    upper.limits <- rep(Inf, p)
  } else if (length(upper.limits) == 1) {
    upper.limits <- rep(upper.limits, p)
  } else if (length(upper.limits) != p) {
    stop("Lower limits dimension mismatch.")
  }
  
  if (!is.null(par)) par <- as.numeric(par)
  
  list(
    X = X,
    y = y,
    n = n,
    p = p,
    family = family,
    intercept = intercept,
    standardize = standardize,
    newx = newx,
    multiplier = multiplier,
    L = L,
    lower.limits = lower.limits,
    upper.limits = upper.limits,
    tol = tol,
    maxit = as.integer(maxit),
    par = par
  )
}

.capnet_standardize_train <- function(spec, idx_train = NULL) {
  if (is.null(idx_train)) idx_train <- seq_len(spec$n)
  idx_train <- as.integer(idx_train)
  
  X_train_raw <- spec$X[idx_train, , drop = FALSE]
  y_train <- spec$y[idx_train]
  
  if (!spec$standardize) {
    X_center <- rep(0, spec$p)
    X_scale <- rep(1, spec$p)
    X_train <- X_train_raw
    lower.limits <- spec$lower.limits
    upper.limits <- spec$upper.limits
  } else {
    X_center <- colMeans(X_train_raw)
    X_scale <- apply(X_train_raw, 2, sd)
    X_scale[X_scale == 0] <- 1
    X_train <- sweep(X_train_raw, 2, X_center, "-")
    X_train <- sweep(X_train, 2, X_scale, "/")
    
    lower.limits <- spec$lower.limits * X_scale
    upper.limits <- spec$upper.limits * X_scale
  }
  scaling_factor <- X_scale
  
  par0 <- spec$par
  if (is.null(par0)) {
    par0 <- rep(0, 1 + spec$p)
  } else {
    if (length(par0) != (1 + spec$p)) {
      stop("par must be length 1 + ncol(X).")
    }
  }
  
  list(
    X = X_train,
    y = y_train,
    X_center = X_center,
    X_scale = X_scale,
    scaling_factor = scaling_factor,
    lower.limits = lower.limits,
    upper.limits = upper.limits,
    family = spec$family,
    intercept = spec$intercept,
    standardize = spec$standardize,
    tol = spec$tol,
    maxit = spec$maxit,
    par = par0
  )
}

.capnet_cap_context <- function(spec, idx_cap = NULL) {
  newx <- spec$newx
  multiplier <- spec$multiplier
  
  if (!is.null(idx_cap)) {
    idx_cap <- as.integer(idx_cap)
    newx <- newx[idx_cap, , drop = FALSE]
    multiplier <- multiplier[idx_cap]
  }
  
  if (ncol(newx) != spec$p) stop("cap newx must have ncol = ncol(X).")
  if (length(multiplier) != nrow(newx)) stop("multiplier must align with cap newx rows.")
  
  list(
    newx = newx,
    multiplier = multiplier,
    L = spec$L
  )
}


# Preparing input for capnet fit
.capnet_prepare <- function(X, y, lambda, alpha, mu, L, newx = NULL,
                            par = NULL, multiplier = 1, 
                            family = NULL, intercept = TRUE, standardize = TRUE,
                            lower.limits = NULL, upper.limits = NULL,
                            tol = 1e-8, maxit = 10000) {
  X_raw <- X
  y_raw <- y
  
  if (is.null(newx)) {
    newx <- X_raw
  } else {
    newx <- as.matrix(newx)
  }

  n <- nrow(X)
  p <- ncol(X)
  m <- nrow(newx)
  
  if (is.null(par)) {
    par <- rep(0, p + 1)
  }
  
  if (is.data.frame(L) || is.matrix(L)) L <- as.numeric(unlist(L))
  if (is.list(L)) L <- unlist(L)
  if (length(L) == 1) L <- rep(L, p)
  if (length(L) != p) stop("L must be length 1 or match the number of covariates (p = ", p, ")")
  
  multiplier = as.numeric(multiplier)
  if (length(multiplier) == 1) multiplier <- rep(multiplier, m)
  if (length(multiplier) != m) stop("multiplier must be length 1 or match the number of observations")
  
  # Validate model family
  family <- normalize_family(family)
  validate_family_supported(family)
  
  # Standardize data
  if (standardize) {
    y <- scale(y_raw)
    y_center <- attr(y, "scaled:center")
    y_scale <- attr(y, "scaled:scale")
    
    X <- scale(X_raw)
    X_center <- attr(X, "scaled:center")
    X_scale <- attr(X, "scaled:scale")
    
    scaling_factor <- X_scale / y_scale
  } else {
    scaling_factor <- rep(1, p)
    y_scale <- 1
    y_center <- 0
    X_center <- rep(0, p)
  }
  
  # Set default lower and upper bounds
  if (is.null(lower.limits)) {
    lower.limits <- rep(-Inf, p)
  } else if (length(lower.limits) == 1) {
    lower.limits <- rep(lower.limits, p)
  } else if (length(lower.limits) != p) {
    stop("Lower limits dimension mismatch.")
  }
  
  if (is.null(upper.limits)) {
    upper.limits <- rep(Inf, p)
  } else if (length(upper.limits) == 1) {
    upper.limits <- rep(upper.limits, p)
  } else if (length(upper.limits) != p) {
    stop("Lower limits dimension mismatch.")
  }
  
  list(
    X = X,
    y = y,
    scaling_factor = scaling_factor,
    y_scale = y_scale,
    X_center = X_center,
    y_center = y_center,
    lambda = lambda,
    alpha = alpha,
    mu = mu,
    L = L,
    newx = newx,
    par = par,
    multiplier = multiplier,
    family = family,
    intercept = intercept,
    standardize = standardize,
    lower.limits = lower.limits,
    upper.limits = upper.limits,
    tol = tol,
    maxit = maxit
  )
}
