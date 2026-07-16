.capnet_spec <- function(
  X, y, L,
  family = "gaussian",
  intercept = TRUE,
  standardize = TRUE,
  z = NULL,
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
  
  z <- z %||% X
  z <- as.matrix(z)
  if (ncol(z) != p) stop("z must have same number of columns as X.")

  m <- nrow(z)
  
  multiplier <- as.numeric(multiplier)
  if (length(multiplier) == 1) {
    multiplier <- rep(multiplier, m)
  } else if (length(multiplier) != m) {
    stop("multiplier must be length 1 or length nrow(z).")
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
    z = z,
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
  z <- spec$z
  multiplier <- spec$multiplier

  if (!is.null(idx_cap)) {
    idx_cap <- as.integer(idx_cap)
    z <- z[idx_cap, , drop = FALSE]
    multiplier <- multiplier[idx_cap]
  }

  if (ncol(z) != spec$p) stop("z must have ncol = ncol(X).")
  if (length(multiplier) != nrow(z)) stop("multiplier must align with z rows.")

  list(
    z = z,
    multiplier = multiplier,
    L = spec$L
  )
}


# Preparing input for capnet fit
.capnet_prepare <- function(X, y, lambda, alpha, gamma, L, z = NULL,
                            par = NULL, multiplier = 1, 
                            family = NULL, intercept = TRUE, standardize = TRUE,
                            lower.limits = NULL, upper.limits = NULL,
                            tol = 1e-8, maxit = 10000) {
  X_raw <- X
  y_raw <- y
  
  if (is.null(z)) {
    z <- X_raw
  } else {
    z <- as.matrix(z)
  }

  n <- nrow(X)
  p <- ncol(X)
  m <- nrow(z)
  
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
    gamma = gamma,
    L = L,
    z = z,
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
