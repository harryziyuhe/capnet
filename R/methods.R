#' @export
coef.capnet <- function(object, ...) {
  coefs <- rbind(object$a0, object$beta)
  rownames(coefs) <- c("(Intercept)", colnames(object$newx))
  colnames(coefs) <- "beta"
  coefs
}

#' @export
coef.capnet_walk <- function(object, index = NULL, ...) {
  coefs <- cbind(object$intercepts, object$betas)
  if (!is.null(index)) {
    coefs <- coefs[index,]
  }
  coefs
}

#' @export
predict.capnet <- function(object, newx  = NULL,...) {
  if (!is.null(newx)) {
    newx <- object$newx
  }
  predictions <- rowSums(sweep(newx, 2, object$beta, "*")) + object$a0
  colnames(predictions) <- "prediction"
  predictions
}

#' @export
predict.capnet_walk <- function(object, ...) {
  object$predictions
}

#' @export
plot.cv_capnet <- function(object, alpha = NULL, lambda = NULL, ...) {
  cv_errors <- object$cv_errors
  mean_errors <- object$mean_errors
  alphas <- object$alpha
  lambdas <- object$lambda
  
  dims <- dim(cv_errors)
  nx <- dims[1]
  ny <- dims[2]
  ns <- dims[3]
  
  melt_matrix <- function(mat) {
    expand.grid(alpha = alphas, lambda = log10(lambdas)) %>% 
      transform(value = as.vector(mat))
  }
  
  if (is.null(alpha) && is.null(lambda)) {
    df <- melt_matrix(mean_errors)
    df <- df %>%
      mutate(
        rank_value = rank(value, ties.method = "average"),
        rank_scaled = (rank_value - min(rank_value)) / 
          (max(rank_value) - min(rank_value))
      )
    
    p <- ggplot(df, aes(x = alpha, y = lambda, fill = rank_scaled)) +
      geom_tile() + 
      scale_fill_viridis_c() +
      labs(title = "Ranked Mean CV Errors",
           x = "α", y = "λ (log)", fill = "Rank") +
      theme_minimal()
    
    print(p)
    return(invisible(object))
  }
  
  if (!is.null(alpha)) {
    if (!(alpha %in% alphas)) {
      stop("alpha value not evaluated during cross validation")
    }
    idx <- which(alphas == alpha)
    errors <- cv_errors[idx, , ]
    errors <- as.matrix(errors)
    lambda.min <- log(lambdas[which(rowSums(errors) == min(rowSums(errors)))])
    
    df <- data.frame(
      lambda = rep(log(lambdas), times = ns),
      fold = rep(1:ns, each = nx),
      error = as.vector(errors)
    )
    
    p <- ggplot(df, aes(x = lambda, y = error)) +
      stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.25, color = "gray60") +
      stat_summary(fun = mean, geom = "point", size = 2, color = "red") +
      geom_vline(xintercept = lambda.min, linetype = 2) +
      labs(
        title = bquote("Cross-validation errors for " * alpha == .(alpha)),
        x = expression(log(lambda)),
        y = "CV error"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()
      )
    
    print(p)
    return(invisible(object))
  }
  
  if (!is.null(lambda)) {
    if (!(lambda %in% lambdas)) {
      stop("lambda value not evaluated during cross validation")
    }
    idx <- which(lambdas == lambda)
    errors <- cv_errors[, idx, ]
    errors <- as.matrix(errors)
    alpha.min <- alphas[which(rowSums(errors) == min(rowSums(errors)))]
    
    df <- data.frame(
      alpha= rep(alphas, times = ns),
      fold = rep(1:ns, each = ny),
      error = as.vector(errors)
    )
    
    p <- ggplot(df, aes(x = alpha, y = error)) +
      stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.25, color = "gray60") +
      stat_summary(fun = mean, geom = "point", size = 2, color = "red") +
      geom_vline(xintercept = alpha.min, linetype = 2) +
      labs(
        title = bquote("Cross-validation errors for " * lambda == .(lambda)),
        x = expression(alpha),
        y = "CV error"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()
      )
    
    print(p)
    return(invisible(object))
  }
}

#' @export
plot.capnet_path <- function(object) {
  param <- names(object)[1]
  betas <- object %>% 
    pivot_longer(!{{param}}, names_to = "feature", values_to = "coef") %>% 
    rename(param = {{param}})
  p <- ggplot(betas, aes(x = param, y = coef, color = feature)) +
    geom_line() +
    labs(
      x = param,
      y = "Coefficients"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
}





















