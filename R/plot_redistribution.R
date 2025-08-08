#' Plot Redistribution of Coefficients and Contributions
#'
#' Generates side-by-side or delta bar plots comparing the regression coefficients and 
#' feature contributions between an uncapped and capped \code{capnet} model.
#'
#' @import ggplot2
#'
#' @param uncapped An object returned by \code{capnet()} without contribution caps.
#' @param capped An object returned by \code{capnet()} with contribution caps applied.
#' @param X Optional matrix of features to use when computing contributions. 
#'   Defaults to \code{capped\$newx}.
#' @param multiplier Optional vector of scaling multipliers for contributions. If not 
#'   provided, uses \code{capped\$multiplier}.
#' @param top_n Integer. Number of top features (by appearance order) to include in the plots. 
#'   Use \code{Inf} to include all features.
#' @param plot One of \code{"raw"} or \code{"delta"}. If \code{"raw"}, shows side-by-side 
#'   bar charts of beta and contribution values. If \code{"delta"}, shows bar charts of 
#'   changes in beta and contributions.
#' @return If \code{plot = "delta"}, returns a combined plot using \pkg{patchwork} or \pkg{gridExtra}.
#'   If \code{plot = "raw"}, the function creates and prints two separate \code{ggplot2} plots.
#' @details Contributions are calculated as the elementwise product of \code{X}, \code{beta},
#'   and \code{multiplier}. If \code{standardize = TRUE} was used during training, \code{multiplier}
#'   typically adjusts for scaling.
#' @import ggplot2
#' @examples
#' \dontrun{
#' fit1 <- capnet(X, y, lambda = 0, mu = 0, L = rep(0.1, ncol(X)))  # uncapped
#' fit2 <- capnet(X, y, lambda = 1, mu = 0.5, L = rep(0.1, ncol(X))) # capped
#' plot_redistribution(fit1, fit2, plot = "delta", top_n = 20)
#' }
#' @export

plot_redistribution <- function(uncapped, capped, X = NULL,
                                multiplier = NULL, top_n = Inf,
                                plot = c("raw", "delta")) {
  plot <- match.arg(plot)
  
  beta_uncapped <- uncapped$beta
  beta_capped <- capped$beta
  
  X <- if (!is.null(X)) X else capped$newx
  if (!is.null(multiplier)) {
    if (length(multiplier) == 1) multiplier <- rep(multiplier, nrow(X))
  } else multiplier <- capped$multiplier
  
  contribution_uncapped <- sweep(sweep(X, 2, beta_uncapped, "*"), 1, multiplier, "*")
  contribution_capped <- sweep(sweep(X, 2, beta_capped, "*"), 1, multiplier, "*")
  
  mean_contribution_uncapped <- colMeans(contribution_uncapped)
  mean_contribution_capped <- colMeans(contribution_capped)
  
  df <- data.frame(
    variable = colnames(X),
    beta_uncapped = beta_uncapped,
    beta_capped = beta_capped,
    contribution_uncapped = mean_contribution_uncapped,
    contribution_capped = mean_contribution_capped
  )
  
  df$delta_beta <- df$beta_capped - df$beta_uncapped
  df$delta_contribution <- df$contribution_capped - df$contribution_uncapped
  
  df_top <- head(df, top_n)
  
  if (plot == "delta") {
    p1 <- ggplot(df_top, aes(x = reorder(variable, delta_beta), y = delta_beta)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      theme_minimal() + 
      labs(
        title = "Change in Coefficients",
        x = NULL, y = NULL
      )
    
    p2 <- ggplot(df_top, aes(x = reorder(variable, delta_contribution), y = delta_contribution)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      theme_minimal() + 
      labs(
        title = "Change in Contributions",
        x = NULL, y = NULL
      )
  } else {
    df_beta <- data.frame(
      variable = rep(df_top$variable, times = 2),
      model = rep(c("Uncapped", "Capped"), each = nrow(df_top)),
      beta = c(df_top$beta_uncapped, df_top$beta_capped)
    )
    p1 <- ggplot(df_beta, aes(x = variable, y = beta, fill = model)) +
      geom_col(position = position_dodge(width = 0.6), width = 0.6) +
      labs(title = "Beta Comparison", y = "Beta", x = NULL, fill = "Model") +
      theme_minimal()
    
    df_contrib <- data.frame(
      variable = rep(df_top$variable, times = 2),
      model = rep(c("Uncapped", "Capped"), each = nrow(df_top)),
      contribution = c(df_top$contribution_uncapped, df_top$contribution_capped)
    )
    p2 <- ggplot(df_contrib, aes(x = variable, y = contribution, fill = model)) +
      geom_col(position = position_dodge(width = 0.6), width = 0.6) +
      labs(title = "Contribution Comparison", y = "Contribution", x = NULL, fill = "Model") +
      theme_minimal()
  }
  list(
    "beta" = p1,
    "contribution" = p2
  )
}