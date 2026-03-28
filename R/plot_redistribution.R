#' Plot redistribution of coefficients and contributions
#' 
#' Compares regression coefficients and per-feature contributions between an
#' uncapped and a capped \code{capnet} fit. Produces either side-by-side
#' ("raw") bar charts of levels or ("delta") bar charts of changes.
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_bar coord_flip labs
#' @importFrom ggplot2 theme_minimal position_dodge
#' @importFrom stats reorder
#' @importFrom utils head
#'
#' @param uncapped An object returned by \code{capnet()} fit without a 
#'  contribution cap(i.e., \code{gamma = 0} or large \code{L}).
#' @param capped An object returned by \code{capnet()} fit with contribution
#'  caps applied.
#' @param X Optional numeric matrix of features used to compute contributions.
#'  If \code{NULL}, the function tries \code{capped$newx}.
#' @param multiplier Optional numeric scalar or length-\eqn{n} vector of 
#'  multipliers to scale contributions. If \code{NULL}, uses 
#'  \code{capped$multiplier} when available; otherwise defaults to 1.
#' @param top_n Integer (\eqn{\ge 1}) or \code{Inf}. Number of features to 
#'  display. Features are ordered by input column order; set \code{top_n = Inf}
#'  to include all.
#' @param plot Character string, one of \code{"raw"} or \code{"delta"}.
#'  If \code{"raw"}, shows side-by-side bar charts of coefficients and
#'  contributions for each fit. If \code{"delta"}, shows bar charts of the
#'  capped minus uncapped changes.
#'
#' @return A list of two \code{ggplot} objects, which are also printed::
#' \describe{
#'   \item{\code{p1}}{If \code{plot = "raw"}: side-by-side bar chart of
#'     coefficients for uncapped vs capped models. If \code{plot = "delta"}:
#'     bar chart of changes in coefficients (\code{capped - uncapped}).}
#'   \item{\code{p2}}{If \code{plot = "raw"}: side-by-side bar chart of
#'     contributions for uncapped vs capped models. If \code{plot = "delta"}:
#'     bar chart of changes in contributions (\code{capped - uncapped}).}
#' }
#' The function prints both plots and invisibly returns the list
#' \code{list("beta" = p1, "contribution" = p2)}
#' 
#' @details
#' Contributions are calculated as column-wise products
#' \eqn{contrib_j=\mathrm{mean}(X_{\cdot j})\times\beta_j\times multiplier_j}.
#' 
#' @seealso [capnet()]
#' 
#' @examples
#' set.seed(1)
#' n <- 40; p <- 6
#' X <- matrix(rnorm(n * p), n, p)
#' beta <- c(1.5, 0.8, 0.2, rep(0, p - 3))
#' y <- as.numeric(X %*% beta + rnorm(n))
#' fit_uncap <- capnet(X, y, lambda = 0.05, alpha = 0.5, gamma = 0, L = 1)
#' fit_cap   <- capnet(X, y, lambda = 0.05, alpha = 0.5, gamma = 1, L = 1)
#' plot_redistribution(fit_uncap, fit_cap, plot = "delta")
#' 
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
  
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("V", seq_len(ncol(X)))
  }
  
  df <- data.frame(
    variable = colnames(X),
    beta_uncapped = beta_uncapped,
    beta_capped = beta_capped,
    contribution_uncapped = mean_contribution_uncapped,
    contribution_capped = mean_contribution_capped
  )
  
  df$delta_beta <- df$beta_capped - df$beta_uncapped
  df$delta_contribution <- df$contribution_capped - df$contribution_uncapped
  
  if (is.finite(top_n)) {
    df_top <- head(df, top_n)
    var_level = unique(colnames(X))[1:top_n]
  } else {
    df_top <- df
    var_level = unique(colnames(X))
  }
  
  if (plot == "delta") {
    p1 <- ggplot(df_top, aes(x = reorder(.data$variable, .data$delta_beta), 
                             y = .data$delta_beta)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      theme_minimal() + 
      labs(
        title = "Change in Coefficients",
        x = NULL, y = NULL
      )
    
    p2 <- ggplot(df_top, aes(x = reorder(.data$variable, .data$delta_contribution), 
                             y = .data$delta_contribution)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      theme_minimal() + 
      labs(
        title = "Change in Contributions",
        x = NULL, y = NULL
      )
  } else {
    df_beta <- data.frame(
      variable = factor(rep(df_top$variable, times = 2), levels = var_level),
      model = rep(c("Uncapped", "Capped"), each = nrow(df_top)),
      beta = c(df_top$beta_uncapped, df_top$beta_capped)
    )
    p1 <- ggplot(df_beta, aes(x = .data$variable, 
                              y = .data$beta, 
                              fill = .data$model)) +
      geom_col(position = position_dodge(width = 0.6), width = 0.6) +
      labs(title = "Beta Comparison", y = "Beta", x = NULL, fill = "Model") +
      theme_minimal()
    
    df_contrib <- data.frame(
      variable = factor(rep(df_top$variable, times = 2), levels = var_level),
      model = rep(c("Uncapped", "Capped"), each = nrow(df_top)),
      contribution = c(df_top$contribution_uncapped, df_top$contribution_capped)
    )
    p2 <- ggplot(df_contrib, aes(x = .data$variable, 
                                 y = .data$contribution, 
                                 fill = .data$model)) +
      geom_col(position = position_dodge(width = 0.6), width = 0.6) +
      labs(title = "Contribution Comparison", y = "Contribution", x = NULL, fill = "Model") +
      theme_minimal()
  }
  list(
    "beta" = p1,
    "contribution" = p2
  )
}
