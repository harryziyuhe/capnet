# Shared, executable experiment definitions. Source after helpers.R.
# The vignettes and tools/run_paper_experiments.R use these same functions.
paper_finance_data <- function() {
  vars <- c("X", "target", "S2", "E2", "E3", "P9", "S1", "S5", "I2", "P8", "P10", "P12", "P13")
  read.csv(system.file("extdata", "data.csv", package = "capnet"))[, vars]
}
paper_wiki_data <- function() {
  w <- read.csv(system.file("extdata", "wiki_traffic.csv", package = "capnet"))
  articles <- setdiff(names(w), c("date", "total"))
  w$date <- as.Date(w$date)
  w <- w[order(w$date), ]
  w$target <- c(tail(w$total, -1), NA_real_)
  w[articles] <- w[articles] / 1000
  w <- w[!is.na(w$target), c("date", "target", articles)]
  w
}
wiki_L <- function(train_use, target_col, feature_cols, y_train) {
  vapply(feature_cols, function(v) {
    uni <- suppressWarnings(stats::glm(train_use[[target_col]] ~ train_use[[v]], family = "poisson"))
    3 * abs(stats::coef(uni)[2])
  }, numeric(1))
}
paper_finance <- function() {
  rolling_enet(paper_finance_data(), alpha = 0.5, window_years = 25,
    method = c("capnet_walk", "capnet", "postcap"), eval_start = "2020-01-01")
}
paper_wiki <- function() {
  # Five-row, week-like conditional batches: later observed covariates in each
  # slice enter estimation. This is not five-day-ahead prediction from day one.
  rolling_enet(paper_wiki_data(), date_col = "date", target_col = "target",
    alpha = 0, window_years = 1, method = c("capnet_walk", "capnet", "postcap"),
    capnet_gamma = 0, walk_gamma = c(10, 100, 1000, 10000, 100000),
    min_obs = 200, eval_start = "2023-01-01", family = "poisson",
    cap_formula = wiki_L, lambda_override = 0, walk = 5, parallel = FALSE)
}
paper_scores <- function(results, date_col, year = NULL) {
  d <- results$predictions
  if (!is.null(year)) d <- d[format(as.Date(d[[date_col]]), "%Y") == as.character(year), ]
  d |> dplyr::group_by(method, gamma) |>
    dplyr::group_modify(~ data.frame(n = nrow(.x), oos_r2 = calc_oos_r2(.x, date_col))) |>
    dplyr::ungroup()
}
paper_snapshot <- function(results, date_col, date, gamma, features, scale = 1) {
  d <- results$contributions
  methods <- c("Uncapped", "Post-Fit Cap", "Walk-forward Capnet")
  out <- data.frame(feature = features)
  for (method in methods) {
    row <- d[d[[date_col]] == date & d$method == method &
      d$gamma == if (method == "Walk-forward Capnet") gamma else 0, , drop = FALSE]
    stopifnot(nrow(row) == 1)
    out[[method]] <- as.numeric(row[1, features]) * scale
  }
  out
}
paper_plots <- function(results, date_col, gamma, window, wiki = FALSE) {
  pr <- results$predictions; ct <- results$contributions
  choose <- function(d) d[d$method %in% c("Uncapped", "Post-Fit Cap") |
    (d$method == "Walk-forward Capnet" & d$gamma == gamma), , drop = FALSE]
  pr <- choose(pr); ct <- choose(ct)
  pieces <- split(pr, pr$method)
  rolled <- dplyr::bind_rows(lapply(pieces, function(d) {
    ans <- rolling_oos_r2(d, window, date_col)
    ans$date <- as.Date(ans[[date_col]]); ans
  }))
  ct$date <- as.Date(ct[[date_col]])
  ct$maximum <- apply(abs(as.matrix(ct[, results$feature_cols, drop = FALSE])), 1, max)
  stopifnot(all(ct$maximum >= 0), all(is.finite(ct$maximum)))
  colors <- c("Uncapped" = "steelblue", "Post-Fit Cap" = "firebrick",
    "Walk-forward Capnet" = "darkolivegreen4")
  styling <- list(ggplot2::scale_color_manual(values = colors), ggplot2::theme_bw(base_size = 13),
    ggplot2::theme(legend.position = "top", legend.title = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()))
  if (wiki) {
    # A log axis compresses the single uncapped spike while retaining all
    # three methods, dates, and contributions on the same numeric scale.
    stopifnot(all(ct$maximum > 0))
    p2 <- ggplot2::ggplot(ct,ggplot2::aes(date,maximum,color=method)) +
      ggplot2::geom_line(linewidth=.65) +
      ggplot2::scale_y_log10(breaks=c(.05,.1,.2,.5,1,2,5),
        labels=c("0.05","0.1","0.2","0.5","1","2","5")) +
      ggplot2::labs(x="Date",y="Maximum absolute contribution (log-mean; log axis)") + styling
    p1 <- NULL
  } else {
    # A 252-observation score first exists after the initial evaluation year.
    visible <- rolled[is.finite(rolled$r2_roll),,drop=FALSE]
    p1 <- ggplot2::ggplot(visible, ggplot2::aes(date, r2_roll, color = method)) +
      ggplot2::geom_line(linewidth = .65) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      ggplot2::scale_x_date(limits=range(visible$date),expand=ggplot2::expansion(mult=c(0,.01))) +
      ggplot2::labs(x = "Date", y = expression(R[OOS]^2)) + styling
    p2 <- ggplot2::ggplot(ct, ggplot2::aes(date, maximum, color = method)) +
      ggplot2::geom_line(linewidth = .65) +
      ggplot2::labs(x = "Date", y = "Maximum absolute contribution") + styling
  }
  list(rolling = p1, contributions = p2, rolling_data = rolled, maximum_data = ct)
}
paper_runtime <- function(repetitions = 3L) {
  # Entire serial backtest, including monthly CV and output assembly. The active
  # penalty is identical for joint and walk fits. Seeds reset on every call.
  args <- list(glmnet = list(method = "glmnet"),
    capnet = list(method = "capnet", capnet_gamma = 0.1),
    walk_capnet = list(method = "capnet_walk", walk_gamma = 0.1))
  rows <- list()
  for (rep in seq_len(repetitions)) for (model in names(args)) {
    message("Timing ", model, " repetition ", rep)
    gc()
    elapsed <- system.time(ans <- do.call(rolling_enet, c(list(df = paper_finance_data(),
      alpha = .5, window_years = 25, eval_start = "2020-01-01"), args[[model]])))[["elapsed"]]
    stopifnot(all(is.finite(as.numeric(ans$predictions$predicted))))
    rows[[length(rows)+1L]] <- data.frame(repetition = rep, model = model, elapsed = elapsed)
  }
  dplyr::bind_rows(rows)
}
