calc_moments <- function(data, na.rm = TRUE) {
  # keep only numeric columns
  data <- data[, sapply(data, is.numeric), drop = FALSE]
  
  res <- lapply(names(data), function(name) {
    x <- data[[name]]
    if (na.rm) x <- x[!is.na(x)]
    
    n <- length(x)
    if (n < 4) {
      return(data.frame(
        variable = name,
        mean = NA_real_,
        variance = NA_real_,
        excess_kurtosis = NA_real_
      ))
    }
    
    m <- mean(x)
    v <- var(x)
    s <- sd(x)
    
    # bias-corrected excess kurtosis
    if (s == 0) {
      k <- NA_real_
    } else {
      z <- (x - m) / s
      k <- (n*(n+1)/((n-1)*(n-2)*(n-3))) * sum(z^4) -
        (3*(n-1)^2/((n-2)*(n-3)))
    }
    
    data.frame(
      variable = name,
      mean = m,
      variance = v,
      excess_kurtosis = k
    )
  })
  
  dplyr::bind_rows(res)
}

moments_latex <- function(data, digits = 3) {
  tab <- calc_moments(data) %>%
    dplyr::mutate(
      dplyr::across(-variable, ~ round(., digits))
    )
  
  xt <- xtable::xtable(tab)
  
  print(
    xt,
    include.rownames = FALSE,
    sanitize.text.function = identity
  )
}

postcap_predict <- function(object, newx, L) {
  beta <- as.numeric(object$beta)
  a0 <- if (!is.null(object$a0)) as.numeric(object$a0) else 0
  
  contrib <- sweep(newx, 2, beta, `*`)
  
  if (length(L) == 1L) {
    contrib_capped <- pmax(pmin(contrib, L), -L)
  } else {
    if (length(L) != ncol(newx)) {
      stop("Length of L must be 1 or equal to ncol(newx).")
    }
    L_mat <- matrix(rep(L, each = nrow(newx)), nrow = nrow(newx))
    contrib_capped <- pmax(pmin(contrib, L_mat), -L_mat)
  }
  
  as.numeric(a0 + rowSums(contrib_capped))
}

rolling_enet <- function(df,
                         date_col = "X",
                         target_col = "target",
                         alpha = 0.5,
                         window_years = 20,
                         nfolds = 10,
                         standardize = TRUE,
                         method = c("capnet_walk", "capnet", "glmnet", "postcap"),
                         capnet_mu = c(0),
                         walk_mu = c(0.1, 1, 10, 100),
                         min_obs = 100,
                         eval_start = "2020-01-01",
                         seed = 42) {
  eval_start <- as.Date(eval_start)
  
  stopifnot(date_col %in% names(df), target_col %in% names(df))
  
  df <- df %>% 
    dplyr::mutate(
      .date = as.Date(.data[[date_col]]),
      .refit_month = lubridate::floor_date(.date, "month")
    ) %>% 
    dplyr::filter(!is.na(.date), !is.na(.data[[target_col]])) %>% 
    dplyr::arrange(.date)
  
  feature_cols <- setdiff(names(df), c(date_col, target_col, ".date", ".refit_month"))
  feature_cols <- feature_cols[sapply(df[feature_cols], is.numeric)]
  
  if (length(feature_cols) == 0) {
    stop("No numeric feature columns available.")
  }
  
  refit_months <- sort(unique(df$.refit_month))
  refit_months <- refit_months[refit_months >= lubridate::floor_date(eval_start, "month")]
  
  pred_list <- vector("list", length(refit_months))
  lambda_vec <- rep(NA_real_, length(refit_months))
  
  for (i in seq_along(refit_months)) {
    pred_month <- refit_months[i]
    train_end <- pred_month - lubridate::days(1)
    train_start <- pred_month %m-% lubridate::years(window_years)
    
    train_use <- df %>% 
      dplyr::filter(.date >= train_start, .date <= train_end) %>% 
      dplyr::select(dplyr::all_of(c(target_col, feature_cols))) %>% 
      dplyr::filter(dplyr::if_all(dplyr::everything(), ~ !is.na(.)))
    
    test_use <- df %>% 
      dplyr::filter(.refit_month == pred_month) %>% 
      dplyr::select(dplyr::all_of(c(date_col, target_col, feature_cols)), .refit_month) %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(feature_cols), ~ !is.na(.)))
    
    if (nrow(train_use) < min_obs || nrow(test_use) == 0) {
      pred_list[[i]] <- NULL
      next
    }
    
    x_train <- as.matrix(train_use[, feature_cols, drop = FALSE])
    y_train <- train_use[[target_col]]
    x_test <- as.matrix(test_use[, feature_cols, drop = FALSE])
    
    L <- sapply(feature_cols, function(v) {
      7 * abs(stats::cor(train_use[[v]], y_train, use = "complete.obs")) / 350
    })
    
    set.seed(seed)
    cvfit <- glmnet::cv.glmnet(
      x = x_train,
      y = y_train,
      alpha = alpha,
      nfolds = nfolds,
      standardize = standardize
    )
    lambda <- cvfit$lambda.min
    lambda_vec[i] <- lambda
    seed <- seed + 1
    
    month_preds <- list()
    idx <- 1L
    
    if ("glmnet" %in% method) {
      preds <- as.numeric(stats::predict(cvfit, newx = x_test, s = lambda))
      
      pred_list[[idx]] <- tibble::tibble(
        !!date_col := test_use[[date_col]],
        refit_month = test_use$.refit_month,
        method = "glmnet",
        mu = NA,
        lambda = lambda,
        actual = test_use[[target_col]],
        predicted = preds
      )
      idx <- idx + 1L
    }
    
    if ("capnet" %in% method || "postcap" %in% method) {
      for (mu_j in capnet_mu) {
        fit_cap <- capnet::capnet(
          X = x_train,
          y = y_train,
          L = L,
          newx = x_test,
          lambda = lambda,
          alpha = alpha,
          mu = mu_j
        )
        
        if ("capnet" %in% method) {
          preds <- as.numeric(stats::predict(fit_cap, newdata = x_test))
          label <- ifelse(mu_j == 0, "Uncapped", "Capnet")
          
          month_preds[[idx]] <- tibble::tibble(
            !!date_col := test_use[[date_col]],
            refit_month = pred_month,
            method = label,
            mu = mu_j,
            lambda = lambda,
            actual = test_use[[target_col]],
            predicted = preds
          )
          idx <- idx + 1L
        }
        
        if ("postcap" %in% method & mu_j == 0) {
          preds_postcap <- postcap_predict(fit_cap, x_test, L)
          
          month_preds[[idx]] <- tibble::tibble(
            !!date_col := test_use[[date_col]],
            refit_month = pred_month,
            method = "Post-Fit Cap",
            mu = mu_j,
            lambda = lambda,
            actual = test_use[[target_col]],
            predicted = preds_postcap
          )
          idx <- idx + 1L
        }
      }
    } 
    
    if ("capnet_walk" %in% method) {
      for (mu_j in walk_mu) {
        fit_walk <- capnet::walk_capnet(
          X = x_train,
          y = y_train,
          L = L,
          newx = x_test,
          lambda = lambda,
          alpha = alpha,
          mu = mu_j
        )
        
        month_preds[[idx]] <- tibble::tibble(
          !!date_col := test_use[[date_col]],
          refit_month = pred_month,
          method = "Walk-forward Capnet",
          mu = mu_j,
          lambda = lambda,
          actual = test_use[[target_col]],
          predicted = fit_walk$predictions
        )
        idx <- idx + 1L
      }
    }
    pred_list[[i]] <- dplyr::bind_rows(month_preds)
  }
  
  predictions <- dplyr::bind_rows(pred_list) %>% 
    dplyr::arrange(.data[[date_col]])
  
  metrics <- predictions %>% 
    dplyr::group_by(method, mu) %>% 
    dplyr::summarise(
      rmse = sqrt(mean((actual - predicted)^2, na.rm = TRUE)),
      mae = mean(abs(actual - predicted), na.rm = TRUE),
      cor = stats::cor(actual, predicted, use = "complete.obs"),
      .groups = "drop"
    )
  
  lambda_table <- tibble::tibble(
    refit_month = refit_months,
    lambda = lambda_vec
  ) %>%
    dplyr::filter(!is.na(lambda))
  
  list(
    predictions = predictions,
    metrics = metrics,
    feature_cols = feature_cols,
    lambdas = lambda_table,
    capnet_mu = capnet_mu,
    walk_mu = walk_mu,
    method = method
  )
}

calc_oos_r2 <- function(pred_df) {
  
  pred_df <- pred_df |>
    arrange(X)
  
  # expanding historical mean benchmark
  pred_df$mean_benchmark <- c(NA, cummean(head(pred_df$actual, -1)))
  
  pred_df <- pred_df |>
    dplyr::filter(!is.na(mean_benchmark))
  
  sse_model <- sum((pred_df$actual - pred_df$predicted)^2)
  sse_bench <- sum((pred_df$actual - pred_df$mean_benchmark)^2)
  
  r2_oos <- 1 - sse_model / sse_bench
  
  return(r2_oos)
}

rolling_oos_r2 <- function(pred_df, window = 20) {
  
  pred_df <- pred_df |>
    arrange(X)
  
  pred_df$mean_benchmark <- c(NA, cummean(head(pred_df$actual, -1)))
  
  pred_df <- pred_df |>
    filter(!is.na(mean_benchmark))
  
  pred_df <- pred_df |>
    mutate(
      se_model = (actual - predicted)^2,
      se_bench = (actual - mean_benchmark)^2
    )
  
  pred_df |>
    mutate(
      r2_roll = 1 - rollapply(se_model, window, sum, fill = NA, align = "right") /
        rollapply(se_bench, window, sum, fill = NA, align = "right")
    )
}