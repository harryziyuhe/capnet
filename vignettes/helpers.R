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

postcap_contrib <- function(object, newx, L) {
  beta <- as.numeric(object$beta)
  
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
  
  contrib_capped
}

make_contrib_tibble <- function(contrib,
                                test_use,
                                date_col,
                                pred_month,
                                method,
                                gamma,
                                lambda,
                                feature_cols) {
  contrib_df <- as.data.frame(contrib)
  names(contrib_df) <- feature_cols
  
  tibble::tibble(
    !!date_col := test_use[[date_col]],
    refit_month = pred_month,
    method = method,
    gamma = gamma,
    lambda = lambda
  ) %>%
    dplyr::bind_cols(tibble::as_tibble(contrib_df))
}

rolling_enet <- function(df,
                         date_col = "X",
                         target_col = "target",
                         alpha = 0.5,
                         window_years = 20,
                         nfolds = 10,
                         standardize = TRUE,
                         method = c("capnet_walk", "capnet", "glmnet", "postcap"),
                         capnet_gamma = c(0),
                         walk_gamma = c(0.1, 1, 10, 100),
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
  contrib_list <- vector("list", length(refit_months))
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
    month_contribs <- list()
    idx <- 1L
    
    if ("glmnet" %in% method) {
      preds <- as.numeric(stats::predict(cvfit, newx = x_test, s = lambda))
      
      beta <- as.numeric(glmnet::coef.glmnet(cvfit$glmnet.fit, s = lambda))[-1]
      contrib <- sweep(x_test, 2, beta, `*`)
      
      month_preds[[idx]] <- tibble::tibble(
        !!date_col := test_use[[date_col]],
        refit_month = test_use$.refit_month,
        method = "glmnet",
        gamma = NA,
        lambda = lambda,
        actual = test_use[[target_col]],
        predicted = preds
      )
      
      month_contribs[[idx]] <- make_contrib_tibble(
        contrib = contrib,
        test_use = test_use,
        date_col = date_col,
        pred_month = pred_month,
        method = "glmnet",
        gamma = NA,
        lambda = lambda,
        feature_cols = feature_cols
      )
      
      idx <- idx + 1L
    }
    
    if ("capnet" %in% method || "postcap" %in% method) {
      for (gamma_j in capnet_gamma) {
        fit_cap <- capnet::capnet(
          X = x_train,
          y = y_train,
          L = L,
          newx = x_test,
          lambda = lambda,
          alpha = alpha,
          gamma = gamma_j
        )
        
        if ("capnet" %in% method) {
          preds <- as.numeric(stats::predict(fit_cap, newdata = x_test))
          label <- ifelse(gamma_j == 0, "Uncapped", "Capnet")
          
          beta <- as.numeric(fit_cap$beta)
          contrib <- sweep(x_test, 2, beta, "*")
          
          month_preds[[idx]] <- tibble::tibble(
            !!date_col := test_use[[date_col]],
            refit_month = pred_month,
            method = label,
            gamma = gamma_j,
            lambda = lambda,
            actual = test_use[[target_col]],
            predicted = preds
          )
          
          month_contribs[[idx]] <- make_contrib_tibble(
            contrib = contrib,
            test_use = test_use,
            date_col = date_col,
            pred_month = pred_month,
            method = label,
            gamma = gamma_j,
            lambda = lambda,
            feature_cols = feature_cols
          )
          
          idx <- idx + 1L
        }
        
        if ("postcap" %in% method & gamma_j == 0) {
          a0 <- if (!is.null(fit_cap$a0)) as.numeric(fit_cap$a0) else 0
          contrib <- postcap_contrib(fit_cap, x_test, L)
          preds_postcap <- a0 + rowSums(contrib)
          
          month_preds[[idx]] <- tibble::tibble(
            !!date_col := test_use[[date_col]],
            refit_month = pred_month,
            method = "Post-Fit Cap",
            gamma = gamma_j,
            lambda = lambda,
            actual = test_use[[target_col]],
            predicted = preds_postcap
          )
          
          month_contribs[[idx]] <- make_contrib_tibble(
            contrib = contrib,
            test_use = test_use,
            date_col = date_col,
            pred_month = pred_month,
            method = "Post-Fit Cap",
            gamma = gamma_j,
            lambda = lambda,
            feature_cols = feature_cols
          )
          
          idx <- idx + 1L
        }
      }
    } 
    
    if ("capnet_walk" %in% method) {
      for (gamma_j in walk_gamma) {
        fit_walk <- capnet::walk_capnet(
          X = x_train,
          y = y_train,
          L = L,
          newx = x_test,
          lambda = lambda,
          alpha = alpha,
          gamma = gamma_j
        )
        
        month_preds[[idx]] <- tibble::tibble(
          !!date_col := test_use[[date_col]],
          refit_month = pred_month,
          method = "Walk-forward Capnet",
          gamma = gamma_j,
          lambda = lambda,
          actual = test_use[[target_col]],
          predicted = fit_walk$predictions
        )
        
        contrib <- fit_walk$feature_contributions
        
        month_contribs[[idx]] <- make_contrib_tibble(
          contrib = contrib,
          test_use = test_use,
          date_col = date_col,
          pred_month = pred_month,
          method = "Walk-forward Capnet",
          gamma = gamma_j,
          lambda = lambda,
          feature_cols = feature_cols
        )
        
        idx <- idx + 1L
      }
    }
    pred_list[[i]] <- dplyr::bind_rows(month_preds)
    contrib_list[[i]] <- dplyr::bind_rows(month_contribs)
  }
  
  predictions <- dplyr::bind_rows(pred_list) %>% 
    dplyr::arrange(.data[[date_col]])
  
  contributions <- dplyr::bind_rows(contrib_list) %>%
    dplyr::arrange(.data[[date_col]])
  
  metrics <- predictions %>% 
    dplyr::group_by(method, gamma) %>% 
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
    contributions = contributions,
    metrics = metrics,
    feature_cols = feature_cols,
    lambdas = lambda_table,
    capnet_gamma = capnet_gamma,
    walk_gamma = walk_gamma,
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