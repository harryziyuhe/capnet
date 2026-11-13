# Predeclared replacement experiment for manuscript Tables 1 and 2.
paper_simulation <- function(repetitions = 10L, n_train = 400L, n_test = 100L,
                            seed = 20261112L, noise_high = 10,
                            high_prob_train = .10, high_prob_test = .20) {
  beta <- c(2.5, 1.9, 1.3, .9, .8, 0, 0, 0)
  p <- length(beta); Sigma <- outer(seq_len(p), seq_len(p), function(i,j) .5^abs(i-j))
  gammas <- c(.1, 1, 10, 100); kappas <- c(.5, 1, 2, 4)
  run_rep <- function(rep) {
    rows <- list(); index <- 0L
    message("Simulation replication ", rep, "/", repetitions)
    set.seed(seed + rep)
    latent <- matrix(rnorm((n_train+n_test)*p), ncol=p) %*% chol(Sigma)
    y <- drop(latent %*% beta + rnorm(n_train+n_test))
    proxy_noise <- matrix(rnorm((n_train+n_test)*p), ncol=p)
    spike_prob <- c(rep(high_prob_train,n_train),rep(high_prob_test,n_test))
    flip <- cbind(rbinom(n_train+n_test,1,spike_prob),
                  rbinom(n_train+n_test,1,spike_prob))
    for (scenario in c("homoskedastic", "mixture")) {
      noise_sd <- matrix(.25, n_train+n_test, p)
      if (scenario == "mixture") noise_sd[,1:2] <- .25 + (noise_high-.25)*flip
      observed <- latent + proxy_noise*noise_sd
      X <- observed[seq_len(n_train),]; z <- observed[n_train+seq_len(n_test),]
      yy <- y[seq_len(n_train)]; test_y <- y[n_train+seq_len(n_test)]
      colnames(X) <- colnames(z) <- paste0("X",seq_len(p))
      # Identical training folds within a replication. No test responses in tuning.
      set.seed(100000L+rep)
      foldid <- sample(rep(seq_len(10), length.out=n_train))
      cv <- glmnet::cv.glmnet(X, yy, alpha=.5, foldid=foldid, standardize=TRUE)
      lambda <- cv$lambda.min
      slopes <- vapply(seq_len(p), function(j) cov(X[,j],yy)/var(X[,j]), numeric(1))
      base_L <- pmax(abs(slopes)*apply(abs(X),2,quantile,probs=.95), .05*sd(yy))
      baseline <- capnet::capnet(X, yy, base_L, alpha=.5, lambda=lambda, gamma=0)
      if (baseline$convergence < 0 && baseline$convergence != -1001) stop("Baseline optimizer failed")
      base_mse <- mean((as.numeric(predict(baseline,z))-test_y)^2)
      add <- function(evaluation,gamma,kappa,pred, effective_gamma, stop_code=NA_integer_) {
        if (any(!is.finite(pred))) stop("Non-finite simulation prediction")
        index <<- index+1L
        rows[[index]] <<- data.frame(replication=rep, scenario=scenario, evaluation=evaluation,
          gamma=gamma,kappa=kappa,lambda=lambda,mse=mean((as.numeric(pred)-test_y)^2),
          baseline_mse=base_mse, reduced_rows=sum(effective_gamma<gamma),
          convergence=stop_code)
      }
      for (evaluation in c("joint","walk")) add(evaluation,0,NA_real_,predict(baseline,z),0,baseline$convergence)
      for (kappa in kappas) for (gamma in gammas) {
        L <- kappa*base_L
        fit <- capnet::capnet(X, yy, L, z=z,alpha=.5,lambda=lambda,gamma=gamma)
        if (fit$convergence < 0 && fit$convergence != -1001) stop("Joint optimizer failed")
        add("joint",gamma,kappa,predict(fit,z),gamma,fit$convergence)
        fitw <- capnet::walk_capnet(X, yy, L, z=z,alpha=.5,lambda=lambda,gamma=gamma,walk=1)
        add("walk",gamma,kappa,as.numeric(fitw$predictions),as.numeric(fitw$gammas))
      }
    }
    dplyr::bind_rows(rows)
  }
  workers <- as.integer(Sys.getenv("PAPER_WORKERS", "4"))
  stopifnot(is.finite(workers), workers >= 1)
  pieces <- if (.Platform$OS.type == "unix" && workers > 1) {
    parallel::mclapply(seq_len(repetitions), run_rep, mc.cores=workers,
      mc.preschedule=FALSE, mc.set.seed=FALSE)
  } else lapply(seq_len(repetitions),run_rep)
  if (any(vapply(pieces,inherits,logical(1),"try-error"))) stop("Simulation worker failed")
  raw <- dplyr::bind_rows(pieces)
  summary <- paper_simulation_summary(raw)
  list(raw=raw,summary=summary, config=list(repetitions=repetitions,n_train=n_train,n_test=n_test,
    beta=beta,Sigma=Sigma,seed=seed,noise_low=.25,noise_high=noise_high,
    high_prob_train=high_prob_train,high_prob_test=high_prob_test,
    unstable=1:2,alpha=.5,kappas=kappas,gammas=gammas))
}

# Compute standard errors before naming the group mean "mse", so dplyr does
# not replace the replication-level mse vector before sd() is evaluated.
paper_simulation_summary <- function(raw) {
  raw |> dplyr::mutate(delta=mse-baseline_mse) |>
    dplyr::group_by(scenario,evaluation,gamma,kappa) |>
    dplyr::summarise(repetitions=dplyr::n(), mcse=sd(mse)/sqrt(dplyr::n()),
      mse=mean(mse), paired_delta=mean(delta), paired_mcse=sd(delta)/sqrt(dplyr::n()),
      reduced_rows=sum(reduced_rows), .groups="drop")
}
