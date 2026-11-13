# Integrity checks on the saved scientific outputs; does not re-fit models.
library(dplyr)
source("vignettes/paper/helpers.R")
source("vignettes/paper/experiments.R")
# Regression fixture: a very large error followed by small errors must not
# cause cumulative-sum cancellation or impossible R-squared above one.
fixture <- data.frame(X=seq_len(12),actual=seq_len(12),predicted=seq_len(12)+.5)
fixture$predicted[2] <- 1e15
rolled <- rolling_oos_r2(fixture,3L)
y <- fixture$actual[-1]; forecast <- fixture$predicted[-1]
b <- cumsum(head(fixture$actual,-1))/seq_len(11)
direct <- vapply(seq_along(y),function(i) {
  if(i<3) return(NA_real_)
  ix <- (i-2):i
  1-sum((y[ix]-forecast[ix])^2)/sum((y[ix]-b[ix])^2)
},numeric(1))
stopifnot(isTRUE(all.equal(rolled$r2_roll,direct,tolerance=1e-12)),
  all(rolled$r2_roll<=1,na.rm=TRUE))
cat("Rolling score: extreme-scale cancellation regression check passed.\n")
out <- "manuscript/results"
sim <- readRDS(file.path(out,"simulation.rds"))
stopifnot(nrow(sim$raw)==30L*2L*34L, all(is.finite(sim$raw$mse)),
  identical(sort(unique(sim$raw$replication)),seq_len(30L)),
  all(sim$summary$repetitions==30L),all(is.finite(sim$summary$mcse)),all(is.finite(sim$summary$paired_mcse)))
stopifnot(all((sim$raw |> group_by(replication,scenario) |>
  summarise(n_lambda=n_distinct(lambda),.groups="drop"))$n_lambda==1L))
baselines <- sim$raw |> filter(gamma==0) |> group_by(replication,scenario) |>
  summarise(n_mse=n_distinct(mse),.groups="drop")
stopifnot(all(baselines$n_mse==1L))
cat("Simulation: 30 paired replications in both scenarios; shared lambda and baseline verified.\n")
for(name in c("finance","wiki")) {
  fit <- readRDS(file.path(out,paste0(name,".rds")))
  date_col <- if(name=="finance") "X" else "date"
  stopifnot(all(is.finite(as.numeric(fit$predictions$predicted))),sum(fit$diagnostics$missing)==0,
    all(is.na(fit$diagnostics$convergence) | fit$diagnostics$convergence>=0 | fit$diagnostics$convergence == -1001))
  reported <- read.csv(file.path(out,paste0(name,"-scores.csv")))
  for(i in seq_len(nrow(reported))) {
    r <- reported[i,]
    pr <- fit$predictions[fit$predictions$method==r$method & fit$predictions$gamma==r$gamma,]
    pr <- pr[order(pr[[date_col]]),]
    y <- pr$actual; forecast <- as.numeric(pr$predicted)
    benchmark <- cumsum(head(y,-1))/seq_len(length(y)-1L)
    direct <- 1-sum((tail(y,-1)-tail(forecast,-1))^2)/sum((tail(y,-1)-benchmark)^2)
    stopifnot(isTRUE(all.equal(direct,r$oos_r2,tolerance=1e-12)))
  }
  for(method in unique(fit$predictions$method)) {
    d <- fit$predictions[fit$predictions$method==method,]
    for(g in unique(d$gamma)) {
      curve <- rolling_oos_r2(d[d$gamma==g,],if(name=="finance")252 else 180,date_col)
      stopifnot(all(curve$r2_roll<=1+1e-10,na.rm=TRUE))
    }
  }
  maxima <- read.csv(file.path(out,paste0(name,"-maximum-contributions.csv")))
  for(i in seq_len(nrow(maxima))) {
    row <- maxima[i,]
    ct <- fit$contributions[as.character(fit$contributions[[date_col]])==row$date &
      fit$contributions$method==row$method & fit$contributions$gamma==row$gamma,fit$feature_cols,drop=FALSE]
    stopifnot(nrow(ct)==1L,abs(max(abs(as.numeric(ct)))-row$maximum)<1e-12,row$maximum>=0)
  }
  # Post-fit forecasts must differ from uncapped forecasts ONLY by the change
  # in summed contributions, on the family's link scale. No intercept refit.
  uc <- fit$contributions[fit$contributions$method=="Uncapped",]
  pc <- fit$contributions[fit$contributions$method=="Post-Fit Cap",]
  up <- fit$predictions[fit$predictions$method=="Uncapped",]
  pp <- fit$predictions[fit$predictions$method=="Post-Fit Cap",]
  stopifnot(identical(as.character(uc[[date_col]]),as.character(pc[[date_col]])),
    identical(as.character(uc[[date_col]]),as.character(up[[date_col]])),
    identical(as.character(pc[[date_col]]),as.character(pp[[date_col]])))
  link <- if(name=="wiki") log else identity
  observed_change <- link(as.numeric(pp$predicted))-link(as.numeric(up$predicted))
  expected_change <- rowSums(pc[,fit$feature_cols])-rowSums(uc[,fit$feature_cols])
  stopifnot(max(abs(observed_change-expected_change))<1e-8)
  cat(name,": finite forecasts, numerical stops, pooled scores, per-model maxima, and post-fit identity verified.\n")
}
runtime <- readRDS(file.path(out,"runtime.rds"))
stopifnot(nrow(runtime)==9L,all(table(runtime$model)==3L),all(is.finite(runtime$elapsed)),all(runtime$elapsed>0))
cat("Runtime: three finite positive repetitions per method.\n")
cat("All saved-result integrity checks passed.\n")
