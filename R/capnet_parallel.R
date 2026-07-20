# Dispatch a list of independent tasks, optionally in parallel via the base
# parallel package. A PSOCK cluster (parallel::makeCluster()/parLapply()) is
# used instead of future/future.apply: the latter incurs a large, fixed
# per-task dispatch overhead (a blocking worker-liveness handshake plus a
# global-environment reset on every single future) that dominates and can
# make parallel execution *slower* than sequential for the many small
# (fold, alpha) / walk-forward-step tasks capnet dispatches. PSOCK clusters
# work cross-platform, including Windows.
.capnet_lapply <- function(X, FUN, parallel = FALSE, workers = NULL) {
  if (!isTRUE(parallel)) return(lapply(X, FUN))

  # FUN is a closure defined inside the caller (cv_capnet()/walk_capnet()),
  # so its enclosing frame may still hold unforced argument promises (e.g.
  # gamma, alpha, lambda) if the caller never touched them before creating
  # FUN. An unforced promise evaluates in the *caller's* environment, which
  # does not exist on a fresh worker process, so it must be forced (turned
  # into a concrete value bound in FUN's own frame) before FUN is serialized.
  .capnet_force_env(environment(FUN))

  workers <- workers %||% max(1, parallel::detectCores() - 1)
  cl <- parallel::makeCluster(workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  parallel::parLapply(cl, X, FUN)
}

.capnet_force_env <- function(envir) {
  nms <- setdiff(ls(envir = envir, all.names = TRUE), "...")
  for (nm in nms) {
    tryCatch(get(nm, envir = envir), error = function(e) NULL)
  }
  invisible(NULL)
}
