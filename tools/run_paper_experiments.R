# Run from the repository root with a source-installed capnet on R_LIBS_USER.
# Usage: Rscript tools/run_paper_experiments.R finance|wiki|simulation|runtime|export|all
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1)
suppressPackageStartupMessages({library(capnet); library(dplyr); library(tidyr); library(lubridate)})
stopifnot(requireNamespace("glmnet", quietly = TRUE))
stopifnot(as.character(packageVersion("capnet")) == read.dcf("DESCRIPTION")[1,"Version"])
source("vignettes/paper/helpers.R")
source("vignettes/paper/experiments.R")
out <- "manuscript/results"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
args <- commandArgs(TRUE); task <- if (length(args)) args[1] else "all"
record <- function(name, value) {
  saveRDS(value, file.path(out, paste0(name, ".rds")))
  writeLines(c(paste("Completed:", Sys.time()), capture.output(sessionInfo()),
    "Thread environment:", capture.output(Sys.getenv(c("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","VECLIB_MAXIMUM_THREADS"))),
    "Input and source MD5:", capture.output(tools::md5sum(c("DESCRIPTION", list.files("R", full.names=TRUE),
      "vignettes/paper/helpers.R", "vignettes/paper/experiments.R", "vignettes/paper/simulation.R", list.files("inst/extdata", full.names=TRUE))))),
    file.path(out, paste0(name, "-session.txt")))
}
for (name in c("finance", "wiki", "simulation", "runtime")) if (task %in% c(name,"all")) {
  message("Starting ", name, " at ", Sys.time())
  if (name == "simulation") source("vignettes/paper/simulation.R")
  fun <- get(paste0("paper_", name))
  value <- fun()
  record(name, value)
  if (name %in% c("finance","wiki")) {
    stopifnot(all(is.finite(as.numeric(value$predictions$predicted))))
    print(paper_scores(value, if(name=="finance") "X" else "date"))
    print(value$diagnostics |> group_by(method, requested_gamma) |>
      summarise(n=n(), missing=sum(missing), reductions=sum(effective_gamma < requested_gamma, na.rm=TRUE),
        bad_stops=sum(!is.na(convergence) & convergence<0 & convergence != -1001), .groups="drop"))
  }
}
if (task %in% c("export","all")) source("vignettes/paper/export_results.R")
