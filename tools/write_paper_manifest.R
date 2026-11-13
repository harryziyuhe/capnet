# Final source/input inventory; per-experiment session files retain run times.
files <- c("DESCRIPTION", list.files("R",full.names=TRUE),
  list.files("vignettes/paper",pattern="[.](R|Rmd)$",full.names=TRUE),
  list.files("inst/extdata",full.names=TRUE),
  list.files("tools",pattern="paper.*[.](R|sh)$",full.names=TRUE))
writeLines(c(paste("Finalized:", Sys.time()), capture.output(sessionInfo()),
  "Source and data MD5 checksums:",capture.output(tools::md5sum(files))),
  "manuscript/results/reproduction-manifest.txt")
