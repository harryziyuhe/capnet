# Render from the latest saved experiment objects; no fitting is repeated.
for (name in c("example_simulation", "example_elan", "example_wiki")) {
  input <- file.path("vignettes/paper", paste0(name,".Rmd"))
  message("Rendering ", input)
  output <- rmarkdown::render(input,quiet=TRUE)
  # Trailing padding in printed R output is not needed in the saved HTML.
  writeLines(sub("[[:blank:]]+$", "", readLines(output,warn=FALSE)),output)
}
