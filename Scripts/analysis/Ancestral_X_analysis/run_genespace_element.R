args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript run_genespace_element.R <working_directory>")
}

Sys.setenv(TZ = "America/New_York")

library(GENESPACE)

wd <- normalizePath(
  args[1],
  mustWork = TRUE
)

path2mcscanx <- "/N/slate/dwtally/MCScanX"

cat("GENESPACE RUN\n")
cat("Working directory:", wd, "\n")
cat("Started:", as.character(Sys.time()), "\n")

gpar <- init_genespace(
  wd = wd,
  path2mcscanx = path2mcscanx
)

out <- run_genespace(gpar)

writeLines(
  capture.output(sessionInfo()),
  file.path(wd, "sessionInfo.txt")
)

cat("GENESPACE run complete\n")
cat("Finished:", as.character(Sys.time()), "\n")
