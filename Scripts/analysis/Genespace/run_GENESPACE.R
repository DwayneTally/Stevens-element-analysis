library(GENESPACE)
library(ggplot2)
wd <- "genespace_run"
path2mcscanx <- "~/MCScanX"

# Initialize run
gpar <- init_genespace(
  wd = wd,
  path2mcscanx = path2mcscanx
)

out <- run_genespace(gpar)

#Save session info
writeLines(capture.output(sessionInfo()), file.path(wd, "sessionInfo.txt"))
cat("GENESPACE run complete.\n")

