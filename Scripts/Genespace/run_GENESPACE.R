library(GENESPACE)
library(ggplot2)
# Directories
wd <- "/N/project/Bracewell_fly/Dwayne/helixer_results/Genes_on_neoX"
path2mcscanx <- "/N/slate/dwtally/MCScanX"

# Initialize run
gpar <- init_genespace(
  wd = wd,
  path2mcscanx = path2mcscanx
)

# Execute run
out <- run_genespace(gpar)

# Save session info
writeLines(capture.output(sessionInfo()), file.path(wd, "sessionInfo.txt"))
cat("✅ GENESPACE run complete.\n")

