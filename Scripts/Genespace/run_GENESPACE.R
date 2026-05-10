Sys.setenv(TZ = "America/New_York")
#.libPaths("/N/scratch/dwtally/R_libs")

library(GENESPACE)
library(ggplot2)
# Directories
#genomeRepo <- "/N/scratch/dwtally/genespace_input"
wd <- "/N/project/Bracewell_fly/Dwayne/helixer_results/Genes_on_neoX"
path2mcscanx <- "/N/slate/dwtally/MCScanX"

# Automatically get species dirs and use them as genomeIDs
#speciesDirs <- list.dirs(genomeRepo, full.names = FALSE, recursive = FALSE)
#genomeIDs <- speciesDirs

# Register and parse genome annotations

#parse_annotations(
#  rawGenomeRepo = "/N/scratch/dwtally/genespace_input",
#  genomeDirs = list.files("/N/scratch/dwtally/genespace_input"),  # auto-collect all species
#  genomeIDs = list.files("/N/scratch/dwtally/genespace_input"),   # use same names for simplicity
#  gffString = "gff",   # matches *.gff
#  faString = "faa",  # matches *.protein.faa
#  genespaceWd = "/N/scratch/dwtally/genespace_run"
#)


# Initialize run
gpar <- init_genespace(
  wd = wd,
  path2mcscanx = path2mcscanx
#  blkSize = 3,
#  blkRadius = 250,
#  synBuff = 200
)

# Execute run
out <- run_genespace(gpar)

# Save session info
writeLines(capture.output(sessionInfo()), file.path(wd, "sessionInfo.txt"))
cat("✅ GENESPACE run complete.\n")

