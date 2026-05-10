#!/usr/bin/env Rscript

Sys.setenv(TZ = "America/New_York")

suppressPackageStartupMessages({
  library(GENESPACE)
  library(ggplot2)
  library(data.table)
})

# === User parameters ===
genespace_dir      <- "rerun_coccinelloidea_tricast"
gsparams_file      <- file.path(genespace_dir, "results", "gsParams.rda")
output_plot        <- file.path(genespace_dir, "Coccinelloidea_custom_order_stevens_tricast.pdf")
custom_labels_file <- file.path(genespace_dir, "Coccinelloidea_output.csv")
chr_order_csv      <- file.path(genespace_dir, "riparian", "Tribolium_order.csv")

reference_genome   <- "Tribolium_castaneum_GCF_000002335.3"

#manually order genomes
genome_order <- c(
  "Tribolium_castaneum_GCF_000002335.3",
  "Dastarcus_helophoroides_GCA_028583605.1",
  "Cartodere_bifasciata_GCA_963920755.1",
  "Endomychus_coccineus_GCA_958510875.1",
  "Serangium_japonicum_GCA_040543525.2",
  "Megalocaria_dilatata_GCA_034642435.1",
  "Halyzia_sedecimguttata_GCA_937662695.2",
  "Myrrha_octodecimguttata_GCA_958510865.1",
  "Pjab_Bracewell_helixer_helixer",
  "Calvia_quattuordecimguttata_GCA_964059645.1",
  "Hippodamia_convergens_GCA_038502775.1",
  "Aphidecta_obliterata_GCA_963966015.1",
  "Harmonia_axyridis_GCF_914767665.1",
  "Coccinella_septempunctata_GCF_907165205.1",
  "Anisosticta_novemdecimpunctata_GCA_964188095.1",
  "Adalia_bipunctata_GCA_910592335.1",
  "Adalia_decempunctata_GCA_951802165.1"
)

if (!file.exists(gsparams_file)) stop(paste("❌ Cannot find gsParams.rda at", gsparams_file))
cat("✅ Loading gsParams from:", gsparams_file, "\n")
load(gsparams_file, verbose = FALSE)

missing_genomes <- setdiff(genome_order, gsParam$genomeIDs)
if (length(missing_genomes) > 0) {
  stop(paste0("❌ The following genomes are missing from gsParams:\n", paste(missing_genomes, collapse = "\n")))
}

#custom chromosome labels
if (!file.exists(custom_labels_file)) stop(paste("Cannot find custom_labels_file at", custom_labels_file))
cat("Loading custom chromosome labels from:", custom_labels_file, "\n")
custom_labels <- fread(custom_labels_file)

#key from the 'chr' field
extract_chr_key <- function(chr_field) {
  sapply(chr_field, function(chr_single) {
    if (grepl("Tribolium_freemani", chr_single)) {
      chr_key <- sub(".*chromosome:_", "", chr_single)
      chr_key <- sub(" .*", "", chr_key)
      return(chr_key)
    }
    if (grepl("Gnatocerus_cornutus", chr_single) && grepl("ctg00000023", chr_single)) {
      return("8_2")
    }
    accession <- regmatches(chr_single, regexpr("([A-Z]{2,3}[0-9A-Z_]*\\.[0-9]+)", chr_single, perl = TRUE))
    if (length(accession) ==0) return(NA_character_)
    return(accession)
  }, USE.NAMES = FALSE)
}

# Build mapping from chrSimple label
custom_labels[, cleanKey := extract_chr_key(chr)]
chr_map <- setNames(custom_labels$chrSimple, custom_labels$cleanKey)

# Display labels; unknown becomes "NA"
chrLabFun <- function(x) {
  x <- as.character(unlist(x))
  cleanX <- extract_chr_key(x)
  mapped <- chr_map[cleanX]
  mapped[is.na(mapped)] <- "NA"
  return(mapped)
}

ggthemes <- ggplot2::theme(panel.background = ggplot2::element_rect(fill = "black"))

canon <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("\\s+", "", x)   # remove whitespace
  x <- gsub("_", "", x)      # remove underscores (NC_007416.3 -> nc007416.3)
  return(x)
}

if (!file.exists(chr_order_csv)) stop(paste("chr_order_csv does not exist:", chr_order_csv))
cat("Loading chromosome order CSV:", chr_order_csv, "\n")
ord_df <- fread(chr_order_csv)

if (!all(c("genome", "chrLab") %in% names(ord_df))) {
  stop("chr_order_csv must have columns: genome, chrLab")
}

ord_df[, chrLab := canon(chrLab)]
ord_df[, rowRank := .I]  # EXACT file row order

dup_chr <- ord_df[duplicated(chrLab), unique(chrLab)]
if (length(dup_chr) > 0) {
  stop(
    "chrLab values are not unique.\n",
    "Duplicates:\n", paste(dup_chr, collapse = "\n"),
    "\nIf this happens, we must encode genome into the key."
  )
}

rank_map <- setNames(ord_df$rowRank, ord_df$chrLab)
norm_token <- function(x) canon(extract_chr_key(as.character(x)))

#If chrLab exists in CSV: use its exact row 
#If chrLab not in CSV: return NA, chromosome will not be plotted
ordFun <- function(chr) {
  tok <- norm_token(chr)
  r <- unname(rank_map[tok])

  miss <- unique(tok[is.na(r)])
  if (length(miss) > 0) {
    message("Dropping chromosomes not present in CSV: ", paste(miss, collapse = ", "))
  }
  return(r)
}

#Riparian plot
ripDat <- plot_riparian(
  gsParam = gsParam,
  genomeIDs = rev(genome_order),
  refGenome = reference_genome,
  refChrOrdFun = ordFun,
  reorderBySynteny = FALSE,
  useOrder = FALSE,
  addThemes = ggthemes,
  chrFill = "lightgrey",
  gapProp = 0.05,
  chrLabFun = chrLabFun
)

# Extract ggplot
if (inherits(ripDat, "ggplot")) {
  ripPlot <- ripDat
} else if (is.list(ripDat) && !is.null(ripDat$plotData$ggplotObj) && inherits(ripDat$plotData$ggplotObj, "ggplot")) {
  ripPlot <- ripDat$plotData$ggplotObj
} else {
  stop("Could not locate a ggplot object in ripDat.")
}

#Relabel genomeIDs for y-axis
b <- ggplot2::ggplot_build(ripPlot)
current_limits <- b$layout$panel_params[[1]]$y$get_labels()
short_genomeIDs <- sub("(_GCA|_GCF)_.*", "", current_limits)

ripPlot <- ripPlot +
  scale_y_discrete(limits = current_limits, labels = short_genomeIDs) +
  labs(y = NULL)

# Save plot
cat("Saving final plot to:", output_plot, "\n")
ggsave(filename = output_plot, plot = ripPlot)

cat("Done: chromosome order follows CSV row order; chromosomes not in CSV are dropped (not plotted).\n")

