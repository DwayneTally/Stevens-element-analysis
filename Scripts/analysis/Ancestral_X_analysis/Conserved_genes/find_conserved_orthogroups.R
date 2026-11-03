#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage:\n",
    "  Rscript scan_conserved_thresholds.R Orthogroups.GeneCount.tsv 0.70 [out_prefix]\n\n",
    "Example:\n",
    "  Rscript scan_conserved_thresholds.R Orthogroups.GeneCount.tsv 0.70 conserved_scan_from_70\n"
  )
}

orthogroups_file <- args[1]
start_threshold <- as.numeric(args[2])
out_prefix <- if (length(args) >= 3) args[3] else "conserved_orthogroup_threshold_scan"

if (is.na(start_threshold) || start_threshold <= 0 || start_threshold > 1) {
  stop("start_threshold must be > 0 and <= 1")
}

step <- 0.05

# Read file
df <- read_tsv(orthogroups_file, show_col_types = FALSE)

species_cols <- setdiff(colnames(df), c("Orthogroup", "Total"))
n_species <- length(species_cols)

if (n_species == 0) {
  stop("No species columns found in input file.")
}

# Presence = at least one gene; duplication allowed
presence_mat <- df[, species_cols] >= 1
presence_count <- rowSums(presence_mat, na.rm = TRUE)

# Thresholds from start to 1.00
thresholds <- seq(start_threshold, 1.0, by = step)
thresholds[length(thresholds)] <- min(thresholds[length(thresholds)], 1.0)

results <- lapply(thresholds, function(threshold) {
  threshold_species <- ceiling(threshold * n_species)
  n_orthogroups <- sum(presence_count >= threshold_species)

  data.frame(
    threshold_fraction = threshold,
    threshold_percent = threshold * 100,
    minimum_species_required = threshold_species,
    orthogroups_passing = n_orthogroups
  )
})

results_df <- bind_rows(results)

# Write results table
write_csv(results_df, paste0(out_prefix, ".csv"))

cat("Total species:", n_species, "\n\n")
print(results_df)
cat("\nWrote table to:", paste0(out_prefix, ".csv"), "\n")

# Plot
p <- ggplot(results_df, aes(x = factor(threshold_percent), y = orthogroups_passing)) +
  geom_col() +
  labs(
    x = "Presence threshold (%)",
    y = "Orthogroups passing",
    title = "Conserved orthogroups across threshold scan"
  ) +
  theme_classic(base_size = 14)

ggsave(paste0(out_prefix, ".pdf"), p, width = 8, height = 5)
ggsave(paste0(out_prefix, ".png"), p, width = 8, height = 5, dpi = 300)

cat("Wrote plot to:", paste0(out_prefix, ".pdf"), "\n")
cat("Wrote plot to:", paste0(out_prefix, ".png"), "\n")
