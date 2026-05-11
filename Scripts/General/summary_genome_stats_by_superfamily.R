#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

infile <- "genome_stats_by_species.csv"
outfile_wide <- "genome_stats_summary_by_superfamily_wide.csv"
outfile_long <- "genome_stats_summary_by_superfamily_long.csv"

df <- read_csv(infile, show_col_types = FALSE)

# Ensure expected columns exist
needed <- c("group", "chromosomes", "genome_size_mb", "gc_percent", "genes")
missing <- setdiff(needed, names(df))
if (length(missing) > 0) {
  stop(paste("Missing required columns in CSV:", paste(missing, collapse = ", ")))
}

# Long format: one row per (group, metric, value)
long <- df %>%
  select(group, chromosomes, genome_size_mb, gc_percent, genes) %>%
  pivot_longer(
    cols = c(chromosomes, genome_size_mb, gc_percent, genes),
    names_to = "metric",
    values_to = "value"
  )

# Summary stats per group/metric
summary_long <- long %>%
  group_by(group, metric) %>%
  summarise(
    n = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    iqr = IQR(value, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(metric, group)

# Wide version (columns like mean_genome_size_mb, median_genes, etc.)
summary_wide <- summary_long %>%
  pivot_wider(
    id_cols = group,
    names_from = metric,
    values_from = c(n, mean, sd, median, iqr, min, max),
    names_glue = "{.value}_{metric}"
  ) %>%
  arrange(group)

# Write outputs
write_csv(summary_long, outfile_long)
write_csv(summary_wide, outfile_wide)

# Print a readable preview to stdout (nice for logs)
print(summary_long, n = 200)

message("Wrote: ", outfile_long)
message("Wrote: ", outfile_wide)

