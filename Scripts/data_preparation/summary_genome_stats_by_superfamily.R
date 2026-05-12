#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

infile <- "genome_stats_by_species.csv"
outfile <- "genome_stats_summary_by_superfamily.csv"

df <- read_csv(infile, show_col_types = FALSE)

#Double check columns exist
needed <- c("group", "chromosomes", "genome_size_mb", "gc_percent", "genes")
missing <- setdiff(needed, names(df))
if (length(missing) > 0) {
  stop(paste("Missing required columns in CSV:", paste(missing, collapse = ", ")))
}

#one row per (group, metric, value)
summ <- df %>%
  select(group, chromosomes, genome_size_mb, gc_percent, genes) %>%
  pivot_longer(
    cols = c(chromosomes, genome_size_mb, gc_percent, genes),
    names_to = "metric",
    values_to = "value"
  )

#Summary stats per group
summary_init <- summ %>%
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

#columns like mean_genome_size_mb, median_genes, etc.
summary <- summary_init %>%
  pivot_wider(
    id_cols = group,
    names_from = metric,
    values_from = c(n, mean, sd, median, iqr, min, max),
    names_glue = "{.value}_{metric}"
  ) %>%
  arrange(group)

#Write outputs
write_csv(summary, outfile)

print(summary_long, n = 200)
message("Wrote: ", outfile)

