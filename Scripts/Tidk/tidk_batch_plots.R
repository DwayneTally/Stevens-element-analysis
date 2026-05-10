#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)

base_dir <- if (length(args) >= 1) args[1] else "/N/scratch/dwtally/tidk_search_results"
window_bp <- if (length(args) >= 2) as.numeric(args[2]) else 10000
min_contig_bp <- if (length(args) >= 3) as.numeric(args[3]) else 250000

message("[info] base_dir      = ", base_dir)
message("[info] window_bp     = ", window_bp)
message("[info] min_contig_bp = ", min_contig_bp)

# find all *telomeric_repeat_windows.tsv anywhere under base_dir
tsv_files <- list.files(
  path = base_dir,
  pattern = "telomeric_repeat_windows\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(tsv_files) == 0) {
  stop("[error] No telomeric_repeat_windows.tsv files found under: ", base_dir)
}

message("[info] Found ", length(tsv_files), " TSV files")

plot_forward <- function(df) {
  ggplot(df, aes(window, forward_repeat_number, color = telomeric_repeat)) +
    geom_line() +
    facet_wrap(~id, scales = "free_x") +
    labs(x = "Window", y = "Forward repeat number", color = "Telomeric repeat") +
    theme_bw()
}

plot_reverse <- function(df) {
  ggplot(df, aes(window, reverse_repeat_number, color = telomeric_repeat)) +
    geom_line() +
    facet_wrap(~id, scales = "free_x") +
    labs(x = "Window", y = "Reverse repeat number", color = "Telomeric repeat") +
    theme_bw()
}

plot_total <- function(df) {
  ggplot(df, aes(window, telo_repeat_total, color = telomeric_repeat)) +
    geom_line() +
    facet_wrap(~id, scales = "free_x") +
    labs(x = "Window", y = "Total repeats (forward + reverse)", color = "Telomeric repeat") +
    theme_bw()
}

for (tsv in tsv_files) {
  genome_dir <- dirname(tsv)
  genome_name <- basename(genome_dir)

  message("\n[info] Processing: ", genome_name)
  message("[info]   TSV: ", tsv)

  # read tidk window table
  df <- read.table(file = tsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

  # sanity check required columns
  required <- c("id", "window", "forward_repeat_number", "reverse_repeat_number", "telomeric_repeat")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    warning("[warn] Skipping ", genome_name, " (missing columns: ", paste(missing, collapse = ", "), ")")
    next
  }

  df <- df %>%
    mutate(
      telo_repeat_total = forward_repeat_number + reverse_repeat_number
    )

  # estimate contig length from max(window) * window_bp
  contig_sizes <- df %>%
    group_by(id) %>%
    summarize(
      max_window = suppressWarnings(max(as.numeric(window), na.rm = TRUE)),
      contig_len_bp_est = max_window * window_bp,
      .groups = "drop"
    )

  df2 <- df %>%
    left_join(contig_sizes, by = "id")

  # filter to contigs >= min_contig_bp
  df_filt <- df2 %>%
    filter(!is.na(contig_len_bp_est) & contig_len_bp_est >= min_contig_bp)

  if (nrow(df_filt) == 0) {
    warning("[warn] No contigs passed filter in ", genome_name, " (min_contig_bp=", min_contig_bp, ")")
    # still write the augmented TSV for completeness
  }

  # write augmented TSV
  out_tsv <- file.path(genome_dir, paste0(genome_name, ".tidk_telomeric_repeat_windows.with_totals.tsv"))
  write.table(df2, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  message("[ok]   Wrote: ", out_tsv)

  # make plots (using filtered set if available, otherwise all)
  df_plot <- if (nrow(df_filt) > 0) df_filt else df2

  out_forward_pdf <- file.path(genome_dir, paste0(genome_name, ".tidk_forward.faceted.pdf"))
  out_reverse_pdf <- file.path(genome_dir, paste0(genome_name, ".tidk_reverse.faceted.pdf"))
  out_total_pdf   <- file.path(genome_dir, paste0(genome_name, ".tidk_total.faceted.pdf"))

  ggsave(out_forward_pdf, plot_forward(df_plot), width = 14, height = 10)
  ggsave(out_reverse_pdf, plot_reverse(df_plot), width = 14, height = 10)
  ggsave(out_total_pdf,   plot_total(df_plot),   width = 14, height = 10)

  message("[ok]   Plots:")
  message("       ", out_forward_pdf)
  message("       ", out_reverse_pdf)
  message("       ", out_total_pdf)
}

message("\n[done] All finished.")

