#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: extract_chr_order.R <input_rda_path> <output_path_or_dir>")
}

rda_path <- args[1]
out_arg  <- args[2]

#Load the RDA
load(rda_path)

#Pull chromosomes table
chr_df <- srcd$sourceData$chromosomes
out_df <- chr_df[, c("genome", "chrLab")]

#If 2nd arg is a directory, write default filename inside it
is_dir_like <- dir.exists(out_arg) || grepl("/$", out_arg)

out_path <- if (is_dir_like) {
  dir.create(out_arg, recursive = TRUE, showWarnings = FALSE)
  file.path(out_arg, "Tribolium_order.csv")
} else {
  # Ensure parent dir exists
  dir.create(dirname(out_arg), recursive = TRUE, showWarnings = FALSE)
  out_arg
}

write.csv(out_df, file = out_path, row.names = FALSE, quote = FALSE)
cat("Wrote:", out_path, "\n")
