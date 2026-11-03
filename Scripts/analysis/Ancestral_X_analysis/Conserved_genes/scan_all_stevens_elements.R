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
    "  Rscript scan_all_stevens_elements.R <base_directory> <start_threshold> [out_prefix]\n",
  )
}

base_dir <- args[1]
start_threshold <- as.numeric(args[2])

out_prefix <- if (length(args) >= 3) {
  args[3]
} else {
  "stevens_conserved_orthogroup_scan"
}

if (
  is.na(start_threshold) ||
  start_threshold <= 0 ||
  start_threshold > 1
) {
  stop("start_threshold must be > 0 and <= 1")
}

if (!dir.exists(base_dir)) {
  stop("Base directory does not exist: ", base_dir)
}

step <- 0.05

elements <- c(
  "A",
  "C",
  "E",
  "G",
  "H",
  "X"
)

stevens_colors <- c(
  "X" = "#c83838ff",
  "A" = "#21405fff",
  "C" = "#fec52dff",
  "E" = "#aca793ff",
  "G" = "#007e87ff",
  "H" = "#a0892cff"
)


find_gene_count_file <- function(element_dir) {

  files <- list.files(
    element_dir,
    pattern = "^Orthogroups\\.GeneCount\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(files) == 0) {
    return(NA_character_)
  }

  if (length(files) > 1) {

    cat(
      "\nWARNING: Multiple Orthogroups.GeneCount.tsv files found in:\n",
      element_dir,
      "\n",
      sep = ""
    )

    for (f in files) {
      cat("  ", f, "\n", sep = "")
    }

    info <- file.info(files)

    files <- files[
      order(
        info$mtime,
        decreasing = TRUE
      )
    ]

    cat(
      "Using newest:\n  ",
      files[1],
      "\n",
      sep = ""
    )
  }

  return(files[1])
}

process_element <- function(element) {

  element_dir <- file.path(
    base_dir,
    element
  )

  cat("\n")
  cat("STEVENS ELEMENT ", element, "\n", sep = "")

  if (!dir.exists(element_dir)) {

    cat(
      "Directory missing: ",
      element_dir,
      "\n",
      sep = ""
    )

    return(NULL)
  }


  orthogroups_file <- find_gene_count_file(
    element_dir
  )

  if (is.na(orthogroups_file)) {

    cat(
      "No Orthogroups.GeneCount.tsv found.\n"
    )

    return(NULL)
  }

  cat(
    "Gene count file:\n  ",
    orthogroups_file,
    "\n",
    sep = ""
  )

  df <- read_tsv(
    orthogroups_file,
    show_col_types = FALSE
  )


  species_cols <- setdiff(
    colnames(df),
    c(
      "Orthogroup",
      "Total"
    )
  )

  n_species <- length(
    species_cols
  )

  if (n_species == 0) {

    stop(
      "No species columns found for element ",
      element
    )
  }


  tcast_candidates <- species_cols[
    grepl(
      "Tribolium_castaneum",
      species_cols,
      fixed = TRUE
    )
  ]

  if (length(tcast_candidates) != 1) {

    stop(
      "Expected exactly one Tribolium castaneum column for element ",
      element,
      ", but found: ",
      paste(
        tcast_candidates,
        collapse = ", "
      )
    )
  }

  tcast_col <- tcast_candidates[1]

  cat(
    "Tribolium castaneum column:\n  ",
    tcast_col,
    "\n",
    sep = ""
  )


  # Total T. castaneum genes

  total_tcast_genes <- sum(
    df[[tcast_col]],
    na.rm = TRUE
  )

  if (total_tcast_genes == 0) {

    stop(
      "No Tribolium castaneum genes found for element ",
      element
    )
  }

  cat(
    "Total T. castaneum genes represented: ",
    total_tcast_genes,
    "\n",
    sep = ""
  )


  # Presence matrix
  # Presence = >=1 gene in a species

  presence_mat <- (
    df[, species_cols] >= 1
  )

  presence_count <- rowSums(
    presence_mat,
    na.rm = TRUE
  )


  thresholds <- seq(
    start_threshold,
    1.0,
    by = step
  )

  thresholds <- thresholds[
    thresholds <= 1.000001
  ]

  if (
    length(thresholds) == 0 ||
    tail(thresholds, 1) < 0.999
  ) {

    thresholds <- c(
      thresholds,
      1.0
    )
  }

  thresholds <- unique(
    round(
      thresholds,
      10
    )
  )


  # Calculate conserved reference genes

  results <- lapply(
    thresholds,
    function(threshold) {

      threshold_species <- ceiling(
        threshold * n_species
      )

      passing <- (
        presence_count >= threshold_species
      )

      # Raw number of orthogroups passing
      n_orthogroups <- sum(
        passing
      )

      # Number of T. castaneum reference genes
      # represented in passing orthogroups
      conserved_reference_genes <- sum(
        df[[tcast_col]][passing],
        na.rm = TRUE
      )

      # Percent of reference genes conserved
      conserved_percent <- (
        conserved_reference_genes /
        total_tcast_genes
      ) * 100


      data.frame(
        element = element,
        n_species = n_species,
        threshold_fraction = threshold,
        threshold_percent = threshold * 100,
        minimum_species_required = threshold_species,
        orthogroups_passing = n_orthogroups,
        total_tcast_genes = total_tcast_genes,
        conserved_reference_genes = conserved_reference_genes,
        conserved_percent = conserved_percent
      )
    }
  )


  results_df <- bind_rows(
    results
  )


  results_df <- results_df %>%
    arrange(
      desc(threshold_percent)
    )


  element_out_dir <- file.path(
    element_dir,
    "conserved_orthogroup_scan"
  )

  dir.create(
    element_out_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )


  element_csv <- file.path(
    element_out_dir,
    paste0(
      element,
      "_conserved_threshold_scan.csv"
    )
  )

  write_csv(
    results_df,
    element_csv
  )

  threshold_levels <- sort(
    unique(
      results_df$threshold_percent
    ),
    decreasing = TRUE
  )


  # Element-specific bar plot

  p <- ggplot(
    results_df,
    aes(
      x = factor(
        threshold_percent,
        levels = threshold_levels
      ),
      y = conserved_percent
    )
  ) +

    geom_col(
      fill = stevens_colors[[element]]
    ) +

    scale_y_continuous(
      limits = c(
        0,
        100
      ),
      breaks = seq(
        0,
        100,
        by = 20
      ),
      expand = expansion(
        mult = c(
          0,
          0.02
        )
      )
    ) +

    labs(
      x = "Presence threshold (%)",
      y = "T. castaneum reference genes conserved (%)",

      title = paste(
        "Stevens element",
        element
      ),

      subtitle = paste0(
        n_species,
        " species; ",
        total_tcast_genes,
        " T. castaneum reference genes"
      )
    ) +

    theme_classic(
      base_size = 14
    )


  ggsave(
    file.path(
      element_out_dir,
      paste0(
        element,
        "_conserved_threshold_scan_normalized.pdf"
      )
    ),
    p,
    width = 8,
    height = 5
  )


  ggsave(
    file.path(
      element_out_dir,
      paste0(
        element,
        "_conserved_threshold_scan_normalized.png"
      )
    ),
    p,
    width = 8,
    height = 5,
    dpi = 300
  )


  cat(
    "Species: ",
    n_species,
    "\n",
    sep = ""
  )

  cat(
    "T. castaneum reference genes: ",
    total_tcast_genes,
    "\n",
    sep = ""
  )

  print(
    results_df
  )

  return(
    results_df
  )
}

all_results <- lapply(
  elements,
  process_element
)


all_results <- all_results[
  !vapply(
    all_results,
    is.null,
    logical(1)
  )
]


if (length(all_results) == 0) {

  stop(
    "No Stevens element datasets were successfully processed."
  )
}


combined_df <- bind_rows(
  all_results
)

combined_df$element <- factor(
  combined_df$element,
  levels = elements
)


threshold_levels <- sort(
  unique(
    combined_df$threshold_percent
  ),
  decreasing = TRUE
)


combined_df$threshold_factor <- factor(
  combined_df$threshold_percent,
  levels = threshold_levels
)

combined_csv <- file.path(
  base_dir,
  paste0(
    out_prefix,
    "_normalized.csv"
  )
)


write_csv(
  combined_df,
  combined_csv
)


denominator_df <- combined_df %>%
  select(
    element,
    n_species,
    total_tcast_genes
  ) %>%
  distinct()


denominator_csv <- file.path(
  base_dir,
  paste0(
    out_prefix,
    "_reference_gene_counts.csv"
  )
)


write_csv(
  denominator_df,
  denominator_csv
)

p_combined <- ggplot(
  combined_df,
  aes(
    x = threshold_factor,
    y = conserved_percent,
    fill = element
  )
) +

  geom_col(
    position = position_dodge(
      width = 0.8
    ),
    width = 0.75
  ) +

  scale_fill_manual(
    values = stevens_colors
  ) +

  scale_y_continuous(
    limits = c(
      0,
      100
    ),
    breaks = seq(
      0,
      100,
      by = 20
    ),
    expand = expansion(
      mult = c(
        0,
        0.02
      )
    )
  ) +

  labs(
    x = "Presence threshold (%)",
    y = "T. castaneum reference genes conserved (%)",
    fill = "Stevens element",
    title = "Gene conservation across Stevens elements"
  ) +

  theme_classic(
    base_size = 14
  )


combined_pdf <- file.path(
  base_dir,
  paste0(
    out_prefix,
    "_normalized_grouped.pdf"
  )
)


combined_png <- file.path(
  base_dir,
  paste0(
    out_prefix,
    "_normalized_grouped.png"
  )
)


ggsave(
  combined_pdf,
  p_combined,
  width = 10,
  height = 6
)


ggsave(
  combined_png,
  p_combined,
  width = 10,
  height = 6,
  dpi = 300
)


p_line <- ggplot(
  combined_df,
  aes(
    x = threshold_percent,
    y = conserved_percent,
    group = element,
    color = element
  )
) +

  geom_line(
    linewidth = 1
  ) +

  geom_point(
    size = 2.5
  ) +

  scale_color_manual(
    values = stevens_colors
  ) +

  scale_x_reverse(
    breaks = threshold_levels
  ) +

  scale_y_continuous(
    limits = c(
      0,
      100
    ),
    breaks = seq(
      0,
      100,
      by = 20
    ),
    expand = expansion(
      mult = c(
        0,
        0.02
      )
    )
  ) +

  labs(
    x = "Presence threshold (%)",
    y = "T. castaneum reference genes conserved (%)",
    color = "Stevens element",
    title = "Gene conservation across Stevens elements"
  ) +

  theme_classic(
    base_size = 14
  )


line_pdf <- file.path(
  base_dir,
  paste0(
    out_prefix,
    "_normalized_line.pdf"
  )
)


line_png <- file.path(
  base_dir,
  paste0(
    out_prefix,
    "_normalized_line.png"
  )
)


ggsave(
  line_pdf,
  p_line,
  width = 9,
  height = 6
)


ggsave(
  line_png,
  p_line,
  width = 9,
  height = 6,
  dpi = 300
)


cat("DONE\n")


cat(
  "\nProcessed elements:\n  ",
  paste(
    unique(
      as.character(
        combined_df$element
      )
    ),
    collapse = ", "
  ),
  "\n",
  sep = ""
)


cat(
  "\nT. castaneum denominator counts:\n"
)


print(
  denominator_df
)


cat(
  "\nCombined normalized CSV:\n  ",
  combined_csv,
  "\n",
  sep = ""
)


cat(
  "\nReference gene-count CSV:\n  ",
  denominator_csv,
  "\n",
  sep = ""
)


cat(
  "\nGrouped normalized plot:\n  ",
  combined_pdf,
  "\n",
  sep = ""
)


cat(
  "\nLine normalized plot:\n  ",
  line_pdf,
  "\n",
  sep = ""
)
