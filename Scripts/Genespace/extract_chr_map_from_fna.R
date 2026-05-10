#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript extract_chr_map_from_fna.R <input_directory> <output.csv>")
}

input_dir <- args[1]
output_csv <- args[2]

# === Discover input files (handle .fna and .fna.gz), keep only no-family names ===
fna_files <- list.files(input_dir, pattern = "\\.fna(\\.gz)?$", full.names = TRUE)

# wanted: Genus_species_GC[AF]_digits.d_genomic.fna(.gz)?
nofam_pat <- "^[A-Za-z]+_[A-Za-z]+_GC[AF]_[0-9]+\\.[0-9]+_genomic\\.fna(\\.gz)?$"
fna_files <- fna_files[grepl(nofam_pat, basename(fna_files), ignore.case = TRUE)]

# deduplicate by real path (avoid symlink + real duplicates)
real_paths <- normalizePath(fna_files, mustWork = TRUE)
fna_files <- fna_files[!duplicated(real_paths)]

cat("Found", length(fna_files),
    "files after filtering for no-family names and dedup by real path.\n")

if (length(fna_files) == 0) {
  stop("❌ No matching .fna files found after filtering (expecting no-family names).")
}


# === OPTIONAL: prefix for contig numbers ("" means plain 1,2,3) ===
CONTIG_PREFIX <- ""   # e.g., set to "ctg" to get ctg1, ctg2, ...

# --- Species-specific scaffold map (force-include) ---
# Keys must match the tolower() base name (without "_genomic.fna")
# e.g., "tribolium_confusum_gca_019155225.1"
scaffold_map_list <- list(
  "tribolium_confusum_gca_019155225.1" = c(
    "JAGFVK010000006.1" = "Chr7_1",
    "JAGFVK010000007.1" = "Chr7_2",
    "JAGFVK010000009.1" = "Chr9_1",
    "JAGFVK010000010.1" = "Chr9_2"
  )
)

# --- helper: order contig-like IDs by numeric part (if any) ---
order_by_numeric_suffix <- function(ids) {
  if (length(ids) == 0) return(integer(0))
  # extract the last run of digits in each ID; fallback to NA
  num <- suppressWarnings(as.integer(sub(".*?(\\d+)(?:\\.[0-9]+)?$", "\\1", ids)))
  # if all NAs, keep original order; otherwise order by num (NAs last)
  if (all(is.na(num))) {
    return(seq_along(ids))
  } else {
    # stable order: numeric ascending, NAs after, tie-breaker = original order
    ord <- order(is.na(num), num, seq_along(ids))
    return(ord)
  }
}

extract_chr_info <- function(header, species_name, species_scaffold_map) {
  header <- sub("^>", "", header)
  parts <- strsplit(header, " ", fixed = TRUE)[[1]]
  seq_id <- parts[1]

  # 1) Force-map specific scaffolds if provided
  if (!is.null(species_scaffold_map) && seq_id %in% names(species_scaffold_map)) {
    chr_label <- species_scaffold_map[[seq_id]]
    return(data.frame(chr = seq_id,
                      chrSimple = chr_label,
                      genome = species_name,
                      stringsAsFactors = FALSE))
  }

  # 2) Parse chromosome / linkage group labels from header, if present
  chr_label <- NA
  if (grepl("chromosome[: ]\\s*[^, ]+", header, ignore.case = TRUE)) {
    chr_label <- sub(".*chromosome[: ]\\s*([^, ]+).*", "\\1", header, ignore.case = TRUE)
  } else if (grepl("linkage group\\s*LG[0-9XY]+", header, ignore.case = TRUE)) {
    chr_label <- sub(".*linkage group\\s*(LG[0-9XY]+).*", "\\1", header, ignore.case = TRUE)
  }

  # 3) Keep T. castaneum only NC_ entries (user rule)
  if (grepl("tribolium_castaneum", species_name) && !startsWith(seq_id, "NC_")) {
    return(NULL)
  }

  # If we found a clear chromosome-like label, standardize and return it.
  if (!is.na(chr_label) && nzchar(chr_label)) {
    chr_label <- toupper(gsub("\\s", "", chr_label))
    return(data.frame(chr = seq_id,
                      chrSimple = chr_label,
                      genome = species_name,
                      stringsAsFactors = FALSE))
  }

  # Otherwise signal "no label" -> let the caller decide how to number contigs later.
  return(data.frame(chr = seq_id,
                    chrSimple = NA_character_,
                    genome = species_name,
                    stringsAsFactors = FALSE))
}

all_chr_maps <- list()

for (fna in fna_files) {
  cat("📂 Processing:", fna, "\n")

  # Read only headers to avoid loading entire genome into memory
  con <- file(fna, open = "r")
  on.exit(close(con), add = TRUE)
  hdrs <- character()
  repeat {
    ln <- readLines(con, n = 10000, warn = FALSE)
    if (length(ln) == 0) break
    hdrs <- c(hdrs, ln[startsWith(ln, ">")])
  }

  #species <- tolower(gsub("_genomic\\.fna$", "", basename(fna)))
  fna_real <- normalizePath(fna, mustWork = TRUE)           # resolves symlinks
  base <- tolower(sub("_genomic\\.fna$", "", basename(fna_real)))
# If you *still* want to be robust to family-in-name, normalize it away:
  base <- sub("^([a-z]+_[a-z]+)_[a-z]+_(gc[af]_[0-9]+\\.[0-9]+)$", "\\1_\\2", base, perl = TRUE)

  species <- base
  species_scaffold_map <- scaffold_map_list[[species]]

  per_header <- lapply(hdrs, extract_chr_info,
                       species_name = species,
                       species_scaffold_map = species_scaffold_map)
  per_species_df <- do.call(rbind, per_header)

  if (is.null(per_species_df) || nrow(per_species_df) == 0) {
    cat("⚠️  No usable entries found in", fna, "\n")
    next
  }

  # Split into already-labeled vs unlabeled (contigs)
  labeled    <- per_species_df[!is.na(per_species_df$chrSimple) & nzchar(per_species_df$chrSimple), , drop = FALSE]
  unlabeled  <- per_species_df[is.na(per_species_df$chrSimple) | !nzchar(per_species_df$chrSimple), , drop = FALSE]

  # For unlabeled, assign chrSimple = "1","2","3",... (optionally with prefix)
  if (nrow(unlabeled) > 0) {
    ids <- unlabeled$chr
    ord <- order_by_numeric_suffix(ids)
    ids_ord <- ids[ord]
    assigned_numbers <- seq_along(ids_ord)
    unlabeled_numbered <- data.frame(
      chr       = ids_ord,
      chrSimple = paste0(CONTIG_PREFIX, assigned_numbers),
      genome    = species,
      stringsAsFactors = FALSE
    )
    # keep labeled rows as-is + add newly numbered contigs
    combined <- rbind(labeled, unlabeled_numbered)
  } else {
    combined <- labeled
  }

  # De-duplicate within species
  combined <- combined[!duplicated(combined[, c("chr", "genome")]), ]
  rownames(combined) <- NULL

  if (nrow(combined) > 0) {
    all_chr_maps[[length(all_chr_maps) + 1]] <- combined
  } else {
    cat("⚠️  No usable entries after processing", fna, "\n")
  }
}

# Combine and write
if (length(all_chr_maps) == 0) {
  writeLines("chr,chrSimple,genome", con = output_csv)
  cat("⚠️  No usable entries across all genomes — wrote empty file:", output_csv, "\n")
} else {
  final_df <- do.call(rbind, all_chr_maps)
  final_df <- unique(final_df)
  # For safety: ensure chrSimple is character
  final_df$chrSimple <- as.character(final_df$chrSimple)
  write.csv(final_df, file = output_csv, row.names = FALSE, quote = FALSE)
  cat("✅ Combined output written to:", output_csv, "\n")
}

