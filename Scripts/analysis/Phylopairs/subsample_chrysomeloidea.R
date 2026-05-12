#!/usr/bin/env Rscript

#Randomly subsamples Chrysomeloidea to match Curculionoidea
args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "."
out_dir  <- if (length(args) >= 2) args[2] else "chrysomeloidea_subsampled"
seed     <- if (length(args) >= 3) as.integer(args[3]) else 42

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

ref_species <- "Gastrophysa_polygoni_GCA_963576655.1"

load_dat <- function(f) {
  dat <- read.csv(f, stringsAsFactors = FALSE)
  colnames(dat)[colnames(dat) == "genome1"] <- "sp1"
  colnames(dat)[colnames(dat) == "genome2"] <- "sp2"
  dat
}

#Get target species count from Curculionoidea 
curcu_a <- load_dat(file.path(data_dir,
                   "Curculionoidea_pairwis_block_stats_with_tree_distance.csv"))
curcu_b <- load_dat(file.path(data_dir,
                   "Curculionidea_O_rusci_distance.csv"))

target_a <- length(unique(c(curcu_a$sp1, curcu_a$sp2)))
target_b <- length(unique(c(curcu_b$sp1, curcu_b$sp2)))

cat("  Curculionoidea (a) species:", target_a, "| rows:", nrow(curcu_a), "\n")
cat("  Curculionoidea (b) species:", target_b, "| rows:", nrow(curcu_b), "\n")

target_species <- target_a

dat_a <- load_dat(file.path(data_dir,
                  "Chrysomeloidea_pairwis_block_stats_with_tree_distance.csv"))
dat_b <- load_dat(file.path(data_dir,
                  "Chrysomeliodea_G_Polygoni_distance.csv"))

cat("  Chrysomeloidea (a) rows:", nrow(dat_a),
    "| Species:", length(unique(c(dat_a$sp1, dat_a$sp2))), "\n")
cat("  Chrysomeloidea (b) rows:", nrow(dat_b),
    "| Species:", length(unique(c(dat_b$sp1, dat_b$sp2))), "\n")

all_sp_a <- unique(c(dat_a$sp1, dat_a$sp2))
all_sp_b <- unique(c(dat_b$sp1, dat_b$sp2))

if (!ref_species %in% all_sp_a)
  stop("Reference species not found in dataset (a): ", ref_species)
if (!ref_species %in% all_sp_b)
  stop("Reference species not found in dataset (b): ", ref_species)

if (length(all_sp_a) < target_species)
  stop("Not enough species in Chrysomeloidea (a) to reach target of ", target_species)

#Ensure ref species, randomly fill remaining slots
cat("\n-- Sampling shared species set (n =", target_species, ") --\n")
remaining      <- setdiff(all_sp_a, ref_species)
shared_species <- c(ref_species, sample(remaining, target_species - 1))

missing_from_b <- setdiff(shared_species, all_sp_b)
if (length(missing_from_b) > 0) {
  cat("\n  [warn] The following sampled species are absent from dataset (b):\n")
  cat(paste0("    ", missing_from_b, collapse = "\n"), "\n")
  cat("  These will not appear in the (b) subsample.\n")
}

dat_a_sub <- dat_a[dat_a$sp1 %in% shared_species & dat_a$sp2 %in% shared_species, ]
dat_b_sub <- dat_b[dat_b$sp1 %in% shared_species & dat_b$sp2 %in% shared_species, ]

cat("\n  Dataset (a) after subsample:", nrow(dat_a_sub), "rows |",
    length(unique(c(dat_a_sub$sp1, dat_a_sub$sp2))), "species\n")
cat("  Dataset (b) after subsample:", nrow(dat_b_sub), "rows |",
    length(unique(c(dat_b_sub$sp1, dat_b_sub$sp2))), "species\n")

#save output
out_a <- file.path(out_dir, "Chrysomeloidea_stevens_subsampled.csv")
out_b <- file.path(out_dir, "Chrysomeloidea_Gpolygoni_subsampled.csv")

write.csv(dat_a_sub, file = out_a, row.names = FALSE)
write.csv(dat_b_sub, file = out_b, row.names = FALSE)

cat("\n  Saved:", out_a, "\n")
cat("  Saved:", out_b, "\n")

#saved shared species list for reference
species_out <- file.path(out_dir, "Chrysomeloidea_subsampled_species.txt")
writeLines(sort(shared_species), species_out)
cat("  Species list saved:", species_out, "\n")

cat("\nDone.\n")
