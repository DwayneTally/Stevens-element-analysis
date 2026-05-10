#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

# =========================================================
# ARGUMENTS
# =========================================================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("
Usage:
  Rscript parse_blk_with_genomes_unionCoverage.R phasedBlks.csv genome_stats_by_species.csv
")
}

blk_csv <- args[1]
genome_stats_csv <- args[2]

cat("Using blkfile: ", blk_csv, "\n")
cat("Using genome stats: ", genome_stats_csv, "\n\n")

# =========================================================
# LOAD BLKFILE
# =========================================================
dt <- fread(blk_csv)

# remove self comparisons
dt <- dt[genome1 != genome2]

# normalize coordinates
dt[, `:=`(
  s1 = pmin(startBp1, endBp1),
  e1 = pmax(startBp1, endBp1),
  s2 = pmin(startBp2, endBp2),
  e2 = pmax(startBp2, endBp2)
)]

# =========================================================
# CANONICALIZE PAIRS
# =========================================================
dt[, genome1_c := pmin(genome1, genome2)]
dt[, genome2_c := pmax(genome1, genome2)]

dt[, `:=`(
  chr_genome1 = ifelse(genome1 == genome1_c, chr1, chr2),
  chr_genome2 = ifelse(genome1 == genome1_c, chr2, chr1),

  s_genome1 = ifelse(genome1 == genome1_c, s1, s2),
  e_genome1 = ifelse(genome1 == genome1_c, e1, e2),

  s_genome2 = ifelse(genome1 == genome1_c, s2, s1),
  e_genome2 = ifelse(genome1 == genome1_c, e2, e1)
)]

dt[, `:=`(genome1 = genome1_c, genome2 = genome2_c)]
dt[, c("genome1_c","genome2_c") := NULL]

# =========================================================
# INTERVAL UNION FUNCTION
# =========================================================
interval_union_bp <- function(start, end) {

  if (length(start) == 0) return(0)

  o <- order(start, end)
  start <- start[o]
  end   <- end[o]

  total <- 0
  cur_s <- start[1]
  cur_e <- end[1]

  if (length(start) >= 2) {
    for (i in 2:length(start)) {
      if (start[i] <= cur_e) {
        cur_e <- max(cur_e, end[i])
      } else {
        total <- total + (cur_e - cur_s + 1)
        cur_s <- start[i]
        cur_e <- end[i]
      }
    }
  }

  total + (cur_e - cur_s + 1)
}

# =========================================================
# UNION COVERAGE (chromosome aware)
# =========================================================
cov1_chr <- dt[, .(
  bp = interval_union_bp(s_genome1, e_genome1)
), by = .(genome2, genome1, chr_genome1)]

cov2_chr <- dt[, .(
  bp = interval_union_bp(s_genome2, e_genome2)
), by = .(genome2, genome1, chr_genome2)]

pair_cov <- merge(
  cov1_chr[, .(bp_cov_genome1 = sum(bp)), by = .(genome2, genome1)],
  cov2_chr[, .(bp_cov_genome2 = sum(bp)), by = .(genome2, genome1)],
  by = c("genome2","genome1"),
  all = TRUE
)

pair_cov[is.na(bp_cov_genome1), bp_cov_genome1 := 0]
pair_cov[is.na(bp_cov_genome2), bp_cov_genome2 := 0]

pair_cov[, genome_avg_bp :=
           (bp_cov_genome1 + bp_cov_genome2) / 2]

# =========================================================
# BLOCK COUNTS
# =========================================================
pair_blocks <- dt[, .(
  n_blocks = .N,
  unique_blkIDs = uniqueN(blkID)
), by = .(genome2, genome1)]

# =========================================================
# MEDIAN BLOCK LENGTH (blkID-aware, chr-aware, overlap-safe)
#   1) For each (pair, blkID, chr) compute union bp
#   2) Sum across chr to get per-(pair, blkID) length
#   3) Median across blkIDs within the pair
# =========================================================

# make sure numeric (prevents integer/double weirdness)
dt[, `:=`(
  s_genome1 = as.integer(s_genome1),
  e_genome1 = as.integer(e_genome1),
  s_genome2 = as.integer(s_genome2),
  e_genome2 = as.integer(e_genome2)
)]

# union bp per blkID per chr (genome1 side)
blk1_chr <- dt[, .(
  bp = interval_union_bp(s_genome1, e_genome1)
), by = .(genome2, genome1, blkID, chr_genome1)]

blk1 <- blk1_chr[, .(
  blk_len_genome1 = sum(bp, na.rm = TRUE)
), by = .(genome2, genome1, blkID)]

# union bp per blkID per chr (genome2 side)
blk2_chr <- dt[, .(
  bp = interval_union_bp(s_genome2, e_genome2)
), by = .(genome2, genome1, blkID, chr_genome2)]

blk2 <- blk2_chr[, .(
  blk_len_genome2 = sum(bp, na.rm = TRUE)
), by = .(genome2, genome1, blkID)]

# merge the two sides (some blkIDs might be missing on one side)
blk_len <- merge(blk1, blk2, by = c("genome2","genome1","blkID"), all = TRUE)
blk_len[is.na(blk_len_genome1), blk_len_genome1 := 0]
blk_len[is.na(blk_len_genome2), blk_len_genome2 := 0]

pair_medians <- blk_len[, .(
  median_block_len_genome1 = median(as.numeric(blk_len_genome1), na.rm = TRUE),
  median_block_len_genome2 = median(as.numeric(blk_len_genome2), na.rm = TRUE)
), by = .(genome2, genome1)]

pair_medians[, median_block_len_pairmean :=
               (median_block_len_genome1 + median_block_len_genome2) / 2]

# =========================================================
# % blkIDs duplicated
# =========================================================
dup_stats <- dt[, {
  id_counts <- table(blkID)
  n_ids <- length(id_counts)
  n_dup_ids <- sum(id_counts > 1)

  .(pct_blkIDs_duplicated =
      (n_dup_ids / n_ids) * 100)
}, by = .(genome2, genome1)]

# =========================================================
# MERGE STATS
# =========================================================
pairwise_stats <- Reduce(
  function(x,y)
    merge(x,y,
          by=c("genome2","genome1"),
          all=TRUE),
  list(pair_blocks,
       pair_cov,
       pair_medians,
       dup_stats)
)

# =========================================================
# GENOME SIZES
# =========================================================
gs <- fread(genome_stats_csv)

gs[, species :=
      sub("\\.gff3$", "", species_file)]

gs[, genome_size_bp :=
      genome_size_mb * 1e6]

gs1 <- gs[, .(
  genome1 = species,
  genome_size1_bp = genome_size_bp
)]

gs2 <- gs[, .(
  genome2 = species,
  genome_size2_bp = genome_size_bp
)]

pairwise_stats <- merge(pairwise_stats, gs1,
                        by="genome1", all.x=TRUE)

pairwise_stats <- merge(pairwise_stats, gs2,
                        by="genome2", all.x=TRUE)

pairwise_stats[, mean_genome_size_bp :=
                 (genome_size1_bp +
                  genome_size2_bp) / 2]

# =========================================================
# COVERAGE METRICS (FINAL)
# =========================================================
pairwise_stats[, pct_cov_genome1 :=
                 bp_cov_genome1 /
                 genome_size1_bp * 100]

pairwise_stats[, pct_cov_genome2 :=
                 bp_cov_genome2 /
                 genome_size2_bp * 100]

pairwise_stats[, pct_cov_pairmean :=
                 (pct_cov_genome1 +
                  pct_cov_genome2) / 2]

# =========================================================
# COLUMN ORDER
# =========================================================
target_cols <- c(
  "genome2","genome1",
  "n_blocks","unique_blkIDs",
  "bp_cov_genome1","bp_cov_genome2","genome_avg_bp",
  "median_block_len_genome1",
  "median_block_len_genome2",
  "median_block_len_pairmean",
  "pct_blkIDs_duplicated",
  "genome_size1_bp","genome_size2_bp",
  "mean_genome_size_bp",
  "pct_cov_genome1","pct_cov_genome2",
  "pct_cov_pairmean"
)

pairwise_stats <- pairwise_stats[, ..target_cols]

# =========================================================
# OUTPUT
# =========================================================
fwrite(pairwise_stats,
       "pairwise_block_stats.csv")

cat("\nWrote: pairwise_block_stats.csv\n")
