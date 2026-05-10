library(data.table)

# ---- load ----
dt <- fread("Tribolium_castaneum_GCF_000002335.3_phasedBlks.csv")

# ---- remove self comparisons ----
dt <- dt[genome1 != genome2]

# ---- block lengths ----
dt[, len1 := endBp1 - startBp1 + 1]
dt[, len2 := endBp2 - startBp2 + 1]

# =========================================================
# 1) SYNTENIC BLOCK NUMBER
# =========================================================

pair_blocks <- dt[, .(
  n_blocks = .N,
  unique_blkIDs = uniqueN(blkID)
), by = .(genome1, genome2)]

# =========================================================
# 2) GENOME AVERAGE BASES COVERED
# =========================================================

pair_cov <- dt[, .(
  bp_cov_genome1 = sum(len1),
  bp_cov_genome2 = sum(len2)
), by = .(genome1, genome2)]

pair_cov[, genome_avg_bp :=
           (bp_cov_genome1 + bp_cov_genome2) / 2]

# =========================================================
# 3) % blkIDs duplicated
# =========================================================

dup_stats <- dt[, {
  id_counts <- table(blkID)

  n_ids <- length(id_counts)
  n_dup_ids <- sum(id_counts > 1)

  .(
    pct_blkIDs_duplicated =
      n_dup_ids / n_ids * 100
  )
}, by = .(genome1, genome2)]

# =========================================================
# MERGE ALL STATS
# =========================================================

pairwise_stats <- Reduce(
  function(x,y) merge(x,y,by=c("genome1","genome2")),
  list(pair_blocks, pair_cov, dup_stats)
)

# =========================================================
# SUPERFAMILY AVERAGES
# =========================================================

superfamily_summary <- pairwise_stats[, .(

  # block counts
  mean_blocks = mean(n_blocks),
  median_blocks = median(n_blocks),

  # coverage
  mean_genome_avg_bp = mean(genome_avg_bp),
  median_genome_avg_bp = median(genome_avg_bp),

  # duplication
  mean_pct_blkID_dup = mean(pct_blkIDs_duplicated),
  median_pct_blkID_dup = median(pct_blkIDs_duplicated),

  n_pairs = .N
)]

print(superfamily_summary)

# ---- outputs ----
fwrite(pairwise_stats,
       "pairwise_block_stats.tsv",
       sep="\t")

fwrite(superfamily_summary,
       "superfamily_summary.tsv",
       sep="\t")
