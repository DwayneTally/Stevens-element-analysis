#!/bin/bash
set -euo pipefail

ASTRAL_TREE="Astral/concord.tre"
ALIGNMENT_DIR="trimmed_alignments_4_original/"
GENE_TREES="ALL.tree"

THREADS=64
BOOTSTRAPS=1000
SCFL_REPS=100

echo "Step 1: ML bootstrap + sCFL on fixed ASTRAL topology"

iqtree3 \
  -te "$ASTRAL_TREE" \
  -p "$ALIGNMENT_DIR" \
  -b "$BOOTSTRAPS" \
  --scfl "$SCFL_REPS" \
  --prefix concord3_boot_scfl \
  -T "$THREADS"

echo "Step 2: Add gene concordance factors"

iqtree3 \
  -te concord3_boot_scfl.cf.tree \
  --gcf "$GENE_TREES" \
  --prefix busco_phylo_boot_gcf_scfl \
  -T "$THREADS"

echo "Finished"

