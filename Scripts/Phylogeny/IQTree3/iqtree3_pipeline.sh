#!/bin/bash
set -euo pipefail

ASTRAL_TREE="../../Data/concord.tre"
ALIGNMENT_DIR="trimmed_alignments_4_original/"
GENE_TREES="ALL.tree"

THREADS=64
BOOTSTRAPS=1000
SCFL_REPS=100

iqtree3 \
  -te "$ASTRAL_TREE" \
  -p "$ALIGNMENT_DIR" \
  -b "$BOOTSTRAPS" \
  --scfl "$SCFL_REPS" \
  --prefix concord_boot_scfl \
  -T "$THREADS"

iqtree3 \
  -te concord_boot_scfl.cf.tree \
  --gcf "$GENE_TREES" \
  --prefix beetle \
  -T "$THREADS"

echo "Finished"

