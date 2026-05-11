#!/bin/bash

set -euo pipefail

module load conda

FNA_LIST="fna_list.txt"
OUTDIR="tidk_search_results"
mkdir -p "$OUTDIR"
FNA=$(sed -n "$((ARRAY_TASK_ID))p" "$FNA_LIST")

BASENAME=$(basename "$FNA" .fna)

echo "[$(date)] Running TIDK on $BASENAME"

tidk search \
  -f "$FNA" \
  -o "$OUTDIR/${BASENAME}_tidk" \
  -t $SLURM_ntasks_per_node 

echo "[$(date)] Finished $BASENAME"

