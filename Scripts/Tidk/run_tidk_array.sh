#!/bin/bash
#SBATCH --job-name=tidk_search
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=60G
#SBATCH --time=72:00:00
#SBATCH --partition=general
#SBATCH --array=0-191%10
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH -A r00259

set -euo pipefail

module load conda
conda activate earlgrey   # or whatever env has tidk

FNA_LIST="fna_list.txt"
OUTDIR="/N/scratch/dwtally/tidk_search_results"
mkdir -p "$OUTDIR"

# ---- pick genome ----
FNA=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$FNA_LIST")

BASENAME=$(basename "$FNA" .fna)

echo "[$(date)] Running TIDK on $BASENAME"

tidk search \
  -f "$FNA" \
  -o "$OUTDIR/${BASENAME}_tidk" \
  -t $SLURM_ntasks_per_node 

echo "[$(date)] Finished $BASENAME"

