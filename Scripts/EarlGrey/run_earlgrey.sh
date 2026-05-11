#!/bin/bash
#SBATCH --job-name=earlgrey_list
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dwtally@iu.edu
#SBATCH --time=96:00:00
#SBATCH --partition=general
#SBATCH -A r00259
#SBATCH --array=1-${N}

set -eo pipefail

: "${GENOME_LIST:?Need GENOME_LIST env var (text file w/ full paths to *_genomic.fna)}"

# -----------------------------
# OUTBASE now fixed to scratch
# -----------------------------
OUTBASE="/N/scratch/dwtally/earlgrey_runs"

module load conda
conda activate earlgrey
module load blast

mkdir -p "$OUTBASE" logs

# 1-indexed array
GENOME="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$GENOME_LIST" | sed 's/\r$//')"
if [[ -z "${GENOME:-}" ]]; then
  echo "No genome found for task ${SLURM_ARRAY_TASK_ID}" >&2
  exit 1
fi
if [[ ! -f "$GENOME" ]]; then
  echo "Genome file not found: $GENOME" >&2
  exit 1
fi

BASENAME="$(basename "$GENOME")"
PREFIX="${BASENAME%_genomic.fna}"   # e.g. Genus_species_GCA_xxx
RUNDIR="$OUTBASE/$PREFIX"

# Genus_species for -s
SPECIES="$(echo "$BASENAME" | cut -d'_' -f1,2)"

mkdir -p "$RUNDIR"

echo -e "task_id\tgenome\toutdir\tstart" >> "$OUTBASE/run_manifest.tsv"
echo -e "${SLURM_ARRAY_TASK_ID}\t${GENOME}\t${RUNDIR}\t$(date -Is)" >> "$OUTBASE/run_manifest.tsv"

cd "$RUNDIR"

earlGrey \
  -g "$GENOME" \
  -t 32 \
  -o "$RUNDIR" \
  -s "$SPECIES" \
  -d yes \
  -e yes \
  > "$RUNDIR/${PREFIX}.earlGrey.stdout.log" 2>&1

echo "[info] earlGrey finished for $PREFIX"

# -----------------------------
# Tar results after completion
# -----------------------------
tar -C "$OUTBASE" -czf "$OUTBASE/${PREFIX}.tar.gz" "$PREFIX"

echo "[info] Created tarball: $OUTBASE/${PREFIX}.tar.gz"

# OPTIONAL: remove uncompressed directory to save scratch space
 rm -rf "$RUNDIR"

echo -e "task_id\tend\ttarball" >> "$OUTBASE/run_manifest.tsv"
echo -e "${SLURM_ARRAY_TASK_ID}\t$(date -Is)\t${OUTBASE}/${PREFIX}.tar.gz" >> "$OUTBASE/run_manifest.tsv"
