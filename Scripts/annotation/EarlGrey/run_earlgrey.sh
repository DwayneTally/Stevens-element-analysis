#!/bin/bash
set -eo pipefail

: "${GENOME_LIST:?Need GENOME_LIST env var (text file w/ full paths to *_genomic.fna)}"
OUTBASE="earlgrey_runs"

module load conda
conda activate earlgrey
module load blast

mkdir -p "$OUTBASE" logs

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
PREFIX="${BASENAME%_genomic.fna}"
RUNDIR="$OUTBASE/$PREFIX"

#Genus species
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
tar -C "$OUTBASE" -czf "$OUTBASE/${PREFIX}.tar.gz" "$PREFIX"

echo "[info] Created tarball: $OUTBASE/${PREFIX}.tar.gz"

#remove uncompressed directory to save file space
 rm -rf "$RUNDIR"

echo -e "task_id\tend\ttarball" >> "$OUTBASE/run_manifest.tsv"
echo -e "${SLURM_ARRAY_TASK_ID}\t$(date -Is)\t${OUTBASE}/${PREFIX}.tar.gz" >> "$OUTBASE/run_manifest.tsv"
