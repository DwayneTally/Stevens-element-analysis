#!/bin/bash
#SBATCH --job-name=Helixer
#SBATCH --output=logs/Helixer_%A_%a.out
#SBATCH --error=logs/Helixer_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dwtally@iu.edu
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --partition=gpu
#SBATCH --mem=128G
#SBATCH -A r00259

set -euo pipefail

module load helixer

# File with one job per line:
# fasta_path<TAB>species_name<TAB>gff_output_path
JOB_LIST="helixer_jobs.txt"

# Get the line for this array task
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOB_LIST")

FASTA=$(echo "$LINE" | cut -f1)
SPECIES=$(echo "$LINE" | cut -f2)
GFF_OUT=$(echo "$LINE" | cut -f3)

echo "Running task ${SLURM_ARRAY_TASK_ID}"
echo "FASTA:   $FASTA"
echo "SPECIES: $SPECIES"
echo "GFF OUT: $GFF_OUT"

Helixer.py \
  --lineage invertebrate \
  --fasta-path "$FASTA" \
  --species "$SPECIES" \
  --gff-output-path "$GFF_OUT"
