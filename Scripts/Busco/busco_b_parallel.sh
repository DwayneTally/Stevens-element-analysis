#!/bin/bash

# Load required modules
module load java
module load gnu
module load mvapich
module load hmmer
module load busco

INPUT_DIR="genomic_fna"  # directory containing .fna files from ncbi
CORES_PER_JOB=10   
MAX_JOBS=10        # max concurrent jobs

# Function to count current background jobs
function wait_for_jobs {
  while (( $(jobs -rp | wc -l) >= MAX_JOBS )); do
    sleep 5
  done
}

# Loop over all .fna files and run BUSCO
for fna_file in "$INPUT_DIR"/*.fna; do
  # Get base name for output
  base_name=$(basename "$fna_file" .fna)
  out_dir="./busco_${base_name}"

  echo "Starting BUSCO on: $fna_file"

  # Run BUSCO in background
  busco -i "$fna_file" \
        -l endopterygota_odb10 \
        -o "$out_dir" \
        -m genome \
        -c "$CORES_PER_JOB" \
        -f &

  # Wait if too many background jobs
  wait_for_jobs
done

# Wait for any remaining background jobs to finish
wait

echo "All BUSCO jobs completed."

