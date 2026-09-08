#!/usr/bin/env bash

# ============================================================
# RUN MULTIPLE ELEMENT-SPECIFIC GENESPACE ANALYSES ON COLO
# ============================================================

module load conda
conda activate genespace_env
set -euo pipefail

BASE="genespace_by_element"

RSCRIPT="run_genespace_element.R"

MAX_JOBS=6

ELEMENTS=(
    A
    C
    E
    G
    H
    X
)

mkdir -p "${BASE}/logs"

echo "ELEMENT-SPECIFIC GENESPACE"
echo "Host: $(hostname)"
echo "Started: $(date)"
echo "Maximum concurrent runs: ${MAX_JOBS}"
echo

run_element() {

    element="$1"

    wd="${BASE}/${element}"
    log="${BASE}/logs/${element}_genespace.log"

    echo "Starting Stevens ${element}: $(date)"

    Rscript "$RSCRIPT" "$wd" \
        > "$log" 2>&1

    echo "Finished Stevens ${element}: $(date)"
}

export -f run_element
export BASE
export RSCRIPT

# Launch up to MAX_JOBS simultaneously.
#
# As soon as one finishes, xargs launches the next one.

printf "%s\n" "${ELEMENTS[@]}" \
    | xargs -P "$MAX_JOBS" -I {} \
        bash -c 'run_element "$@"' _ {}

echo "ALL ELEMENT RUNS FINISHED"
echo "Finished: $(date)"
