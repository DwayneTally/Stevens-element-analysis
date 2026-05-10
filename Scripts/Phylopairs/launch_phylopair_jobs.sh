#!/bin/bash
# ============================================================
# launch_jobs.sh
# Launches one background R job per superfamily on colo38.
# Each job uses 4 cores (Stan chains). With 10 superfamilies
# this requests up to 40 cores simultaneously — adjust
# MAX_PARALLEL below if the node has fewer cores available.
#
# Usage:
#   bash launch_jobs.sh <data_dir> <tree_file>
#
# Example:
#   bash launch_jobs.sh /path/to/data full_tree_superfamily_calibrated.newick
# ============================================================

DATA_DIR=${1:?  "ERROR: provide data directory as argument 1"}
TREE_FILE=${2:? "ERROR: provide tree file path as argument 2"}
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
R_SCRIPT="${SCRIPT_DIR}/run_phylopairs_superfamily_updated_beta_final_product.R"

# Maximum number of superfamily jobs running at the same time.
# Each job uses 4 cores, so MAX_PARALLEL=8 uses up to 32 cores.
MAX_PARALLEL=10

TOTAL_SF=10   # number of entries in file_sets
LOG_DIR="${SCRIPT_DIR}/logs/cai"
mkdir -p "${LOG_DIR}"

echo "============================================"
echo " Synteny phylopairs — parallel job launcher"
echo " Data dir  : ${DATA_DIR}"
echo " Tree file : ${TREE_FILE}"
echo " Script    : ${R_SCRIPT}"
echo " Max parallel jobs: ${MAX_PARALLEL}"
echo " Log dir   : ${LOG_DIR}"
echo "============================================"
echo ""

running=0

for i in $(seq 1 ${TOTAL_SF}); do
    # Wait if we have hit the parallel job limit
    while [ "${running}" -ge "${MAX_PARALLEL}" ]; do
        sleep 30
        # Recount still-running jobs
        running=$(jobs -r | wc -l)
    done

    LOG_FILE="${LOG_DIR}/sf_${i}.log"
    echo "Launching superfamily index ${i} → log: ${LOG_FILE}"

    nohup Rscript "${R_SCRIPT}" "${i}" "${DATA_DIR}" "${TREE_FILE}" \
        > "${LOG_FILE}" 2>&1 &

    running=$(jobs -r | wc -l)
    sleep 2   # small stagger so Stan compilation doesn't collide
done

echo ""
echo "All ${TOTAL_SF} jobs submitted. Monitor progress with:"
echo "  tail -f ${LOG_DIR}/sf_<N>.log"
echo ""
echo "Check for completion with:"
echo "  grep 'Done:' ${LOG_DIR}/*.log"
echo ""
echo "Results will appear in: ${SCRIPT_DIR}/results/<SuperfamilyName>/"
