#!/bin/bash
# =============================================================================
# step1_run_simulation.sh  (Section 5.4: TGNQ vs QGNAR/gqnar)
# -----------------------------------------------------------------------------
# Every (case, additive, block, seed) is its own shell job
# (8000 total = 4 cases x 2 DGPs x 2 networks x 500 seeds). Each job is light,
# so MAX_JOBS is kept reasonably large.
#
# For each combination:
#     Rscript simulation.R  case additive seed block
#
# Logs: ${LOG_DIR}/gqnar_case<case>_add<add>_seed<seed>_block<block>.log
# =============================================================================

# ===================== Auto-locate project root =====================
# Script lives in: <DIR>/code/simulation/5.4_compare_with_gqnar/
# so go up three levels to reach <DIR>.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DIR="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

R_SCRIPT="${DIR}/code/simulation/5.4_compare_with_gqnar/simulation.R"
LOG_DIR="${DIR}/log"

MAX_JOBS=128

mkdir -p "${LOG_DIR}"
cd "${DIR}"

running=0    # number of launched-but-not-yet-reaped jobs

# ===================== Parameter loops =====================
# seed     : Monte Carlo replicate (1..500)
# case     : design index (1..4); see README for (N, T, Nblock)
# additive : 1 = additive DGP, 0 = multiplicative DGP
# block    : 1 = SBM network, 0 = power-law network
for seed in $(seq 1 500); do
  for case in 1 2 3 4; do
    for additive in 0 1; do
      for block in 0 1; do

        # -------- Submit one (case, additive, seed, block) job --------
        # Arguments to simulation.R:
        #   $1 = case     (design index)
        #   $2 = additive (1 = additive DGP)
        #   $3 = seed     (Monte Carlo replicate)
        #   $4 = block    (1 = SBM network)
        nohup Rscript "${R_SCRIPT}" \
          "${case}" "${additive}" "${seed}" "${block}" \
          > "${LOG_DIR}/gqnar_case${case}_add${additive}_seed${seed}_block${block}.log" 2>&1 &

        echo "Started: case=${case} additive=${additive} seed=${seed} block=${block}"

        # -------- Throttle: when MAX_JOBS launched, wait for the batch --------
        running=$((running + 1))
        if [ "$running" -ge "$MAX_JOBS" ]; then
          wait
          running=0
        fi

      done
    done
  done
done

wait
echo "All jobs done."