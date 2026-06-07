#!/bin/bash
# =============================================================================
# step1_run_simulation.sh  (Section 5.2 candidate (G,H) experiments)
# -----------------------------------------------------------------------------
# For each (case, additive, block, seed):
#     Rscript simulation.R  case additive seed block
# Each call runs ONE replication that fits seven candidate (G,H) models
# and writes .rda files under output/res_selectGH_<block|power>_<add|mul>/.
#
# Concurrency capped at MAX_JOBS parallel R processes; per-job logs at
#   ${LOG_DIR}/selectGH_case*_add*_block*_seed*.log
# =============================================================================

# ===================== Auto-locate project root =====================
# Script lives in: <DIR>/code/simulation/5.2_select_GH/
# so go up three levels to reach <DIR>.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DIR="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

R_SCRIPT="${DIR}/code/simulation/5.2_select_GH/simulation.R"
LOG_DIR="${DIR}/log"

MAX_JOBS=128   # maximum number of concurrent R jobs

mkdir -p "${LOG_DIR}"
cd "${DIR}"

running=0      # number of launched-but-not-yet-reaped jobs

# ===================== Parameter loops =====================
# seed     : Monte Carlo replication index, 1..500
# case     : design index 1..5
# additive : 1 = additive DGP, 0 = multiplicative / general DGP
# block    : 1 = block (SBM) network, 0 = power-law network
for seed in $(seq 1 500); do
  for case in 1 2 3 4 5; do
    for additive in 0 1; do
      for block in 0 1; do

        # -------- Submit one Monte Carlo replication --------
        # Arguments to simulation.R:
        #   $1 = case      (design index)
        #   $2 = additive  (1 = additive DGP)
        #   $3 = seed      (replication index)
        #   $4 = block     (1 = SBM, 0 = power-law)
        nohup Rscript "${R_SCRIPT}" \
          "${case}" "${additive}" "${seed}" "${block}" \
          > "${LOG_DIR}/selectGH_case${case}_add${additive}_block${block}_seed${seed}.log" 2>&1 &

        echo "Started: case=${case} additive=${additive} block=${block} seed=${seed}"

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