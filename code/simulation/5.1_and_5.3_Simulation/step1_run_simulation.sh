#!/bin/bash
# =============================================================================
# step1_run_simulation.sh  (Sections 5.1 and 5.3 simulation experiments)
# -----------------------------------------------------------------------------
# A single run of this driver produces the Monte Carlo results for BOTH
# Section 5.1 (Estimation and inference) and Section 5.3 (Misspecified models):
# every replication fits the oracle, general, special and misspecified
# estimators at once.
#
#   Pass 1: G0 = 2, case in {1,2,3,4}
#   Pass 2: G0 = 3, case in {2,3,4,5}
#
# For each (case, additive, block): Rscript simulation.R G0 case additive block ncores
# Each R job iterates seeds 1..500 in parallel via doParallel (NCORES cores),
# so the shell runs ONE heavy R job at a time (MAX_JOBS=1).
# =============================================================================

# ===================== Auto-locate project root =====================
# Script lives in: <DIR>/code/simulation/5.1_and_5.3_Simulation/
# so go up three levels to reach <DIR>.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DIR="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

R_SCRIPT="${DIR}/code/simulation/5.1_and_5.3_Simulation/simulation.R"
LOG_DIR="${DIR}/log"

MAX_JOBS=1   # one heavy R job at a time (parallelism is inside R)
NCORES=128     # cores used by doParallel inside each R job

mkdir -p "${LOG_DIR}"
cd "${DIR}"

running=0    # number of launched-but-not-yet-reaped jobs

# =============================================================================
# Pass 1: G0 = 2, case in {1,2,3,4}
# =============================================================================
for case in 1 2 3 4; do
  for additive in 0 1; do
    for block in 0 1; do

      # Arguments: G0 case additive block ncores
      nohup Rscript "${R_SCRIPT}" \
        2 "${case}" "${additive}" "${block}" "${NCORES}" \
        > "${LOG_DIR}/misspec_G0_2_case${case}_add${additive}_block${block}.log" 2>&1 &

      echo "Started: G0=2 case=${case} additive=${additive} block=${block} ncores=${NCORES}"

      running=$((running + 1))
      if [ "$running" -ge "$MAX_JOBS" ]; then
        wait          # wait for the current batch to finish
        running=0
      fi

    done
  done
done

# =============================================================================
# Pass 2: G0 = 3, case in {2,3,4,5}
# =============================================================================
for case in 2 3 4 5; do
  for additive in 0 1; do
    for block in 0 1; do

      nohup Rscript "${R_SCRIPT}" \
        3 "${case}" "${additive}" "${block}" "${NCORES}" \
        > "${LOG_DIR}/misspec_G0_3_case${case}_add${additive}_block${block}.log" 2>&1 &

      echo "Started: G0=3 case=${case} additive=${additive} block=${block} ncores=${NCORES}"

      running=$((running + 1))
      if [ "$running" -ge "$MAX_JOBS" ]; then
        wait
        running=0
      fi

    done
  done
done

wait
echo "All jobs done."