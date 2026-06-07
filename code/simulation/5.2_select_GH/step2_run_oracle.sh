#!/bin/bash
# =============================================================================
# step2_run_oracle.sh  (Section 5.2 ORACLE experiments)
# -----------------------------------------------------------------------------
# For each (case, additive, block, seed):
#     Rscript simulation_or.R  case additive seed block
# Fits general + special TGNQ at the TRUE memberships (G0, H0) and saves
# .rda under output/<block|power>_select_or_<add|mul>/.
# Needed by step3_output_tex.R for the "Oracle" rows.
# =============================================================================

# ===================== Auto-locate project root =====================
# Script lives in: <DIR>/code/simulation/5.2_select_GH/
# so go up three levels to reach <DIR>.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DIR="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

R_SCRIPT="${DIR}/code/simulation/5.2_select_GH/simulation_or.R"
LOG_DIR="${DIR}/log"

MAX_JOBS=128   # maximum number of concurrent R jobs

mkdir -p "${LOG_DIR}"
cd "${DIR}"

running=0    # number of launched-but-not-yet-reaped jobs

# ===================== Parameter loops =====================
#       For the full run use e.g.  for seed in $(seq 1 500);  and  for case in 1 2 3 4 5;
for seed in $(seq 1 500); do
  for case in 1 2 3 4 5; do
    for additive in 0 1; do
      for block in 0 1; do

        # -------- Submit one oracle replication --------
        nohup Rscript "${R_SCRIPT}" \
          "${case}" "${additive}" "${seed}" "${block}" \
          > "${LOG_DIR}/selectGH_or_case${case}_add${additive}_block${block}_seed${seed}.log" 2>&1 &

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