#!/bin/bash
# =============================================================================
# step1_selectGH.sh
#
# Fit the TGNQ model on a 5x5 grid of (G, H) with T_train = 70 by launching
# background Rscript jobs via nohup, throttled to MAX_JOBS at a time.
#
# Each job:
#   - calls real_tgnq.R with arguments  T_train  G  H  rq_lambda
#   - writes:   output/real_res/tgnq_70-<G>-<H>.rda
#   - logs to:  output/log/<G><H>.log
# =============================================================================

# ===================== Auto-locate project root =====================
# Script lives in: <DIR>/code/real_data/
# so go up two levels to reach <DIR>.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Common settings
T_TRAIN=70
RQ_LAMBDA=0.01
SCRIPT="${DIR}/code/real_data/real_tgnq.R"
LOGDIR="${DIR}/output/log"

MAX_JOBS=1     # maximum number of concurrent R jobs

# Create the log directory (-p: no error if it already exists)
mkdir -p "${LOGDIR}"
cd "${DIR}"

running=0      # number of launched-but-not-yet-reaped jobs

# Loop over the 5x5 grid of (G, H)
for G in 1 2 3 4 5; do
  for H in 1 2 3 4 5; do

    nohup Rscript "$SCRIPT" "$T_TRAIN" "$G" "$H" "$RQ_LAMBDA" \
      > "$LOGDIR/${G}${H}.log" 2>&1 &

    echo "Started: G=${G} H=${H}"

    # -------- Throttle: when MAX_JOBS launched, wait for the batch --------
    running=$((running + 1))
    if [ "$running" -ge "$MAX_JOBS" ]; then
      wait
      running=0
    fi

  done
done

wait   # wait for the last (partial) batch to finish
echo "All TGNQ (G,H) grid jobs finished."