# Code

This folder contains all R code used in the paper

"Two-way Homogeneity Pursuit for Quantile Network Vector Autoregression".

The code is organized into four subfolders, **each with its own dedicated
`README.md`**. This top-level README only describes the overall layout and how
the pieces fit together; for per-experiment design grids, argument orders, and
`.rda` schemas, see the README inside each subfolder.

All scripts assume the project root is the `TGNQ/` directory.

---

## 1. Folder structure

```
code/
├── generator/      # Data-generating processes (DGPs) for the simulations
│   ├── generate_data.R          # baseline / multiplicative TGNQ DGP
│   ├── additive_parameter.R     # additive parameter overrides
│   └── README.md
│
├── simulation/     # Monte Carlo experiments for Sections 5.1–5.4
│   ├── 5.1_and_5.3_Simulation/  # Tables 1 & 2 (estimation/inference + misspec.)
│   ├── 5.2_select_GH/           # Table 3 (QIC selection of G, H)
│   ├── 5.4_compare_with_gqnar/  # Table 4 (TGNQ vs one-way QGNAR)
│   └── README.md
│
├── real_data/      # Empirical application to Sina Weibo (Section 6)
│   ├── real_tgnq.R              # core TGNQ fitter (one (G,H) at a time)
│   ├── estimator_gqnar.R        # QGNAR benchmark estimator
│   ├── step1_selectGH.sh ... step9_inference.R
│   └── README.md
│
└── utils/          # Auxiliary helpers used across simulations
    ├── assess.R                 # RMSE, CI coverage, clustering error,
    │                            # and LaTeX-table reshaping utilities
    └── README.md
```

---

## 2. Mapping to the paper

| Paper section | Code location |
|---|---|
| Section 5.1 (estimation & inference) | `simulation/5.1_and_5.3_Simulation/` |
| Section 5.2 (selection of `G`, `H`)  | `simulation/5.2_select_GH/` |
| Section 5.3 (model misspecification) | `simulation/5.1_and_5.3_Simulation/` |
| Section 5.4 (vs. QGNAR)              | `simulation/5.4_compare_with_gqnar/` |
| Section 6 (Sina Weibo application)   | `real_data/` |

> **Note.** Sections 5.1 and 5.3 are produced by a **single run** inside
> `5.1_and_5.3_Simulation/` (every replication fits all estimators at once), so
> the `simulation/` folder has only **three** subfolders rather than four.

---

## 3. What each subfolder does

- **`generator/`** — defines the DGPs (network generators, group-parameter
  functions, and the panel simulator) that every Section-5 simulation sources.
  Not meant to be run directly.
- **`simulation/`** — the Monte Carlo suite. Every experiment follows the same
  two-stage pattern: a Bash driver loops over designs/seeds and calls
  `simulation.R` → `simulator.R`, then a `step*_output_*_tex.R` script
  aggregates the saved `.rda` files into the LaTeX tables.
- **`real_data/`** — the nine-step Weibo pipeline, from `(G, H)` selection
  (QIC + rolling out-of-sample) through the QGNAR comparison to final inference
  and figures.
- **`utils/`** — shared evaluation helpers (`assess.R`) called by the simulation
  post-processing scripts.

Refer to each subfolder's README for the exact scripts, command-line arguments,
case→`(N, T, Nblock)` mappings, and output naming.

---

## 4. R package dependencies

| Purpose | Packages |
|---|---|
| Quantile estimation                  | `twmq` (proposed estimator), `quantreg` |
| Discrete optimization (group search) | `CEoptim` |
| Numerics                             | `MASS`, `abind`, `compiler` |
| Parallelism                          | `parallel`, `foreach`, `doParallel`, `doSNOW` |
| Plotting                             | `ggplot2`, `reshape2`, `gridExtra`, `patchwork` |
| Tables / paths                       | `xtable`, `kableExtra`, `dplyr`, `here` |

The TGNQ estimator itself is provided by the package **`twmq`**, which must be
installed separately (it implements `twmq.estimate.auto.parallel`,
`update_NARG_twmq_parallel`, `Refine_G_parallel`, `Refine_H_parallel`,
`twmq.estimate_thetaGH.member.iterate`, `twmq_ci`, and their serial variants).

A multi-core machine is recommended for both the simulations and the empirical
application; scripts use `numCores = 16` or `20` by default.

---

## 5. Suggested reproduction workflow

From the project root `TGNQ/`:

```bash
# ---- Section 5: simulations ------------------------------------------------
# Each subfolder works the same way: launch the batch, then post-process.
bash    code/simulation/5.1_and_5.3_Simulation/step1_run_simulation.sh
Rscript code/simulation/5.1_and_5.3_Simulation/step2_output_tex.R

bash    code/simulation/5.2_select_GH/step1_run_simulation.sh
bash    code/simulation/5.2_select_GH/step2_run_oracle.sh
Rscript code/simulation/5.2_select_GH/step3_output_tex.R

bash    code/simulation/5.4_compare_with_gqnar/step1_run_simulation.sh
Rscript code/simulation/5.4_compare_with_gqnar/step2_output_tex.R

# ---- Section 6: empirical application --------------------------------------
bash    code/real_data/step1_selectGH.sh
Rscript code/real_data/step2_QIC_plot.R
Rscript code/real_data/step3_cf_selectGH.R
Rscript code/real_data/step4_fit_tgnq_alldata.R 77 4 3 0.01
Rscript code/real_data/step5_fit_qnarg.R
Rscript code/real_data/step6_predict_qnarg.R
Rscript code/real_data/step7_fit_compare.R
Rscript code/real_data/step8_pred_plot.R
Rscript code/real_data/step9_inference.R
```

The simulation batch scripts and Step 1 of the empirical pipeline submit many
parallel `Rscript` jobs; on smaller machines, lower `MAX_JOBS` (set at the top
of each driver) or split them into batches.

---

## 6. Output organization

All scripts write to the project-level `output/` directory, with per-job logs in
`log/`:

```
output/
├── block/  power/                         # 5.1 & 5.3 per-seed snapshots
├── res_block/  res_power/                 # 5.1 & 5.3 aggregated results
├── res_selectGH_*/  *_select_or_*/        # 5.2 candidate (G,H) and oracle
├── res_gqnar/res_*/                       # 5.4 TGNQ vs QGNAR
├── real_res/                             # Section 6 fitted models & predictions
└── figs/                                 # Section 6 PDF figures
log/                                       # one log file per Monte Carlo / real-data job
```

File-naming conventions are documented in the per-subfolder READMEs.