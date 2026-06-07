# TGNQ: Two-way Homogeneity Pursuit for Quantile Network Vector Autoregression

This repository contains all materials needed to **reproduce the results** in the paper

"Two-way Homogeneity Pursuit for Quantile Network Vector Autoregression".

It includes the Monte Carlo simulations of **Section 5**, the **Sina Weibo**
empirical application of **Section 6**, the data-generating code, and the
post-processing scripts that build the tables and figures of the paper.

> **The TGNQ estimator itself lives in a separate package repository.**
> To keep the package cleanly installable, the estimation package `twmq` has
> been split into its own GitHub repository and can be installed directly with
> ```r
> remotes::install_github("https://github.com/tgnq-repro/twmq-package")
> ```
> *This* repository contains only the reproduction code (simulations, real-data
> analysis, generators, and utilities). See §3 below.

---

## 1. Repository structure

```
TGNQ/                         # project root (all scripts assume this as cwd)
├── code/
│   ├── generator/            # data-generating processes (DGPs)
│   ├── simulation/           # Monte Carlo experiments (Sections 5.1–5.4)
│   │   ├── 5.1_and_5.3_Simulation/   # Tables 1 & 2 (single run)
│   │   ├── 5.2_select_GH/            # Table 3 (QIC selection of G, H)
│   │   └── 5.4_compare_with_gqnar/   # Table 4 (TGNQ vs QGNAR)
│   ├── real_data/            # Sina Weibo application (Section 6)
│   ├── utils/                # shared evaluation helpers (assess.R)
│   └── README.md             # overview of the code folder
│
├── data/
│   ├── weibo.rda             # Sina Weibo dataset (Section 6 input)
│   └── README.md             # data dictionary
│
├── output/                   # all generated results, tables, and figures
├── log/                      # per-job stdout/stderr from batch scripts
├── manuscript/               # paper sources
└── README.md                 # (this file)
```

Each subfolder under `code/` has its **own `README.md`** documenting the exact
scripts, command-line arguments, design grids, and `.rda` schemas. This
top-level README only covers global setup, requirements, and reproduction order.

---

## 2. Mapping to the paper

| Paper section | Code location | Output |
|---|---|---|
| Section 5.1 (estimation & inference) | `code/simulation/5.1_and_5.3_Simulation/` | Table 1 |
| Section 5.2 (selection of `G`, `H`)  | `code/simulation/5.2_select_GH/`          | Table 3 |
| Section 5.3 (model misspecification) | `code/simulation/5.1_and_5.3_Simulation/` | Table 2 |
| Section 5.4 (vs. QGNAR)              | `code/simulation/5.4_compare_with_gqnar/` | Table 4 |
| Section 6 (Sina Weibo application)   | `code/real_data/`                         | Figures + comparison tables |

Sections 5.1 and 5.3 are produced by a **single run** inside
`5.1_and_5.3_Simulation/` (every replication fits all estimators at once).

---

## 3. Software requirements

### 3.1 R and platform

- **R version:** ≥ 4.2.0 (tested on R 4.3.x).
- **Platform:** Linux x86-64 (tested on AMD EPYC 7763). A multi-core machine is
  strongly recommended (see §5).
- **Paths:** all scripts locate files **relative to the project root** via the
  [`here`](https://here.r-lib.org/) package, so no `setwd()` editing is needed.

### 3.2 The `twmq` package (proposed estimator)

The TGNQ estimator is provided by the standalone package **`twmq`**, hosted in a
separate repository and installed with:

```r
install.packages("remotes")
remotes::install_github("https://github.com/tgnq-repro/twmq-package")
```

It supplies `twmq.estimate.auto(.parallel)`, `update_NARG_twmq(.parallel)`,
`Refine_G_parallel`, `Refine_H_parallel`,
`twmq.estimate_thetaGH.member.iterate`, `twmq_ci`, and `twmq.label.switch`.
Version numbers and detailed requirements are listed in that repository's
own README.

### 3.3 CRAN dependencies

```r
install.packages(c(
  "here",                                  # relative paths
  "quantreg", "CEoptim", "MASS", "compiler",  # estimation / optimization
  "parallel", "foreach", "doParallel", "doSNOW", # parallelism
  "xtable", "kableExtra", "dplyr",         # LaTeX tables
  "ggplot2", "reshape2", "gridExtra", "patchwork" # figures
))
```

The specific versions used to produce the results in the paper are:

| Package | Version | Package | Version |
|---|---|---|---|
| `here`       | 1.0.1  | `doParallel` | 1.0.17 |
| `quantreg`   | 5.97   | `doSNOW`     | 1.0.20 |
| `CEoptim`    | 1.3    | `xtable`     | 1.8-4  |
| `MASS`       | 7.3-60 | `kableExtra` | 1.3.4  |
| `foreach`    | 1.5.2  | `dplyr`      | 1.1.4  |
| `ggplot2`    | 3.5.1  | `reshape2`   | 1.4.4  |
| `gridExtra`  | 2.3    | `patchwork`  | 1.2.0  |



### 3.4 Which packages each part needs

| Part | Packages |
|---|---|
| Simulations 5.1 / 5.3 | `twmq`, `CEoptim`, `doParallel`, `foreach`, `here`, `kableExtra`, `dplyr` |
| Simulation 5.2        | `twmq`, `CEoptim`, `here`, `xtable` |
| Simulation 5.4        | `twmq`, `CEoptim`, `quantreg`, `here`, `xtable` |
| Real data (Section 6) | `twmq`, `quantreg`, `CEoptim`, `MASS`, `compiler`, `parallel`, `foreach`, `doParallel`, `doSNOW`, `ggplot2`, `reshape2`, `gridExtra`, `patchwork`, `here` |

---

## 4. How to reproduce

From the project root `TGNQ/`:

```bash
# ---- Section 5: simulations ------------------------------------------------
# Each subfolder: launch the batch, then run the post-processing script.
bash    code/simulation/5.1_and_5.3_Simulation/step1_run_simulation.sh
Rscript code/simulation/5.1_and_5.3_Simulation/step2_output_tex.R   # Tables 1 & 2

bash    code/simulation/5.2_select_GH/step1_run_simulation.sh
bash    code/simulation/5.2_select_GH/step2_run_oracle.sh
Rscript code/simulation/5.2_select_GH/step3_output_tex.R            # Table 3

bash    code/simulation/5.4_compare_with_gqnar/step1_run_simulation.sh
Rscript code/simulation/5.4_compare_with_gqnar/step2_output_tex.R   # Table 4

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

Per-experiment details (argument order, `case → (N, T, Nblock)` grids, output
naming) are in the README inside each subfolder. The batch drivers cap
concurrency via `MAX_JOBS` (set at the top of each `.sh`); lower it on smaller
machines.

---

## 5. Expected runtime

Reproducing the **full set of simulations and the empirical application in the
main text** (excluding the appendix) takes roughly **two weeks of wall-clock
time** on **two servers**, each equipped with an **AMD EPYC 7763 64-Core
Processor (128 logical cores)**. Individual components are much faster; the cost
is dominated by the 500 Monte Carlo replications per design across all sections.
Runtime scales roughly linearly with the number of available cores.

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
├── real_res/                              # Section 6 fitted models & predictions
└── figs/                                  # Section 6 PDF figures
log/                                       # one log file per simulation / real-data job
```

File-naming conventions are documented in the per-subfolder READMEs.

---

## 7. Where to look for more detail

- **Code overview:** `code/README.md`
- **DGPs:** `code/generator/README.md`
- **Simulations:** `code/simulation/README.md` and each experiment subfolder
- **Real-data pipeline:** `code/real_data/README.md`
- **Evaluation helpers:** `code/utils/README.md`
- **Dataset dictionary:** `data/README.md`
- **The estimator package:** the standalone `twmq` repository.