# Simulation code for Sections 5.1 and 5.3

This folder contains the Monte Carlo code for
"Two-way Homogeneity Pursuit for Quantile Network Vector Autoregression".

**A single run of the simulations in this folder reproduces both**
**Table 1 (Section 5.1, Estimation and inference) and Table 2 (Section 5.3,Misspecified models) of the main paper,**
**together with their supplementary counterparts.**
For every replication, `simulator.R` fits *all* estimators at once
— the oracle, the TGNQ-general, the TGNQ-special, and the TGNQ-misspecified estimators (each with and without the refinement step).
Because all estimators are produced in one pass, a single post-processing script
(`step2_output_tex.R`) builds the LaTeX for both tables:

* **Table 1 (Section 5.1):** oracle estimator and TGNQ-with-refinement results,
  reported as `RMSE (AEcp)`;
* **Table 2 (Section 5.3):** head-to-head General vs. Additive vs.
  Multiplicative (MIS) comparison, reported as RMSE only.

---

## 1. Pipeline overview

```
step1_run_simulation.sh   -- loops over (G0, case, additive, block) -->
        |
        v
simulation.R              -- parses args, sources DGP + simulator,
        |                    runs simulation(seed) in parallel (seeds 1..500)
        v
simulator.R::simulation() -- one replication: fits oracle / general /
        |                    special / misspecified estimators
        v
   ┌───────────────────────────────────────────────────────────────┐
   │ output/<block|power>/<...>-seed_<seed>.rda    (per-seed fits) │──┐
   │ output/res_<block|power>/res-G_<...>.rda      (error vectors) │──┤
   └───────────────────────────────────────────────────────────────┘  │
                                                                      v
                                                            step2_output_tex.R --> Table 1 (Section 5.1) + Table 2 (Section 5.3)
```

---

## 2. How to reproduce

From the project root (e.g. `TGNQ/`), run the following in order.

1. **Run all Monte Carlo simulations.**

   ```bash
   bash code/simulation/5.1_and_5.3_Simulation/step1_run_simulation.sh
   ```

   The shell driver runs **two passes** of the loop:

   * Pass 1: `G0 = 2`, `case ∈ {1, 2, 3, 4}`,
   * Pass 2: `G0 = 3`, `case ∈ {2, 3, 4, 5}`,

   each combined with `additive ∈ {0, 1}` and `block ∈ {0, 1}`. For every
   combination it calls

   ```
   Rscript simulation.R  G0  case  additive  block  ncores
   ```

   and writes a log to `log/`. Each R job iterates seeds `1..500` in parallel
   via `doParallel` (`NCORES` cores), so only one heavy R job is launched at a
   time (`MAX_JOBS = 1`).

2. **Print the Table 1 and Table 2 LaTeX.**

   In R:

   ```r
   source("code/simulation/5.1_and_5.3_Simulation/step2_output_tex.R")
   ```

   This single script aggregates the saved outputs and prints all tables to the
   console; paste them directly into the manuscript.

---

## 3. Command-line arguments and `case` mapping

`simulation.R` takes five arguments:

| arg        | meaning                                                        |
|------------|---------------------------------------------------------------|
| `G0`       | true number of row groups (`2` or `3`)                         |
| `case`     | design index (1–5), maps to `(N, T, Nblock)` below            |
| `additive` | `1` = additive DGP, `0` = multiplicative DGP                   |
| `block`    | `1` = SBM network, `0` = power-law network                    |
| `ncores`   | cores used by `doParallel` inside the R job                    |

Mapping `case` → `(N, T, Nblock)`:

| case | N   | T   | Nblock |
|------|-----|-----|--------|
| 1    | 100 | 50  | 5      |
| 2    | 100 | 100 | 5      |
| 3    | 100 | 200 | 5      |
| 4    | 200 | 200 | 10     |
| 5    | 200 | 400 | 10     |

---

## 4. Estimators fitted per replication

For each seed, `simulator.R::simulation()` generates the (fixed) network and the TGNQ data, then fits:

* **Oracle** estimators (true memberships): general and correctly specified
  (`ci_ge_or`, `ci_or`).
* **TGNQ-general** estimator (no structural assumption), with refinement
  (`ci_ge_gtnq`, `ci_ge_refine`).
* **TGNQ-special** estimator (correct structural assumption), with refinement
  (`ci_gtnq`, `ci_refine`).
* **TGNQ-misspecified** estimator (the *wrong* structural form: additive ↔
  multiplicative is swapped) (`ci_gtnq_mis`).

The correct structural form is selected automatically from the DGP:
`method = "additive"` (and `method_mis = "multiplicative"`) when `additive = 1`,
and vice versa.

`simulation(seed)` returns a 4-element list of error vectors:

1. `res_vec`        — special (correctly specified) model,
2. `res_vec_ge`     — general TGNQ model,
3. `res_vec_mis`    — misspecified TGNQ model,
4. `res_vec_mis_T`  — misspecified model at the correctly specified memberships.

These four vectors are aggregated over the 500 seeds into `vec_list`. In
`step2_output_tex.R`, streams 1–2 build **Table 1 (Section 5.1)** (the oracle and
TGNQ-with-refinement columns), while streams 2–4 build **Table 2 (Section 5.3)**
(the General / Additive / Multiplicative-MIS comparison).

---

## 5. Output script: tables produced and configuration switches

`step2_output_tex.R` assembles three data frames and emits each as LaTeX:

| object    | paper output             | content                                       | cell format        |
|-----------|--------------------------|-----------------------------------------------|--------------------|
| `res_spe` | Table 1 (Section 5.1)    | oracle / special-model + refined estimators   | `RMSE (AEcp)`      |
| `res_ge`  | Table 1 (Section 5.1)    | TGNQ-general + refined estimators             | `RMSE (AEcp)`      |
| `res_mis` | Table 2 (Section 5.3)    | General vs. Additive vs. Multiplicative (MIS) | RMSE only          |

Two settings are edited at the top of the script:

* **Network type** — which aggregated folder to summarize:

  ```r
  folder_path = "output/res_block/"   # SBM network
  # folder_path = "output/res_power/" # power-law network
  ```

* **DGP type** — by default the additive-DGP files (names ending in `a_TRUE`)
  are selected; switch the pattern to `"FA"` for the multiplicative-DGP files
  (`a_FALSE`):

  ```r
  file_me = r_script_files[grepl("TR", r_script_files)] %>% sort()
  ```

Files are processed in sorted name order, so the `G0 = 2` designs (cases 1–4)
fill rows 1–20 and the `G0 = 3` designs (cases 2–5) fill rows 21–40. 
Each design contributes 5 consecutive rows, one per quantile
`τ ∈ {0.1, 0.3, 0.5, 0.7, 0.9}`, for **40 rows** total.

---

## 6. Files in this folder

### 6.1 `step1_run_simulation.sh`

Shell driver that auto-locates the project root, then launches all
`(G0, case, additive, block)` jobs in two passes (see §2). Redirects each job's
stdout/stderr to `log/`. Concurrency is `MAX_JOBS = 1` (parallelism lives inside
R); `NCORES` controls the within-R core count.

### 6.2 `simulation.R`

R driver for one `(G0, case, additive, block)` design. It parses the five
command-line arguments, sources the DGP from `code/generator/`
(`generate_data.R`, plus `additive_parameter.R` when `additive = 1`), selects
the correct (`method`) and misspecified (`method_mis`) structural forms, ensures
the output folders exist, sources `simulator.R`, runs `simulation(seed)` in
parallel for `seed = 1..500`, and saves the aggregated `vec_list` to
`output/res_<block|power>/res-G_<G0>-c_<case>-a_<additive>.rda`.

### 6.3 `simulator.R`

Defines `simulation(seed)` — one full Monte Carlo replication shared by Table 1
and Table 2 (see §4). In addition to returning the four error vectors, it writes
a **per-seed snapshot** `.rda` to `output/<block|power>/` containing the fitted
objects and the simulated data.

### 6.4 `step2_output_tex.R`

The single post-processing script. It reads the aggregated `vec_list` files,
averages over seeds, unpacks them with `cal_res()` (`code/utils/assess.R`), and
prints the LaTeX for **Table 1 (Section 5.1)** (`res_spe`, `res_ge`) and
**Table 2 (Section 5.3)** (`res_mis`) via `kableExtra::kable`.
R packages required: `kableExtra`, `dplyr`, `here`.

---

## 7. Output directory layout

```
output/
├── block/            # per-seed snapshots, SBM network        -> read by step2_output_tex.R
│   └── <G0>-c_<case>-a_<additive>-seed_<seed>.rda
├── power/            # per-seed snapshots, power-law network
│   └── <G0>-c_<case>-a_<additive>-seed_<seed>.rda
├── res_block/        # aggregated vec_list per design, SBM     -> read by step2_output_tex.R
│   └── res-G_<G0>-c_<case>-a_<additive>.rda
└── res_power/        # aggregated vec_list per design, power
    └── res-G_<G0>-c_<case>-a_<additive>.rda
log/                  # one log per (G0, case, additive, block) job
```

Each per-seed snapshot stores:

| Object                                     | Description                                  |
|--------------------------------------------|----------------------------------------------|
| `X_tensor`, `Ymat`                         | simulated covariates / response panel        |
| `W`, `Amat`                                | row-normalized / raw adjacency network       |
| `member_G0`, `member_H0`                   | true row / column group memberships          |
| `ci_ge_or`                                 | oracle general estimator + CI                |
| `res_ge_update1`                           | TGNQ-general estimator (no refinement)       |
| `ci_ge_refine`                             | refined TGNQ-general estimator + CI          |
| `member_G_refine_ge`, `member_H_refine_ge` | refined row / column memberships (general)   |