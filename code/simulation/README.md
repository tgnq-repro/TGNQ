# Simulation code

This folder contains all Monte Carlo simulation code used in **Section 5** of

"Two-way Homogeneity Pursuit for Quantile Network Vector Autoregression".

The four experiments of Section 5 are organized into **three** subfolders
(Sections 5.1 and 5.3 share one run), and together they produce **Tables 1–4**
of the main text and the corresponding supplementary tables.

This top-level README only explains the **shared** design and how to navigate
the suite; for the exact scripts, per-experiment design grid, argument order,
and `.rda` schema of each experiment, see the **`README.md` inside each
subfolder**.

---

## 1. Mapping to the paper

| Subfolder | Paper section(s) | Produces |
|-----------|------------------|----------|
| `5.1_and_5.3_Simulation/` | 5.1 and 5.3 | **Table 1** (estimation & inference) and **Table 2** (misspecified models), from a single run. |
| `5.2_select_GH/`          | 5.2 | **Table 3**: QIC-based selection of `(G, H)` vs. candidates and the oracle. |
| `5.4_compare_with_gqnar/` | 5.4 | **Table 4**: RMSE of the implied AR matrix `B(τ)`, TGNQ vs. one-way QGNAR. |

---

## 2. Shared pipeline

Every subfolder follows the same two-stage structure (simulate → aggregate):

```
step1_run_simulation.sh            (5.2 also has step2_run_oracle.sh)
        │  Bash driver: loops over the design switches and seeds,
        │  calls Rscript on simulation.R (5.2 also simulation_or.R).
        ▼
simulation.R  ->  simulator.R::simulation(seed)
        │  Generates the network and panel data, fits the estimators
        │  required by that experiment, saves .rda to output/.
        ▼
stepN_output_*_tex.R
        │  Aggregates over seeds/designs and prints the LaTeX tables.
        ▼
LaTeX tables for Section 5.
```

The post-processing entry point differs per experiment:

| Experiment | Post-processing script |
|------------|------------------------|
| 5.1 / 5.3 | `5.1_and_5.3_Simulation/step2_output_tex.R` |
| 5.2       | `5.2_select_GH/step3_output_tex.R` |
| 5.4       | `5.4_compare_with_gqnar/step2_output_tex.R` |

> **Command-line arguments differ across experiments** (see each subfolder's
> README). In brief:
> * `5.1_and_5.3_Simulation`: `Rscript simulation.R  G0 case additive block ncores`
>   — the 500 seeds are looped **inside** R via `doParallel`.
> * `5.2_select_GH` and `5.4_compare_with_gqnar`: `Rscript simulation.R  case additive seed block`
>   — one shell job per seed (`G0` is fixed internally: 3 in 5.2, 2 in 5.4).

---

## 3. Common conventions

- **Replications.** 500 Monte Carlo replications per design.
- **Quantile levels.** `τ ∈ {0.1, 0.3, 0.5, 0.7, 0.9}`.
- **Networks.** Block (SBM) via `getBlockW()` and power-law via
  `getPowerLawW()`, both in `code/generator/generate_data.R`.
- **DGPs.** Additive `θ_{gh}(τ)=α_g(τ)+β_h(τ)`
  (`code/generator/additive_parameter.R`) and multiplicative
  `θ_{gh}(τ)=α_g(τ)·β_h(τ)` (`code/generator/generate_data.R`).
- **Estimators.** Depending on the experiment: `oracle`, general (`ge`),
  special (`spe`), misspecified (`mis`, used in 5.3), and `QGNAR` (used in 5.4);
  TGNQ estimators are reported both without and with the Section-4.1 refinement.
- **Group-number selection (5.2).** QIC penalty
  `λ = N^0.1 · log(T-1) / ((T-1) · c · κ)` with `c = 10` and
  `κ ∈ {8.34 (N=100), 7.44 (N=200)}`; see `5.2_select_GH/README.md`.

The `case → (N, T, Nblock)` mapping is **not identical** across experiments
(5.2 uses a different ordering); the exact grid actually used is given in each
subfolder's README.

---

## 4. How to reproduce

From the project root `TGNQ/`, run the shell driver of each experiment, e.g.

```bash
bash "code/simulation/5.1_and_5.3_Simulation/step1_run_simulation.sh"
bash "code/simulation/5.2_select_GH/step1_run_simulation.sh"
bash "code/simulation/5.2_select_GH/step2_run_oracle.sh"
bash "code/simulation/5.4_compare_with_gqnar/step1_run_simulation.sh"
```

Each driver caps concurrency via `MAX_JOBS` (set at the top of the script;
note 5.1/5.3 uses `MAX_JOBS = 1` because parallelism lives inside R). Per-job
logs are written to `log/` at the project root. Then source the post-processing
script listed in §2 to print the LaTeX tables.

---

## 5. Output organization

All output is under `output/`; logs are in the sibling folder `log/`.

```
output/
├── block/   power/                       # 5.1 & 5.3 per-seed snapshots
├── res_block/   res_power/               # 5.1 & 5.3 aggregated results
│
├── res_selectGH_{block,power}_{add,mul}/ # 5.2 candidate (G, H)
├── {block,power}_select_or_{add,mul}/    # 5.2 oracle
│
└── res_gqnar/res_{block,power}_{add,mul}/ # 5.4 TGNQ vs QGNAR
log/                                       # one log file per Monte Carlo job
```

File-name conventions within each subfolder are documented in that subfolder's
README.

---

## 6. Dependencies

R packages (loaded inside the per-experiment `simulation.R`):

- `twmq` — the TGNQ estimation package for this paper;
- `CEoptim` — discrete optimization in the enhanced algorithm;
- `quantreg` — used by the QGNAR estimator (5.4);
- `doParallel`, `foreach` — within-R parallelism (5.1/5.3);
- `xtable`, `kableExtra` — LaTeX output in the post-processing scripts.

Helpers for RMSE, CI coverage, and clustering error rates live in
`code/utils/assess.R`; data-generating code in `code/generator/`.

---

## 7. Where to look for more detail

- Per-experiment scripts, design grids, argument order and `.rda` schema:
  the `README.md` inside `5.1_and_5.3_Simulation/`, `5.2_select_GH/`,
  `5.4_compare_with_gqnar/`.
- Data-generating code: `code/generator/README.md`.
- Post-processing helpers: `code/utils/README.md`.