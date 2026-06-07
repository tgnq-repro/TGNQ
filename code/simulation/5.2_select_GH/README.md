# Simulation code for Section 5.2 (Selection of G and H)

This folder contains the code for the experiments in **Section 5.2 (“select G and H”)** of the paper

"Two-way Homogeneity Pursuit for Quantile Network Vector Autoregression",

where we study how a QIC-type information criterion selects the number of row groups `G` and column groups `H`, and compare the resulting estimator against:

* the **oracle** estimator (using the true `(G0, H0)` and the true memberships),
* candidate **mis-specified** models with `(G, H)` in  `{ (2,2), (2,3), (3,2), (3,3), (3,4), (4,3), (4,4) }`,
* both the **TGNQ general (`ge`)** model and the **special (`spe`)** model (additive or multiplicative) used in the paper.

The true design throughout Section 5.2 uses `G0 = H0 = 3`.


---

## 1. Pipeline overview

```
step1_run_simulation.sh  -- loops over (case, additive, block, seed) -->
        |
        v
simulation.R             -- parses args, sources DGP code -->
        |
        v
simulator.R::simulation()  -- fits seven candidate (G, H) models per replication
        |
        v
output/res_selectGH_<block|power>_<add|mul>/<case>-<G>-<H>-<seed>.rda


step2_run_oracle.sh      -- loops over (case, additive, block, seed) -->
        |
        v
simulation_or.R          -- parses args, sources DGP code -->
        |
        v
simulator.R::simulation_or() -- fits oracle estimators at (G0, H0)
        |
        v
output/<block|power>_select_or_<add|mul>/<case>-<seed>.rda


step3_output_tex.R       -- aggregates everything, applies QIC, prints LaTeX tables for Section 5.2.
```

---

## 2. How to reproduce Section 5.2

From the project root (`TGNQ/`), run the following in order:

1. **Candidate (G, H) simulations (TGNQ).**

   ```bash
   bash "code/simulation/5.2_select_GH/step1_run_simulation.sh"
   ```

   Result files are written to:

   * `output/res_selectGH_block_add/`, `output/res_selectGH_block_mul/`
   * `output/res_selectGH_power_add/`, `output/res_selectGH_power_mul/`

2. **Oracle simulations (true `(G0, H0)` only).**

   ```bash
   bash "code/simulation/5.2_select_GH/step2_run_oracle.sh"
   ```

   Result files are written to:

   * `output/block_select_or_add/`, `output/block_select_or_mul/`
   * `output/power_select_or_add/`, `output/power_select_or_mul/`

3. **Post-process and print LaTeX tables.**

   In R:

   ```r
   source("code/simulation/5.2_select_GH/step3_output_tex.R")
   ```

Each step is fully parallel and logs each Monte Carlo job to `log/`.

---

## 3. Case → design mapping

The `case` argument indexes the (N, T, Nblock) design:

| case | N   | T   | Nblock |
|------|-----|-----|--------|
| 1    | 100 | 100 | 5      |
| 2    | 100 | 200 | 5      |
| 3    | 200 | 200 | 10     |
| 4    | 200 | 400 | 10     |
| 5    | 100 | 50  | 5      |

---

## 4. Files in this folder

### 4.1 `step1_run_simulation.sh`

Shell driver for the **candidate (G, H)** grid. Loops over seeds 1..500
and the four design switches `(case, additive, block)`, submitting

```
Rscript simulation.R  case additive seed block
```

with at most `MAX_JOBS = 128` concurrent R processes; each job's stdout/stderr is logged in `log/selectGH_case*_add*_block*_seed*.log`.

### 4.2 `simulation.R`

R driver for **one** replication of the candidate-grid simulation:

* parses the four command-line arguments `(case, additive, seed, block)`,
* sets `G0 = H0 = 3`,
* picks `method = "additive"` or `"multiplicative"` based on `additive`,
* makes sure `output/res_selectGH_*/` exists,
* sources `simulator.R` and runs `simulation(seed)`.

### 4.3 `simulator.R`

Defines two functions:

* **`simulation(seed)`** — performs one replication for **all seven** candidate `(G, H)` pairs. 
  For each pair it runs `twmq.estimate.auto`  (with `ntrial = 100` or `ntrial = 1000` depending on `(G, H)`),
  selects the three initializations with the smallest loss, refines them with `update_NARG_twmq`, 
  takes the best one, and finally runs one more `update_NARG_twmq` to obtain the special-model estimator.
  Output files have the form  `output/res_selectGH_*/<case>-<G>-<H>-<seed>.rda`.
* **`simulation_or(seed)`** — performs one replication of the **oracle** experiment: 
  it generates data the same way but only fits two estimators at the true memberships, 
  namely the general and the special TGNQ models (`twmq.estimate_thetaGH.member.iterate`).  
  Output files have the form  `output/<block|power>_select_or_<add|mul>/<case>-<seed>.rda`.

### 4.4 `simulation_or.R`

R driver for **one** oracle replication: parses the arguments, creates the target folder if needed, sources `simulator.R`, and calls `simulation_or(seed)`.

### 4.5 `step2_run_oracle.sh`

Shell driver for the **oracle** simulations; identical structure to `step1_run_simulation.sh` but calls `simulation_or.R`.

### 4.6 `step3_output_tex.R`

R post-processing script. For each design `(case, additive, block)`:

* loads every `.rda` file in `output/res_selectGH_*/` whose seed produced results for all seven `(G, H)` pairs,
* for each replication, computes the QIC value of every candidate `(G, H)` and selects the minimizer,
* aggregates the **selection frequencies** for the QIC criterion under both the general and the special model,
* computes the **RMSE** of the selected models (and of every candidate individually) for `α/β`, `αβ` (combined network effect), `ν`, and `γ`,
* computes the **clustering error rates** `ρ_G`, `ρ_H` using `err_rate_mapping`,
* finally, it overlays the **oracle** RMSEs from `output/<block|power>_select_or_<add|mul>/*.rda` and prints the consolidated LaTeX table for Section 5.2.

The selection penalty in the QIC is

```
λ = N^0.1 / (T - 1) / c / κ * log(T - 1)
```
with `c = 10` and `κ ∈ {8.34 (N = 100), 7.44 (N = 200)}`, matching the choices used in the paper.


---

## 5. Output structure

```
output/
├── res_selectGH_block_add/    # candidate (G,H) results, block + additive
├── res_selectGH_block_mul/    # candidate (G,H) results, block + multiplicative
├── res_selectGH_power_add/    # candidate (G,H) results, power-law + additive
├── res_selectGH_power_mul/    # candidate (G,H) results, power-law + multiplicative
├── block_select_or_add/       # oracle results, block + additive
├── block_select_or_mul/       # oracle results, block + multiplicative
├── power_select_or_add/       # oracle results, power-law + additive
└── power_select_or_mul/       # oracle results, power-law + multiplicative
log/                           # one log file per Monte Carlo job
```

Each `res_selectGH_*` file stores the two estimators `res_ge_update`, `res_update` for one `(case, G, H, seed)`, together with the true memberships `member_G0`, `member_H0`.

Each `*_select_or_*` file stores the two oracle estimators `or_ge`, `or_spe` together with `member_G0`, `member_H0`.


