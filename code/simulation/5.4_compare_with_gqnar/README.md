# Simulation code for Section 5.4 (Comparison with QGNAR)

This folder contains the code for the experiments in **Section 5.4 (“compare with gqnar”)** of the paper

"Two-way Homogeneity Pursuit for Quantile Network Vector Autoregression",

where the proposed **TGNQ** estimator is compared head-to-head with the
**grouped quantile network autoregressive (QGNAR / gqnar)** model of
Zhu et al., in terms of how accurately each method recovers the
**autoregressive matrix** $B(\tau)$ of the network VAR.

The two methods differ in their assumptions:

* **TGNQ** uses a *two-way* group structure.
* **QGNAR** uses a *one-way* group structure.

Both yield an implied autoregressive matrix $B(\tau) \in \mathbb{R}^{N\times N}$ and the comparison metric is the RMSE of the nonzero entries of $B(\tau)$ against the truth $B_0(\tau)$, averaged over Monte Carlo seeds.


---

## 1. Pipeline overview

```
step1_run_simulation.sh   -- loops over (case, additive, block, seed) -->
        |
        v
simulation.R              -- parses args; sources DGP and simulator   -->
        |
        v
simulator.R::simulation() -- one seed: fits TGNQ + QGNAR, saves .rda
        |
        v
output/res_gqnar/res_<block|power>_<add|mul>/<G>-<case>-<seed>.rda

step2_output_tex.R        -- rebuilds B_0(τ), reads TGNQ + QGNAR fits,
                             computes RMSE of B(τ), prints LaTeX tables.
```

---

## 2. How to reproduce Section 5.4

From the project root (`TGNQ/`), run in order:

1. **Run all simulation experiments.**

   ```bash
   bash "code/simulation/5.4_compare_with_gqnar/step1_run_simulation.sh"
   ```

   This loops over

   * `case`     ∈ {1, 2, 3, 4} — design index (see table below),
   * `additive` ∈ {0, 1}       — multiplicative vs additive DGP,
   * `block`    ∈ {0, 1}       — power-law vs SBM network,
   * `seed`     ∈ {1, …, 500}  — Monte Carlo replicate index.

   For each combination, it calls `simulation.R` and writes one `.rda` snapshot to:

   ```
   output/res_gqnar/res_block_add/   output/res_gqnar/res_block_mul/
   output/res_gqnar/res_power_add/   output/res_gqnar/res_power_mul/
   ```

   Logs are written to `log/gqnar_case<case>_add<add>_seed<seed>_block<block>.log`.

2. **Aggregate and print LaTeX tables.**

   In R:

   ```r
   source("code/simulation/5.4_compare_with_gqnar/step2_output_tex.R")
   ```

   This reads all snapshot `.rda` files, reconstructs the implied $B(\tau)$ matrices for TGNQ and QGNAR, computes the RMSE of the nonzero entries against the true $B_0(\tau)$, and prints LaTeX tables via `xtable`.

---

## 3. Case → design mapping

In Section 5.4, only `G0 = H0 = 2` is used. The case index maps to
`(N, T, Nblock)` as:

| case | N   | T   | Nblock |
|------|-----|-----|--------|
| 1    | 100 | 50  | 5      |
| 2    | 100 | 100 | 5      |
| 3    | 100 | 200 | 5      |
| 4    | 200 | 200 | 10    |

True row-group ratios are `(0.5, 0.5)` and true column-group ratios
are `(0.4, 0.6)`.


---

## 4. Output schema

For each seed `simulation(seed)` saves an `.rda` containing:

* `member_G0`, `member_H0` — the true row- and column-group memberships,
* `res_ge_update1`         — the TGNQ-general fit (after label switching), whose `theta_GH$alphabeta_GHs[[k]]` and `theta_GH$theta_Gs[[k]]` give the off-diagonal and diagonal blocks of $B_{TGNQ}(\tau_k)$,
* `narg`                   — the QGNAR fit returned by `gqnar(...)`, whose `narg$psi[[k]]` is a $G \times (G+1+P)$ matrix containing the follower group effects, the own-lag coefficient, and the covariate coefficients for quantile $\tau_k$.

`step2_output_tex.R` reconstructs both `Btgnq` and `Bnarg` per τ and compares them to the true `B0`.

---

## 5. Files in this folder

### 5.1 `step1_run_simulation.sh`

Shell driver that submits all 8000 `(case, additive, block, seed)` jobs, each running `simulation.R` once. Uses `MAX_JOBS = 128` for shell-level concurrency and writes per-job logs to `log/`.

### 5.2 `simulation.R`

R driver for ONE `(case, additive, block, seed)`:

* parses the four command-line arguments,
* sets the working directory and loads `twmq`, `CEoptim`,
* sources:
  * `code/generator/generate_data.R`,
  * `code/generator/additive_parameter.R` (if `additive`),
  * `estimator_gqnar.R`         — the QGNAR implementation,
  * `simulator.R`               — defines `simulation(seed)`,
* ensures the appropriate `output/res_gqnar/res_*` folder exists, calls `simulation(seed)`.

### 5.3 `simulator.R`

Defines `simulation(seed)` which, for a single replication:

* sets `(N, T, Nblock)` from `case`,
* builds the network `W`, the friend lists `FriendW`, `FriendW2`,
* draws the true memberships and simulates the panel `Ymat`,
* fits **TGNQ** via `twmq.estimate.auto` → `update_NARG_twmq` → `twmq.label.switch` (so estimated group labels match the truth),
* fits **QGNAR** via `gqnar(...)` from `estimator_gqnar.R`,
* saves both fits plus the true memberships to `output/res_gqnar/res_<block|power>_<add|mul>/<G>-<case>-<seed>.rda`.

### 5.4 `estimator_gqnar.R`

Implementation of the QGNAR estimator used in the comparison:

* `est.NARG.init(...)`         — k-means based initial-membership search, built on a per-individual lasso quantile regression (`quantreg::rq(..., method = "lasso")`),
* `est.NARG(...)`              — given a membership, fits the group-level parameters $(\beta, \nu, \gamma)$ by joint quantile regression for all quantile levels in `taus`,
* `est.member(...)`            — updates the membership of each individual by minimizing the relevant local check loss; supports an "Initial" mode (only the individual's own loss) and a full mode (their loss + that of all their followers),
* `est.NARG.member(...)`       — two-stage iterative scheme: first local-loss updates, then full-loss updates, until convergence,
* `gqnar(...)`                 — the main wrapper used by `simulator.R`: runs `est.NARG.init` for multiple k-means types and trials, picks the initialization with the lowest loss, and runs `est.NARG.member` on it.


### 5.5 `step2_output_tex.R`

Post-processing script for Section 5.4:

* iterates over `block ∈ {TRUE, FALSE}`, `additive ∈ {FALSE, TRUE}`, `case ∈ {1, 2, 3, 4}`,
* for each design, rebuilds the true $B_0(\tau)$ matrix from the oracle group parameters, then over all seed files, builds $B_{\text{TGNQ}}(\tau)$ and $B_{\text{QGNAR}}(\tau)$, and records the RMSE of the nonzero entries (i.e. those at positions where `W != 0`),
* assembles a wide LaTeX table (one row per `(N, T, τ)`, two competing methods per (additive, network) block) via `xtable`.

---

## 6. Output directory layout

```
output/
└── res_gqnar/
    ├── res_block_add/    # block network, additive DGP
    ├── res_block_mul/    # block network, multiplicative DGP
    ├── res_power_add/    # power-law network, additive DGP
    └── res_power_mul/    # power-law network, multiplicative DGP
log/                      # one log per (case, additive, seed, block) job
```

Each file name encodes `(G, case, seed)` as `<G>-<case>-<seed>.rda`.


---

