# Real-data analysis (Sina Weibo) — Section 6

This folder reproduces the **Sina Weibo application** in **Section 6** of the paper

"Two-way Homogeneity Pursuit for Quantile Network Vector Autoregression".

All scripts assume the project root is `TGNQ/` and the input data file is `data/weibo.rda`. Outputs are written under `output/real_res/`, `output/log/`, and `output/figs/`.

---

## 1. Pipeline at a glance

The analysis is decomposed into nine sequential steps. Each step is a self-contained script; the only coupling between steps is through the `.rda` files written under `output/real_res/`.

| Step | Script(s) | Purpose | Main output |
|------|-----------|---------|-------------|
| 1 | `step1_selectGH.sh` + `real_tgnq.R` | Fit TGNQ on a `(G, H) ∈ {1,…,5}²` grid with `T_train = 70`. | `output/real_res/tgnq_70-<G>-<H>.rda` |
| 2 | `step2_QIC_plot.R` | Compute the QIC criterion across `(G, H)` and produce the QIC figure. | `output/figs/QIC_plot.pdf` |
| 3 | `step3_cf_selectGH.R` | Out-of-sample (rolling) check of `H` at `G = 4` over test days 71–77. | per-day `pred_train70-test<T>-G4-H<H>.rda` |
| 4 | `step4_fit_tgnq_alldata.R` | Refit TGNQ on the full training set with the chosen `(G, H)` (e.g. `T_train = 77, G = 4, H = 3`). | `output/real_res/tgnq_77-4-3.rda` |
| 5 | `step5_fit_qnarg.R` (uses `estimator_gqnar.R`) | Fit the QGNAR benchmark for `T_train ∈ {70, 77}` and `G ∈ {3, 4}`. | `output/real_res/qnarg_G<G>-<T_train>.rda` |
| 6 | `step6_predict_qnarg.R` | One-step-ahead QGNAR predictions for test days 71–77. | `output/real_res/narg_pred_train70-test<T>.rda` |
| 7 | `step7_fit_compare.R` | In-sample goodness-of-fit comparison between TGNQ and QGNAR. | console output |
| 8 | `step8_pred_plot.R` | Aggregate prediction losses over test days; plot loss vs. day by τ. | `output/figs/prediction.pdf` |
| 9 | `step9_inference.R` | Membership refinement, CI computation, and final figures. | `output/real_res/r43-4-3.rda`, `output/figs/{network_effect, covar_effect, activity}.pdf` |

The dependency graph between steps is essentially linear:

```
step1  →  step2 ──┐
   │              ├──► step4 ──► step5 ──► step6 ──► step7
   └────► step3 ──┘                        │
                                           └──► step8
                                step4  ──► step9
```

---

## 2. Data file: `data/weibo.rda`

The `.rda` file contains:

- **`Ymat`** — `N × T` numeric matrix of responses. Column `t` is the response at time `t`.
- **`Amat`** — `N × N` (0/1) directed adjacency matrix of the follow / interaction network. The row-normalized matrix $W = D^{-1} A$ is built inside each script.
- **`Xi`** — `N × 8` matrix of time-invariant node-level covariates.
- **`Xt`** — `T × 2` matrix of node-invariant time-level covariates.

Every script assembles a 3-D covariate tensor `X_train` of shape `N × T_train × 10` in which slices `1:8` are `Xi` (replicated across time) and slices `9:10` are `Xt` (replicated across nodes), with the training/response alignment `Xt[2:(T_train+1), ]` ↔ `Y_{t+1}`.

---

## 3. Estimators

Two estimators are compared.

### 3.1 TGNQ (proposed)

The two-way grouped network quantile model. Implemented in the package `twmq` and called through:

- `twmq.estimate.auto.parallel(...)` — multi-start initialization;
- `update_NARG_twmq_parallel(...)` — local refinement of memberships and parameters;
- `Refine_G_parallel(...)`, `Refine_H_parallel(...)` — Section 4.1 membership refinement;
- `twmq.estimate_thetaGH.member.iterate(...)` — final parameter fit at refined memberships;
- `twmq_ci(...)` — pointwise confidence intervals at the five quantiles.

### 3.2 QGNAR (benchmark)

A grouped quantile network autoregression with a single, common grouping. The full implementation is provided in `estimator_gqnar.R`, with these key functions:

- `est.NARG.init(...)` — produces a pool of K-means initializations from a node-level lasso quantile pre-fit;
- `est.NARG(...)` — given memberships, fits group-level coefficients via L1-regularized quantile regression for all τ jointly (`quantreg::rq` with `method = "lasso"`);
- `est.member(...)` — updates each individual's membership by minimizing its quantile loss given the current parameters;
- `est.NARG.member(...)` — alternates between the two updates with a two-stage criterion (local then global loss);
- `est.NARG.member.parallel(...)` — runs `est.NARG.member` from many K-means starts in parallel and keeps the best.

---

## 4. Required R packages

Estimation: `twmq`, `quantreg`, `CEoptim`, `MASS`, `compiler`.
Parallelism: `parallel`, `foreach`, `doParallel`, `doSNOW`.
Plotting: `ggplot2`, `reshape2`, `gridExtra`, `patchwork`.

A multi-core machine is recommended.

---

## 5. Output organization

```
output/
├── real_res/
│   ├── tgnq_<T_train>-<G>-<H>.rda            # Step 1, 4
│   ├── pred_train70-test<T>-G4-H<H>.rda      # Step 3
│   ├── qnarg_G<G>-<T_train>.rda              # Step 5
│   ├── narg_pred_train70-test<T>.rda         # Step 6
│   ├── rG43.rda, rH43.rda                    # Step 9 (refined memberships)
│   └── r43-4-3.rda                           # Step 9 (final fit + CI)
├── log/                                       # per-job stdout/stderr
└── figs/
    ├── QIC_plot.pdf                          # Step 2
    ├── prediction.pdf                        # Step 8
    ├── network_effect.pdf                    # Step 9
    ├── covar_effect.pdf                      # Step 9
    └── activity.pdf                          # Step 9
```

---

## 6. How to reproduce (recommended order)

From the project root `TGNQ/`:

```bash
# Step 1: fit the (G, H) grid in parallel (background jobs)
bash code/real_data/step1_selectGH.sh

# Step 2: QIC plot (run interactively in R)
Rscript code/real_data/step2_QIC_plot.R

# Step 3: rolling out-of-sample selection of H
Rscript code/real_data/step3_cf_selectGH.R

# Step 4: refit on the full training window
Rscript code/real_data/step4_fit_tgnq_alldata.R 77 4 3 0.01

# Step 5–6: QGNAR benchmark
Rscript code/real_data/step5_fit_qnarg.R
Rscript code/real_data/step6_predict_qnarg.R

# Step 7–8: comparison and prediction plot
Rscript code/real_data/step7_fit_compare.R
Rscript code/real_data/step8_pred_plot.R

# Step 9: inference + final figures
Rscript code/real_data/step9_inference.R
```

Step 1 launches 25 background `Rscript` jobs via `nohup`. Adjust the parallelism (or split into batches) according to your machine's capacity.

---

## 7. Key conventions

- **Quantile grid.** All steps use `taus = c(0.1, 0.3, 0.5, 0.7, 0.9)`.
- **Check function.** `check.func(u, τ) = u (τ − 1{u < 0})`.
- **Row-normalization of W.** `W = Amat / rowSums(Amat)` for the first `N` nodes; nodes with zero out-degree should be filtered out beforehand to avoid division by zero.
- **Training/test alignment.** `Y_train = Ymat[, 1:(T_train + 1)]` (so Y_train uses `t = 0,…,T_train`), and `Xt = Xt[2:(T_train + 1), ]` so that time covariates at time `t` are aligned with response `Y_{t+1}`.
- **QIC criterion** (Step 2):

  $$
    \mathrm{QIC}(G, H) = \log \overline{\ell}_{G,H}
      + \lambda_{NT} \cdot \big( G H + (P+1) G \big),
    \qquad
    \lambda_{NT} = \frac{N^{0.1} \log T}{10\, T \min(\bar d, 10)},
  $$

  with `P = 8` (covariates excluding the intercept) and $\bar d = \mathrm{mean degree} = \sum A / N$.

- **QGNAR `psi`.** In the QGNAR output, `narg_mem_final$psi[[k]]` is a `G × (G + 1 + P)` matrix whose first `G` columns are network-effect coefficients, column `G+1` is the autoregressive coefficient, and the remaining `P` columns are covariate coefficients.

---
