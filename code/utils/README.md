# Utils

This folder contains **auxiliary functions** used by the simulation code in `code/simulation/` for evaluating and reshaping results. They are sourced by other scripts and are not meant to be run directly.

---

## Files

### `assess.R`

Functions for **measuring estimation performance** and **organizing outputs**:

- `cal_para(...)`  
  For the structured (additive / multiplicative) TGNQ model: computes RMSE of
  `alpha`, `beta`, `nu`, `gamma` and CI coverage indicators.

- `cal_para_general(...)`  
  For the general TGNQ model: computes RMSE of `alphabeta_GHs`, `nu`, `gamma`
  and CI coverage indicators.

- `cal_para_mis(...)`, `cal_para_mis_T(...)`, `cal_para_mis_ge(...)`  
  Variants used in **misspecification experiments** to compare correct vs
  misspecified models.

- `cal_res(res_vec)`  
  Reshapes a long result vector (RMSE + coverage + clustering errors across τ)
  into a matrix and extracts clustering error rates, for direct use in
  LaTeX tables.

These functions are mainly called by `step*_output_*_tex.R` scripts to build
the tables reported in the paper.