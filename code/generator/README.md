# Data Generator

This folder contains the **data-generating code (DGP)** used in the paper

*Two-way Homogeneity Pursuit for Quantile Network Vector Autoregression*.

These scripts are sourced by the simulation programs in `code/simulation/` and are **not meant to be run directly**.

---

## Files

### `generate_data.R`

Provides the **baseline DGP** that is **always used** to generate data, regardless of whether the model is additive or multiplicative. It includes:

- Generation of true group memberships.
- Group-specific parameter functions:
  - `alpha.func`, `beta.func` (network effects),
  - `nu.func` (autoregressive effects),
  - `gamma.func` (covariate effects).
- Network generators:
  - `getBlockW()` for block networks,
  - `getPowerLawW()` for power-law networks.
- `simu.Y()` to generate panel responses from the TGNQ model given the network, groups, and covariates.

The parameter functions defined in this file correspond to the **multiplicative** specification
(`theta_{g h}(tau) = alpha_g(tau) * beta_h(tau)`).
When an **additive** DGP is required, these functions are overridden by sourcing
`additive_parameter.R` (see below); all other components of this file
(membership generation, network generators, and `simu.Y()`) remain in use.

---

### `additive_parameter.R`

Defines **alternative parameter functions** for the **additive TGNQ model**:

- `alpha.func`, `beta.func` for additive network effects,
- `nu.func`, `gamma.func` consistent with the additive specification.

When `additive = TRUE` in the simulation scripts, this file is sourced **after**
`generate_data.R` to **override** the corresponding parameter functions, so that the
data are generated under the additive structure
(`theta_{g h}(tau) = alpha_g(tau) + beta_h(tau)`).

---

## Usage

- These scripts are sourced by `simulation.R` files in Sections 5.1–5.4.
- Switching between **additive** and **general/multiplicative** DGPs is controlled by the `additive` flag in the simulation code:
  - `additive = FALSE` (default): use the multiplicative parameters from `generate_data.R`.
  - `additive = TRUE`: source `additive_parameter.R` to override them with the additive parameters.