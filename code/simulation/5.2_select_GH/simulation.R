# =============================================================================
# simulation.R
# -----------------------------------------------------------------------------
# R driver for ONE candidate-grid replication used in Section 5.2 of the TGNQ
# paper. It is invoked from step1_run_simulation.sh via Rscript:
#
#     Rscript simulation.R  case additive seed block
#
# where
#   case     - design index (1..5) mapping to (N, T, Nblock); see README
#   additive - 1 = additive DGP, 0 = multiplicative / general DGP
#   seed     - Monte Carlo replication index
#   block    - 1 = block (SBM) network, 0 = power-law network
#
# Throughout Section 5.2 we fix G0 = H0 = 3 (true number of row/column groups).
# This driver only sets up the environment; the actual simulation is run by
# simulator.R::simulation(seed), which fits seven candidate (G, H) models per
# replication and saves them to output/res_selectGH_<block|power>_<add|mul>/.
# =============================================================================
rm(list=ls())

# --- Parse command-line arguments from step1_run_simulation.sh --------------
args <- commandArgs(TRUE)
case <- as.integer(args[1])
additive <- (args[2]==1)
seed = as.integer(args[3])
block = (args[4]==1)

# --- Fixed true design (used by the generator inside simulator.R) -----------
G0 = H0 = 3

# --- Working directory and packages -----------------------------------------
library(here)
dir <- here()
setwd(dir)


library(twmq)
library(CEoptim)

# --- Source the data-generating process (DGP) -------------------------------
source("code/generator/generate_data.R")

# Override generator functions with the additive versions if needed; 
# this also fixes `method` to the corresponding special-model code that simulator.R expects.
if(additive){
  source("code/generator/additive_parameter.R")
  method = "additive"
}else{
  method = "multiplicative"
}

# --- Make sure the output folder exists -------------------------------------
if (!file.exists(paste0("output/res_selectGH_", ifelse(block, "block", "power"),ifelse(additive, "_add", "_mul"),"/"))) {
  dir.create(paste0("output/res_selectGH_", ifelse(block, "block", "power"),ifelse(additive, "_add", "_mul"),"/"))
}

# --- Source the per-replication routine and run one simulation --------------
# simulator.R defines simulation(seed), which fits the seven candidate (G, H)
# models and saves one .rda file per (G, H) pair.
source("code/simulation/5.2_select_GH/simulator.R")
simulation(seed)
