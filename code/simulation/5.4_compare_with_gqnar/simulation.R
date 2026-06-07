# =============================================================================
# simulation.R
# -----------------------------------------------------------------------------
# R driver for ONE replicate of the Section 5.4 (TGNQ vs QGNAR) experiment.
# Invoked from step1_run_simulation.sh via Rscript:
#
#     Rscript simulation.R  case additive seed block
#
# Arguments:
#   case     - design index (1..4); see README for (N, T, Nblock)
#   additive - 1 = additive DGP, 0 = multiplicative DGP
#   seed     - Monte Carlo replicate id (1..500)
#   block    - 1 = SBM network, 0 = power-law network
#
# Workflow:
#   * source the data generator (and the additive parameter file, if needed),
#   * source estimator_gqnar.R (QGNAR implementation),
#   * make sure the appropriate output/res_gqnar/res_*/ folder exists,
#   * source simulator.R and call simulation(seed),
# which itself fits TGNQ and QGNAR and writes an .rda snapshot.
# =============================================================================
rm(list=ls())

# --- Parse command-line arguments from step1_run_simulation.sh ---------------
args <- commandArgs(TRUE)
case <- as.integer(args[1])
additive <- (args[2]==1)
seed = as.integer(args[3])
block = (args[4]==1)


# --- Working directory and packages ------------------------------------------
library(here)
dir <- here()
setwd(dir)


library(twmq)
library(CEoptim)
source("code/generator/generate_data.R")
source("code/simulation/5.4_compare_with_gqnar/estimator_gqnar.R")

# --- Method based on the DGP -------------------------------------------------
# Section 5.4 fits the TGNQ-general estimator (no special structure assumed),
# but the DGP-side `additive` flag determines how the truth is generated.
if(additive){
  source("code/generator/additive_parameter.R")
  method = "additive"
}else{
  method = "multiplicative"
}

# --- Make sure the output folder exists --------------------------------------
if (!file.exists("output/res_gqnar")) {
  dir.create("output/res_gqnar")
}

if (!file.exists(paste0("output/res_gqnar/res_", ifelse(block,"block_","power_"), ifelse(additive,"add/","mul/")))) {
  dir.create(paste0("output/res_gqnar/res_", ifelse(block,"block_","power_"), ifelse(additive,"add/","mul/")))
}



# --- Source the per-seed routine and run it ----------------------------------
source("code/simulation/5.4_compare_with_gqnar/simulator.R")
simulation(seed)
