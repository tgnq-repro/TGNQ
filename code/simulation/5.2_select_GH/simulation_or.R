# =============================================================================
# simulation_or.R
# -----------------------------------------------------------------------------
# R driver for ONE oracle replication used in Section 5.2 of the TGNQ paper.
# It is invoked from step2_run_oracle.sh via Rscript:
#
#     Rscript simulation_or.R  case additive seed block
#
# Arguments are identical to those of simulation.R:
#   case     - design index (1..5)
#   additive - 1 = additive DGP, 0 = multiplicative DGP
#   seed     - replication index
#   block    - 1 = SBM network, 0 = power-law network
#
# This driver only sets up the environment; the actual simulation is run by
# simulator.R::simulation_or(seed), which fits the GENERAL and SPECIAL TGNQ
# estimators at the TRUE memberships (G0, H0) and saves the result file under
#   output/<block|power>_select_or_<add|mul>/<case>-<seed>.rda
# =============================================================================

rm(list=ls())

# --- Parse command-line arguments from step2_run_oracle.sh ------------------
args <- commandArgs(TRUE)
case <- as.integer(args[1])
additive <- (args[2]==1)
seed = as.integer(args[3])
block = (args[4]==1)

# --- Fixed true design -------------------------------------------------------
G0 = 3

# --- Working directory and packages ------------------------------------------
library(here)
dir <- here()
setwd(dir)

library(twmq)
source("code/generator/generate_data.R")
source("code/utils/assess.R")

# --- Pick the right DGP / special-model code --------------------------------
if(additive){
  source("code/generator/additive_parameter.R")
  method = "additive"
}else{
  method = "multiplicative"
}

# --- Make sure the output folder exists -------------------------------------
if(!file.exists(paste0("output/", ifelse(block, "block", "power"),"_select_or",ifelse(additive, "_add", "_mul"),"/"))){
  dir.create(paste0("output/", ifelse(block, "block", "power"),"_select_or",ifelse(additive, "_add", "_mul"),"/"))
}

# --- Source the per-replication routine and run one oracle simulation ------
source("code/simulation/5.2_select_GH/simulator.R")

simulation_or(seed)