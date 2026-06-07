# =============================================================================
# simulation.R
# -----------------------------------------------------------------------------
# R driver for ONE simulation design of Sections 5.1 and 5.3 of the TGNQ paper.
# A single run fits all estimators (oracle / general / special / misspecified),
# so it feeds BOTH the Section 5.1 tables and the Section 5.3 table.
# It is invoked from step1_run_simulation.sh via Rscript:
#
#     Rscript simulation.R  G0 case additive block ncores
#
# Arguments:
#   G0       - true number of row groups (2 or 3 in Section 5.3)
#   case     - design index (1..5), maps to (N, T, Nblock); see README
#   additive - 1 = additive DGP, 0 = multiplicative DGP
#   block    - 1 = SBM network, 0 = power-law network
#   ncores   - number of cores used by doParallel inside R
#
# Workflow:
#   * pick the CORRECT method ("additive" / "multiplicative") and the MIS-
#     specified one (the opposite),
#   * make sure the output folders exist,
#   * source simulator.R and call simulation(seed) in parallel for seed = 1..500,
#   * save the list of per-seed result vectors to
#       output/res_<block|power>_mis/res-G_<G0>-c_<case>-a_<additive>.rda
# =============================================================================
rm(list=ls())

# --- Parse command-line arguments from step1_run_simulation.sh ---------------
args <- commandArgs(TRUE)
G0 <- as.integer(args[1])
case <- as.integer(args[2])
additive <- (args[3]==1)
block = (args[4]==1)
ncores = as.integer(args[5])


library(here)
dir <- here()
setwd(dir)



library(twmq)
library(CEoptim)
source("code/generator/generate_data.R")
source("code/utils/assess.R")


# --- Pick correct vs misspecified method based on the DGP --------------------
# `method`     = correct structural form (matches the DGP).
# `method_mis` = the WRONG structural form (additive <-> multiplicative).
if(additive){
  source("code/generator/additive_parameter.R")
  method = "additive"
  method_mis  = "multiplicative"
}else{
  method = "multiplicative"
  method_mis  =  "additive"
}


# --- Make sure the per-seed snapshot folder exists ---------------------------
# `<block|power>/` stores one .rda per (G0, case, additive, seed); this is used
# by simulator.R::simulation() to dump per-seed fits/diagnostics, which
# step2_output_tex.R reads back for the Section 5.1 tables.
if(block){
  if (!file.exists("output/block")) {
    dir.create("output/block")
  }
}else{
  if (!file.exists("output/power")) {
    dir.create("output/power")
  }
}


# --- Make sure the aggregated-results folder exists --------------------------
# `res_<block|power>/` stores ONE file per (G0, case, additive), holding the
# full vec_list returned by foreach() over all 500 seeds (consumed by
# step2_output_tex.R for both the Section 5.1 and Section 5.3 tables).
if(!file.exists(paste0("output/res_",ifelse(block,"block", "power")))){
  dir.create(paste0("output/res_",ifelse(block,"block", "power")))
}

# --- Source the per-seed routine and register the parallel backend -----------
source("code/simulation/5.1_and_5.3_Simulation/simulator.R")


library(doParallel)
library(foreach)
library(parallel)


registerDoParallel(cores = ncores)

# --- Run simulation(seed) in parallel for seeds 1..500 -----------------------
# Each simulation(seed) returns a 4-element list:
#   [[1]] res_vec       - special (correctly specified) model
#   [[2]] res_vec_ge    - general TGNQ model
#   [[3]] res_vec_mis   - misspecified TGNQ model
#   [[4]] res_vec_mis_T - misspecified TGNQ at the correctly-specified memberships
vec_list = foreach(seed = 1:500) %dopar% {
  simulation(seed)
}

# --- Save aggregated results -------------------------------------------------
save(vec_list, file = paste0(paste0("output/res_",ifelse(block,"block", "power")),"/res-G_",G0,"-c_",case,"-a_",additive,".rda"))

stopImplicitCluster()
