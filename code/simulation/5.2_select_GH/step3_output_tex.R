# =============================================================================
# step3_output_tex.R
# -----------------------------------------------------------------------------
# Post-processing of the Section 5.2 ("select G and H") simulations.
#
# For each design (case, additive, block) the script:
#   * loads the seven candidate-(G, H) result files per seed from
#     output/res_selectGH_<block|power>_<add|mul>/,
#   * applies the QIC criterion
#         QIC(G, H) = log(L_bar) + lambda * (G*H + 3*G)
#     to pick the best (G, H) under both the general (ge) and the special
#     (spe) TGNQ models,
#   * aggregates SELECTION FREQUENCIES of (G, H) over seeds,
#   * computes RMSEs (alpha, beta, alphabeta, nu, gamma) for each candidate
#     and for the selected model, and clustering error rates rho_G, rho_H,
#   * overlays the ORACLE RMSEs from output/<block|power>_select_or_<add|mul>/
#     into the table,
#   * prints the consolidated LaTeX table for Section 5.2.
#
# Switch between designs at the top of the script via:
#   power    = FALSE -> block network ; TRUE -> power-law network
#   additive = TRUE  -> additive DGP  ; FALSE -> multiplicative DGP
# =============================================================================

library(here)
dir <- here()
setwd(dir)

library(dplyr)
library(twmq)

# ---- Choose which design to summarize --------------------------------------
power    <- FALSE   # FALSE = block network, TRUE = power-law network
additive <- TRUE    # TRUE  = additive DGP,  FALSE = multiplicative DGP
# Load the generator functions matching the DGP (for true-parameter values).
if(additive){
  source("code/generator/additive_parameter.R")
}else{
  source("code/generator/generate_data.R")
}

# ---- True-parameter setup --------------------------------------------------
c=10
taus = c(0.1,0.3,0.5,0.7,0.9)
G0 = H0 = 3

# True alpha, beta, nu, gamma at each quantile 
alpha0 = sapply(taus, function(tau)
  sapply(1:G0, function(g) alpha.func(qnorm(tau), g)))
beta0 = sapply(taus, function(tau)
  sapply(1:H0, function(h) beta.func(qnorm(tau), h)))
nu0 = sapply(taus, function(tau)
  sapply(1:G0, function(g) nu.func(qnorm(tau), g, G0)))
gamma0 = lapply(taus, function(tau)
  t(sapply(1:G0, function(g) c(
                               gamma.func(qnorm(tau), g, 2),
                               gamma.func(qnorm(tau), g, 3)))))
# True combined "alphabeta" matrix at each quantile, additive or multiplicative
alphabeta0 <- list()
if(additive){
  for (i in 1:5) {
    alphabeta0[[i]] <- outer(alpha0[, i],rep(1,3))   +  outer(rep(1,3),beta0[, i])
  }
}else{
  for (i in 1:5) {
  alphabeta0[[i]] <- outer(alpha0[, i],beta0[, i])
  }
}

# ---- Pick the right candidate-results folder -------------------------------
if(power){
  if(additive){
    folder_path <- "output/res_selectGH_power_add/"
  }else{
    folder_path <- "output/res_selectGH_power_mul/"
  }
}else{
  if(additive){
    folder_path <- "output/res_selectGH_block_add/"
  }else{
    folder_path <- "output/res_selectGH_block_mul/"
  }
}





file_names <- list.files(path = folder_path, full.names = TRUE)


# -----------------------------------------------------------------------------
# result_ge / result_spe: for one `seed`, gather (seed, GH, # parameters,
# scaled loss) for each of the seven candidate (G, H) fits. Used to evaluate
# the QIC criterion in selectGH().
#
# Columns:
#   1) seed
#   2) integer encoding of (G, H), e.g. (3, 4) -> 34
#   3) effective parameter count  (G*H + 3*G  -- one per (g,h) block plus
#      three group-level parameters per row group)
#   4) average per-cell loss      (sum(loss) / N / T / 5)
# -----------------------------------------------------------------------------
result_ge = function(seed){
  seed_file = matching_files[seed_vec == seed]
  res = list()
  i = 1
  for (f in seed_file) {
    load(f)
    G = max(res_ge_update$member_G)
    H = max(res_ge_update$member_H)
    loss = res_ge_update$theta_GH$loss
    res[[i]] = c(as.integer(seed),as.integer(paste0(G,H)),G*H+3*G,sum(loss)/N/Time/5)
    i = i + 1
  }
  res = do.call(rbind, res)
  return(res)
}


result_spe = function(seed){
  seed_file = matching_files[seed_vec == seed]
  res = list()
  i = 1
  for (f in seed_file) {
    load(f)
    G = max(res_update$member_G)
    H = max(res_update$member_H)
    loss = res_update$theta_GH$loss
    res[[i]] = c(as.integer(seed),as.integer(paste0(G,H)),G*H+3*G,sum(loss)/N/Time/5)
    i = i + 1
  }
  res = do.call(rbind, res)
  return(res)
}


# -----------------------------------------------------------------------------
# selectGH(df, lambda): returns the (G, H) (encoded as integer) that
# minimizes QIC = log(loss) + lambda * (# params).
# -----------------------------------------------------------------------------
selectGH = function(df,lambda){
  QIC = log(df[,4]) + lambda * df[,3]
  GH = df[,2][which.min(QIC)]
}

# -----------------------------------------------------------------------------
# cal_RMSE(seed, GH_ge, GH_spe): for one seed, compute per-quantile RMSE-like
# averages of |estimate - truth| for the seven candidate (G, H) fits, then
# stack the RMSEs of the QIC-selected models at the end.
#
# Returned vector layout (length 7 * 7 + 4 + 3 = 56):
#   * 7 candidate (G, H) fits in the order corresponding to `matching_files`,
#     each contributing a length-7 RMSE vector:
#       [1] alpha  (special)
#       [2] beta   (special)   -- NB: kept as 0 in this script's vector layout
#       [3] nu     (special)
#       [4] gamma  (special)
#       [5] alphabeta (general)
#       [6] nu        (general)
#       [7] gamma     (general)
#   * the special-model RMSEs of the QIC-selected (G, H) under the special
#     model (4 numbers),
#   * the general-model RMSEs of the QIC-selected (G, H) under the general
#     model (3 numbers).
# -----------------------------------------------------------------------------
cal_RMSE = function(seed,GH_ge, GH_spe){
  seed_file = matching_files[seed_vec == seed]
  rmse = numeric(7)
  res = c()
  gh = c(22,23,32,33,34,43,44)
  i_ge = which(gh == GH_ge)
  i_spe = which(gh == GH_spe)
  
  i=0
  for (f in seed_file) {
    load(f)
    print(f)
    # Average per-quantile |est - true|; divide by 5 (quantiles) at the end.
    for (q in 1:5) {
      
      # ---- General-model errors -------------------------------------------
      theta_all_ge = res_ge_update$theta_GH$alphabeta_GHs[[q]][res_ge_update$member_G,res_ge_update$member_H]
      theta_all0 = alphabeta0[[q]][member_G0, member_H0]
      rmse[5] = rmse[5] + mean(abs(theta_all_ge - theta_all0))
      
      columns <- lapply(res_ge_update$theta_GH$theta_Gs, function(mat) mat[, 1])
      nu <- do.call(cbind,columns)
      nu_all_ge = nu[,q][res_ge_update$member_G]
      nu_all0 = nu0[,q][member_G0]
      rmse[6] = rmse[6] + mean(abs(nu_all_ge - nu_all0))
      
      gamma = lapply(res_ge_update$theta_GH$theta_Gs, function(mat) mat[, 2:3])
      gamma_all_ge = gamma[[q]][res_ge_update$member_G,]
      gamma_all0 = gamma0[[q]][member_G0,]
      r = sqrt(rowSums((gamma_all_ge - gamma_all0)^2))
      rmse[7] = rmse[7] + mean(r)
      
      # ---- Special-model errors (additive or multiplicative) ---------------
      alpha = res_update$theta_GH$alpha
      beta = res_update$theta_GH$beta
      alpha_all = alpha[,q][res_update$member_G]
      beta_all = beta[,q][res_update$member_H]
      if(additive){
        alphabeta_all = outer(alpha_all, rep(1, length(alpha_all))) + outer(rep(1,length(beta_all)), beta_all)
      }else{
        alphabeta_all = outer(alpha_all, rep(1, length(alpha_all))) * outer(rep(1,length(beta_all)), beta_all)
      }
      
      rmse[1] = rmse[1] +mean(abs(alphabeta_all - theta_all0))

      columns <- lapply(res_update$theta_GH$theta_Gs, function(mat) mat[, 2])
      nu <- do.call(cbind,columns)
      nu_all_ge = nu[,q][res_update$member_G]
      nu_all0 = nu0[,q][member_G0]
      rmse[3] = rmse[3] + mean(abs(nu_all_ge - nu_all0))
      
      gamma = lapply(res_update$theta_GH$theta_Gs, function(mat) mat[, 3:4])
      gamma_all_ge = gamma[[q]][res_update$member_G,]
      gamma_all0 = gamma0[[q]][member_G0,]
      r = sqrt(rowSums((gamma_all_ge - gamma_all0)^2))
      rmse[4] = rmse[4] + mean(r)
    }
    rmse = rmse / 5 # average over the 5 quantile levels
    res= c(res,rmse)
    
    i=i+1
    if(i == i_ge){
      temp_ge = rmse[5:7]# RMSE block of the QIC pick (ge)
    }
    if(i == i_spe){
      temp_spe = rmse[1:4]# RMSE block of the QIC pick (spe)
    }
    rmse = numeric(7)# reset for next candidate
  }
  res = c(res, temp_spe,temp_ge)
  return(res)
}


# -----------------------------------------------------------------------------
# cal_rho(seed, GH_ge, GH_spe): for one seed, compute clustering error rates
# rho_g, rho_h for each of the seven candidate (G, H) fits (both ge and spe)
# and append the value corresponding to the QIC-selected (G, H).
# Leading 0 is a placeholder used to align with the table's "Oracle" row.
# -----------------------------------------------------------------------------
cal_rho = function(seed, GH_ge, GH_spe){
  gh = c(22,23,32,33,34,43,44)
  seed_file = matching_files[seed_vec == seed]
  rho_g_ge = c()
  rho_h_ge = c()
  rho_g_spe = c()
  rho_h_spe = c()
  for (f in seed_file) {
    load(f)
    rho_g_ge = c(rho_g_ge, err_rate_mapping(res_ge_update$member_G,member0 = member_G0))
    rho_g_spe = c(rho_g_spe, err_rate_mapping(res_update$member_G,member0 = member_G0))
    rho_h_ge = c(rho_h_ge, err_rate_mapping(res_ge_update$member_H,member0 = member_H0))
    rho_h_spe = c(rho_h_spe, err_rate_mapping(res_update$member_H,member0 = member_H0))
  }
  rho_g_ge = c(0,rho_g_ge, rho_g_ge[which(gh == GH_ge)])
  rho_g_spe = c(0,rho_g_spe, rho_g_spe[which(gh == GH_spe)])
  rho_h_ge = c(0,rho_h_ge, rho_h_ge[which(gh == GH_ge)])
  rho_h_spe = c(0,rho_h_spe, rho_h_spe[which(gh == GH_spe)])
  res = list()
  res$rho_g_ge = rho_g_ge
  res$rho_g_spe = rho_g_spe
  res$rho_h_ge = rho_h_ge
  res$rho_h_spe = rho_h_spe
  return(res)
}

# -----------------------------------------------------------------------------
# Allocate the table that will eventually be printed as LaTeX. 5 designs x 9
# rows per design (Oracle + 7 candidates + the QIC pick) = 45 rows. Columns
# split into:
#   1-3:  design info (N, T, method label)
#   4-10: special-model results  (alpha, beta, nu, gamma, MSR, rho_g, rho_h)
#   11-16: general-model results (theta, nu, gamma, MSR, rho_g, rho_h)
# (Column 5 — beta — is later dropped before printing.)
# -----------------------------------------------------------------------------
res_over = matrix(NA, 45, 16)
res_over = as.data.frame(res_over)
names(res_over) = c("N", "T", "method", "alpha", "beta", "nu", "gamma", "MSR", "rho_g","rho_h",
                    "theta", "nu", "gamma", "MSR", "rho_g","rho_h")
res_over[c(1,10,19,28,37),1] = c(100,100,200,200,100)
res_over[c(1,10,19,28,37),2] = c(100,200,200,400,50)
method = c("Oracle", "GTNQ-22","GTNQ-23","GTNQ-32","GTNQ-33","GTNQ-34","GTNQ-43","GTNQ-44","GTNQ-wh G wh H")
res_over[,3] = rep(method,5)


# -----------------------------------------------------------------------------
# Main loop: for each design `case`, aggregate candidate results into res_over.
# -----------------------------------------------------------------------------
for(case in 1:5){
  
  # --- Design parameters and corresponding output rows ----------------------
  if(case==1) {
    N=100
    Time=100
    Nblock = 5
  }
  if(case==2) {
    N=100
    Time=200
    Nblock = 5
  }
  if(case==3) {
    N=200
    Time=200
    Nblock = 10
  }
  if(case==4) {
    N=200
    Time=400
    Nblock = 10
  }
  if(case == 5){
    N=100
    Time=50
    Nblock = 5
  }

  # --- Match files by leading digit (case) and pick QIC penalty lambda ------
  # lambda = N^0.1 / (T - 1) / c / kappa * log(T - 1), with kappa depending
  # on N (8.34 for N = 100, 7.44 for N = 200).
  if(case == 1){
    matching_files <- list.files(path = folder_path, pattern = "^1", full.names = TRUE)
    matching_files = sort(matching_files)
    lambda =  N^0.1/(Time-1)/c/8.34*log(Time-1)
    rows = 1:9
  }
  if(case == 2){
    matching_files <- list.files(path = folder_path, pattern = "^2", full.names = TRUE)
    matching_files = sort(matching_files)
    lambda =  N^0.1/(Time-1)/c/8.34*log(Time-1)
    rows = 10:18
  }
  if(case == 3){
    matching_files <- list.files(path = folder_path, pattern = "^3", full.names = TRUE)
    matching_files = sort(matching_files)
    lambda =  N^0.1/(Time-1)/c/7.44*log(Time-1)
    rows = 19:27
  }
  if(case == 4){
    matching_files <- list.files(path = folder_path, pattern = "^4", full.names = TRUE)
    matching_files = sort(matching_files)
    lambda =  N^0.1/(Time-1)/c/7.44*log(Time-1)
    rows = 28:36
  }
  if(case == 5){
    matching_files <- list.files(path = folder_path, pattern = "^5", full.names = TRUE)
    matching_files = sort(matching_files)
    lambda =  N^0.1/(Time-1)/c/8.34*log(Time-1)
    rows = 37:45
  }
  
  # --- Identify seeds for which all 7 candidate files exist -----------------
  seed_vec <- sub(".+-(\\d+).+", "\\1", matching_files)
  tb = table(seed_vec)
  seed_success = names(tb[tb == 7])
  file_vec = matching_files[seed_vec %in% seed_success]
  
  # --- General-model QIC selection and selection-frequency table ------------
  res_ge = lapply(seed_success, function(i)result_ge(i))
  GH_vec_ge = unlist(lapply(res_ge, function(r)selectGH(r,lambda)))
  # Append the candidate-grid labels plus 0 and 5 to ensure all 9 categories
  # show up in the frequency table; subtract 1 to undo this padding.
  temp = c(GH_vec_ge, 22,23,32,33,34,43,44,0,5)
  tb = table(temp)-1
  tb = tb / sum(tb)
  tb = tb[order(names(tb))]
  res_over[rows,14] = tb
  
  
  # --- Special-model QIC selection and selection-frequency table ------------
  res_spe = lapply(seed_success, function(i)result_spe(i))
  GH_vec_spe = unlist(lapply(res_spe, function(r)selectGH(r,lambda)))
  temp = c(GH_vec_spe, 22,23,32,33,34,43,44,0,5)
  tb = table(temp)-1
  tb = tb / sum(tb)
  tb = tb[order(names(tb))]
  res_over[rows,8] = tb
  
  # --- Clustering error rates (averages over seeds) -------------------------
  rho_list = lapply(1:length(seed_success), function(i)cal_rho(seed_success[i], GH_vec_ge[i], GH_vec_spe[i]))
  
  rho_g_ge = lapply(rho_list,function(a)a$rho_g_ge)
  rho_g_ge = do.call(rbind,rho_g_ge)
  rho_g_ge = colMeans(rho_g_ge)
  res_over[rows,15] = rho_g_ge
  
  rho_h_ge = lapply(rho_list,function(a)a$rho_h_ge)
  rho_h_ge = do.call(rbind,rho_h_ge)
  rho_h_ge = colMeans(rho_h_ge)
  res_over[rows,16] = rho_h_ge
  
  rho_g_spe = lapply(rho_list,function(a)a$rho_g_spe)
  rho_g_spe = do.call(rbind,rho_g_spe)
  rho_g_spe = colMeans(rho_g_spe)
  res_over[rows,9] = rho_g_spe
  
  rho_h_spe = lapply(rho_list,function(a)a$rho_h_spe)
  rho_h_spe = do.call(rbind,rho_h_spe)
  rho_h_spe = colMeans(rho_h_spe)
  res_over[rows,10] = rho_h_spe
  
  # --- RMSEs (averages over seeds) for candidates + QIC pick ---------------
  rmse_list = lapply(1:length(seed_success), function(i)cal_RMSE(seed_success[i], GH_vec_ge[i], GH_vec_spe[i]))
  rmse_list = do.call(rbind, rmse_list)
  rmse = colMeans(rmse_list)
  rmse = t(matrix(rmse,7,8)) # 8 rows x 7 cols (cand + QIC pick)
  res_over[rows[-1], c(4:7,11:13)] = rmse
}

# ---- Format the table: percentages with one decimal place, blank NAs ------
res_over[,4:16] = format(round(res_over[,4:16]*100,1),nsmall = 1)
res_over[is.na(res_over)]=""
res_over[res_over=="NA"]=""



# =============================================================================
# Oracle results: fill the "Oracle" row of each design block.
# =============================================================================
if(power){
  if(additive){
    folder_path <- "output/power_select_or_add/"
  }else{
    folder_path <- "output/power_select_or_mul/"
  }
}else{
  if(additive){
    folder_path <- "output/block_select_or_add/"
  }else{
    folder_path <- "output/block_select_or_mul/"
  }
}


file_names <- list.files(path = folder_path, full.names = TRUE)

# -----------------------------------------------------------------------------
# cal_RMSE() for one oracle file. Returns a length-7 vector of average
# |est - true| over the 5 quantile levels (same layout as in the candidate
# version above):
#   [1] alpha           (special, alpha+beta combined into alphabeta)
#   [2] (unused)
#   [3] nu              (special)
#   [4] gamma           (special)
#   [5] alphabeta       (general)
#   [6] nu              (general)
#   [7] gamma           (general)
# -----------------------------------------------------------------------------
cal_RMSE = function(f){
  load(f)
  rmse = numeric(7)
  for (q in 1:5) {
    # General-model oracle errors
    theta_all_ge = or_ge$alphabeta_GHs[[q]][member_G0,member_H0]
    theta_all0 = alphabeta0[[q]][member_G0, member_H0]
    rmse[5] = rmse[5] + mean(abs(theta_all_ge - theta_all0))
    
    columns <- lapply(or_ge$theta_Gs, function(mat) mat[, 1])
    nu <- do.call(cbind,columns)
    nu_all_ge = nu[,q][member_G0]
    nu_all0 = nu0[,q][member_G0]
    rmse[6] = rmse[6] + mean(abs(nu_all_ge - nu_all0))
    
    gamma = lapply(or_ge$theta_Gs, function(mat) mat[, 2:3])
    gamma_all_ge = gamma[[q]][member_G0,]
    gamma_all0 = gamma0[[q]][member_G0,]
    r = sqrt(rowSums((gamma_all_ge - gamma_all0)^2))
    rmse[7] = rmse[7] + mean(r)
    
    # Special-model oracle errors
    alpha = or_spe$alpha
    beta = or_spe$beta
    alpha_all = alpha[,q][member_G0]
    beta_all = beta[,q][member_H0]
    if(additive){
      alphabeta_all = outer(alpha_all, rep(1, length(alpha_all))) + outer(rep(1,length(beta_all)), beta_all)
    }else{
      alphabeta_all = outer(alpha_all, rep(1, length(alpha_all))) * outer(rep(1,length(beta_all)), beta_all)
    }
    rmse[1] = rmse[1] +mean(abs(alphabeta_all - theta_all0))

    
    columns <- lapply(or_spe$theta_Gs, function(mat) mat[, 2])
    nu <- do.call(cbind,columns)
    nu_all_ge = nu[,q][member_G0]
    nu_all0 = nu0[,q][member_G0]
    rmse[3] = rmse[3] + mean(abs(nu_all_ge - nu_all0))
    
    gamma = lapply(or_spe$theta_Gs, function(mat) mat[, 3:4])
    gamma_all_ge = gamma[[q]][member_G0,]
    gamma_all0 = gamma0[[q]][member_G0,]
    r = sqrt(rowSums((gamma_all_ge - gamma_all0)^2))
    rmse[4] = rmse[4] + mean(r)
  }
  rmse = rmse /5
  return(rmse)
}

# Fill the Oracle row for each design.
for(case in 1:5){
  if(case == 1){
    matching_files <- list.files(path = folder_path, pattern = "^1", full.names = TRUE)
    row = 1
  }
  if(case == 2){
    matching_files <- list.files(path = folder_path, pattern = "^2", full.names = TRUE)
    row = 10
  }
  if(case == 3){
    matching_files <- list.files(path = folder_path, pattern = "^3", full.names = TRUE)
    row = 19
  }
  if(case == 4){
    matching_files <- list.files(path = folder_path, pattern = "^4", full.names = TRUE)
    row = 28
  }
  if(case == 5){
    matching_files <- list.files(path = folder_path, pattern = "^5", full.names = TRUE)
    row = 37
  }
  res = lapply(matching_files, function(f) cal_RMSE(f))
  res = do.call(rbind,res)
  res = colMeans(res)
  res_over[row,c(4:7,11:13)] = format(round(res*100,1),nsmall = 1)
  res_over[row,c(8:10,14:16)] = "-"
  res_over[row+8,c(8,14)] = "-"
  
}

# Some cells are repeated (defensive). Make sure they stay as "-".
for(case in 1:5){
  if(case == 1){
    row = 1
  }
  if(case == 2){
    row = 10
  }
  if(case == 3){
    row = 19
  }
  if(case == 4){
    row = 28
  }
  if(case == 5){
    row = 37
  }
  res_over[row,c(8:10,14:16)] = "-"
  res_over[row+8,c(8,14)] = "-"
  
}

# Reorder: move case = 5 (N=100, T=50) to the top of the table, then drop the
# unused "beta" column (column 5).
res_over=res_over[c(37:45,1:36),-5]

print(res_over)




# -----------------------------------------------------------------------------
# Print the LaTeX table.
# -----------------------------------------------------------------------------
library(xtable)
df = as.data.frame(res_over)
latex_table <- xtable(df)
print(latex_table, include.rownames = FALSE)


