# =============================================================================
# assess.R
# -----------------------------------------------------------------------------
# Utility functions for evaluating estimation quality in the TGNQ simulations
# and for reshaping result vectors into tidy matrices used to build LaTeX
# tables.
#
# Two groups of functions are provided:
#
#   (A) Per-replication evaluation
#       cal_para            -- additive / multiplicative TGNQ model
#       cal_para_general    -- general TGNQ model
#       cal_para_mis*       -- variants used in misspecification experiments
#
#       Each returns a long numeric vector that concatenates
#         (i)  RMSE of each parameter block at five quantile levels, and
#         (ii) 0/1 indicators of whether the truth lies inside the 95% CI.
#
#   (B) Table formatting
#       cal_res             -- reshapes the long result vector into a matrix
#                              and extracts the clustering-error rates.
# =============================================================================


# -----------------------------------------------------------------------------
# cal_para
# Evaluation routine for the structured (additive / multiplicative) TGNQ model.
#
# Arguments
#   res     : estimated parameters from one replication, with fields
#               $alpha, $beta, $theta_Gs (list with columns [, 1]=theta,
#                                                          [, 2]=nu,
#                                                          [, 3:4]=gamma)
#   ci      : list with $up and $low, each having the same structure as `res`,
#             representing upper / lower bounds of the 95% CI.
#   alpha0, beta0, nu0, gamma0 : true parameter values.
#
# Returns
#   A numeric vector concatenating
#     RMSE(alpha), RMSE(beta), RMSE(nu), RMSE(gamma),
#     CI coverage indicators for alpha, beta, nu, gamma.
# -----------------------------------------------------------------------------

cal_para = function(res, ci, alpha0, beta0, nu0, gamma0){
  # ---- extract estimates from `res` ----
  alpha = res$alpha
  beta = res$beta
  second_columns <- lapply(res$theta_Gs, function(mat) mat[, 2])
  nu <- do.call(cbind, second_columns) # G x K matrix
  gamma = lapply(res$theta_Gs, function(mat) mat[, 3:4]) # one G x p block per quantile
  
  # ---- RMSE for each parameter block ----
  rmse_alpha = sqrt(colSums((alpha - alpha0)^2))
  rmse_beta = sqrt(colSums((beta - beta0)^2))
  rmse_nu = sqrt(colSums((nu-nu0)^2))
  rmse_gamma = sqrt(sapply(1:5, function(i) sum((gamma[[i]] - gamma0[[i]])^2)))
  
  # ---- CI upper bounds ----
  alpha_up = ci$up$alpha
  beta_up = ci$up$beta
  second_columns_up <- lapply(ci$up$theta_Gs, function(mat) mat[, 2])
  nu_up <- do.call(cbind, second_columns_up)
  gamma_up= lapply(ci$up$theta_Gs, function(mat) mat[, 3:4])
  
  # ---- CI lower bounds ----
  alpha_low = ci$low$alpha
  beta_low = ci$low$beta
  second_columns_low <- lapply(ci$low$theta_Gs, function(mat) mat[, 2])
  nu_low <- do.call(cbind, second_columns_low)
  gamma_low = lapply(ci$low$theta_Gs, function(mat) mat[, 3:4])
  
  # ---- coverage indicators: 1 if true value lies within [low, up] ----
  ci_t_alpha = as.numeric((alpha0 >= alpha_low) & (alpha0 <= alpha_up))
  ci_t_beta = as.numeric((beta0 >= beta_low) & (beta0 <= beta_up))
  ci_t_nu = as.numeric((nu0 >= nu_low) & (nu0 <= nu_up))
  ci_t_gamma = lapply(1:5, function(i){as.numeric(gamma0[[i]] >= (gamma_low[[i]]) & ((gamma0[[i]] <=(gamma_up[[i]]))))})
  ci_t_gamma = unlist(ci_t_gamma)
  
  return(c(rmse_alpha, rmse_beta, rmse_nu, rmse_gamma, ci_t_alpha, ci_t_beta, ci_t_nu, ci_t_gamma))
}

# -----------------------------------------------------------------------------
# cal_para_general
# Evaluation routine for the GENERAL TGNQ model.
#
# Here the interaction parameters are stored as full G x H matrices per
# quantile (in `res$alphabeta_GHs`). The "true" interaction matrices are
# reconstructed from (alpha0, beta0) depending on whether the DGP is additive.
# -----------------------------------------------------------------------------
cal_para_general = function(res, ci, alpha0, beta0, nu0, gamma0, additive){
  
  # ---- extract estimates ----
  alphabeta = res$alphabeta_GHs # list of G x H matrices
  columns <- lapply(res$theta_Gs, function(mat) mat[, 1])
  nu <- do.call(cbind,columns) # G x K
  gamma = lapply(res$theta_Gs, function(mat) mat[, 2:3])
  
  # ---- build true interaction matrices ----
  if(additive){
    # theta_{gh} = alpha_g + beta_h
    alphabeta0 <- list()
    for (i in 1:5) {
      alphabeta0[[i]] <- outer(alpha0[, i],rep(1,nrow(alphabeta[[1]])))   +  outer(rep(1,nrow(alphabeta[[1]])),beta0[, i])
    }
  }else{
    # theta_{gh} = alpha_g * beta_h
    alphabeta0 <- list()
    for (i in 1:5) {
      alphabeta0[[i]] <- outer(alpha0[, i],beta0[, i])
    }
  }
  
  # ---- RMSEs ----
  rmse_alphabeta = sqrt(sapply(1:5, function(i) sum((alphabeta[[i]] - alphabeta0[[i]])^2)))
  rmse_nu = sqrt(colSums((nu-nu0)^2))
  rmse_gamma = sqrt(sapply(1:5, function(i) sum((gamma[[i]] - gamma0[[i]])^2)))
  
  # ---- CI bounds for the interaction matrices ----
  alphabeta_up = ci$up$alphabeta_GHs
  columns_up <- lapply(ci$up$theta_Gs, function(mat) mat[, 1])
  nu_up <- do.call(cbind, columns_up)
  gamma_up= lapply(ci$up$theta_Gs, function(mat) mat[, 2:3])
  
  alphabeta_low = ci$low$alphabeta_GHs
  columns_low <- lapply(ci$low$theta_Gs, function(mat) mat[, 1])
  nu_low <- do.call(cbind, columns_low)
  gamma_low = lapply(ci$low$theta_Gs, function(mat) mat[, 2:3])
  
  # ---- coverage indicators ----
  ci_t_alphabeta = lapply(1:5, function(i){as.numeric((alphabeta0[[i]] >= alphabeta_low[[i]]) & (alphabeta0[[i]] <= alphabeta_up[[i]]))})
  ci_t_alphabeta = unlist(ci_t_alphabeta)
  ci_t_nu = as.numeric(t((nu0 >= nu_low) & (nu0 <= nu_up)))
  ci_t_gamma = lapply(1:5, function(i){as.numeric((gamma0[[i]] >= gamma_low[[i]]) & (gamma0[[i]] <= gamma_up[[i]]))})
  ci_t_gamma = unlist(ci_t_gamma)
  
  return(c(rmse_alphabeta, rmse_nu, rmse_gamma, ci_t_alphabeta, ci_t_nu, ci_t_gamma))
}

# -----------------------------------------------------------------------------
# cal_para_mis
# Variant used in the misspecification experiments (Section 5.3) when the
# WORKING model is the structured (alpha/beta) form but the TRUE DGP may use
# the OPPOSITE interaction form.
#
# Specifically:
#   * if additive == TRUE,  the truth is additive   and the working model is
#                            misspecified as multiplicative;
#   * if additive == FALSE, the truth is multiplicative and the working model
#                            is misspecified as additive.
#
# Returns RMSE of the implied G x H interaction matrices, plus RMSE of nu and
# gamma. No CI information is returned because it is not used in the
# misspecification tables.
# -----------------------------------------------------------------------------
cal_para_mis = function(res, ci, alpha0, beta0, nu0, gamma0, additive){
  # ---- extract working-model estimates ----
  alpha = res$alpha
  beta = res$beta
  second_columns <- lapply(res$theta_Gs, function(mat) mat[, 2])
  nu <- do.call(cbind, second_columns)
  gamma = lapply(res$theta_Gs, function(mat) mat[, 3:4])
  
  # ---- build the implied interaction matrices under the OPPOSITE structure
  # NOTE: the branches below match the misspecified working model used in 5.3
  if(additive){
    alphabeta <- list()
    for (i in 1:5) {
      alphabeta[[i]] <- outer(alpha[, i],beta[, i])
    }
  }else{
    alphabeta <- list()
    for (i in 1:5) {
      alphabeta[[i]] <- outer(alpha[, i],rep(1,nrow(alpha)))   +  outer(rep(1,nrow(alpha)),beta[, i])
    }
  }
  
  
  # ---- true interaction matrices ----
  if(additive){
    alphabeta0 <- list()
    for (i in 1:5) {
      alphabeta0[[i]] <- outer(alpha0[, i],rep(1,nrow(alphabeta[[1]])))   +  outer(rep(1,nrow(alphabeta[[1]])),beta0[, i])
    }
  }else{
    alphabeta0 <- list()
    for (i in 1:5) {
      alphabeta0[[i]] <- outer(alpha0[, i],beta0[, i])
    }
  }
  rmse_alphabeta = sqrt(sapply(1:5, function(i) sum((alphabeta[[i]] - alphabeta0[[i]])^2)))
  rmse_nu = sqrt(colSums((nu-nu0)^2))
  rmse_gamma = sqrt(sapply(1:5, function(i) sum((gamma[[i]] - gamma0[[i]])^2)))
  

  return(c(rmse_alphabeta, rmse_nu, rmse_gamma))
}



# -----------------------------------------------------------------------------
# cal_para_mis_T
# Same as cal_para_mis() but with the additive / multiplicative branches
# swapped. Used when the WORKING model is specified to match the TRUTH while
# the comparison target is the OTHER structure (for Table 2 in Section 5.3).
# -----------------------------------------------------------------------------
cal_para_mis_T = function(res, ci, alpha0, beta0, nu0, gamma0, additive){
  alpha = res$alpha
  beta = res$beta
  second_columns <- lapply(res$theta_Gs, function(mat) mat[, 2])
  nu <- do.call(cbind, second_columns)
  gamma = lapply(res$theta_Gs, function(mat) mat[, 3:4])
  
  # branches swapped relative to cal_para_mis
  if(!additive){
    alphabeta <- list()
    for (i in 1:5) {
      alphabeta[[i]] <- outer(alpha[, i],beta[, i])
    }
  }else{
    alphabeta <- list()
    for (i in 1:5) {
      alphabeta[[i]] <- outer(alpha[, i],rep(1,nrow(alpha)))   +  outer(rep(1,nrow(alpha)),beta[, i])
    }
  }
  
  if(additive){
    alphabeta0 <- list()
    for (i in 1:5) {
      alphabeta0[[i]] <- outer(alpha0[, i],rep(1,nrow(alphabeta[[1]])))   +  outer(rep(1,nrow(alphabeta[[1]])),beta0[, i])
    }
  }else{
    alphabeta0 <- list()
    for (i in 1:5) {
      alphabeta0[[i]] <- outer(alpha0[, i],beta0[, i])
    }
  }
  rmse_alphabeta = sqrt(sapply(1:5, function(i) sum((alphabeta[[i]] - alphabeta0[[i]])^2)))
  rmse_nu = sqrt(colSums((nu-nu0)^2))
  rmse_gamma = sqrt(sapply(1:5, function(i) sum((gamma[[i]] - gamma0[[i]])^2)))
  
  
  return(c(rmse_alphabeta, rmse_nu, rmse_gamma))
}



# -----------------------------------------------------------------------------
# cal_para_mis_ge
# Variant of cal_para_mis() used when the WORKING model is the GENERAL TGNQ
# model. The estimated interaction matrix is taken directly from
# res$alphabeta_GHs; the true matrices are still rebuilt from (alpha0, beta0).
# -----------------------------------------------------------------------------
cal_para_mis_ge = function(res, ci, alpha0, beta0, nu0, gamma0, additive){
  alphabeta = res$alphabeta_GHs
  columns <- lapply(res$theta_Gs, function(mat) mat[, 1])
  nu <- do.call(cbind,columns)
  gamma = lapply(res$theta_Gs, function(mat) mat[, 2:3])
  
  if(additive){
    alphabeta0 <- list()
    for (i in 1:5) {
      alphabeta0[[i]] <- outer(alpha0[, i],rep(1,nrow(alphabeta[[1]])))   +  outer(rep(1,nrow(alphabeta[[1]])),beta0[, i])
    }
  }else{
    alphabeta0 <- list()
    for (i in 1:5) {
      alphabeta0[[i]] <- outer(alpha0[, i],beta0[, i])
    }
  }
  rmse_alphabeta = sqrt(sapply(1:5, function(i) sum((alphabeta[[i]] - alphabeta0[[i]])^2)))
  rmse_nu = sqrt(colSums((nu-nu0)^2))
  rmse_gamma = sqrt(sapply(1:5, function(i) sum((gamma[[i]] - gamma0[[i]])^2)))
  
  return(c(rmse_alphabeta, rmse_nu, rmse_gamma))
}



# -----------------------------------------------------------------------------
# cal_res
# Reshape the long Monte Carlo result vector `res_vec` (concatenation of RMSEs,
# CI coverage proportions, and clustering errors across five quantile levels)
# into a tidy matrix that the `step*_output_*_tex.R` scripts can consume.
#
# The function recognizes the result-vector length to dispatch among the
# different simulation settings (G=2 vs G=3, with vs without CI, etc.) and
# performs the corresponding bookkeeping:
#   * the first columns of the returned matrix store RMSE blocks at each
#     quantile,
#   * the remaining columns store the *average* coverage error |coverage - 0.95|
#     for each parameter block,
#   * the clustering-error rates (rho values) are returned separately in
#     `res_list$rho`.
#
# Returns
#   A list with components
#     $rho     : clustering-error rates for row and column memberships,
#     $res_mat : the assembled result matrix.
# -----------------------------------------------------------------------------
cal_res = function(res_vec){
  
  # Each `if (length(res_vec) == ...)` branch below corresponds to a particular
  # simulation configuration. The numerical constants encode the exact ordering
  # of blocks produced by the simulation drivers; do NOT modify them unless the
  # corresponding simulation script changes its output layout.
  if(length(res_vec) == 244){
    res_mat = matrix(NA, 5,24)
    rho_vec = res_vec[c(161:162,243:244)]
    temp_vec =  res_vec[-c(161:162,243:244)]
    col = 1
    while (length(temp_vec)>0) {
      # 4 RMSE columns (one per parameter block) at the current quantile
      res_mat[1:5, col:(col+3)] = temp_vec[1:20]
      temp_vec = temp_vec[-(1:20)]
      col = col + 4
      # average CI-coverage error for three "scalar" blocks
      for (i in 1:3) {
        res_mat[1:5, col] = colMeans(matrix(abs(temp_vec[1:10]-0.95),ncol = 5))
        temp_vec = temp_vec[-(1:10)]
        col=col+1
      }
      # average CI-coverage error for the "vector" block (gamma)
      res_mat[1:5, col] = colMeans(matrix(abs(temp_vec[1:30]-0.95),ncol = 5))
      temp_vec = temp_vec[-(1:30)]
      col=col+1
    }
  }
  
  if(length(res_vec) == 334){
    res_mat = matrix(NA, 5,24)
    rho_vec = res_vec[c(221:222,333:334)]
    temp_vec =  res_vec[-c(221:222,333:334)]
    col = 1
    while (length(temp_vec)>0) {
      res_mat[1:5, col:(col+3)] = temp_vec[1:20]
      temp_vec = temp_vec[-(1:20)]
      col = col + 4
      for (i in 1:3) {
        res_mat[1:5, col] = colMeans(matrix(abs(temp_vec[1:15]-0.95),ncol = 5))
        temp_vec = temp_vec[-(1:15)]
        col=col+1
      }
      res_mat[1:5, col] = colMeans(matrix(abs(temp_vec[1:45]-0.95),ncol = 5))
      temp_vec = temp_vec[-(1:45)]
      col=col+1
    }
  }
  
  if(length(res_vec) == 229){
    res_mat = matrix(NA, 5,18)
    rho_vec = res_vec[c(151:152,228:229)]
    temp_vec =  res_vec[-c(151:152,228:229)]
    col = 1
    while (length(temp_vec)>0) {
      res_mat[1:5, col:(col+2)] = temp_vec[1:15]
      temp_vec = temp_vec[-(1:15)]
      col = col + 3
      
      res_mat[1:5, col] = colMeans(matrix(abs(temp_vec[1:20]-0.95),ncol = 5))
      temp_vec = temp_vec[-(1:20)]
      col=col+1
      
      res_mat[1:5, col] = colMeans(matrix(abs(temp_vec[1:10]-0.95),ncol = 5))
      temp_vec = temp_vec[-(1:10)]
      col=col+1
      
      res_mat[1:5, col] = colMeans(matrix(abs(temp_vec[1:30]-0.95),ncol = 5))
      temp_vec = temp_vec[-(1:30)]
      col=col+1
    }
  }
  
  if(length(res_vec) == 364){
    res_mat = matrix(NA, 5,18)
    rho_vec = res_vec[c(241:242,363:364)]
    temp_vec =  res_vec[-c(241:242,363:364)]
    col = 1
    while (length(temp_vec)>0) {
      res_mat[1:5, col:(col+2)] = temp_vec[1:15]
      temp_vec = temp_vec[-(1:15)]
      col = col + 3
      
      res_mat[1:5, col] = colMeans(matrix(abs(temp_vec[1:45]-0.95),ncol = 5))
      temp_vec = temp_vec[-(1:45)]
      col=col+1
      
      res_mat[1:5, col] = colMeans(matrix(abs(temp_vec[1:15]-0.95),ncol = 5))
      temp_vec = temp_vec[-(1:15)]
      col=col+1
      
      res_mat[1:5, col] = colMeans(matrix(abs(temp_vec[1:45]-0.95),ncol = 5))
      temp_vec = temp_vec[-(1:45)]
      col=col+1
    }
  }
  
  
  if(length(res_vec) ==214){
    res_mat = matrix(NA, 5,24)
    rho_vec = res_vec[c(141:142,213:214)]
    temp_vec =  res_vec[-c(141:142,213:214)]
    col = 1
    while (length(temp_vec)>0) {
      res_mat[1:5, col:(col+3)] = temp_vec[1:20]
      temp_vec = temp_vec[-(1:20)]
      col = col + 4
      for (i in 1:3) {
        res_mat[1:5, col] = colMeans(matrix(abs(temp_vec[1:10]-0.95),ncol = 5))
        temp_vec = temp_vec[-(1:10)]
        col=col+1
      }
      res_mat[1:5, col] = colMeans(matrix(abs(temp_vec[1:20]-0.95),ncol = 5))
      temp_vec = temp_vec[-(1:20)]
      col=col+1
    }
  }
  
  if(length(res_vec) ==289){
    res_mat = matrix(NA, 5,24)
    rho_vec = res_vec[c(191:192,288:289)]
    temp_vec =  res_vec[-c(191:192,288:289)]
    col = 1
    while (length(temp_vec)>0) {
      res_mat[1:5, col:(col+3)] = temp_vec[1:20]
      temp_vec = temp_vec[-(1:20)]
      col = col + 4
      for (i in 1:3) {
        res_mat[1:5, col] = colMeans(matrix(abs(temp_vec[1:15]-0.95),ncol = 5))
        temp_vec = temp_vec[-(1:15)]
        col=col+1
      }
      res_mat[1:5, col] = colMeans(matrix(abs(temp_vec[1:30]-0.95),ncol = 5))
      temp_vec = temp_vec[-(1:30)]
      col=col+1
    }
  }
  
  if(length(res_vec) == 319){
    res_mat = matrix(NA, 5,18)
    rho_vec = res_vec[c(211:212,318:319)]
    temp_vec =  res_vec[-c(211:212,318:319)]
    col = 1
    while (length(temp_vec)>0) {
      res_mat[1:5, col:(col+2)] = temp_vec[1:15]
      temp_vec = temp_vec[-(1:15)]
      col = col + 3
      
      res_mat[1:5, col] = colMeans(matrix(abs(temp_vec[1:45]-0.95),ncol = 5))
      temp_vec = temp_vec[-(1:45)]
      col=col+1
      
      res_mat[1:5, col] = colMeans(matrix(abs(temp_vec[1:15]-0.95),ncol = 5))
      temp_vec = temp_vec[-(1:15)]
      col=col+1
      
      res_mat[1:5, col] = colMeans(matrix(abs(temp_vec[1:30]-0.95),ncol = 5))
      temp_vec = temp_vec[-(1:30)]
      col=col+1
    }
  }
  
  if(length(res_vec) == 199){
    res_mat = matrix(NA, 5,18)
    rho_vec = res_vec[c(131:132,198:199)]
    temp_vec =  res_vec[-c(131:132,198:199)]
    col = 1
    while (length(temp_vec)>0) {
      res_mat[1:5, col:(col+2)] = temp_vec[1:15]
      temp_vec = temp_vec[-(1:15)]
      col = col + 3
      
      res_mat[1:5, col] = colMeans(matrix(abs(temp_vec[1:20]-0.95),ncol = 5))
      temp_vec = temp_vec[-(1:20)]
      col=col+1
      
      res_mat[1:5, col] = colMeans(matrix(abs(temp_vec[1:10]-0.95),ncol = 5))
      temp_vec = temp_vec[-(1:10)]
      col=col+1
      
      res_mat[1:5, col] = colMeans(matrix(abs(temp_vec[1:20]-0.95),ncol = 5))
      temp_vec = temp_vec[-(1:20)]
      col=col+1
    }
  }
  
  res_list = list()
  res_list$rho = rho_vec
  res_list$res_mat = res_mat
  return(res_list)
}












