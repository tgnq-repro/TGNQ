# =============================================================================
# simulator.R
# -----------------------------------------------------------------------------
# Defines simulation(seed), one Monte Carlo replication shared by Sections 5.1
# and 5.3. The function:
#
#   1. Sets (N, T, Nblock) from `case`.
#   2. Builds a (fixed across seeds) block or power-law network W.
#   3. Generates true memberships and the response panel Y from the TGNQ DGP.
#   4. Fits FIVE estimators and records their RMSE / coverage / clustering
#      errors against the truth:
#        - oracle GE / SPE   (true memberships, "general" / correct method)
#        - TGNQ-GE / SPE     (estimated memberships, both methods)
#        - TGNQ-MIS          (estimated memberships, WRONG method)
#        - refined GE / SPE  (one Refine_G + Refine_H pass on the above)
#        - "transferred" MIS errors at the correctly-specified memberships.
#   5. Saves a per-seed snapshot .rda and returns a 4-element list of error
#      vectors for downstream aggregation.
#
# Expected globals (set in simulation.R): case, additive, block, G0,
# method, method_mis.
# =============================================================================

simulation = function(seed){
  # ---- (1) Set (N, T, Nblock) according to `case` --------------------------
  H0=G0
  if(case==1) {
    N=100
    Time=50
    Nblock = 5
  }
  if(case==2) {
    N=100
    Time=100
    Nblock = 5
  }
  if(case==3) {
    N=100
    Time=200
    Nblock = 5
  }
  if(case==4) {
    N=200
    Time=200
    Nblock = 10
  }
  if(case==5) {
    N=200
    Time=400
    Nblock = 10
  }
  
  taus = c(0.1, 0.3, 0.5, 0.7, 0.9) # quantile levels evaluated
  n_taus = length(taus) 
  
  
  # ---- (2) Build the (fixed) network ---------------------------------------
  set.seed(123)
  p = 2
  # 设置网络
  if(block){
    Amat = getBlockW(N,Nblock, normalize = F)
  }else{
    Amat = getPowerLawW(N)
  }
  
  W = Amat[1:N,1:N]
  rw = rowSums(W)
  W[rw==0, sample(1:N, 3)] = 1
  W = W/rowSums(W)
  
  FriendW <-FriendW2 <- vector("list", N)
  for(i in 1:N){
    FriendW[[i]] <- as.vector(which(W[,i]!=0))
    for (j in  FriendW[[i]]) {
      FriendW2[[i]] <- unique(c(FriendW2[[i]],as.vector(which(W[j,]!=0))))
    }
  }
  
  
  # ---- (3) True group memberships ------------------------------------------
  # G0 = 2 -> ratios (0.5, 0.5) for row groups, (0.4, 0.6) for column groups.
  # G0 = 3 -> ratios (0.3, 0.3, 0.4) and (0.4, 0.3, 0.3).
  if (G0==2){
    member_G0 = simu.member0(N, ratios = c(0.5, 0.5))
    member_H0 = simu.member0(N, ratios = c(0.4, 0.6))
  }else{
    member_G0 = simu.member0(N, ratios = c(0.3, 0.3, 0.4))
    member_H0 = simu.member0(N, ratios = c(0.4, 0.3, 0.3))
  }
  
  
  # ---- (4) True parameter matrices (one column per quantile) ---------------
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
  
  # ---- (5) Pre-allocate output vectors -------------------------------------
  # Layout:
  #   res_vec_ge = [oracle | gtnq | refine] x (15 + G0*H0*5 + G0*3*5) entries
  #              + (rho_G, rho_H) for gtnq and refine = 4 trailing entries.
  #   res_vec    = [oracle | gtnq | refine] x (20 + G0*5*5) entries
  #              + (rho_G, rho_H) for gtnq and refine = 4 trailing entries.
  sim_times = 0
  res_vec_ge = numeric((15 + (G0*H0*5 + G0*3*5))*3 + 4)
  res_vec = numeric((20 + G0 * 5 * 5 )*3 + 4)
  
  # ---- (6) Run ONE replication; retry on numerical failure -----------------
  while (sim_times < 1) {
    tryCatch({
      # Reseed and bump by 500 so any retry uses a fresh draw.
      set.seed(seed)
      seed = seed + 500

      # ------ (6.1) Simulate covariates and response panel ------------------
      X_tensor = array(rnorm(N*p*(Time+10)), dim = c(N, p, (Time+10)))
      X_tensor = abs(X_tensor)
      G = max(member_G0); H = max(member_H0)
      set.seed(seed)
      Ymat = simu.Y(W, member_G0, member_H0, X_tensor, additive = additive, verbose = F)
      X_tensor = X_tensor[,,-(1:11)]
      
      
      # ============================================================
      # (6.2) ORACLE estimators (true memberships, general + special)
      # ============================================================
      ci_ge_or = twmq_ci(Ymat, X_tensor, W,
                                                         member_G0, member_H0,
                                                         taus, conquer = F,
                                                         h_conquer = 0.05,
                                                         method = "general", n_iter.max = 100, verbose = F)
      ci_or = twmq_ci(Ymat, X_tensor, W,
                                                      member_G0, member_H0,
                                                      taus, conquer = F,
                                                      h_conquer = 0.05,
                                                      method = method, n_iter.max = 100, verbose = F)
      
      # RMSE/CP/AE for each oracle estimator
      res_vec_ge[1:(15+(G0*H0*5+G0*3*5))] = cal_para_general(ci_ge_or$point,ci_ge_or,  alpha0, beta0, nu0, gamma0, additive = additive)
      res_vec[1:(20 + G0 * 5 * 5 )] = cal_para(ci_or$point,ci_or, alpha0, beta0, nu0, gamma0)
      
      
      # ============================================================
      # (6.3) TGNQ-GENERAL estimator (no refinement)
      # ============================================================
      ## Initial estimator of the model
      res_ge_ini = twmq.estimate.auto(Ymat, X_tensor, W,G,H, taus, method = "general", verbose = F,conquer = F, h_conquer=0.05, ntrial = 100)
      ###Updated the parameter algorithm further
      try <- min(3,length(res_ge_ini$loss)) ##Update the first "try" esimated model with lowest loss
      res_tmp <- list()
      idx <- sort(res_ge_ini$loss,index=TRUE)$ix
      try_loss <- numeric(try)
      for(i in 1:try){
        res <- res_ge_ini$res_all[[idx[i]]]
        ##############################
        res_tmp[[i]] <-   update_NARG_twmq(Ymat, W,X_tensor, res$member_G,res$member_H,taus,res$theta,method ="general",G,H,conquer=F,h_conquer=NULL,Iter=10,frac=0.5,MaxOutIter=100,Maxit=5)
        try_loss[i] <- sum(res_tmp[[i]]$theta_GH$loss)
      }
      res_ge_update <- res_tmp[[which.min(try_loss)]]
      ###when G=G0, H=H0, one needs to do label switching for estimated parameters.
      res_ge_update1 <-  twmq.label.switch(res_ge_update, member_G0,member_H0, method="general")
      ci_ge_gtnq = twmq_ci(Ymat, X_tensor, W,
                                                           res_ge_update1$member_G, res_ge_update1$member_H,
                                                           taus, conquer = F,
                                                           h_conquer = 0.05,
                                                           method ="general", n_iter.max = 100, verbose = F)
      
      start_ge = (15+(G0*H0*5+G0*3*5)) + 1
      end_ge = (15+(G0*H0*5+G0*3*5)) *2
      res_vec_ge[start_ge:end_ge] = cal_para_general(ci_ge_gtnq$point,ci_ge_gtnq, alpha0, beta0, nu0, gamma0, additive = additive)
      res_vec_ge[(end_ge+1):(end_ge+2)] = c(err_rate_mapping(res_ge_update1$member_G, member_G0),
                                            err_rate_mapping(res_ge_update1$member_H, member_H0))
      
      # ============================================================
      # (6.4) TGNQ-SPECIAL estimator (no refinement)
      # ------------------------------------------------------------
      # Initialize from the general-model fit, then refine the special
      # structure via one update_NARG_twmq pass with `method` (the correct
      # special form).
      # ============================================================
      res =  twmq_ci(Ymat, X_tensor, W,
                                                     res_ge_update1$member_G, res_ge_update1$member_H,
                                                     taus, conquer = F,
                                                     h_conquer = 0.05,
                                                     method =method, n_iter.max = 100, verbose = F)
      res_update = update_NARG_twmq(Ymat, W,X_tensor, res_ge_update1$member_G,res_ge_update1$member_H,taus,res$point,method,G,H,conquer=F,h_conquer=NULL,Iter=10,frac=0.5,MaxOutIter=100,Maxit=5)
      res_update1 <-  twmq.label.switch(res_update, member_G0,member_H0, method)
      ci_gtnq = twmq_ci(Ymat, X_tensor, W,
                                                        res_update1$member_G, res_update1$member_H,
                                                        taus, conquer = F,
                                                        h_conquer = 0.05,
                                                        method =method, n_iter.max = 100, verbose = F)
      start = (20+G0*5*5) + 1
      end = (20+G0*5*5) *2
      res_vec[start:end] = cal_para(ci_gtnq$point, ci_gtnq,alpha0, beta0, nu0, gamma0)
      res_vec[(end+1):(end+2)] = c(err_rate_mapping(res_update1$member_G, member_G0),
                                   err_rate_mapping(res_update1$member_H, member_H0))

      
      # ============================================================
      # (6.5) MISSPECIFIED TGNQ estimator
      # ------------------------------------------------------------
      # Same procedure as in (6.3) but with `method_mis` instead of "general".
      # We then also produce a "transferred" variant: misspecified parameters
      # evaluated AT the correctly specified memberships (`res_vec_mis_T`).
      # ============================================================
      res_ini_mis = twmq.estimate.auto(Ymat, X_tensor, W,G,H, taus, method = method_mis, verbose = F,conquer = F, h_conquer=0.05, ntrial = 100)
      ###Updated the parameter algorithm further
      try <- min(3,length(res_ini_mis$loss)) ##Update the first "try" esimated model with lowest loss
      res_tmp <- list()
      idx <- sort(res_ini_mis$loss,index=TRUE)$ix
      try_loss <- numeric(try)
      for(i in 1:try){
        res_mis <- res_ini_mis$res_all[[idx[i]]]
        ##############################
        res_tmp[[i]] <-   update_NARG_twmq(Ymat, W,X_tensor, res_mis$member_G,res_mis$member_H,taus,res_mis$theta,method =method_mis,G,H,conquer=F,h_conquer=NULL,Iter=10,frac=0.5,MaxOutIter=100,Maxit=5)
        try_loss[i] <- sum(res_tmp[[i]]$theta_GH$loss)
      }
      res_update_mis <- res_tmp[[which.min(try_loss)]]
      ###when G=G0, H=H0, one needs to do label switching for estimated parameters.
      res_update1_mis <-  twmq.label.switch(res_update_mis, member_G0,member_H0, method=method_mis)
      ci_gtnq_mis = twmq_ci(Ymat, X_tensor, W,
                                                            res_update1_mis$member_G, res_update1_mis$member_H,
                                                            taus, conquer = F,
                                                            h_conquer = 0.05,
                                                            method =method_mis, n_iter.max = 100, verbose = F)
      
      res_vec_mis = cal_para_mis(ci_gtnq_mis$point, ci_gtnq_mis,alpha0, beta0, nu0, gamma0,additive = additive)
      res_vec_mis = c(res_vec_mis, err_rate_mapping(res_update1_mis$member_G, member_G0),
                      err_rate_mapping(res_update1_mis$member_H, member_H0))
      res_vec_mis_T = cal_para_mis_T(ci_gtnq$point, ci_gtnq,alpha0, beta0, nu0, gamma0,additive = additive)
      res_vec_mis_T = c(res_vec_mis_T, res_vec[(end+1):(end+2)])

      
      # ============================================================
      # (6.6) Refinement step on the TGNQ-GE estimator
      # ------------------------------------------------------------
      # Run Refine_G then Refine_H, only accepting reassignments whose
      # relative loss improvement exceeds 1/sqrt(T).
      # ============================================================
      ##Refine membership G
      res_ge_update1$theta_GH$Loss.fun = res_ge_update1$theta_GH$Loss_fun
      member_G_refine_ge <- res_ge_update1$member_G
      ##############################
      obj_G_refine_ge <-  Refine_G(Ymat, X_tensor, W, res_ge_update1$member_G, res_ge_update1$member_H, taus, res_ge_update1$theta, method ="general", G, H)
      ##Check how many H memberships needs to be refined
      idx_G <- which((obj_G_refine_ge$loss_old-obj_G_refine_ge$loss_new)/obj_G_refine_ge$loss_old>1/sqrt(Time))
      if(length(idx_G)>0) member_G_refine_ge[idx_G]=obj_G_refine_ge$g_r[idx_G]
      ##Refine membership H
      member_H_refine_ge <- res_ge_update1$member_H
      ##############################
      obj_H_refine_ge <-  Refine_H(Ymat, X_tensor, W, member_G_refine_ge, res_ge_update1$member_H, taus, res_ge_update1$theta, method ="general", G, H, FriendW, FriendW2)
      ##Check how many H memberships needs to be refined
      idx_H <- which((obj_H_refine_ge$loss_old-obj_H_refine_ge$loss_new)/obj_H_refine_ge$loss_old>1/sqrt(Time))
      if(length(idx_H)>0) member_H_refine_ge[idx_H]=obj_H_refine_ge$h_r[idx_H]
      ci_ge_refine = twmq_ci(Ymat, X_tensor, W,
                                                             member_G_refine_ge, member_H_refine_ge,
                                                             taus, conquer = F,
                                                             h_conquer = 0.05,
                                                             method ="general", n_iter.max = 100, verbose = F)
      start_ge = end_ge + 3
      end_ge = length(res_vec_ge)-2
      res_vec_ge[start_ge:end_ge] = cal_para_general(ci_ge_refine$point,ci_ge_refine, alpha0, beta0, nu0, gamma0, additive = additive)
      res_vec_ge[(end_ge+1):(end_ge+2)] = c(err_rate_mapping(member_G_refine_ge, member_G0),err_rate_mapping(member_H_refine_ge, member_H0))

      # ============================================================
      # (6.7) Refinement step on the TGNQ-SPECIAL estimator
      # ============================================================
      res_update1$theta_GH$Loss.fun = res_update1$theta_GH$Loss_fun
      member_G_refine <- res_update1$member_G
      obj_G_refine <-  Refine_G(Ymat, X_tensor, W, res_update1$member_G, res_update1$member_H, taus, res_update1$theta, method =method, G, H)
      ##Check how many H memberships needs to be refined
      idx_G <- which((obj_G_refine$loss_old-obj_G_refine$loss_new)/obj_G_refine$loss_old>1/sqrt(Time))
      if(length(idx_G)>0) member_G_refine[idx_G]=obj_G_refine$g_r[idx_G]
      ##Refine membership H
      member_H_refine <- res_update1$member_H
      obj_H_refine <-  Refine_H(Ymat, X_tensor, W, member_G_refine, res_update1$member_H, taus, res_update1$theta, method =method, G, H, FriendW, FriendW2)
      ##Check how many H memberships needs to be refined
      idx_H <- which((obj_H_refine$loss_old-obj_H_refine$loss_new)/obj_H_refine$loss_old>1/sqrt(Time))
      if(length(idx_H)>0) member_H_refine[idx_H]=obj_H_refine$h_r[idx_H]
      ci_refine = twmq_ci(Ymat, X_tensor, W,
                                                          member_G_refine, member_H_refine,
                                                          taus, conquer = F,
                                                          h_conquer = 0.05,
                                                          method =  method, n_iter.max = 100, verbose = F)
      start = end + 3
      end = length(res_vec)-2
      res_vec[start:end] = cal_para(ci_refine$point, ci_refine,alpha0, beta0, nu0, gamma0)
      res_vec[(end+1):(end+2)] = c(err_rate_mapping(member_G_refine, member_G0),
                                   err_rate_mapping(member_H_refine, member_H0))
      
      # ============================================================
      # (6.8) Save the full per-seed snapshot and return error vectors
      # ============================================================
      sim_times = sim_times + 1

      fn = paste0("output/",ifelse(block, "block/", "power/"),G0,"-c_",case,"-a_",additive,"-seed_",seed-500,".rda")

      save(ci_ge_or, ci_ge_refine, 
           res_ge_update1,
           member_G_refine_ge, member_H_refine_ge,
           X_tensor, Ymat, W, Amat,
           member_G0, member_H0, file = fn)
      return(list(res_vec,res_vec_ge,res_vec_mis,res_vec_mis_T))
    }, error = function(e){
      cat("Error in iteration", seed - 500, ":", conditionMessage(e), "\n")
    })
  }
}