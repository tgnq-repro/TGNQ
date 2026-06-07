# =============================================================================
# simulator.R
# -----------------------------------------------------------------------------
# Defines simulation(seed), ONE Monte Carlo replicate of the Section 5.4
# (TGNQ vs QGNAR) experiment. The function:
#
#   1. Sets (N, T, Nblock) from `case`.
#   2. Builds the (fixed across seeds) block or power-law network W.
#   3. Generates true memberships and the response panel Y from the TGNQ DGP.
#   4. Fits two estimators:
#        - TGNQ-general (twmq.estimate.auto + update_NARG_twmq + label switch),
#        - QGNAR via gqnar(...) from estimator_gqnar.R.
#   5. Saves the true memberships and both fits to
#        output/res_gqnar/res_<block|power>_<add|mul>/<G>-<case>-<seed>.rda
#      for downstream aggregation by step2_output_tex.R.
#
# Expected globals (set in simulation.R): case, additive, block.
# =============================================================================
simulation = function(seed){
  sim_times = 0
  
  # ---- One retry loop: re-draw on numerical failure -----------------------
  while (sim_times < 1) {
    tryCatch({
      set.seed(seed)
      seed = seed + 500 # offset to avoid colliding with other seeds
      
      # ---- (1) Design parameters from `case` -----------------------------
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
      
      
      G = G0 = 2
      H = H0 = 2
      taus = c(0.1, 0.3, 0.5, 0.7, 0.9)
      n_taus = length(taus) 
      p=2
      
      # ---- (2) Build the network (fixed across seeds via set.seed(123)) --
      set.seed(123)
      p = 2
      
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
      
      # ---- (3) True group memberships -----------------------------------
      # In Section 5.4 only G0 = 2 is used; the if/else covers G0 = 3 as well
      # in case the script is reused with a different G0.
      if (G0==2){
        member_G0 = simu.member0(N, ratios = c(0.5, 0.5))
        member_H0 = simu.member0(N, ratios = c(0.4, 0.6))
      }else{
        member_G0 = simu.member0(N, ratios = c(0.3, 0.3, 0.4))
        member_H0 = simu.member0(N, ratios = c(0.4, 0.3, 0.3))
      }
      G=G0
      
      
      
      # ---- (4) Simulate covariates and the response panel ----------------
      X_tensor = array(rnorm(N*p*(Time+10)), dim = c(N, p, (Time+10)))
      X_tensor = abs(X_tensor)
      Ymat = simu.Y(W, member_G0, member_H0, X_tensor, additive = additive, verbose = F)
      X_tensor = X_tensor[,,-(1:11)]
      
      # ============================================================
      # (5) Fit TGNQ-general
      # ------------------------------------------------------------
      # Initialization: many random starts via twmq.estimate.auto, then
      # take the top `try = min(3, .)` losses and refine each one with
      # update_NARG_twmq. Keep the refinement with the lowest loss.
      # Finally label-switch so estimated group labels match the truth.
      # ============================================================
      res_ge_ini = twmq.estimate.auto(Ymat, X_tensor, W,G,H, taus, method = "general", verbose = F,conquer = F, h_conquer=0.05, ntrial = 100)
      ### Updated the parameter algorithm further
      try <- min(3,length(res_ge_ini$loss)) ##Update the first "try" esimated model with lowest loss
      res_tmp <- list()
      idx <- sort(res_ge_ini$loss,index=TRUE)$ix
      try_loss <- numeric(try)
      for(i in 1:try){
        res <- res_ge_ini$res_all[[idx[i]]]
        res_tmp[[i]] <-   update_NARG_twmq(Ymat, W,X_tensor, res$member_G,res$member_H,taus,res$theta,method ="general",G,H,conquer=F,h_conquer=NULL,Iter=10,frac=0.5,MaxOutIter=100,Maxit=5)
        try_loss[i] <- sum(res_tmp[[i]]$theta_GH$loss)
      }
      res_ge_update <- res_tmp[[which.min(try_loss)]]
      res_ge_update1 <-  twmq.label.switch(res_ge_update, member_G0,member_H0, method="general")
      
      # ============================================================
      # (6) Fit QGNAR (gqnar)
      # ------------------------------------------------------------
      # gqnar expects X_tensor in (N, Time, P) layout, whereas twmq uses
      # (N, P, Time). We permute axes before calling gqnar.
      # ============================================================
      X_tensor = aperm(X_tensor,c(1,3,2))
      narg = gqnar(Ymat, W, X_tensor, G,taus)
      sim_times = 1
      
      # ---- (7) Save snapshot -------------------------------------------
      # File name: <G>-<case>-<original seed>.rda
      save(member_G0, member_H0,res_ge_update1,narg, 
           file = paste0("output/res_gqnar/res_", ifelse(block,"block_","power_"), ifelse(additive,"add/","mul/"), G,"-", case, "-", seed-500, ".rda"))
    }, error = function(e){
      cat("Error in iteration", seed - 500, ":", conditionMessage(e), "\n")
    })
  }
}