# =============================================================================
# simulator.R
# -----------------------------------------------------------------------------
# Defines TWO functions used in the Section 5.2 (select G, H) experiments:
#
#   simulation(seed)     -- one candidate-grid replication; fits seven TGNQ
#                            models with (G, H) in
#                            {(2,2),(2,3),(3,2),(3,3),(3,4),(4,3),(4,4)},
#                            saves results to output/res_selectGH_*.
#
#   simulation_or(seed)  -- one oracle replication; fits the general and the
#                            special TGNQ models at the TRUE (G0, H0)
#                            memberships, saves to output/*_select_or_*.
#
# Both functions expect the following objects to exist in the calling scope
# (they are created in simulation.R / simulation_or.R):
#   G0       - true row-group number (3 in Section 5.2)
#   case     - design index (1..5)
#   additive - logical, additive vs. multiplicative DGP
#   block    - logical, SBM vs. power-law network
#   method   - "additive" or "multiplicative"
#   seed     - Monte Carlo replication index
# =============================================================================


# -----------------------------------------------------------------------------
# simulation(seed): one candidate-grid replication.
# -----------------------------------------------------------------------------
simulation = function(seed){
  # ---- (1) Set the (N, T, Nblock) design according to `case` ---------------
  H0=G0
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
  taus = c(0.1, 0.3, 0.5, 0.7, 0.9) # quantile levels evaluated
  n_taus = length(taus) 
  
  # ---- (2) Generate the (fixed) network for this design --------------------
  # Held fixed across replications so that randomness comes only from data.
  set.seed(123)
  p = 2
  if(block){
    Amat = getBlockW(N,Nblock, normalize = F)
  }else{
    Amat = getPowerLawW(N)
  }
  W = Amat[1:N,1:N]
  rw = rowSums(W)
  W[rw==0, sample(1:N, 3)] = 1   # patch zero-degree rows
  W = W/rowSums(W)  # row-normalize
  
  # First- and second-order follower sets for each node.
  FriendW <-FriendW2 <- vector("list", N)
  for(i in 1:N){
    FriendW[[i]] <- as.vector(which(W[,i]!=0))
    for (j in  FriendW[[i]]) {
      FriendW2[[i]] <- unique(c(FriendW2[[i]],as.vector(which(W[j,]!=0))))
    }
  }
  
  
  # ---- (3) True group memberships ------------------------------------------
  # Throughout Section 5.2 we fix G0 = H0 = 3 with ratios (0.3, 0.3, 0.4) and
  # (0.4, 0.3, 0.3) respectively.
  member_G0 = simu.member0(N, ratios = c(0.3, 0.3, 0.4))
  member_H0 = simu.member0(N, ratios = c(0.4, 0.3, 0.3))
  
  # ---- (4) Run ONE replication; retry on rare numerical failures -----------
  sim_times = 0
  
  while (sim_times < 1) {
    tryCatch({
      # Reseed and bump by 500 so any retry uses a fresh draw.
      set.seed(seed)
      seed = seed + 500
      
      # ---- (4.1) Simulate covariates and the response panel -----------------
      X_tensor = array(rnorm(N*p*(Time+10)), dim = c(N, p, (Time+10)))
      X_tensor = abs(X_tensor)
      set.seed(seed)
      Ymat = simu.Y(W, member_G0, member_H0, X_tensor, additive = additive, verbose = F)
      X_tensor = X_tensor[,,-(1:11)]
      
      
      # ---- (4.2) Auto-initialization for each candidate (G, H) --------------
      # Off-diagonal (G != H) pairs ((2,3),(3,2),(4,3),(3,4)) are harder to initialize,
      # so they use ntrial = 1000; the remaining pairs ((2,2),(3,3),(4,4)) use ntrial = 100.
      int5 = twmq.estimate.auto(Ymat, X_tensor, W,G=4,H=3, taus, method = "general", verbose = F,conquer = F, h_conquer=0.05, ntrial = 1000)
      print(paste0("Initialization:", 43))
      int6 = twmq.estimate.auto(Ymat, X_tensor, W,G=3,H=4, taus, method = "general", verbose = F,conquer = F, h_conquer=0.05, ntrial = 1000)
      print(paste0("Initialization: ", 34))
      int2 = twmq.estimate.auto(Ymat, X_tensor, W,G=2,H=3, taus, method = "general", verbose = F,conquer = F, h_conquer=0.05, ntrial = 1000)
      print(paste0("Initialization: ", 23))
      int3 = twmq.estimate.auto(Ymat, X_tensor, W,G=3,H=2, taus, method = "general", verbose = F,conquer = F, h_conquer=0.05, ntrial = 1000)
      print(paste0("Initialization: ", 32))
      int7 = twmq.estimate.auto(Ymat, X_tensor, W,G=4,H=4, taus, method = "general", verbose = F,conquer = F, h_conquer=0.05, ntrial = 100)
      print(paste0("Initialization: ", 44))
      int1 = twmq.estimate.auto(Ymat, X_tensor, W,G=2,H=2, taus, method = "general", verbose = F,conquer = F, h_conquer=0.05, ntrial = 100)
      print(paste0("Initialization: ", 22))
      int4 = twmq.estimate.auto(Ymat, X_tensor, W,G=3,H=3, taus, method = "general", verbose = F,conquer = F, h_conquer=0.05, ntrial = 100)
      print(paste0("Initialization: ", 33))
      print(paste0("Initialization complete: ", seed-500))
      
      # ---- (4.3) Loop over the seven candidate (G, H) pairs -----------------
      # For each pair: pick the 3 initializations with the smallest loss,
      # refine each via update_NARG_twmq, keep the best one (res_ge_update),
      # then run ANOTHER update_NARG_twmq from the chosen point to obtain
      # the special-model estimator (res_update).
      for (gh in c(5,6,2,3,7,1,4)) {
        if(gh==1){
          G=2
          H=2
          res_ge_ini = int1
        }
        if(gh==2){
          G=2
          H=3
          res_ge_ini = int2
        }
        if(gh==3){
          G=3
          H=2
          res_ge_ini = int3
        }
        if(gh==4){
          G=3
          H=3
          res_ge_ini = int4
        }
        if(gh==5){
          G=4
          H=3
          res_ge_ini = int5
        }
        if(gh==6){
          G=3
          H=4
          res_ge_ini = int6
        }
        if(gh==7){
          G=4
          H=4
          res_ge_ini = int7
        }
        
        
        # Update the top-`try` initializations and keep the best one (lowest post-update loss). This is the same pattern as in Section 5.1.
        try <- min(3,length(res_ge_ini$loss))
        res_tmp <- list()
        idx <- sort(res_ge_ini$loss,index=TRUE)$ix
        try_loss <- numeric(try)
        for(i in 1:try){
          res <- res_ge_ini$res_all[[idx[i]]]
          res_tmp[[i]] <-   update_NARG_twmq(Ymat, W,X_tensor, res$member_G,res$member_H,taus,res$theta,method ="general",G,H,conquer=F,h_conquer=NULL,Iter=10,frac=0.5,MaxOutIter=100,Maxit=5)
          try_loss[i] <- sum(res_tmp[[i]]$theta_GH$loss)
        }
        res_ge_update <- res_tmp[[which.min(try_loss)]]
        
        # Special-model estimator: starts from the general-model fit and uses `method` ("additive" or "multiplicative") to enforce the corresponding parameter structure.
        res_update = update_NARG_twmq(Ymat, W,X_tensor, res_ge_update$member_G,res_ge_update$member_H,taus,res_ge_update$theta_GH,method,G,H,conquer=F,h_conquer=NULL,Iter=10,frac=0.5,MaxOutIter=100,Maxit=5)
        
        # ---- (4.4) Save the two estimators for this (G, H) ------------------
        fn = paste0("output/res_selectGH_",ifelse(block, "block", "power"),ifelse(additive, "_add", "_mul"),"/",case, "-",G, "-", H, "-", seed-500, ".rda")
        
        save(res_ge_update,res_update,member_G0,member_H0, file = fn)
      }
      sim_times = sim_times + 1
      return(0)
    }, error = function(e){
      cat("Error in iteration", seed - 500, ":", conditionMessage(e), "\n")
    })
  }
}



# -----------------------------------------------------------------------------
# simulation_or(seed): one oracle replication, used to compute the Oracle row
# in the Section 5.2 tables.
# -----------------------------------------------------------------------------
simulation_or = function(seed){
  # ---- (1) Set (N, T, Nblock) according to `case` --------------------------
  H0=G0
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
  if(case==5) {
    N=100
    Time=50
    Nblock = 5
  }
  
  taus = c(0.1, 0.3, 0.5, 0.7, 0.9) 
  n_taus = length(taus) 
  
  # ---- (2) Build the (fixed) network ---------------------------------------
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
  
  # ---- (3) True memberships ------------------------------------------------
  member_G0 = simu.member0(N, ratios = c(0.3, 0.3, 0.4))
  member_H0 = simu.member0(N, ratios = c(0.4, 0.3, 0.3))
  
  # ---- (4) True parameters (used only by step3_output_tex.R) ----------------
  # alpha0[g, k]: true row-group effect, beta0[h, k]: col-group effect,
  # nu0[g, k]: true autoregressive effect, gamma0[[k]]: list of G x p
  # covariate-effect matrices, one per quantile.
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
  
  # ---- (5) One replication; retry on numerical failure ---------------------
  sim_times = 0
  
  while (sim_times < 1) {
    tryCatch({
      set.seed(seed)
      seed = seed + 500
      
      # Simulate data the same way as in simulation()
      X_tensor = array(rnorm(N*p*(Time+10)), dim = c(N, p, (Time+10)))
      X_tensor = abs(X_tensor)
      G = max(member_G0); H = max(member_H0)
      set.seed(seed)
      Ymat = simu.Y(W, member_G0, member_H0, X_tensor, additive = additive, verbose = F)
      X_tensor = X_tensor[,,-(1:11)]
      
      
      # ---- (5.1) Oracle estimators ------------------------------------------
      # `or_ge`: general TGNQ at true memberships.
      # `or_spe`: special (additive / multiplicative) TGNQ at true memberships.
      or_ge = twmq.estimate_thetaGH.member.iterate(Ymat, X_tensor, W,
                                                   member_G0, member_H0,
                                                   taus, conquer = F,
                                                   h_conquer = 0.05,
                                                   method = "general", n_iter.max = 100, verbose = F)
      or_spe = twmq.estimate_thetaGH.member.iterate(Ymat, X_tensor, W,
                                                    member_G0, member_H0,
                                                    taus, conquer = F,
                                                    h_conquer = 0.05,
                                                    method = method, n_iter.max = 100, verbose = F)
      # ---- (5.2) Save the oracle outputs ------------------------------------
      save(member_G0,member_H0,or_ge,or_spe, 
           file = paste0("output/", ifelse(block, "block", "power"),"_select_or_",ifelse(additive, "add", "mul"),"/",case,"-",seed-500, ".rda"))
      return(seed - 500)
    }, error = function(e){
      cat("Error in iteration", seed - 500, ":", conditionMessage(e), "\n")
    })
  }
}
