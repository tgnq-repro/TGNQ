# =============================================================================
# step5_fit_qnarg.R
#
# Fit the QGNAR (one-way grouped quantile network autoregression) benchmark
# on the Weibo data, for two training windows (T_train in {70, 77}) and two
# group counts (G in {3, 4}).
#
# Output: output/real_res/qnarg_G<G>-<T_train>.rda  (object: narg_mem_final)
# =============================================================================

rm(list=ls())

for(T_train in c(70,77)){
  for (G in 3:4) {
    library(abind)
    
    
    library(here)
    dir <- here()
    setwd(dir)
    
    gqnar_script <- here("code", "real_data", "estimator_gqnar.R")
    source(gqnar_script)
    

    # ---- Load data and assemble training inputs ----
    load("data/weibo.rda")
    N = nrow(Ymat)
    taus=c(0.1,0.3,0.5,0.7,0.9)
    

    X_tensor <- aperm(replicate(T_train, Xi, simplify = "array"), c(1,3,2))
    X_tensor_new <- array(0, dim = c(N, T_train, 10))
    X_tensor_new[,,1:8] <- X_tensor
    Xt_all = Xt
    Xt = Xt[2:(T_train+1),]
    for (i in 1:N) {
      X_tensor_new[i,,9:10] = Xt
    }
    p <- 10
    
    X_train = X_tensor_new
    Y_train = Ymat[,1:(T_train+1)]
    
    
    W = Amat[1:N,1:N]
    W = as.matrix(W)
    W = W/rowSums(W)

    # ---- Fit QGNAR with multi-start parallel K-means initializations -------
    numCores=20
    taus=c(0.1,0.3,0.5,0.7,0.9)
    
    p = dim(X_tensor)[3]
    N = nrow(Ymat)
    
    narg_mem_final = est.NARG.member.parallel(Y_train, W, X_train, G, taus,numCores = numCores, n.initial = 100,member.int=NULL, script_path = gqnar_script)
    print(narg_mem_final)
    save(narg_mem_final, file = paste0("output/real_res/qnarg_G", G, "-", T_train,".rda"))
    
  }
  
}


