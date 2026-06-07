# =============================================================================
# real_tgnq.R
#
# Fit the proposed TGNQ (Two-way Grouped Network Quantile) model on the Weibo
# data for ONE configuration (T_train, G, H, rq_lambda).
#
# Called from: step1_selectGH.sh (in parallel, once per (G, H) cell).
#
# Command-line arguments:
#   args[1] = T_train     # length of the training window
#   args[2] = G           # number of row groups
#   args[3] = H           # number of column groups
#   args[4] = rq_lambda   # L1 penalty for the initial rq fits
#
# Output: output/real_res/tgnq_<T_train>-<G>-<H>.rda  (object: res_ge_update)
# =============================================================================
rm(list=ls())

set.seed(123)

# ---- 1. Parse command-line arguments ----------------------------------------
args <- commandArgs(TRUE)
T_train <- as.integer(args[1])
G <- as.integer(args[2])
H <- as.integer(args[3])
rq_lambda = as.numeric(args[4])

library(twmq)


library(here)
dir <- here()
setwd(dir)

source("code/real_data/estimator_gqnar.R")


# ---- 2. Load data and set up output directory --------------------------------
load("data/weibo.rda")
N = nrow(Ymat)
taus=c(0.1,0.3,0.5,0.7,0.9)

numCores = 16


if (!file.exists("output/real_res")) {
  dir.create("output/real_res")
}

# ---- 3. Build the 3-D covariate tensor X_train (N x T_train x 10) ------------
#   Slices 1:8 = node-level covariates Xi (constant in time)
#   Slices 9:10 = time-level covariates Xt (constant across nodes)
X_tensor <- aperm(replicate(T_train, Xi, simplify = "array"), c(1,3,2))
X_tensor_new <- array(0, dim = c(N, T_train, 10))
X_tensor_new[,,1:8] <- X_tensor

# Align Xt[t] with response Y_{t+1}
Xt = Xt[2:(T_train+1),]
for (i in 1:N) {
  X_tensor_new[i,,9:10] = Xt
}
p <- 10

X_train = X_tensor_new
Y_train = Ymat[,1:(T_train+1)]

W = Amat[1:N,1:N]


# ---- 4. Build the row-normalized network weight matrix W --------------------
W = as.matrix(W)
W = W/rowSums(W)

# ---- 5. Estimate TGNQ -------------------------------------------------------
if(G==1 & H==1){
  res_ge_update = twmq.estimate_thetaGH.member.general(Y_train, aperm(X_train, c(1,3,2)), W, 
                                                       member_G=rep(1,N), member_H=rep(1,N), 
                                                       taus=taus,conquer=F,h_conquer=0.05)
  loss = lapply(1:5, function(k)check.func(res_ge_update$resi[,k], taus[[k]]))
  loss = lapply(loss, sum)
  res_ge_update$theta_GH$loss = unlist(loss)
  res_ge_update$member_G = rep(1,N)
  res_ge_update$member_H = rep(1,N)
}else{
  # ---- 5a. Multi-start initial estimator ------------------------------------
  ## Initial estimator of the model
  res_ge_ini = twmq.estimate.auto.parallel(Y_train, aperm(X_train, c(1,3,2)), W,G=G,H=H,
                                           taus=taus, method = "general", 
                                           verbose = F,conquer = F, h_conquer=0.05, ntrial = 100,numCores = numCores,rq_lambda=rq_lambda)
  # ---- 5b. Refine the top-3 starts with the iterative algorithm -------------
  ### Updated the parameter algorithm further
  try <- min(3,length(res_ge_ini$loss)) ##Update the first "try" esimated model with lowest loss
  res_tmp <- list()
  idx <- sort(res_ge_ini$loss,index=TRUE)$ix
  try_loss <- numeric(try)
  for(i in 1:try){
    res <- res_ge_ini$res_all[[idx[i]]]
    res_tmp[[i]] <- update_NARG_twmq_parallel(Y_train, W, aperm(X_train, c(1,3,2)), 
                                              res$member_G,res$member_H, taus=c(0.1,0.3,0.5,0.7,0.9),res$theta,
                                              method ="general",G=G,H=H,conquer=F,h_conquer=0.05,
                                              Iter=10,frac=0.5,numCores=numCores,MaxOutIter=100,Maxit=5,rq_lambda=rq_lambda)
    try_loss[i] <- sum(res_tmp[[i]]$theta_GH$loss)
  }
  # Keep the start that achieves the smallest total loss across τ
  res_ge_update <- res_tmp[[which.min(try_loss)]]
}

# ---- 6. Save -----------------------------------------------------------------
save(res_ge_update, file = paste0("output/real_res/tgnq_",T_train,"-",G,"-",H,".rda"))
