# =============================================================================
# step7_fit_compare.R
#
# In-sample goodness-of-fit comparison between TGNQ and QGNAR (G = 3 and 4)
# at T_train = 70, using a quantile pseudo-R^2:
#
#       1 - V_TGNQ / V_QGNAR
#
# where V is the sum of pinball losses for the corresponding model.
# Computed and printed for each tau in {0.1, 0.3, 0.5, 0.7, 0.9}.
# =============================================================================

rm(list=ls())
library(here)
dir <- here()
setwd(dir)



T_train=70
# ---- 1. Load data and assemble training tensor ------------------------------
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

# ---- 2. Prediction helper (full panel) --------------------------------------
pred = function(theta_Gs, alphabeta_GHs, Y1, W, X, member_G, member_H){
  G = nrow(alphabeta_GHs)
  H = ncol(alphabeta_GHs)
  N = length(Y1)
  WY1 <- W%*%Y1
  WY1 = matrix(WY1,ncol=1)
  Yh <- list()
  for (h in 1:H) {
    Yh[[h]] = W%*%(Y1*(member_H==h))
  }
  
  Y_pred <- numeric(N)

  for (g in 1:G)
  {
    ind_g = which(member_G == g)
    Y1_g=Y1[ind_g,,drop=F]
    Z1 <- as.vector(Y1_g)
    Zh <- NULL
    for (h in 1:H) {
      Zh <- cbind(Zh,as.vector((Yh[[h]])[ind_g,,drop=F]))
    }
    X1 = X[ind_g,]
    Z_all <- cbind(Z1,X1[,,drop=T],Zh)
    
    theta = cbind(theta_Gs,alphabeta_GHs)
    
    Y_pred[ind_g] = Z_all%*%theta[g,]
  }
  return(Y_pred)
}

# ---- 3. TGNQ in-sample predictions over all training time points -----------
load("output/real_res/tgnq_70-4-3.rda")
TGNQ.Y_pred_taus = lapply(1:5, function(k){
  pred_mat = lapply(1:T_train, function(t){
    pred(res_ge_update$theta_GH$theta_Gs[[k]], res_ge_update$theta_GH$alphabeta_GHs[[k]],
         matrix(Y_train[,t],ncol=1), W, X_train[,t,],
         res_ge_update$member_G, res_ge_update$member_H)
  })
  pred_mat = do.call(cbind, pred_mat)
  pred_mat
})

# ---- 4. QGNAR in-sample predictions for G = 3 and G = 4 --------------------
load("output/real_res/qnarg_G3-70.rda")
G=3
GNAR3.Y_pred_taus = lapply(1:5, function(k){
  pred_mat = lapply(1:T_train, function(t){
    pred(narg_mem_final$psi[[k]][,-(1:G)], narg_mem_final$psi[[k]][,1:G],
         matrix(Y_train[,t],ncol=1), W, X_train[,t,],
         narg_mem_final$member, narg_mem_final$member)
  })
  pred_mat = do.call(cbind, pred_mat)
  pred_mat
})

load("output/real_res/qnarg_G4-70.rda")
G=4
GNAR4.Y_pred_taus = lapply(1:5, function(k){
  pred_mat = lapply(1:T_train, function(t){
    pred(narg_mem_final$psi[[k]][,-(1:G)], narg_mem_final$psi[[k]][,1:G],
         matrix(Y_train[,t],ncol=1), W, X_train[,t,],
         narg_mem_final$member, narg_mem_final$member)
  })
  pred_mat = do.call(cbind, pred_mat)
  pred_mat
})






# ---- 5. Quantile pseudo-R^2 (Koenker & Machado, 1999) ----------------------
goodfit <- function(resid, resid_nl,resid_full, tau){
  # minimum sum of deviations
  V1 <- resid * (tau - (resid < 0))
  V1 <- sum(V1, na.rm = T) 
  
  # null sum of deviations
  V0 <- resid_nl * (tau - (resid_nl < 0))
  V0 <- sum(V0, na.rm = T) 
  
  # null sum of deviations
  Vf <- resid_full * (tau - (resid_full < 0))
  Vf <- sum(Vf, na.rm = T) 
  
  # explained deviance
  out <- 1 - (V1-Vf)/(V0-Vf)
  
  # exceptions for output
  if(any(c(Inf, -Inf) %in% out)) out <- NA
  if(V1 > V0) out <- NA
  
  return(out)
  
}






# Compare TGNQ to QGNAR(G = 3)
taus = seq(0.1,0.9,0.2)


for (k in 1:5) {
  resid_nl = Y_train[,-1] - GNAR3.Y_pred_taus[[k]]
  # resid_full = node.resi_taus[[k]]
  resid = Y_train[,-1] - TGNQ.Y_pred_taus[[k]]
  print(paste0("tau=",taus[k],":", round(goodfit(resid, resid_nl,0,  taus[k]),4)))
}




# Compare TGNQ to QGNAR(G = 4)
taus = seq(0.1,0.9,0.2)

for (k in 1:5) {
  resid_nl = Y_train[,-1] - GNAR4.Y_pred_taus[[k]]
  # resid_full = node.resi_taus[[k]]
  resid = Y_train[,-1] - TGNQ.Y_pred_taus[[k]]
  print(paste0("tau=",taus[k],":", round(goodfit(resid, resid_nl,0,  taus[k]),4)))
}



