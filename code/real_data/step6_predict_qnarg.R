# =============================================================================
# step6_predict_qnarg.R
#
# One-step-ahead QGNAR predictions on test days T_test = 71,...,77, using the
# QGNAR fits from Step 5 (T_train = 70, G in {3, 4}).
#
# Output: output/real_res/narg_pred_train70-test<T>.rda
#         (objects: Y_test, GNAR3.Y_pred_taus, GNAR4.Y_pred_taus)
# =============================================================================

rm(list = ls())
library(here)
dir <- here()
setwd(dir)


T_train = 70
G = 4
check.func <- function(u, tau)
{
  u * (tau - (u < 0) * 1)
}

for (T_test in 71:77) {
  # ---- Reload data and reset Xt every iteration --------------------------
  load("data/weibo.rda")
  N = nrow(Ymat)
  taus = c(0.1, 0.3, 0.5, 0.7, 0.9)
  
  X_tensor <- aperm(replicate(T_train, Xi, simplify = "array"), c(1, 3, 2))
  X_tensor_new <- array(0, dim = c(N, T_train, 10))
  X_tensor_new[, , 1:8] <- X_tensor
  Xt_all = Xt
  Xt = Xt[2:(T_train + 1), ]
  for (i in 1:N) {
    X_tensor_new[i, , 9:10] = Xt
  }
  p <- 10
  
  X_train = X_tensor_new
  Y_train = Ymat[, 1:(T_train + 1)]
  Y1_test = Ymat[, T_test]
  Y_test = Ymat[, T_test + 1]
  X_test = cbind(Xi, matrix(0, N, 2))
  X_test[, 9] = Xt_all[T_test + 1, 1]
  X_test[, 10] = Xt_all[T_test + 1, 2]
  
  
  W = Amat[1:N, 1:N]
  W = as.matrix(W)
  W = W / rowSums(W)
 
  # Same prediction helper as in step3_cf_selectGH.R, but specialized to a
  # single time point so member_G == member_H (one-way grouping in QGNAR).
  pred = function(theta_Gs,
                  alphabeta_GHs,
                  Y1,
                  W,
                  X,
                  member_G,
                  member_H) {
    G = nrow(alphabeta_GHs)
    H = ncol(alphabeta_GHs)
    N = length(Y1)
    WY1 <- W %*% Y1
    WY1 = matrix(WY1, ncol = 1)
    Yh <- list()
    for (h in 1:H) {
      Yh[[h]] = W %*% (Y1 * (member_H == h))
    }
    
    Y_pred <- numeric(N)
    
    for (g in 1:G)
    {
      ind_g = which(member_G == g)
      Y1_g = Y1[ind_g]
      Z1 <- as.vector(Y1_g)
      Zh <- NULL
      for (h in 1:H) {
        Zh <- cbind(Zh, as.vector((Yh[[h]])[ind_g, , drop = F]))
      }
      X1 = X[ind_g, ]
      Z_all <- cbind(Z1, X1[, , drop = T], Zh)
      
      theta = cbind(theta_Gs, alphabeta_GHs)
      
      Y_pred[ind_g] = Z_all %*% theta[g, ]
    }
    return(Y_pred)
  }
  
  # ---- QGNAR with G = 3 ---------------------------------------------------
  load(paste0("output/real_res/qnarg_G3-", T_train, ".rda"))
  G = 3
  GNAR3.Y_pred_taus = lapply(1:5, function(k) {
    pred_vec = pred(
      narg_mem_final$psi[[k]][, -(1:G)],# covariate + AR + intercept
      narg_mem_final$psi[[k]][, 1:G],# network-effect block (G columns)
      matrix(Y1_test, ncol = 1),
      W,
      X_test,
      narg_mem_final$member,
      narg_mem_final$member
    )
  })
  
  # ---- QGNAR with G = 4 ---------------------------------------------------
  load(paste0("output/real_res/qnarg_G4-", T_train, ".rda"))
  G = 4
  GNAR4.Y_pred_taus = lapply(1:5, function(k) {
    pred_vec = pred(
      narg_mem_final$psi[[k]][, -(1:G)],
      narg_mem_final$psi[[k]][, 1:G],
      matrix(Y1_test, ncol = 1),
      W,
      X_test,
      narg_mem_final$member,
      narg_mem_final$member
    )
  })
  
  
  save(
    Y_test,
    GNAR3.Y_pred_taus,
    GNAR4.Y_pred_taus,
    file = paste0(
      "output/real_res/narg_pred_train",
      T_train,
      "-test",
      T_test,
      ".rda"
    )
  )
  
  
}
