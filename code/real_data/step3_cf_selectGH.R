# =============================================================================
# step3_cf_selectGH.R
#
# Rolling out-of-sample selection of H given G = 4 and T_train = 70.
#
# For H in {2,3,4,5} and test days T_test in {71,...,77}:
#   1. load the TGNQ fit from Step 1,
#   2. produce a one-step-ahead forecast for day (T_test+1),
#   3. compute the average quantile loss across the 5 quantiles,
#   4. save it to output/real_res/pred_train70-test<T>-G4-H<H>.rda.
#
# Finally, average the loss over test days for each H and print the best H.
# =============================================================================
rm(list=ls())
library(here)
dir <- here()
setwd(dir)

T_train = 70
G=4

# Quantile check (loss) function
check.func<-function(u, tau)
{
  u*(tau-(u<0)*1)
}

# ---- Loop over (T_test, H) --------------------------------------------------
for(T_test in 71:77){
  for (H in 2:5) {
    # Reload data inside the loop because some objects (Xt) are overwritten below
    load("data/weibo.rda")
    N = nrow(Ymat)
    taus=c(0.1,0.3,0.5,0.7,0.9)
    
    # Build training tensor
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
    
    # One-step-ahead test inputs at day T_test
    Y1_test = Ymat[,T_test] # predictor: response at time T_test
    Y_test = Ymat[,T_test+1] # truth:    response at time T_test+1
    X_test = cbind(Xi,matrix(0,N,2))
    X_test[,9] = Xt_all[T_test+1,1]
    X_test[,10] = Xt_all[T_test+1,2]
    
    # Row-normalized W
    W = Amat[1:N,1:N]
    W = as.matrix(W)
    W = W/rowSums(W)

    # ---- Helper: one-step-ahead prediction at one quantile k --------------------
    # Given group-level coefficients (theta_Gs and alphabeta_GHs) and memberships,
    # build the design matrix for each row group and return Yhat.
    pred = function(theta_Gs, alphabeta_GHs, Y1, W, X, member_G, member_H){
      G = nrow(alphabeta_GHs)
      H = ncol(alphabeta_GHs)
      N = length(Y1)
      WY1 <- W%*%Y1
      WY1 = matrix(WY1,ncol=1)
      Yh <- list()
      # Yh[[h]] = network-weighted lag restricted to column-group h
      for (h in 1:H) {
        Yh[[h]] = W%*%(Y1*(member_H==h))
      }
      
      Y_pred <- numeric(N)
      
      for (g in 1:G)
      {
        ind_g = which(member_G == g)
        Y1_g=Y1[ind_g]
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
    
    # Load the TGNQ fit for the current (G, H) and predict
    load(paste0("output/real_res/tgnq_",T_train,"-",G,"-",H,".rda"))
    t = T_train+1
    TGNQ.Y_pred_taus = lapply(1:5, function(k){
      pred_vec = pred(res_ge_update$theta_GH$theta_Gs[[k]], res_ge_update$theta_GH$alphabeta_GHs[[k]],
                      Y1_test, W, X_test,
                      res_ge_update$member_G, res_ge_update$member_H)
    })
    
    # Average quantile loss across τ (and across nodes, days)
    loss = lapply(1:5, function(k){
      check.func(Y_test - TGNQ.Y_pred_taus[[k]],tau = taus[k])
    })
    loss = mean(unlist(loss))
    save(loss,Y_test,TGNQ.Y_pred_taus,file=paste0("output/real_res/pred_train",T_train,"-test",T_test,"-G",G,"-H",H,".rda"))
    
  }
}

# ---- Aggregate loss over test days for each H and pick the best H ----------
loss_H = lapply(2:5, function(H){
  loss = lapply(71:77, function(T_test){
    load(paste0("output/real_res/pred_train",T_train,"-test",T_test,"-G",G,"-H",H,".rda"))
    loss
  })
  c(H,mean(unlist(loss)))
})
loss_H = do.call(rbind,loss_H)
loss_H
print(paste0("H=", c(2:5)[which.min(loss_H[,2])]))
