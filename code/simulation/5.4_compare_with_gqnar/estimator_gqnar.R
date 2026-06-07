# =============================================================================
# estimator_gqnar.R
# -----------------------------------------------------------------------------
# QGNAR (gqnar) implementation used as the competitor in Section 5.4.
#
# The QGNAR model is the group quantile network autoregression with a free
# G x G "follower" effect matrix beta_{g,h}, an own-lag coefficient nu_g, and
# a covariate effect gamma_g. For each quantile tau, the response of unit i
# at time t+1 is modeled as
#
#   Q_tau(Y_{i,t+1} | F_t) = sum_{h=1}^G beta_{g_i, h}(tau) * (W_{i,.} Y_{.,t})_h
#                            + nu_{g_i}(tau) * Y_{i,t}
#                            + X_{i,.,t} * gamma_{g_i}(tau).
#
# Functions in this file:
#   - check.func           : check loss for quantile regression
#   - est.NARG.init        : multi-start k-means initialization of memberships
#   - est.NARG             : group-parameter update via quantreg::rq
#   - est.member           : individual-membership update by minimizing loss
#   - est.NARG.member      : two-stage iterative refinement (local + full loss)
#   - gqnar                : main entry point used by simulator.R
# =============================================================================
library("CEoptim")
library(MASS)
library(parallel)
library(foreach)
library(doParallel)
library(doSNOW)
library(compiler)

# Standard check-loss function for quantile regression at level `tau`.
check.func<-function(u, tau)
{
  u*(tau-(u<0)*1)
}

# -----------------------------------------------------------------------------
# est.NARG.init
# -----------------------------------------------------------------------------
# Multi-start k-means initialization of unit memberships for QGNAR.
#
# For each unit i, runs a per-individual L1-penalized median quantile
# regression of Y_{i,t} on (Y_{i,t-1}, W_{i,.} Y_{.,t-1}), then collects
# three kinds of per-unit features:
#   (type 1) the own-lag coefficient,
#   (type 2) the intercept (fixed effect),
#   (type 3) a (1 + G^2)-dim feature vector combining the own-lag coefficient
#            with G^2 averaged follower coefficients (clustered via kmeans).
# For each feature type and each of `ntrial` k-means restarts, runs
# kmeans(centers = G) to produce a candidate membership and stores the
# unique ones (deduplicated by total within-cluster SS).
#
# Returns a matrix whose rows are candidate initial memberships.
# -----------------------------------------------------------------------------
est.NARG.init <-function(Ymat, W, G, taus,X_tensor, lambda=0.01,ntrial=100)
{
  N = nrow(Ymat)
  Time = ncol(Ymat)
  Ymat_lag = Ymat[,-Time]
  Ymat1 = Ymat[,-1]
  
  beta_lag=beta_fix = matrix(0, N,1)
  fri_lag = matrix(0, N, N)

  # ---- Per-unit lasso quantile regression at the median (tau = 0.5) ------
  for (i in 1:N)
  {
    fri_i = which(W[i,]!=0) # i的朋友
    di = length(fri_i)
    if (di>1)
      X = cbind(Ymat_lag[i,], t(Ymat_lag[fri_i,])/di)  else
        X = cbind(Ymat_lag[i,], (Ymat_lag[fri_i,])/di)
    
    resrq = quantreg::rq(Ymat_lag[i,] ~ X, tau=0.5,  method = "lasso", 
                         lambda = lambda*sqrt(mean(diag(t(X)%*%X))))
    beta_lag[i, 1] <- resrq$coefficients[2]
    beta_fix[i, 1] <- resrq$coefficients[1]
    fri_lag[i, fri_i] = resrq$coefficients[-(1:2)]
  }
  
  # ---- Build a (G^2)-clustered representation of the follower effects ----
  fri_vec = fri_lag[fri_lag!=0]
  km_fri = kmeans(fri_vec, centers = G^2)
  fri_cl = matrix(0, N, N)
  fri_cl[fri_lag!=0] = km_fri$cluster
  fri_G2 = sapply(1:N, function(i){
    sapply(1:(G^2), function(g) {
      indg = which(fri_cl[i,]==g)
      if (length(indg)==0)
        return(0)
      else{
        return(mean(fri_lag[i,indg]))
      }
    })
  })
  fri_G2 = t(fri_G2)
  km_mat = cbind(beta_lag, fri_G2)
  
  # ---- Multi-start k-means over three feature types ----------------------
  Member.int <- NULL
  totss <-NULL
  for(type in c(1,2,3)){
    for(i in 1:ntrial){
      if(type==1)  km = kmeans(beta_lag, centers = G,iter.max = 100,nstart = 10) else
        if(type==2)  km = kmeans(beta_fix, centers = G,iter.max = 100,nstart = 10) else
          if(type==3) {
            km_fri = kmeans(fri_vec, centers = G^2,nstart=10)
            fri_cl = matrix(0, N, N)
            fri_cl[fri_lag!=0] = km_fri$cluster
            fri_G2 = sapply(1:N, function(i){
              sapply(1:(G^2), function(g) {
                indg = which(fri_cl[i,]==g)
                if (length(indg)==0)
                  return(0)
                else{
                  return(mean(fri_lag[i,indg]))
                }
              })
            })
            fri_G2 = t(fri_G2)
            km_mat = cbind(beta_lag, fri_G2)
            km = kmeans(km_mat, centers = G,nstart=10)
          }
        Member.int <- rbind(Member.int,km$cluster)
        totss <- rbind(totss,c(km$totss,km$tot.withinss))
    }
  }
  
  # Deduplicate candidate memberships by their (totss, within-ss) summary.
  totss <- round(totss,6)
  member.int_set <- Member.int[!duplicated(totss),,drop=F]
  member.int_set
}


# -----------------------------------------------------------------------------
# est.NARG
# -----------------------------------------------------------------------------
# Given a membership vector, fit the QGNAR group-level parameters jointly
# across all quantile levels in `taus`.
#
# For each group g:
#   - aggregate Y_{j,t-1} over followees within each group h via W, yielding
#     a G-column design (Y_{i,t-1}'s contribution from group-h neighbors),
#   - append the own-lag Y_{i,t-1} and the covariate slice X_{i,.,t-1},
#   - fit quantreg::rq(...) jointly at all `taus`.
#
# Returns:
#   psi[[k]] : G x (G + 1 + P) parameter matrix at quantile tau_k,
#              columns labeled (fri1, ..., friG, Y_lag, X1, ..., XP),
#   loss[k]  : total check loss at tau_k.
# -----------------------------------------------------------------------------

est.NARG <- function(Ymat, W, X_tensor, G, member, taus,lambda = 0, verbose = F, residual = F) {
  
  N <- nrow(Ymat)
  Time <- ncol(Ymat)
  Y_lag <- Ymat[, -Time]  # drop last time point
  Y <- Ymat[, -1]         # drop first time point
  
  P <- dim(X_tensor)[3]   # number of covariates
  psi <- matrix(0, nrow = G + P + 1, ncol = G)
  resi <- matrix(0, nrow = N, ncol = Time - 1)
  
  # Build the per-group design matrix and fit quantreg::rq jointly.
  ft <- function(g) {
    indg <- which(member == g)
    if (length(indg) == 0) {
      return(c(rep(0, 2 * (G + P + 1) + 1)))
    }
    Ylag_g <- c(Y_lag[indg, ])
    Y_g <- c(Y[indg, ])
    
    Xall_g <- NULL
    W_g <- W[indg, , drop = F]
    for (h in 1:G) {
      indh <- which(member == h)
      Wh <- W_g[, indh, drop = F]
      # Sum (W . Y_lag) restricted to group-h followees.
      Xall_g <- cbind(Xall_g, as.vector(Wh %*% Y_lag[indh, , drop = F]))
    }
    
    Xall_g <- cbind(Xall_g, Ylag_g)
    
    # ---- Joint estimation across all `taus` --------------------------------
    Xg_tensor <- X_tensor[indg, , , drop = F] 
    XT_g <- do.call(rbind, lapply(1:(Time - 1), function(t) Xg_tensor[,t, , drop = T])) 
    Yall_g <- c(Y_g)
    
    Xmat <- cbind(Xall_g, XT_g)
    
    resrq = quantreg::rq(Yall_g ~ Xmat-1, tau=taus)
    psi_g = resrq$coefficients
    loss_g = unlist(lapply(1:length(taus),function(k)sum(check.func(resrq$residuals[,k], taus[k]))))
    list(psi_g=psi_g, loss_g=loss_g)
  }
  
  tmp <- lapply(1:G, ft)
  psi = lapply(1:length(taus), function(k){
    psi_g = lapply(1:G, function(g){
      vec = tmp[[g]]$psi_g[,k]
      names(vec) = c(paste0("fri",1:G), "Y_lag", paste0("X", 1:P))
      vec
    })
    do.call(rbind, psi_g)
  })
  
  loss <- colSums(do.call(rbind,lapply(1:G, function(g){tmp[[g]]$loss_g})))
  
  return(list(psi = psi, loss = loss))

}



# -----------------------------------------------------------------------------
# est.member
# -----------------------------------------------------------------------------
# Update the membership of each individual i. For each candidate group g,
# compute the check-loss contribution of i (and, if Initial = FALSE, of all
# units that follow i) under the assignment member[i] = g, summed over all
# quantile levels. Pick g* = argmin loss.
#
# Initial = TRUE  : use only i's own loss   (fast; used in the first stage).
# Initial = FALSE : add the loss of every j that follows i (j with W[j,i] != 0).
#
# After the pass, any group that became empty is reseeded by splitting the
# currently largest group in half.
# -----------------------------------------------------------------------------
est.member <- function(Ymat, W, X_tensor, psi, member,taus, verbose = F, Initial = FALSE) {
  G <- nrow(psi[[1]])
  n_taus = length(psi)
  beta <- lapply(1:n_taus, function(k)psi[[k]][,1:G])
  nu <- lapply(1:n_taus, function(k)psi[[k]][,G+1])
  gamma <- lapply(1:n_taus, function(k)psi[[k]][,-(1:(1 + G))])
  
  Time <- ncol(Ymat)
  N <- nrow(Ymat)
  
  for (i in 1:N) {
    # In "Initial" mode, only consider unit i's own loss; otherwise add the
    # losses of every unit j that follows i (W[j, i] != 0).
    if (Initial) fri_i_relate <- c(i) else fri_i_relate <- c(i, which(W[, i] != 0))# 哪些个体关注了i
    
    loss_i <- rep(0, G)
    
    for (j in fri_i_relate) {
      # j's followees:
      fri_j <- which(W[j, ] != 0)
      dj <- length(fri_j)
      
      
      if (dj > 1) {
        Ymat_lag_frj <- t(Ymat[fri_j, -Time])
      } else {
        Ymat_lag_frj <- matrix(Ymat[fri_j, -Time], ncol = 1)
      }
      
      # Try each candidate group g for unit i.
      for (g in 1:G) {
        member[i] <- g
        g_j <- member[j]
        
        X_tensor_j <- X_tensor[j, , , drop = T]
        Yj_pred <- lapply(1:n_taus, function(k){
          Ymat_lag_frj %*% beta[[k]][g_j, member[fri_j]] / dj +
            Ymat[j, -Time] * nu[[k]][g_j] + X_tensor_j%*%gamma[[k]][g_j,]
        })
          
          
        Yj <- Ymat[j, -1]
        
        loss = lapply(1:n_taus, function(k){
          sum(check.func(Yj - Yj_pred[[k]], taus[k]))
        })
        loss_i[g] <- loss_i[g] + sum(unlist(loss))
      }
    }
    member[i] <- which.min(loss_i)
  }
  
  # Re-seed any empty groups by splitting the largest group in half.
  ID <- unique(names(table(member)))
  while (length(ID) < G) {
    g_max <- (names(table(member)))[which.max(table(member))]
    g_new <- which(!((1:G) %in% ID))[1]
    idx <- which(member == g_max)
    member[sample(idx, size = length(idx) / 2)] <- g_new
    ID <- unique(member)
  }
  
  return(member)
}


# -----------------------------------------------------------------------------
# gqnar
# -----------------------------------------------------------------------------
# Main entry point used by simulator.R.
#
# Strategy:
#   1. Get a pool of candidate initial memberships via est.NARG.init.
#   2. For each candidate (up to 100), run est.NARG.member and record loss.
#   3. Re-run est.NARG.member from the best-scoring candidate.
# -----------------------------------------------------------------------------
gqnar = function(Ymat, W, X_tensor, G,taus, verbose = F, Maxit=100,tol=10^{-4}){
  member.int_set = est.NARG.init(Ymat, W, G, taus,X_tensor, lambda=0.01,ntrial=100)
  m <- min(dim(member.int_set), 100)
  loss = NULL
  for (L in 1:m) {
    loss = c(loss,sum(est.NARG.member(Ymat, W, X_tensor, G,taus, 
                                      init_member=member.int_set[L,], verbose = verbose, Maxit=Maxit, tol=tol)$loss))
  }
  member = member.int_set[which.min(loss),]
  res = est.NARG.member(Ymat, W, X_tensor, G,taus,  init_member=member,verbose = verbose, Maxit=Maxit, tol=tol)
  res
}


# -----------------------------------------------------------------------------
# est.NARG.member
# -----------------------------------------------------------------------------
# Two-stage iterative refinement.
#   First stage  (Initial = TRUE) : update membership using each individual's
#                                   OWN local loss only (fast and stable).
#   Second stage (Initial = FALSE): update membership using the FULL loss
#                                   (individual's loss + their followers').
# Each stage iterates until the maximum relative change in `psi` falls below
# `tol` (averaged across quantile levels) or `Maxit` iterations are reached.
# -----------------------------------------------------------------------------
est.NARG.member<-function(Ymat, W, X_tensor, G,taus, verbose = F, init_member = NULL,Maxit=100,tol=10^{-4})
{
  # Degenerate case: G = 1 (no clustering to do).
  if (G == 1)
  {
    member = rep(1, nrow(Ymat))
    res = est.NARG(Ymat, W, X,G, member,taus, verbose = F)
    psi = res$psi
    loss = res$loss
    return(list(psi = psi, member = member,
                psi_sd = res$psi_sd,
                psi_loc = psi,
                member_loc=member, 
                loss = loss))
  }
  

  member = init_member

  
  # Initial parameter fit at `init_member`.
  res = est.NARG(Ymat, W, X_tensor,G, member,taus, verbose = F)
  psi = res$psi
  loss = res$loss
  

  # ----- First stage: local-loss membership updates -------------------------
  del = 1
  iter = 1
  while(del>tol & iter <= Maxit)
  {
    member = est.member(Ymat, W,  X_tensor, psi, member,taus, Initial = T)
    res = est.NARG(Ymat, W, X_tensor,G, member,taus, verbose = F)
    
    psi_new = res$psi
    loss = res$loss
    
    del = lapply(1:length(taus), function(k){
      max(abs(psi[[k]] - psi_new[[k]])/sqrt(mean(psi[[k]]^2)))
    })
    del = mean(unlist(del))
    psi = psi_new
    iter = iter + 1
    if (verbose)
    {
      cat("First Round  ", iter, "\t", loss, "\t", del,  "\n")
    }
    
  }
  
  if (verbose)
  {
    cat("\n\n")
  }
  
  # Snapshot of the local-loss-only solution.
  member_loc <- member
  psi_loc <- psi
  
  # ----- Second stage: full-loss membership updates -------------------------
  del = 1
  iter = 1
  while(del>tol & iter <= Maxit)
  {
    member = est.member(Ymat, W,  X_tensor, psi, member,taus, Initial = F)
    res = est.NARG(Ymat, W, X_tensor,G, member,taus, verbose = F)
    
    psi_new = res$psi
    loss = res$loss
    
    del = lapply(1:length(taus), function(k){
      max(abs(psi[[k]] - psi_new[[k]])/sqrt(mean(psi[[k]]^2)))
    })
    del = mean(unlist(del))
    psi = psi_new
    iter = iter + 1
    
    if (verbose)
    {
      cat("Second Round  ", iter, "\t", loss, "\t", del, "\n")
    }
  }
  return(list(psi = psi, member = member,
              psi_loc = psi_loc,
              member_loc=member_loc, 
              loss = loss))
}


