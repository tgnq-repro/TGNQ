# =============================================================================
# estimator_gqnar.R
#
# Implementation of the QGNAR (Grouped Quantile Network AutoRegression)
# benchmark used in Section 6 to compare against TGNQ.
#
# Public functions:
#
#   est.NARG.init(Ymat, W, G, taus, X_tensor, lambda, ntrial)
#       Build a pool of K-means initializations from per-node lasso quantile
#       pre-fits (uses three feature representations: AR coefficient, intercept,
#       network-effect signature).
#
#   est.NARG(Ymat, W, X_tensor, G, member, taus, ...)
#       Given memberships, fit group-level coefficients via L1-regularized
#       quantile regression (`quantreg::rq` with method = "lasso") for all taus
#       jointly. Returns psi (a list of G x (G+P+1) matrices, one per tau)
#       and the per-tau loss.
#
#   est.member(Ymat, W, X_tensor, psi, member, taus, ...)
#       Update each individual's membership by minimizing its quantile loss
#       given the current parameters. Two modes (Initial = TRUE / FALSE):
#       Initial uses only the node's own loss; otherwise the loss is summed
#       over all nodes that have node i as an in-neighbor.
#
#   est.NARG.member(Ymat, W, X_tensor, G, taus, ...)
#       Alternate between est.NARG and est.member until convergence, in a
#       two-stage scheme: a "local" stage (Initial = TRUE) followed by a
#       "global" stage (Initial = FALSE).
#
#   est.NARG.member.parallel(Ymat, W, X_tensor, G, taus, numCores, n.initial, ...)
#       Run est.NARG.member from up to n.initial K-means starts in parallel
#       (via doSNOW) and return the run with the smallest total loss.
#
#   gqnar(Ymat, W, X_tensor, G, taus, ...)
#       Convenience wrapper: serial multi-start + final refit.
#
# Helper:
#
#   check.func(u, tau) = u * (tau - 1{u < 0})
#       Standard quantile (pinball) loss.
# ============================================================================
library("CEoptim")
library(MASS)
library(parallel)
library(foreach)
library(doParallel)
library(doSNOW)
library(compiler)


check.func<-function(u, tau)
{
  u*(tau-(u<0)*1)
}

est.NARG.init <-function(Ymat, W, G, taus,X_tensor, lambda=0.01,ntrial=100)
{
  N = nrow(Ymat)
  Time = ncol(Ymat)
  Ymat_lag = Ymat[,-Time]
  Ymat1 = Ymat[,-1]
  
  beta_lag=beta_fix = matrix(0, N,1)
  fri_lag = matrix(0, N, N)

  
  for (i in 1:N)
  {
    fri_i = which(W[i,]!=0) # 
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
  member.int_set <- Member.int[!duplicated(totss),,drop=F]
  member.int_set

}




est.NARG <- function(Ymat, W, X_tensor, G, member, taus,lambda = 0, verbose = F, residual = F) {
  
  N <- nrow(Ymat)
  Time <- ncol(Ymat)
  Y_lag <- Ymat[, -Time]  
  Y <- Ymat[, -1]         
  
  P <- dim(X_tensor)[3]  
  psi <- matrix(0, nrow = G + P + 1, ncol = G)
  resi <- matrix(0, nrow = N, ncol = Time - 1)
  
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
      Xall_g <- cbind(Xall_g, as.vector(Wh %*% Y_lag[indh, , drop = F]))
    }
    
    Xall_g <- cbind(Xall_g, Ylag_g)
    
    #### Joint estimation
    Xg_tensor <- X_tensor[indg, , , drop = F]  
    XT_g <- do.call(rbind, lapply(1:(Time - 1), function(t) Xg_tensor[,t, , drop = T])) 
    Yall_g <- c(Y_g)
    
    Xmat <- cbind(Xall_g, XT_g)
    
    resrq = quantreg::rq(Yall_g ~ Xmat-1, tau=taus,method="lasso", lambda =0)
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




est.member <- function(Ymat, W, X_tensor, psi, member,taus, verbose = F, Initial = FALSE) {
  G <- nrow(psi[[1]])
  n_taus = length(psi)
  beta <- lapply(1:n_taus, function(k)psi[[k]][,1:G])
  nu <- lapply(1:n_taus, function(k)psi[[k]][,G+1])
  gamma <- lapply(1:n_taus, function(k)psi[[k]][,-(1:(1 + G))])
  
  Time <- ncol(Ymat)
  N <- nrow(Ymat)
  
  for (i in 1:N) {
    if (Initial) fri_i_relate <- c(i) else fri_i_relate <- c(i, which(W[, i] != 0))# 哪些个体关注了i
    
    loss_i <- rep(0, G)
    
    for (j in fri_i_relate) {
      fri_j <- which(W[j, ] != 0)
      dj <- length(fri_j)
      
      
      if (dj > 1) {
        Ymat_lag_frj <- t(Ymat[fri_j, -Time])
      } else {
        Ymat_lag_frj <- matrix(Ymat[fri_j, -Time], ncol = 1)
      }
      
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

est.NARG.member<-function(Ymat, W, X_tensor, G,taus, verbose = F, init_member = NULL,Maxit=100,tol=10^{-4})
{
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

  res = est.NARG(Ymat, W, X_tensor,G, member,taus, verbose = F)
  psi = res$psi
  loss = res$loss
  

  ###Update membership using individual loss
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
  
  member_loc <- member
  psi_loc <- psi
  ##Further update membership with full loss function
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




est.NARG.member.parallel<-function(Ymat, W, X_tensor, G, taus,numCores = 10, n.initial = 10,member.int=NULL,script_path = NULL)
{
  p = dim(X_tensor)[3]
  N = nrow(Ymat)

    narg_mem_list <- list()
    loss <- NULL
    member.int_set <- est.NARG.init(Ymat, W, G,taus, X_tensor, lambda=0.01,ntrial=100)
    
    ###The maximum number of initial values to try is 10, one can increase this number
    m <- min(dim(member.int_set), n.initial)
    
    numCores <- min(numCores,m)
    cl <- makeCluster(numCores)
    registerDoSNOW(cl)
    ###Create progress report
    pb <- txtProgressBar(max = m, style = 3)
    progress <- function(count) setTxtProgressBar(pb, count)
    opts <- list(progress = progress)
    
    narg_mem_list <- foreach(L = 1:m,.options.snow = opts,.combine = c, .verbose = F)%dopar%{
      source(script_path) 
      narg_mem <- tryCatch({
        est.NARG.member(Ymat, W, X_tensor, G, taus, verbose = F, init_member =  member.int_set[L,],
                        Maxit=100,tol=10^(-4))
      }, error = function(err)
      {
        print(c(paste("Est ERROR1:  ",err)))
        list(psi = matrix(0, nrow = G+p+1, ncol = G),
             member = rep(0, N),
             loss = 100000000000,
             psi_sd = matrix(0, nrow = G+p+1, ncol = G)) #  set loss to be a large number
      })
      

      list(narg_mem)
    }
    close(pb)
    #Stop the cluster
    stopCluster(cl)
    
    loss <- lapply(narg_mem_list, function(res)sum(res$loss))
    loss=unlist(loss)
    idx <- which.min(loss)
    print(min(loss))
    narg_mem_final = narg_mem_list[[idx]]
    print(narg_mem_final)

  return(narg_mem_final)
  
}

