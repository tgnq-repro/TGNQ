# =============================================================================
# generate_data.R
# -----------------------------------------------------------------------------
# Baseline data-generating process (DGP) for the TGNQ model.
#
# This file is sourced by the simulation drivers in code/simulation/.
# It provides:
#   * random generation of group memberships,
#   * group-specific parameter functions alpha, beta, nu, gamma,
#   * two network generators (SBM and power-law),
#   * a panel response simulator simu.Y() that draws Y_{it} from the TGNQ model.
#
# The flag `additive` inside simu.Y() switches between
#   theta_{g h}(tau) = alpha_g(tau) + beta_h(tau)   (additive)
#   theta_{g h}(tau) = alpha_g(tau) * beta_h(tau)   (multiplicative)
# =============================================================================

library(poweRlaw)
library(Matrix)

# -----------------------------------------------------------------------------
# simu.member0
# Randomly assign a discrete group label to each of N nodes according to
# the given probability vector `ratios`.
# -----------------------------------------------------------------------------
simu.member0<-function(N, ratios)
{
  G = length(ratios)
  member0 = sample(1:G, N, prob = ratios, replace = T)
  return(member0)
}

# -----------------------------------------------------------------------------
# logit
# Auxiliary monotone link function (centered logistic transform), kept here
# in case it is needed when redefining parameter functions.
# -----------------------------------------------------------------------------
logit<-function(u)
{
  return(0.5* exp(u)/(1+exp(u))-0.5)
}

# -----------------------------------------------------------------------------
# alpha.func / beta.func / nu.func / gamma.func
#
# Group-specific parameter functions of a uniform-like input `u`.
# Each returns the parameter value at the quantile level implied by `u`,
# evaluated for group `g` (row groups) or `h` (column groups).
#
# In the simulation, `u` is drawn as standard normal and then transformed
# through pnorm / pt / punif so that the resulting parameter is a monotone
# function of u -- ensuring the conditional quantile representation in
# equation (2.1) of the paper is well-defined.
# -----------------------------------------------------------------------------
alpha.func<-function(u, g){
  if (g == 1)
    return(u*0+1)
  if (g == 2)
    return(exp(pt(u,df=4)*0.1)+0.2)
  if (g == 3)
    return(exp(pt(u+1,df=4)*0.1)+0.3)
}

beta.func<-function(u, h){
  if (h == 1)
    return(exp(punif(u,-6,6)-3))
  if (h == 2)
    return(exp(punif(u,-6,6)-1.5)*0.9+0.1)
  if (h==3)
    return(exp(punif(u,-6,6)-1)*0.9+0.1)
}

nu.func<-function(u, g, G){
  if (g == 1)
    return(pnorm(u-0.5)*0.2)
  if (g == 2)
    return(pnorm(u)*0.2)
  if (g == 3)
    return(pnorm(u+0.5)*0.2)
}

# Covariate effect for the p-th covariate, group g
gamma.func<-function(u, g, p){
  if (g == 1 & p==2)
    return(pnorm(u)*0.2)
  if (g == 2& p==2)
    return(pnorm(u+0.5)*0.2)
  if (g == 3&p==2)
    return(pnorm(u-0.5)*0.2)
  
  if (g == 1 & p==3)
    return(pt(u,df=4)*0.4)
  if (g == 2& p==3)
    return(pt(u+0.5,df=4)*0.4)
  if (g == 3&p==3)
    return(pt(u-0.5,df=4)*0.4)
}


# -----------------------------------------------------------------------------
# getPowerLawW
# Generate a directed power-law network with N nodes and tail index `alpha`.
#
# For each node i, draw an in-degree k_i from a discrete power-law(1, alpha),
# then sample k_i nodes uniformly at random to follow i. Nodes with zero
# in-degree are patched by randomly selecting 3 followees so that the resulting
# matrix has no zero rows. If `normalize = TRUE`, return the row-normalized
# weighting matrix W; otherwise return the raw adjacency A.
# -----------------------------------------------------------------------------
### get power-law network W
getPowerLawW <- function(N, alpha=2.5, normalize = F) {
  Nfollowers = rpldis(N, 1, alpha)*4  ### generate N random numbers following power-law(1, alpha): k1-kN
  A = sapply(Nfollowers, function(n) {
    ### for node i, randomly select ki nodes to follow it
    vec = rep(0, N)
    vec[sample(1:N, min(n, N))] = 1
    return(vec)
  })
  diag(A) = 0
  ind = which(rowSums(A) == 0)  ### in case some row sums are zero
  for (i in ind) {
    A[i, sample(setdiff(1:N, i), 3)] = 1  ### for those node, randomly select 3 followees
  }
  #A = as(A, "dgCMatrix")
  if (!normalize) 
    return(A)
  W = A/rowSums(A)
  return(W)
} 

# -----------------------------------------------------------------------------
# getBlockW
# Generate an undirected stochastic block model (SBM) network with N nodes
# divided into Nblock approximately equal blocks.
#
# Within-block edges occur with probability ~ 6 log(N) / N (high density);
# between-block edges occur with probability ~ log(N) / N (low density).
# Nodes with zero degree are patched by adding 3 random links.
# When `normalize = TRUE`, the function returns the row-normalized W.
# -----------------------------------------------------------------------------

getBlockW<-function(N, Nblock, normalize = T)                                                          ### get block network
{
  if (N%%Nblock==0){                                                                                   ### if N mod Nblock is integer
    isDiagList = rep(list(matrix(1, nrow = N/Nblock, ncol = N/Nblock)), Nblock)                        ### obtain the diagnal block list
    mList = rep(list(matrix(rbinom((N/Nblock)^2, size = 1, prob = 6*log(N)/N),                       ### generate following relations within the blocks
                            nrow = N/Nblock, ncol = N/Nblock)), Nblock)
  }
  else
  {
    isDiagList = rep(list(matrix(1, nrow = floor(N/Nblock), ncol = floor(N/Nblock))), Nblock-1)        ### if N mod Nblock is not integer
    isDiagList[[length(Nblock)]] = matrix(1, nrow = N%%Nblock, ncol = N%%Nblock)
    
    mList = rep(list(matrix(rbinom(floor(N/Nblock)^2, size = 1, prob = 2*log(N)/N),                  ### generate following relations within the blocks
                            nrow = floor(N/Nblock), ncol = floor(N/Nblock))), Nblock-1)
    mList[[Nblock]] = matrix(rbinom(floor(N/Nblock)^2, size = 1, prob = 2*log(N)/N),                 ### generate following relations within the blocks
                             nrow = floor(N/Nblock), ncol = floor(N/Nblock))
  }
  isDiag = bdiag(isDiagList)                                                                           ### combine the blocks in matrix
  offDiag = which(isDiag == 0, arr.ind = T)                                                            ### to calculate the index of the off digonal indexes
  mList = lapply(mList, function(M){
    ind = which(rowSums(M)==0)
    if (length(ind)>0)
      M[cbind(ind, sample(1:nrow(M), length(ind)))] = 1
    return(M)
  })
  bA = bdiag(mList)
  bA[offDiag] = rbinom(nrow(offDiag), size = 1, prob = log(N)/N)                                          ### people between blocks have 0.3 prob to follow
  bA = as.matrix(bA)
  upperInd = which(upper.tri(bA), arr.ind = T)
  
  ################ transform bA to be a symmetric matrix ##############################################
  bA[upperInd[,2:1]] = bA[upper.tri(bA)]
  diag(bA) = 0
  
  
  ind = which(rowSums(bA)==0)                                                                          ### in case some row sums are zero
  for (i in ind)
  {
    bA[i, sample(setdiff(1:N,i), 3)] = 1                                                               ### for those node, randomly select 3 followees
  }
  
  if (!normalize)
    return(bA)
  W = bA/rowSums(bA)                                                                                   ### row normalize bA
  W = as(W, "dgCMatrix")
  return(W)
}



# -----------------------------------------------------------------------------
# simu.Y
# Simulate the response panel Y from the TGNQ model.
#
# Arguments
#   W         : N x N row-normalized network weighting matrix.
#   member_G  : length-N vector of true row-group memberships in {1,...,G0}.
#   member_H  : length-N vector of true column-group memberships in {1,...,H0}.
#   X_tensor  : N x p x (T + 10) array of exogenous covariates; the first 10
#               periods are used as burn-in.
#   additive  : if TRUE, network effects are additive
#                 theta_{g h} = alpha_g + beta_h;
#               otherwise they are multiplicative
#                 theta_{g h} = alpha_g * beta_h.
#   verbose   : currently unused; reserved for future diagnostic printing.
#
# Returns
#   An N x T matrix of simulated responses (after dropping the 10-period burn-in).
# -----------------------------------------------------------------------------
simu.Y<-function(W, member_G, member_H, X_tensor, additive = F, verbose = F)
{
  N = nrow(W)
  p = dim(X_tensor)[2]
  Time1 = dim(X_tensor)[3]
  Time = Time1 - 10
  # generate U matrix
  U = matrix(rnorm(nrow(W)*(Time1), 0, 1), nrow = nrow(W))
  G = max(member_G); H = max(member_H)
  
  Ymat = matrix(0, N, Time1 + 1)
  for (j in 1:(Time1))
  {
    alpha_all = rep(0, N); beta_all = rep(0, N)
    nu_all = rep(0, N); gamma_all = matrix(0, N, p)
    # obtain parameters
    for (g in 1:G){
      ind_g = which(member_G==g)
      alpha_all[ind_g] = alpha.func(U[ind_g,j], g)
      nu_all[ind_g] = nu.func(U[ind_g,j], g, G)
      gamma_all[ind_g,] = cbind(gamma.func(U[ind_g,j], g , 2),
                                gamma.func(U[ind_g,j], g , 3))
    }
    
    beta_mat = matrix(0, N, N) # beta should be a matrix
    
    for (h in 1:H){
      ind_h = which(member_H==h)
      beta_mat[,ind_h] = matrix(rep(beta.func(U[,j], h), 
                                    length(ind_h)), nrow = N)
    }
    
    if (additive){
      WG = diag(alpha_all)%*%W+W*beta_mat
    }else{
      WG = diag(alpha_all)%*%(W*(beta_mat))
    }
    
    WGY = WG%*%Ymat[,j]
    Ymat[,j+1] = as.vector(WGY) + nu_all*Ymat[,j]+rowSums(X_tensor[,,j]*gamma_all)
  }
  
  
  return(Ymat[,-(1:11)])
}


