# =============================================================================
# step2_output_tex.R
# -----------------------------------------------------------------------------
# Post-processing for the Section 5.4 (TGNQ vs QGNAR) experiment.
#
# Iterates over (block, additive, case) and, within each design:
#   * rebuilds the TRUE autoregressive matrix B0(tau) (one per quantile tau)
#     from the oracle group parameters,
#   * for each seed file in output/res_gqnar/res_<block|power>_<add|mul>/,
#     constructs the TGNQ-implied and QGNAR-implied autoregressive matrices
#     B_TGNQ(tau) and B_QGNAR(tau) from the saved estimates,
#   * computes the RMSE of the NONZERO entries (i.e. positions where W != 0)
#     for each method and each tau,
#   * averages across seeds,
#   * assembles a wide LaTeX table and prints it via xtable.
# =============================================================================

rm(list=ls())


library(here)
dir <- here()
setwd(dir)


library(twmq)

# -----------------------------------------------------------------------------
# Aggregate results across (block, additive, case)
# -----------------------------------------------------------------------------
# `res_all` is a nested list:
#   res_all[[block]]   ($block in {TRUE, FALSE})
#     [[add]]          ($add   in {FALSE, TRUE})
#       a 4 x (n_taus * 2) matrix:
#         rows = case (1..4),
#         columns = (TGNQ_tau1, QGNAR_tau1, TGNQ_tau2, QGNAR_tau2, ...).
res_all = lapply(c(T,F), function(block){
  res = lapply(c(F,T), function(add){
    temp = lapply(1:4, function(case){
      
      # Per-design setup -----------------------------------------------------
      source("code/generator/generate_data.R")
      if(add){source("code/generator/additive_parameter.R")}
      
      
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
      set.seed(123)
      p = 2
      
      # Network (same construction as in simulator.R) -----------------------
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
      
      
      # True group-level parameters -----------------------------------------
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
      
      # Locate all seed files for this design ------------------------------
      files = list.files(paste0("output/res_gqnar/res_", ifelse(block,"block_","power_"), ifelse(add,"add/","mul/")), 
                         pattern = paste0("-", case, "-"),
                         full.names = T)
      # Load one snapshot to obtain member_G0 / member_H0 used to build B0.
      load(files[1])
      
      # ---- Build the TRUE B0(tau) (per tau), keeping only W != 0 entries -
      B0_list = lapply(1:5, function(k){
        B0_mat = matrix(0,N,N)
        for (i in 1:N) {
          for(j in 1:N){
            if(i==j){
              B0_mat[i,j] = nu0[member_G0[i],k]
            }else{
              if(add){
                B0_mat[i,j] = (alpha0[member_G0[i],k] + beta0[member_H0[j],k])*W[i,j]
              }else{
                B0_mat[i,j] = (alpha0[member_G0[i],k]*beta0[member_H0[j],k])*W[i,j]
              }
              
            }
          }
        }
        # Vectorize keeping only entries where W != 0.
        B0_mat = c(B0_mat[W!=0])
        B0_mat
      })
      
      # ---- Per-seed RMSE of B_TGNQ(tau) and B_QGNAR(tau) -----------------
      res=lapply(files, function(f){
        load(f)
        
        # Implied TGNQ autoregressive matrix.
        # Diagonal: theta_Gs[[k]][g_i, 1] = nu_{g_i}(tau_k).
        # Off-diagonal: alphabeta_GHs[[k]][g_i, h_j] * W_{ij}.
        Btgnq_list = lapply(1:5, function(k){
          B_mat = matrix(0,N,N)
          for (i in 1:N) {
            for(j in 1:N){
              if(i==j){
                B_mat[i,j] = res_ge_update1$theta_GH$theta_Gs[[k]][res_ge_update1$member_G[i],1]
              }else{
                B_mat[i,j] = W[i,j]*res_ge_update1$theta_GH$alphabeta_GHs[[k]][res_ge_update1$member_G[i],res_ge_update1$member_H[j]]
              }
            }
          }
          B_mat = c(B_mat[W!=0])
          B_mat
          
        })
        
        # Implied QGNAR autoregressive matrix.
        # Per QGNAR's psi layout: psi[[k]][g, 1:G] = beta_{g,.},
        #                          psi[[k]][g, G+1] = nu_g,
        # so for i = j we use the (G+1)-th column (the Y_lag coefficient).
        # NOTE: this code uses psi[[k]][..., 3] explicitly which corresponds
        # to column G+1 = 3 when G = 2; if G changes this needs updating.
        Bnarg_list = lapply(1:5, function(k){
          B_mat = matrix(0,N,N)
          for (i in 1:N) {
            for(j in 1:N){
              if(i==j){
                B_mat[i,j] = narg$psi[[k]][narg$member[i],3]
              }else{
                B_mat[i,j] = W[i,j]*narg$psi[[k]][narg$member[i],narg$member[j]]
              }
            }
          }
          B_mat = c(B_mat[W!=0])
          B_mat
          
        })
        
        # Per-tau RMSEs: c(TGNQ_tau, QGNAR_tau).
        rmes_tau = lapply(1:5, function(k){
          c(sqrt(mean((B0_list[[k]]-Btgnq_list[[k]])^2)),sqrt(mean((B0_list[[k]]-Bnarg_list[[k]])^2)))
        })
        unlist(rmes_tau)
      })
      
      # Average over seeds for this (block, add, case).
      colMeans(do.call(rbind, res))
      
    })
    
    temp = do.call(rbind,temp)
    temp
  })
  res
})



# =============================================================================
# Assemble the final wide table and print LaTeX
# =============================================================================

# (N, T, tau) header columns (20 rows total = 4 cases x 5 taus).
NTt = matrix(NA,20,3)
NTt[seq(1,20,5),1]=c(100,100,100,200)
NTt[seq(1,20,5),2]=c(50,100,200,200)
NTt[,3] = seq(0.1,0.9,0.2)

# Block-network half: take additive (res_all[[1]][[2]]) and multiplicative
# (res_all[[1]][[1]]) results, reshape each to a (20 x 2) (TGNQ, QGNAR) frame,
# and stack horizontally with the (N, T, tau) header.
df_block = cbind(cbind(NTt,t(matrix(t(res_all[[1]][[2]]),nrow=2))),cbind(NTt,t(matrix(t(res_all[[1]][[1]]),nrow=2))))
df_power = cbind(cbind(NTt,t(matrix(t(res_all[[2]][[2]]),nrow=2))),cbind(NTt,t(matrix(t(res_all[[2]][[1]]),nrow=2))))
df = as.data.frame(cbind(df_block,df_power))
names(df) = c("N","T","tau", "TGNQ", "NARG", "N","T","tau", "TGNQ.1", "NARG.1", "N","T","tau", "TGNQ.2", "NARG.2", "N","T","tau", "TGNQ.3", "NARG.3")

# Scale RMSEs to percent for readability in the table.
df$TGNQ=df$TGNQ*100
df$NARG=df$NARG*100
df$TGNQ.1=df$TGNQ.1*100
df$NARG.1=df$NARG.1*100
df$TGNQ.2=df$TGNQ.2*100
df$NARG.2=df$NARG.2*100
df$TGNQ.3=df$TGNQ.3*100
df$NARG.3=df$NARG.3*100

cols_to_format <- c("TGNQ", "NARG", "TGNQ.1", "NARG.1", "TGNQ.2", "NARG.2", "TGNQ.3", "NARG.3")

# Format RMSE cells to two decimals; NAs are kept as NA for now.
for (col in cols_to_format) {
  df[[col]] <- ifelse(
    is.na(df[[col]]), 
    NA, 
    formatC(as.numeric(df[[col]]), format = "f", digits = 2)
  )
}

# Coerce everything to character and replace NA by empty strings, so the LaTeX
# table shows blank cells where (N, T) headers repeat.
df_char <- as.data.frame(lapply(df, function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}), stringsAsFactors = FALSE)


print(df_char, row.names = FALSE)

library(xtable)
# Drop the duplicated (N, T, tau) header columns from the second/third/fourth
# blocks (we only need them once on the left).
print(xtable(df_char[,-c(6,7,8,11,12,13,16,17,18)]), include.rownames = FALSE, NA.string = "")
