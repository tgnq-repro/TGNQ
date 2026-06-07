# =============================================================================
# step9_inference.R
#
# Final inference + visualization for the chosen TGNQ model (T_train = 77,
# G = 4, H = 3):
#
#   1. Refine row-group memberships (Refine_G_parallel) and column-group
#      memberships (Refine_H_parallel), keeping only changes whose relative
#      loss improvement exceeds 1/sqrt(T - 1).
#   2. Re-estimate parameters at the refined memberships
#      (twmq.estimate_thetaGH.member.iterate).
#   3. Compute pointwise CIs for the parameters (twmq_ci).
#   4. Build three figures:
#        - output/figs/network_effect.pdf
#        - output/figs/covar_effect.pdf
#        - output/figs/activity.pdf
#
# Inputs : data/weibo.rda
#          output/real_res/tgnq_77-4-3.rda     (initial fit from Step 4)
# Outputs: output/real_res/rG43.rda  (refined member_G_refine_ge)
#          output/real_res/rH43.rda  (refined member_H_refine_ge)
#          output/real_res/r43-4-3.rda  (ge_refine + ci)
# =============================================================================
rm(list=ls())

library(twmq)
library(here)
dir <- here()
setwd(dir)

# ---- 1. Load data and assemble training tensor ------------------------------
load("data/weibo.rda")
N = nrow(Ymat)

T_train = 77
X_tensor <- aperm(replicate(T_train, Xi, simplify = "array"), c(1,3,2))
X_tensor_new <- array(0, dim = c(N, T_train, 10))
X_tensor_new[,,1:8] <- X_tensor
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

# ---- 2. Load the initial TGNQ fit -------------------------------------------
load("output/real_res/tgnq_77-4-3.rda")
check.func<-function(u, tau)
{
  u*(tau-(u<0)*1)
}

##Refine membership G
res_ge_update$theta_GH$Loss.fun = check.func

member_G_refine_ge <- res_ge_update$member_G

numCores=20
taus=seq(0.1,0.9,0.2)

# ---- 3. Build neighbor lists used by Refine_H_parallel ----------------------
#   FriendW[[i]]  = first-order in-neighbors of i (j with W[j, i] != 0)
#   FriendW2[[i]] = second-order in-neighbors of i
FriendW <-FriendW2 <- vector("list", N)
for(i in 1:N){
  FriendW[[i]] <- as.vector(which(W[,i]!=0))
  for (j in  FriendW[[i]]) {
    FriendW2[[i]] <- unique(c(FriendW2[[i]],as.vector(which(W[j,]!=0))))
  }
}


# ---- 4. Refine row-group memberships (G) -----------------------------------
obj_G_refine_ge <-  Refine_G_parallel(Y_train, aperm(X_train, c(1,3,2)), W, res_ge_update$member_G, res_ge_update$member_H, taus, res_ge_update$theta, method ="general", G=4, H=3,numCores = numCores)
##Check how many H memberships needs to be refined
idx_G <- which((obj_G_refine_ge$loss_old-obj_G_refine_ge$loss_new)/obj_G_refine_ge$loss_old>1/sqrt(T_train))
if(length(idx_G)>0) member_G_refine_ge[idx_G]=obj_G_refine_ge$g_r[idx_G]
save(member_G_refine_ge,file="output/real_res/rG43.rda")


# ---- 5. Refine column-group memberships (H), conditional on refined G ------
member_H_refine_ge <- res_ge_update$member_H
obj_H_refine_ge <-  Refine_H_parallel(Y_train, aperm(X_train, c(1,3,2)), W, member_G_refine_ge, res_ge_update$member_H, taus, res_ge_update$theta, method ="general", G=4, H=3, FriendW, FriendW2,numCores = numCores)
##Check how many H memberships needs to be refined
idx_H <- which((obj_H_refine_ge$loss_old-obj_H_refine_ge$loss_new)/obj_H_refine_ge$loss_old>1/sqrt(T_train))
if(length(idx_H)>0) member_H_refine_ge[idx_H]=obj_H_refine_ge$h_r[idx_H]
save(member_H_refine_ge,file="output/real_res/rH43.rda")

# ---- 6. Re-estimate parameters at refined memberships + CIs -----------------
ge_refine = twmq.estimate_thetaGH.member.iterate(Y_train,  aperm(X_train, c(1,3,2)), W,
                                                 member_G_refine_ge, member_H_refine_ge,
                                                 taus, conquer = F,
                                                 h_conquer = 0.05,
                                                 method ="general", n_iter.max = 100, verbose = F)
ci = twmq_ci(Y_train, aperm(X_train, c(1,3,2)), W,
                                             member_G_refine_ge, member_H_refine_ge,
                                             taus, conquer = F,
                                             h_conquer = 0.05,
                                             method ="general", n_iter.max = 100, verbose = F)
f  = paste0("output/real_res/r43-", 4, "-", 3, ".rda")
save(ge_refine,ci, file = f)





# =============================================================================
# 7. Visualization
#
#    Three figures are produced, in this order:
#
#    (a) Network-effect figure (output/figs/network_effect.pdf)
#        - 4 row groups x 3 column groups grid
#        - For each (g, h) cell, plots theta_{g,h}(tau) with 95% CIs vs tau
#
#    (b) Covariate-effect figure (output/figs/covar_effect.pdf)
#        - 4 row groups x 3 selected covariates grid
#        - Covariates shown: lagged response Y_{i,t-1}, "Public", "Weekend"
#
#    (c) Activity figure (output/figs/activity.pdf)
#        - Bar chart of average response by row group / column group
#
#    The three blocks below each follow the same template:
#      i)  load r43-4-3.rda
#      ii) flatten ci$point / ci$up / ci$low into a wide matrix
#      iii) build one ggplot per panel, store them in p_list
#      iv) assemble with patchwork and ggsave().
# =============================================================================
format_with_decimal <- function(x) {
  sprintf("%.4f", x)
}
library(here)
dir <- here()
setwd(dir)

load("output/real_res//r43-4-3.rda")
para = matrix(NA,20,16)
para = as.data.frame(para)
names(para) = c("g", "tau","theta_g1","theta_g2","theta_g3",
                "nu_g",
                "intercept","gamma_g1","gamma_g2","gamma_g3","gamma_g4","gamma_g5","gamma_g6","gamma_g7","gamma_t1","gamma_t2")
para[seq(1,20,5),1] = 1:4
para$tau = rep(seq(0.1,0.9,0.2), 4)

para_up = para
para_low = para

for (k in 1:5) {
  para[seq(k,20,5),3:5] = ci$point$alphabeta_GHs[[k]]
  para[seq(k,20,5),6:16] = ci$point$theta_Gs[[k]]
}

for (k in 1:5) {
  para_up[seq(k,20,5),3:5] = ci$up$alphabeta_GHs[[k]]
  para_up[seq(k,20,5),6:16] = ci$up$theta_Gs[[k]]
}

for (k in 1:5) {
  para_low[seq(k,20,5),3:5] = ci$low$alphabeta_GHs[[k]]
  para_low[seq(k,20,5),6:16] = ci$low$theta_Gs[[k]]
}


thetas = para[,c(3:5)]
thetas = cbind(thetas[1:5,],thetas[6:10,],thetas[11:15,],thetas[16:20,])
thetas = thetas
thetas = as.data.frame(thetas)

thetas_up = para_up[,c(3:5)]
thetas_up = cbind(thetas_up[1:5,],thetas_up[6:10,],thetas_up[11:15,],thetas_up[16:20,])
thetas_up = thetas_up
thetas_up = as.data.frame(thetas_up)


thetas_low = para_low[,c(3:5)]
thetas_low = cbind(thetas_low[1:5,],thetas_low[6:10,],thetas_low[11:15,],thetas_low[16:20,])
thetas_low = thetas_low
thetas_low = as.data.frame(thetas_low)

par(mfrow = c(4,3))
ylim = c(min(thetas_low),max(thetas_up))
for (i in 1:12) {
  plot(seq(0.1,0.9,0.2),thetas[,i], col = "#2C5284", type="l",ylim=ylim)
  lines(seq(0.1,0.9,0.2),thetas_up[,i], col = "grey", type="l")
  lines(seq(0.1,0.9,0.2),thetas_low[,i], col = "grey", type="l")
  abline(h=0,col="#DD4B44")
}


library(ggplot2)
ylim <- c(min(thetas_low), max(thetas_up))
p_list=list()

i =  2
df <- data.frame(x = seq(0.1, 0.9, 0.2),
                 y = thetas[,i],
                 upper = thetas_up[,i],
                 lower = thetas_low[,i])

p1=ggplot(df, aes(x = x, y = y)) +
  geom_line(color = "#2C5284",size=1.5)  +
  geom_line(aes(y = upper),linetype="dashed", color = "grey60",size=1) +
  geom_line(aes(y = lower),linetype="dashed", color = "grey60",size=1) +
  geom_hline(yintercept = 0, linetype="dashed",color = "#DD4B44") +
  ylim(ylim[1], ylim[2]) +
  ylab("Row-Group 1") +  # Change y-axis label to "G"
  xlab("")+
  theme_minimal()+
  theme(axis.text = element_blank(),  # Remove axis text
        axis.title.x = element_text(size = 40),  # Remove x-axis title
        axis.title.y = element_text(size = 20),
        axis.ticks = element_blank())  + # Remove minor gridlines
  ggtitle("Column-Group 1")+
  theme(plot.title = element_text(hjust = 0.5,size=20))
p1





i =  1
df <- data.frame(x = seq(0.1, 0.9, 0.2),
                 y = thetas[,i],
                 upper = thetas_up[,i],
                 lower = thetas_low[,i])

p2=ggplot(df, aes(x = x, y = y)) +
  geom_line(color = "#2C5284",size=1.5)  +
  geom_line(aes(y = upper),linetype="dashed", color = "grey60",size=1) +
  geom_line(aes(y = lower),linetype="dashed", color = "grey60",size=1) +
  geom_hline(yintercept = 0, linetype="dashed",color = "#DD4B44") +
  ylim(ylim[1], ylim[2]) +
  theme_minimal()+
  xlab("")+
  theme(axis.text = element_blank(),  # Remove axis text
        axis.title.x = element_text(size = 40),  # Remove x-axis title
        axis.title.y = element_blank(),
        axis.ticks = element_blank())  + # Remove minor gridlines
  ggtitle("Column-Group 2")+
  theme(plot.title = element_text(hjust = 0.5,size=20))
p2


i =  3
df <- data.frame(x = seq(0.1, 0.9, 0.2),
                 y = thetas[,i],
                 upper = thetas_up[,i],
                 lower = thetas_low[,i])

p3=ggplot(df, aes(x = x, y = y)) +
  geom_line(color = "#2C5284",size=1.5)  +
  geom_line(aes(y = upper),linetype="dashed", color = "grey60",size=1) +
  geom_line(aes(y = lower),linetype="dashed", color = "grey60",size=1) +
  geom_hline(yintercept = 0, linetype="dashed",color = "#DD4B44") +
  # ylim(ylim[1], ylim[2])+
  ylab("Network Effet") + 
  xlab("")+
  scale_y_continuous(limits = ylim,position = "right")+
  theme_minimal()+
  theme(axis.text.x= element_blank(),  # Remove axis text
        axis.title.x = element_text(size = 40),  # Remove x-axis title
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 20),  # Adjust size of y-axis title
        axis.text.y = element_text(size = 10))+  # Adjust size of y-axis labels)  + # Remove minor gridlines
  ggtitle("Column-Group 3")+
  theme(plot.title = element_text(hjust = 0.5,size=20))
p3

i=11
df <- data.frame(x = seq(0.1, 0.9, 0.2),
                 y = thetas[,i],
                 upper = thetas_up[,i],
                 lower = thetas_low[,i])
p4=ggplot(df, aes(x = x, y = y)) +
  geom_line(color = "#2C5284",size=1.5)  +
  geom_line(aes(y = upper),linetype="dashed", color = "grey60",size=1) +
  geom_line(aes(y = lower),linetype="dashed", color = "grey60",size=1) +
  geom_hline(yintercept = 0, linetype="dashed",color = "#DD4B44") +
  ylim(ylim[1], ylim[2]) +
  ylab("Row-Group 2") +  # Change y-axis label to "G"
  xlab("")+
  theme_minimal()+
  theme(axis.title.y = element_text(size = 20),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 30),  # Adjust size of y-axis title
        axis.text.x =  element_blank(),
        axis.text.y = element_blank())  # Remove minor gridlines
p4


i = 10
df <- data.frame(x = seq(0.1, 0.9, 0.2),
                 y = thetas[,i],
                 upper = thetas_up[,i],
                 lower = thetas_low[,i])

p5=ggplot(df, aes(x = x, y = y)) +
  geom_line(color = "#2C5284",size=1.5)  +
  geom_line(aes(y = upper),linetype="dashed", color = "grey60",size=1) +
  geom_line(aes(y = lower),linetype="dashed", color = "grey60",size=1) +
  geom_hline(yintercept = 0, linetype="dashed",color = "#DD4B44") +
  ylim(ylim[1], ylim[2]) +# Change y-axis label to "G"
  xlab("")+
  theme_minimal()+
  theme(axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 30),  # Adjust size of y-axis title
        axis.text.x =  element_blank(),
        axis.text.y = element_blank())  # Remove minor gridlines
p5

i =  12
df <- data.frame(x = seq(0.1, 0.9, 0.2),
                 y = thetas[,i],
                 upper = thetas_up[,i],
                 lower = thetas_low[,i])

p6=ggplot(df, aes(x = x, y = y)) +
  geom_line(color = "#2C5284",size=1.5) +
  geom_line(aes(y = upper),linetype="dashed", color = "grey60",size=1) +
  geom_line(aes(y = lower),linetype="dashed", color = "grey60",size=1) +
  geom_hline(yintercept = 0, linetype="dashed",color = "#DD4B44") +
  # ylim(ylim[1], ylim[2])+
  ylab("Network Effet") + 
  xlab("")+
  scale_y_continuous(limits = ylim,position = "right")+
  theme_minimal()+
  theme(axis.text.x=  element_blank(),  # Remove axis text
        axis.title.x = element_text(size = 30),  # Remove x-axis title
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 20),  # Adjust size of y-axis title
        axis.text.y = element_text(size = 10))  # Adjust size of y-axis labels)  + # Remove minor gridlines
p6




i=5
df <- data.frame(x = seq(0.1, 0.9, 0.2),
                 y = thetas[,i],
                 upper = thetas_up[,i],
                 lower = thetas_low[,i])
p7=ggplot(df, aes(x = x, y = y)) +
  geom_line(color = "#2C5284",size=1.5)  +
  geom_line(aes(y = upper),linetype="dashed", color = "grey60",size=1) +
  geom_line(aes(y = lower),linetype="dashed", color = "grey60",size=1) +
  geom_hline(yintercept = 0, linetype="dashed",color = "#DD4B44") +
  ylim(ylim[1], ylim[2]) +
  ylab("Row-Group 3") +  # Change y-axis label to "G"
  xlab("")+
  theme_minimal()+
  theme(axis.title.y = element_text(size = 20),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 30),  # Adjust size of y-axis title
        axis.text.x =  element_blank(),
        axis.text.y = element_blank())  # Remove minor gridlines
p7


i = 4
df <- data.frame(x = seq(0.1, 0.9, 0.2),
                 y = thetas[,i],
                 upper = thetas_up[,i],
                 lower = thetas_low[,i])

p8=ggplot(df, aes(x = x, y = y)) +
  geom_line(color = "#2C5284",size=1.5)  +
  geom_line(aes(y = upper),linetype="dashed", color = "grey60",size=1) +
  geom_line(aes(y = lower),linetype="dashed", color = "grey60",size=1) +
  geom_hline(yintercept = 0, linetype="dashed",color = "#DD4B44") +
  ylim(ylim[1], ylim[2]) +# Change y-axis label to "G"
  xlab("")+
  theme_minimal()+
  theme(axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 30),  # Adjust size of y-axis title
        axis.text.x =  element_blank(),
        axis.text.y = element_blank())  # Remove minor gridlines
p8

i = 6
df <- data.frame(x = seq(0.1, 0.9, 0.2),
                 y = thetas[,i],
                 upper = thetas_up[,i],
                 lower = thetas_low[,i])

p9=ggplot(df, aes(x = x, y = y)) +
  geom_line(color = "#2C5284",size=1.5) +
  geom_line(aes(y = upper),linetype="dashed", color = "grey60",size=1) +
  geom_line(aes(y = lower),linetype="dashed", color = "grey60",size=1) +
  geom_hline(yintercept = 0, linetype="dashed",color = "#DD4B44") +
  ylab("Network Effet") + 
  xlab("")+
  scale_y_continuous(limits = ylim,position = "right")+
  theme_minimal()+
  theme(axis.text.x=  element_blank(),  # Remove axis text
        axis.title.x = element_text(size = 25),  # Remove x-axis title
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 20),  # Adjust size of y-axis title
        axis.text.y = element_text(size = 10))  # Adjust size of y-axis labels)  + # Remove minor gridlines
p9





i =  8
df <- data.frame(x = seq(0.1, 0.9, 0.2),
                 y = thetas[,i],
                 upper = thetas_up[,i],
                 lower = thetas_low[,i])

p10=ggplot(df, aes(x = x, y = y)) +
  geom_line(color = "#2C5284",size=1.5)  +
  geom_line(aes(y = upper),linetype="dashed", color = "grey60",size=1) +
  geom_line(aes(y = lower),linetype="dashed", color = "grey60",size=1) +
  geom_hline(yintercept = 0, linetype="dashed",color = "#DD4B44") +
  ylim(ylim[1], ylim[2]) +
  ylab("Row-Group 4") +  # Change y-axis label to "G"
  xlab(expression(tau))+
  theme_minimal()+
  theme(axis.title.y = element_text(size = 20),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 25),  # Adjust size of y-axis title
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank())  # Remove minor gridlines
p10

i = 7
df <- data.frame(x = seq(0.1, 0.9, 0.2),
                 y = thetas[,i],
                 upper = thetas_up[,i],
                 lower = thetas_low[,i])

p11=ggplot(df, aes(x = x, y = y)) +
  geom_line(color = "#2C5284",size=1.5)  +
  geom_line(aes(y = upper),linetype="dashed", color = "grey60",size=1) +
  geom_line(aes(y = lower),linetype="dashed", color = "grey60",size=1) +
  geom_hline(yintercept = 0, linetype="dashed",color = "#DD4B44") +
  ylim(ylim[1], ylim[2]) +# Change y-axis label to "G"
  xlab(expression(tau))+
  theme_minimal()+
  theme(axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size = 25),  # Adjust size of y-axis title
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank())  # Remove minor gridlines
p11

i =  9
df <- data.frame(x = seq(0.1, 0.9, 0.2),
                 y = thetas[,i],
                 upper = thetas_up[,i],
                 lower = thetas_low[,i])

p12=ggplot(df, aes(x = x, y = y)) +
  geom_line(color = "#2C5284",size=1.5) +
  geom_line(aes(y = upper),linetype="dashed", color = "grey60",size=1) +
  geom_line(aes(y = lower),linetype="dashed", color = "grey60",size=1) +
  geom_hline(yintercept = 0, linetype="dashed",color = "#DD4B44") +
  ylab("Network Effet") + 
  xlab(expression(tau))+
  scale_y_continuous(limits = ylim,position = "right")+
  theme_minimal()+
  theme(axis.text.x= element_text(size = 10),  # Remove axis text
        axis.title.x = element_text(size = 25),  # Remove x-axis title
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 20),  # Adjust size of y-axis title
        axis.text.y = element_text(size = 10))  # Adjust size of y-axis labels)  + # Remove minor gridlines
p12


gridExtra::grid.arrange(grobs = list(p1,p2,p3,
                                     p4,p5,p6,
                                     p7,p8,p9,
                                     p10,p11,p12), ncol = 3)


# 2 1 3
# 11 10 12
# 5 4 6
# 8 7 9



library(patchwork)
p = (p1 + p2 + p3) /  (p4 + p5 + p6)/ (p7 + p8 + p9) /(p10 + p11 + p12)


ggsave("output/figs/network_effect.pdf", plot = p, width = 18, height = 24)




format_with_decimal <- function(x) {
  sprintf("%.4f", x)
}
library(here)
dir <- here()
setwd(dir)

load("output/real_res//r43-4-3.rda")
para = matrix(NA,20,16)
para = as.data.frame(para)
names(para) = c("g", "tau","theta_g1","theta_g2","theta_g3",
                "nu_g",
                "intercept","gamma_g1","gamma_g2","gamma_g3","gamma_g4","gamma_g5","gamma_g6","gamma_g7","gamma_t1","gamma_t2")
para[seq(1,20,5),1] = 1:4
para$tau = rep(seq(0.1,0.9,0.2), 4)

para_up = para
para_low = para

for (k in 1:5) {
  para[seq(k,20,5),3:5] = ci$point$alphabeta_GHs[[k]]
  para[seq(k,20,5),6:16] = ci$point$theta_Gs[[k]]
}

for (k in 1:5) {
  para_up[seq(k,20,5),3:5] = ci$up$alphabeta_GHs[[k]]
  para_up[seq(k,20,5),6:16] = ci$up$theta_Gs[[k]]
}

for (k in 1:5) {
  para_low[seq(k,20,5),3:5] = ci$low$alphabeta_GHs[[k]]
  para_low[seq(k,20,5),6:16] = ci$low$theta_Gs[[k]]
}


thetas = para[,6:16]
thetas = cbind(thetas[1:5,],thetas[6:10,],thetas[11:15,],thetas[16:20,])
thetas = thetas
thetas = as.data.frame(thetas)

thetas_up = para_up[,6:16]
thetas_up = cbind(thetas_up[1:5,],thetas_up[6:10,],thetas_up[11:15,],thetas_up[16:20,])
thetas_up = thetas_up
thetas_up = as.data.frame(thetas_up)


thetas_low = para_low[,6:16]
thetas_low = cbind(thetas_low[1:5,],thetas_low[6:10,],thetas_low[11:15,],thetas_low[16:20,])
thetas_low = thetas_low
thetas_low = as.data.frame(thetas_low)


library(ggplot2)

p_list=list()

i =  1
k=0
df <- data.frame(x = seq(0.1, 0.9, 0.2),
                 y = thetas[,i],
                 upper = thetas_up[,i],
                 lower = thetas_low[,i])
ylim <- c(min(thetas_low[,seq(i,44,11)]), max(thetas_up[,seq(i,44,11)]))
k=k+1;p_list[[k]]=ggplot(df, aes(x = x, y = y)) +
  geom_line(color = "#2C5284",size=1.5)  +
  geom_line(aes(y = upper),linetype="dashed", color = "grey60",size=1) +
  geom_line(aes(y = lower),linetype="dashed", color = "grey60",size=1) +
  geom_hline(yintercept = 0, linetype="dashed",color = "#DD4B44") +
  ylim(ylim[1], ylim[2]) +
  ylab("") +  # Change y-axis label to "G"
  ggtitle(expression(Y[i*","~t-1]))+
  theme_minimal()+
  theme(axis.text.x= element_blank(),  # Remove axis text
        axis.text.y= element_text(size = 8),
        axis.title.x = element_blank(),  # Remove x-axis title
        axis.title.y = element_text(size = 12),
        axis.ticks = element_blank())  + # Remove minor gridlines
  theme(plot.title = element_text(hjust = 0.5,size=20))
p_list[[k]]


i =  i+33
df <- data.frame(x = seq(0.1, 0.9, 0.2),
                 y = thetas[,i],
                 upper = thetas_up[,i],
                 lower = thetas_low[,i])

k=k+1;p_list[[k]]=ggplot(df, aes(x = x, y = y)) +
  geom_line(color = "#2C5284",size=1.5)  +
  geom_line(aes(y = upper),linetype="dashed", color = "grey60",size=1) +
  geom_line(aes(y = lower),linetype="dashed", color = "grey60",size=1) +
  geom_hline(yintercept = 0, linetype="dashed",color = "#DD4B44") +
  ylim(ylim[1], ylim[2]) +
  ylab("") +  # Change y-axis label to "G"
  theme_minimal()+
  theme(axis.text.x= element_blank(),  # Remove axis text
        axis.text.y= element_text(size = 8),
        axis.title.x = element_blank(),  # Remove x-axis title
        axis.title.y = element_text(size = 12),
        axis.ticks = element_blank())  + # Remove minor gridlines
  theme(plot.title = element_text(hjust = 0.5,size=20))
p_list[[k]]


i =  i-22
df <- data.frame(x = seq(0.1, 0.9, 0.2),
                 y = thetas[,i],
                 upper = thetas_up[,i],
                 lower = thetas_low[,i])

k=k+1;p_list[[k]]=ggplot(df, aes(x = x, y = y)) +
  geom_line(color = "#2C5284",size=1.5)  +
  geom_line(aes(y = upper),linetype="dashed", color = "grey60",size=1) +
  geom_line(aes(y = lower),linetype="dashed", color = "grey60",size=1) +
  geom_hline(yintercept = 0, linetype="dashed",color = "#DD4B44") +
  ylim(ylim[1], ylim[2]) +
  ylab("") +  # Change y-axis label to "G"
  theme_minimal()+
  theme(axis.text.x= element_blank(),  # Remove axis text
        axis.text.y= element_text(size = 8),
        axis.title.x = element_blank(),  # Remove x-axis title
        axis.title.y = element_text(size = 12),
        axis.ticks = element_blank())  + # Remove minor gridlines
  theme(plot.title = element_text(hjust = 0.5,size=20))
p_list[[k]]


i =  i+11
df <- data.frame(x = seq(0.1, 0.9, 0.2),
                 y = thetas[,i],
                 upper = thetas_up[,i],
                 lower = thetas_low[,i])

k=k+1;p_list[[k]]=ggplot(df, aes(x = x, y = y)) +
  geom_line(color = "#2C5284",size=1.5)  +
  geom_line(aes(y = upper),linetype="dashed", color = "grey60",size=1) +
  geom_line(aes(y = lower),linetype="dashed", color = "grey60",size=1) +
  geom_hline(yintercept = 0, linetype="dashed",color = "#DD4B44") +
  ylim(ylim[1], ylim[2]) +
  ylab("") +  # Change y-axis label to "G"
  xlab(expression(tau))+
  theme_minimal()+
  theme(axis.text.x= element_text(size = 10),  # Remove axis text
        axis.text.y= element_text(size = 8),
        axis.title.x = element_text(size = 25),  # Remove x-axis title
        axis.title.y = element_text(size = 12),
        axis.ticks = element_blank())  + # Remove minor gridlines
  theme(plot.title = element_text(hjust = 0.5,size=20))
p_list[[k]]



i =  9

df <- data.frame(x = seq(0.1, 0.9, 0.2),
                 y = thetas[,i],
                 upper = thetas_up[,i],
                 lower = thetas_low[,i])
ylim <- c(min(thetas_low[,seq(i,44,11)]), max(thetas_up[,seq(i,44,11)]))
k=k+1;p_list[[k]]=ggplot(df, aes(x = x, y = y)) +
  geom_line(color = "#2C5284",size=1.5)  +
  geom_line(aes(y = upper),linetype="dashed", color = "grey60",size=1) +
  geom_line(aes(y = lower),linetype="dashed", color = "grey60",size=1) +
  geom_hline(yintercept = 0, linetype="dashed",color = "#DD4B44") +
  ylim(ylim[1], ylim[2]) +
  ylab("") +  # Change y-axis label to "G"
  xlab("")+
  ggtitle("Public")+
  theme_minimal()+
  theme( # Remove axis text
    axis.text=  element_blank(),
    axis.title.x =  element_blank(),  # Remove x-axis title
    axis.title.y = element_text(size = 12),
    axis.text.y=element_text(size = 8),
    axis.ticks = element_blank())+   # Remove minor gridlines
  theme(plot.title = element_text(hjust = 0.5,size=20))
p_list[[k]]



i =  i+33
df <- data.frame(x = seq(0.1, 0.9, 0.2),
                 y = thetas[,i],
                 upper = thetas_up[,i],
                 lower = thetas_low[,i])

k=k+1;p_list[[k]]=ggplot(df, aes(x = x, y = y)) +
  geom_line(color = "#2C5284",size=1.5)  +
  geom_line(aes(y = upper),linetype="dashed", color = "grey60",size=1) +
  geom_line(aes(y = lower),linetype="dashed", color = "grey60",size=1) +
  geom_hline(yintercept = 0, linetype="dashed",color = "#DD4B44") +
  ylim(ylim[1], ylim[2]) +
  ylab("") +  # Change y-axis label to "G"
  xlab("")+
  theme_minimal()+
  theme( # Remove axis text
    axis.text=  element_blank(),
    axis.title.x =  element_blank(),  # Remove x-axis title
    axis.title.y = element_text(size = 12),
    axis.text.y=element_text(size = 8),
    axis.ticks = element_blank())   # Remove minor gridlines
p_list[[k]]

i =  i-22
df <- data.frame(x = seq(0.1, 0.9, 0.2),
                 y = thetas[,i],
                 upper = thetas_up[,i],
                 lower = thetas_low[,i])

k=k+1;p_list[[k]]=ggplot(df, aes(x = x, y = y)) +
  geom_line(color = "#2C5284",size=1.5)  +
  geom_line(aes(y = upper),linetype="dashed", color = "grey60",size=1) +
  geom_line(aes(y = lower),linetype="dashed", color = "grey60",size=1) +
  geom_hline(yintercept = 0, linetype="dashed",color = "#DD4B44") +
  ylim(ylim[1], ylim[2]) +
  ylab("") +  # Change y-axis label to "G"
  xlab("")+
  theme_minimal()+
  theme( # Remove axis text
    axis.text= element_blank(),
    axis.title.x =  element_blank(),  # Remove x-axis title
    axis.title.y = element_text(size = 12),
    axis.text.y=element_text(size = 8),
    axis.ticks = element_blank())   # Remove minor gridlines
p_list[[k]]

i =  i+11
df <- data.frame(x = seq(0.1, 0.9, 0.2),
                 y = thetas[,i],
                 upper = thetas_up[,i],
                 lower = thetas_low[,i])

k=k+1;p_list[[k]]=ggplot(df, aes(x = x, y = y)) +
  geom_line(color = "#2C5284",size=1.5)  +
  geom_line(aes(y = upper),linetype="dashed", color = "grey60",size=1) +
  geom_line(aes(y = lower),linetype="dashed", color = "grey60",size=1) +
  geom_hline(yintercept = 0, linetype="dashed",color = "#DD4B44") +
  ylim(ylim[1], ylim[2]) +
  ylab("") +  # Change y-axis label to "G"
  xlab(expression(tau))+
  theme_minimal()+
  theme( # Remove axis text
    axis.text= element_text(size = 10),
    axis.title.x =  element_text(size = 25),  # Remove x-axis title
    axis.title.y = element_text(size = 12),
    axis.text.y=element_text(size = 8),
    axis.ticks = element_blank())   # Remove minor gridlines
p_list[[k]]










i =  10

df <- data.frame(x = seq(0.1, 0.9, 0.2),
                 y = thetas[,i],
                 upper = thetas_up[,i],
                 lower = thetas_low[,i])
ylim <- c(min(thetas_low[,seq(i,44,11)]), max(thetas_up[,seq(i,44,11)])+0.01)
k=k+1;p_list[[k]]=ggplot(df, aes(x = x, y = y)) +
  geom_line(color = "#2C5284",size=1.5)  +
  geom_line(aes(y = upper),linetype="dashed", color = "grey60",size=1) +
  geom_line(aes(y = lower),linetype="dashed", color = "grey60",size=1) +
  geom_hline(yintercept = 0, linetype="dashed",color = "#DD4B44") +
  ylim(ylim[1], ylim[2]) +
  ylab("") +  # Change y-axis label to "G"
  xlab("")+
  ggtitle("Weekend")+
  theme_minimal()+
  theme( # Remove axis text
    axis.text=  element_blank(),
    axis.title.x =  element_blank(),  # Remove x-axis title
    axis.title.y = element_text(size = 12),
    axis.text.y=element_text(size = 8),
    axis.ticks = element_blank())+   # Remove minor gridlines
  theme(plot.title = element_text(hjust = 0.5,size=20))
p_list[[k]]



i =  i+33
df <- data.frame(x = seq(0.1, 0.9, 0.2),
                 y = thetas[,i],
                 upper = thetas_up[,i],
                 lower = thetas_low[,i])

k=k+1;p_list[[k]]=ggplot(df, aes(x = x, y = y)) +
  geom_line(color = "#2C5284",size=1.5)  +
  geom_line(aes(y = upper),linetype="dashed", color = "grey60",size=1) +
  geom_line(aes(y = lower),linetype="dashed", color = "grey60",size=1) +
  geom_hline(yintercept = 0, linetype="dashed",color = "#DD4B44") +
  ylim(ylim[1], ylim[2]) +
  ylab("") +  # Change y-axis label to "G"
  xlab("")+
  theme_minimal()+
  theme( # Remove axis text
    axis.text=  element_blank(),
    axis.title.x =  element_blank(),  # Remove x-axis title
    axis.title.y = element_text(size = 12),
    axis.text.y=element_text(size = 8),
    axis.ticks = element_blank())   # Remove minor gridlines
p_list[[k]]

i =  i-22
df <- data.frame(x = seq(0.1, 0.9, 0.2),
                 y = thetas[,i],
                 upper = thetas_up[,i],
                 lower = thetas_low[,i])

k=k+1;p_list[[k]]=ggplot(df, aes(x = x, y = y)) +
  geom_line(color = "#2C5284",size=1.5)  +
  geom_line(aes(y = upper),linetype="dashed", color = "grey60",size=1) +
  geom_line(aes(y = lower),linetype="dashed", color = "grey60",size=1) +
  geom_hline(yintercept = 0, linetype="dashed",color = "#DD4B44") +
  ylim(ylim[1], ylim[2]) +
  ylab("") +  # Change y-axis label to "G"
  xlab("")+
  theme_minimal()+
  theme( # Remove axis text
    axis.text= element_blank(),
    axis.title.x =  element_blank(),  # Remove x-axis title
    axis.title.y = element_text(size = 12),
    axis.text.y=element_text(size = 8),
    axis.ticks = element_blank())   # Remove minor gridlines
p_list[[k]]

i =  i+11
df <- data.frame(x = seq(0.1, 0.9, 0.2),
                 y = thetas[,i],
                 upper = thetas_up[,i],
                 lower = thetas_low[,i])

k=k+1;p_list[[k]]=ggplot(df, aes(x = x, y = y)) +
  geom_line(color = "#2C5284",size=1.5)  +
  geom_line(aes(y = upper),linetype="dashed", color = "grey60",size=1) +
  geom_line(aes(y = lower),linetype="dashed", color = "grey60",size=1) +
  geom_hline(yintercept = 0, linetype="dashed",color = "#DD4B44") +
  ylim(ylim[1], ylim[2]) +
  ylab("") +  # Change y-axis label to "G"
  xlab(expression(tau))+
  theme_minimal()+
  theme( # Remove axis text
    axis.text= element_text(size = 10),
    axis.title.x =  element_text(size = 25),  # Remove x-axis title
    axis.title.y = element_text(size = 12),
    axis.text.y=element_text(size = 8),
    axis.ticks = element_blank())   # Remove minor gridlines
p_list[[k]]


id = c(seq(1,12,4),seq(2,12,4),seq(3,12,4),seq(4,12,4))
gridExtra::grid.arrange(grobs = p_list[id], ncol = 3)



library(patchwork)

p=(p_list[[1]] | p_list[[5]] | p_list[[9]]) /
  (p_list[[2]] | p_list[[6]] | p_list[[10]]) /
  (p_list[[3]] | p_list[[7]] | p_list[[11]])/
  (p_list[[4]] | p_list[[8]] | p_list[[12]])


ggsave("output/figs/covar_effect.pdf", plot = p, width = 18, height = 24)






library(here)
dir <- here()
setwd(dir)
load("data/weibo.rda")
N = nrow(Ymat)

T_train = 77
# 不标准化###############################################################
X_tensor <- aperm(replicate(T_train, Xi, simplify = "array"), c(1,3,2))
X_tensor_new <- array(0, dim = c(N, T_train, 10))
X_tensor_new[,,1:8] <- X_tensor
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
###############################################################
X_tensor =X_train

load("output/real_res/rG43.rda")
load("output/real_res/rH43.rda")

s = lapply(c(1,4,2,3,2+4,1+4,3+4), function(i){
  if(i<5){
    X = X_tensor[member_G_refine_ge == i,,]
    s_g = apply(X,MARGIN = 3,mean)
  }else{
    X = X_tensor[member_H_refine_ge == i-4,,]
    s_g = apply(X,MARGIN = 3,mean)
  }
})


color = c("#DD4B44","#EF8A47", "#FDC95C","#FFE6B7", "#9CD6D8","#62AFCB","#427C9B")[7:1]

library(ggplot2)
library(gridExtra)

plots <- list()


y_g = unlist(lapply(c(1,4,2,3), function(g)mean(Ymat[member_G_refine_ge==g,])))
y_g = c(y_g,unlist(lapply(c(2,1,3), function(h)mean(Ymat[member_H_refine_ge==h,]))))



group_labels <- c(
  "Row-Group 1",
  "Row-Group 2",
  "Row-Group 3",
  "Row-Group 4",
  "Column-Group 1",
  "Column-Group 2",
  "Column-Group 3"
)


df_response <- data.frame(
  x = factor(group_labels, levels = rev(group_labels)),  # 因子水平反转
  values = y_g
)

p=ggplot(df_response, aes(x = x, y = values, fill = x)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color) +
  labs(title = "Activity", x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(size = 15),
    plot.title = element_text(
      hjust = 0.5,
      size = 20
    ),
    axis.text.x = element_text(size = 12),    # 显示并放大分组标签
    axis.text.y = element_text(color="black",size = 18), 
    axis.ticks.x = element_blank()
  ) +
  coord_flip()


ggsave("output/figs/activity.pdf", plot = p, width = 18, height = 24)

