# =============================================================================
# step2_QIC_plot.R
#
# Compute the QIC criterion across the (G, H) grid produced by Step 1, and
# plot QIC against G, faceted by H (and vice versa).
#
# Input : output/real_res/tgnq_70-<G>-<H>.rda for (G, H) in {1,...,5}^2
# Output: output/figs/QIC_plot.pdf
#
# QIC formula (matching Section 6 of the paper):
#   QIC(G, H) = log(loss_bar(G, H)) + lambda_{NT} * df(G, H)
#   lambda_NT = N^{0.1} / Time / min(mean_degree, 10) / 10 * log(Time)
#   df(G, H)  = G*H + (8 + 1) * G       (8 covariates + 1 intercept)
# =============================================================================

rm(list=ls())

library(here)
dir <- here()
setwd(dir)

T_train = 70

# ---- 1. Load data (only used to compute N and the average degree) -----------
load("data/weibo.rda")
N = nrow(Ymat)
taus=c(0.1,0.3,0.5,0.7,0.9)


# ---- 2. Build training response and tensor (kept here for completeness; the
#         QIC computation only needs the saved losses) -----------------------
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

# ---- 3. Loop over the (G, H) grid and read each saved loss ------------------
G_vec = rep(1:5,each=5)
H_vec = rep(1:5,times=5)

Time = T_train
N = nrow(Amat)

loss_vec = lapply(1:length(G_vec), function(case){
  G = G_vec[case]
  H = H_vec[case]
  load(paste0("output/real_res/tgnq_",T_train,"-",G,"-", H, ".rda" ))
  loss = mean(res_ge_update$theta_GH$loss)/N/Time
})
loss_vec = unlist(loss_vec)

# ---- 4. Compute QIC ---------------------------------------------------------
lambda=N^0.1/Time/min(sum(Amat)/N,10)/10*log(Time)
QIC = log(loss_vec) + lambda * (G_vec*H_vec+(8+1)*G_vec)



# ---- 5. Plot QIC vs H, faceted by G -----------------------------------------
library(ggplot2)
data = data.frame(QIC = c(QIC), H = H_vec,G = G_vec)

label_G <- function(variable, value) {
  return(paste0("G=", value))
}
names(data) = c("QIC","G","group")

p=ggplot(data, aes(x = factor(G), y = QIC, group = factor(group))) +
  geom_line(aes(color = factor(group)), size = 2) +
  geom_point(aes(color = factor(group)), size = 3) +
  geom_hline(yintercept = min(data$QIC,na.rm = T), linetype = "dashed", color = "#DD4B44", size = 1.5) +
  facet_wrap(~ group, nrow = 1, labeller = label_G) + 
  labs(x = "H", y = "QIC") +
  guides(color = FALSE)+
  scale_color_manual(values = c("#427C9B","#427C9B","#427C9B" ,"#FDC95C" ,"#427C9B")) +
  theme_minimal()+
  theme(text = element_text(size = 27))


ggsave("output/figs/QIC_plot.pdf", plot = p, width = 18, height = 6)

