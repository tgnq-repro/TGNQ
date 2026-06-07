# =============================================================================
# step8_pred_plot.R
#
# Aggregate one-step-ahead prediction losses for TGNQ and QGNAR over the test
# window (days 71-77 = Thursday-Wednesday) and plot them by quantile.
#
# Output: output/figs/prediction.pdf
# =============================================================================
rm(list = ls())
library(here)
dir <- here()
setwd(dir)

# ---- 1. TGNQ pinball loss per (tau, T_test) ---------------------------------
all_res = c()

for (T_test in 71:77) {
    G = 4
    H = 3
    load(paste0("output/real_res/pred_train",70,"-test",T_test,"-G",G,"-H",H,".rda"))
    taus=c(0.1,0.3,0.5,0.7,0.9)
    res=lapply(1:5, function(k){
      resid = Y_test  - TGNQ.Y_pred_taus[[k]]
      mean(resid * (taus[k] - (resid < 0)))
    })
    unlist(res)
    all_res=c(all_res,unlist(res))
}
tgnq_pred = matrix(all_res, nrow = 5)# 5 quantiles x 7 days

# ---- 2. QGNAR(3) and QGNAR(4) pinball loss per (tau, T_test) ----------------
narg_all3 = c()
narg_all4 = c()
for (T_test in 71:77) {
  load(paste0("output/real_res//narg_pred_train",70,"-test",T_test,".rda"))
  taus=c(0.1,0.3,0.5,0.7,0.9)
  
  res=lapply(1:5, function(k){
    resid = Y_test  - GNAR3.Y_pred_taus[[k]]
    mean(resid * (taus[k] - (resid < 0)))
  })
  res3=unlist(res)
  
  res=lapply(1:5, function(k){
    resid = Y_test  - GNAR4.Y_pred_taus[[k]]
    mean(resid * (taus[k] - (resid < 0)))
  })
  res4=unlist(res)

  narg_all3=c(narg_all3,res3)
  narg_all4=c(narg_all4,res4)
}
narg_pred3 = matrix(narg_all3, nrow = 5)
narg_pred4 = matrix(narg_all4, nrow = 5)




# ---- 3. Reshape to long format ----------------------------------------------
library(ggplot2)
library(reshape2)

quantiles <- c(0.1, 0.3, 0.5, 0.7, 0.9)
quantile_labels <- paste0("tau==", quantiles)

to_long <- function(mat, model_name) {
  df <- as.data.frame(mat)
  df$Quantile <- quantile_labels
  long <- melt(df, id.vars = "Quantile")
  long$Time <- as.numeric(gsub("[^0-9]", "", long$variable)) + 71  
  long$Matrix <- model_name
  long
}

df1 <- to_long(tgnq_pred, "TGNQ")
df2 <- to_long(narg_pred3, "QGNAR(3)")
df3 <- to_long(narg_pred4, "QGNAR(4)")

all_df <- rbind(df1, df2, df3)
colnames(all_df)[which(names(all_df) == "value")] <- "Value"

all_df$Matrix <- factor(all_df$Matrix, levels = c("TGNQ", "QGNAR(3)", "QGNAR(4)"))
line_colors <- c("TGNQ" = "#DD4B44", "QGNAR(3)" = "#9CD6D8", "QGNAR(4)" = "#427C9B")
line_types <- c("TGNQ" = "solid", "QGNAR(3)" = "dashed", "QGNAR(4)" = "dashed")

library(ggplot2)
library(reshape2)

quantiles <- c(0.1, 0.3, 0.5, 0.7, 0.9)
quantile_labels <- paste0("tau==", quantiles)

to_long <- function(mat, model_name) {
  df <- as.data.frame(mat)
  df$Quantile <- quantile_labels
  long <- melt(df, id.vars = "Quantile")
  long$Time <- as.numeric(gsub("[^0-9]", "", long$variable)) + 71
  long$Matrix <- model_name
  long
}

df1 <- to_long(tgnq_pred, "TGNQ")
df2 <- to_long(narg_pred3, "QGNAR(3)")
df3 <- to_long(narg_pred4, "QGNAR(4)")

all_df <- rbind(df1, df2, df3)
colnames(all_df)[which(names(all_df) == "value")] <- "Value"

all_df$Matrix <- factor(all_df$Matrix, levels = c("TGNQ", "QGNAR(3)", "QGNAR(4)"))

# ---- 4. Plot ----------------------------------------------------------------
line_colors <- c("TGNQ" = "#DD4B44", "QGNAR(3)" = "#9CD6D8", "QGNAR(4)" = "#427C9B")
line_types <- c("TGNQ" = "solid", "QGNAR(3)" = "dashed", "QGNAR(4)" = "dashed")

weekday_labels <- c("Thursday", "Friday", "Saturday", "Sunday", "Monday", "Tuesday", "Wednesday")

p=ggplot(all_df, aes(x = Time, y = Value, color = Matrix, linetype = Matrix)) +
  geom_line(size = 1) +
  facet_wrap(~ Quantile, nrow = 1, labeller = label_parsed, scales = "free_y") +
  scale_color_manual(values = line_colors) +
  scale_linetype_manual(values = line_types) +
  scale_x_continuous(breaks = 72:78, 
                     labels = weekday_labels,
                     expand = c(0.02, 0)) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_blank(),    
    strip.text = element_text(face = "bold", size = 20),
    legend.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 16), 
    panel.border = element_blank()
  ) +
  labs(x = "", y = "Predicted Loss")



ggsave("output/figs/prediction.pdf", plot = p, width = 18, height = 6)
