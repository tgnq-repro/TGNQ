# =============================================================================
# step2_output_tex.R
# -----------------------------------------------------------------------------
# Post-processing for BOTH Section 5.1 (Estimation and inference) and
# Section 5.3 (Misspecified models).
#
# A single Monte Carlo run (see step1_run_simulation.sh / simulation.R) fits the
# oracle, TGNQ-general, TGNQ-special and TGNQ-misspecified estimators in one
# pass, so this one script reproduces all of:
#   * the Section 5.1 tables -> res_spe (oracle/special) and res_ge (general),
#                               with "RMSE (coverage)" cells, and
#   * the Section 5.3 table  -> res_mis (GE vs SPE vs MIS comparison),
#                               with RMSE-only cells.
#
# This version is configured for the design grid in which the available
# G0 = 3 designs are cases {2,3,4,5}.
#
# Table layout (40 rows = 8 designs x 5 quantiles):
#   * rows 1-20  = G0 = 2 designs (cases 1..4),
#   * rows 21-40 = G0 = 3 designs (cases 2..5).
#
# For each aggregated result file in output/res_block/ (or output/res_power/),
# this script:
#   * loads the list of per-seed error vectors `vec_list`,
#   * averages them over seeds,
#   * uses cal_res() from code/utils/assess.R to unpack the averages into
#     the (RMSE, CP, rho) blocks expected by the tables,
#   * builds three formatted data frames:
#       res_spe -- correctly specified (oracle/special-model) results [Sec. 5.1],
#       res_ge  -- TGNQ-general results                               [Sec. 5.1],
#       res_mis -- head-to-head GE vs SPE vs MIS comparison           [Sec. 5.3],
#   * prints LaTeX with kableExtra::kable.
#
# Cells in res_spe / res_ge contain "RMSE (coverage)" pairs; res_mis contains
# RMSE values only.
# =============================================================================

library(kableExtra)
library(dplyr)
library(here)

dir <- here()
setwd(dir)

folder_path = "output/res_block/"
# folder_path = "output/res_power/"

source("code/utils/assess.R")

# All aggregated result files. The names follow
#   res-G_<G0>-c_<case>-a_<additive>.rda
# but historically were prefixed with "TR" in earlier runs; we grep that prefix.
files <- list.files(folder_path)
r_script_files <- files[grep("\\.rda$", files)]

# Aggregated result files follow the naming
#   res-G_<G0>-c_<case>-a_<additive>.rda
# where <additive> is TRUE / FALSE. The "TR" pattern keeps the additive-DGP
# files (a_TRUE); switch the pattern to "FA" for the multiplicative-DGP files
# (a_FALSE).
file_me = r_script_files[grepl("TR", r_script_files)] %>% sort()

# -----------------------------------------------------------------------------
# Pre-allocate the SPECIAL-model table (res_spe). 40 rows:
#   * rows 1-20  = G0 = 2 designs (cases 1..4),
#   * rows 21-40 = G0 = 3 designs (cases 2..5).
# Each design spans 5 consecutive rows (one per quantile).
# -----------------------------------------------------------------------------
res_spe = data.frame()
res_spe[1:40,] = NA
res_spe[,1:20] = NA
res_spe[c(1,21),1] = c(2,3)

# G0 = 2 (rows 1..20): cases 1..4
res_spe[seq(1,20,5),2] = c(100,100,100,200)
res_spe[seq(1,20,5),3] = c(50,100,200,200)

# G0 = 3 (rows 21..40): cases 2..5
res_spe[seq(21,40,5),2] = c(100,100,200,200)
res_spe[seq(21,40,5),3] = c(100,200,200,400)

res_spe[,4] = rep(seq(0.1,0.9,0.2), 8)

names(res_spe) = c("G","N","T","tau","alpha^o", "beta^o", "nu^o", "gamma^o",
                   "alpha", "beta", "nu", "gamma","rho_g", "rho_h",
                   "alpha^r", "beta^r", "nu^r", "gamma^r","rho_g^r", "rho_h^r")

# -----------------------------------------------------------------------------
# Pre-allocate the GENERAL-model table (res_ge). Same row layout.
# -----------------------------------------------------------------------------
res_ge = data.frame()
res_ge[1:40,] = NA
res_ge[,1:17] = NA
res_ge[c(1,21),1] = c(2,3)

# G0 = 2 (rows 1..20): cases 1..4
res_ge[seq(1,20,5),2] = c(100,100,100,200)
res_ge[seq(1,20,5),3] = c(50,100,200,200)

# G0 = 3 (rows 21..40): cases 2..5
res_ge[seq(21,40,5),2] = c(100,100,200,200)
res_ge[seq(21,40,5),3] = c(100,200,200,400)

res_ge[,4] = rep(seq(0.1,0.9,0.2), 8)

names(res_ge) = c("G","N","T","tau","theta^o", "nu^o", "gamma^o",
                  "theta", "nu", "gamma","rho_g", "rho_h",
                  "theta^r", "nu^r", "gamma^r","rho_g^r", "rho_h^r")

# -----------------------------------------------------------------------------
# Pre-allocate the MIS comparison table (res_mis). Same row layout.
# -----------------------------------------------------------------------------
res_mis = data.frame()
res_mis[1:40,] = NA
res_mis[,1:19] = NA
res_mis[c(1,21),1] = c(2,3)

# G0 = 2 (rows 1..20): cases 1..4
res_mis[seq(1,20,5),2] = c(100,100,100,200)
res_mis[seq(1,20,5),3] = c(50,100,200,200)

# G0 = 3 (rows 21..40): cases 2..5
res_mis[seq(21,40,5),2] = c(100,100,200,200)
res_mis[seq(21,40,5),3] = c(100,200,200,400)

res_mis[,4] = rep(seq(0.1,0.9,0.2), 8)

names(res_mis) = c("G","N","T","tau","theta#ge", "nu#ge", "gamma#ge","rho_g#ge", "rho_h#ge",
                   "theta#spe", "nu#spe", "gamma#spe","rho_g#spe", "rho_h#spe",
                   "theta#mis", "nu#mis", "gamma#mis","rho_g#mis", "rho_h#mis")

# -----------------------------------------------------------------------------
# Main loop: for each aggregated file, fill the next block of 5 rows (taus).
# Files are sorted so that the G0 = 2 designs come first (rows 1..20) and the
# G0 = 3 designs follow (rows 21..40), filled continuously.
# -----------------------------------------------------------------------------
row_rho = 1
row = 1

for (file in file_me) {
  
  load(paste0(dir,"/",folder_path,file))
  
  # ---- (1) SPECIAL-model summary -------------------------------------------
  spe_vectors <- sapply(vec_list, function(x) x[[1]])
  spe_vec = rowMeans(spe_vectors)
  spe_list = cal_res(spe_vec)
  
  # ---- (2) GENERAL-model summary -------------------------------------------
  ge_vectors <- sapply(vec_list, function(x) x[[2]])
  ge_vec = rowMeans(ge_vectors)
  ge_list = cal_res(ge_vec)
  
  # Clustering errors (rho_g, rho_h) for SPE (gtnq + refined) and GE.
  res_spe[row_rho, c(13,14,19,20)] = round(spe_list$rho * 100, 3)
  res_ge[row_rho,  c(11,12,16,17)] = round(ge_list$rho * 100, 3)
  
  # ---- (3) MIS results -----------------------------------------------------
  mis_vectors <- sapply(vec_list, function(x) x[[3]])
  mis_vec = round(rowMeans(mis_vectors) * 100, 1)
  
  mis_vectors2 <- sapply(vec_list, function(x) x[[4]])
  mis_vec2 = round(rowMeans(mis_vectors2) * 100, 1)
  
  res_mis[row_rho, c(8,9)]    = round(ge_list$rho * 100, 3)[1:2]
  res_mis[row_rho, c(13,14)]  = round(mis_vec2[16:17], 3)
  res_mis[row_rho, c(18,19)]  = round(mis_vec[16:17], 3)
  
  row_rho = row_rho + 5
  
  # ---- (4) Build "RMSE (CP)" strings for SPE -------------------------------
  mat_spe = spe_list$res_mat
  mat_spe = round(mat_spe * 100, 1)
  mat_spe = format(mat_spe, nsmall = 1)
  mat_spe = as.data.frame(mat_spe)
  
  df_spe <- transform(mat_spe,
                      V1  = paste0(V1,  " (", V5,  ")"),
                      V2  = paste0(V2,  " (", V6,  ")"),
                      V3  = paste0(V3,  " (", V7,  ")"),
                      V4  = paste0(V4,  " (", V8,  ")"),
                      V9  = paste0(V9,  " (", V13, ")"),
                      V10 = paste0(V10, " (", V14, ")"),
                      V11 = paste0(V11, " (", V15, ")"),
                      V12 = paste0(V12, " (", V16, ")"),
                      V17 = paste0(V17, " (", V21, ")"),
                      V18 = paste0(V18, " (", V22, ")"),
                      V19 = paste0(V19, " (", V23, ")"),
                      V20 = paste0(V20, " (", V24, ")"))
  
  res_spe[row:(row+4), c(5:12,15:18)] = df_spe[, c(1:4,9:12,17:20)]
  
  # ---- (5) Build "RMSE (CP)" strings for GE --------------------------------
  mat_ge = ge_list$res_mat
  mat_ge = round(mat_ge * 100, 1)
  mat_ge = format(mat_ge, nsmall = 1)
  mat_ge = as.data.frame(mat_ge)
  
  df_ge <- transform(mat_ge,
                     V1  = paste0(V1,  " (", V4,  ")"),
                     V2  = paste0(V2,  " (", V5,  ")"),
                     V3  = paste0(V3,  " (", V6,  ")"),
                     V7  = paste0(V7,  " (", V10, ")"),
                     V8  = paste0(V8,  " (", V11, ")"),
                     V9  = paste0(V9,  " (", V12, ")"),
                     V13 = paste0(V13, " (", V16, ")"),
                     V14 = paste0(V14, " (", V17, ")"),
                     V15 = paste0(V15, " (", V18, ")"))
  
  res_ge[row:(row+4), c(5:10,13:15)] = df_ge[, c(1:3,7:9,13:15)]
  
  # ---- (6) Build res_mis rows ----------------------------------------------
  res_mis[row:(row+4), c(5:7)]   = mat_ge[, c(7:9)]
  res_mis[row:(row+4), c(10:12)] = format(mis_vec2[1:15], nsmall = 1)
  res_mis[row:(row+4), c(15:17)] = format(mis_vec[1:15],  nsmall = 1)
  
  row = row + 5
}

# =============================================================================
# Final formatting and LaTeX export
# =============================================================================

# G0 = 2 occupies rows 1..20 (left half), G0 = 3 occupies rows 21..40
# (right half of the wide two-column layout).
right_rows <- 21:40

# ---- General-model tables (Section 5.1) -------------------------------------
res_ge[is.na(res_ge)] = ""
res_ge = res_ge[, -c(1,8:12)]  # keep theta/nu/gamma and the refined block only

# Two-column wide layout
res_ge1 = cbind(res_ge[1:20, ], res_ge[right_rows, ])
latex_table <- kable(res_ge1, format = "latex", booktabs = TRUE) %>%
  kable_styling(full_width = FALSE)
cat(latex_table)

# One-column version
latex_table <- kable(res_ge, format = "latex", booktabs = TRUE) %>%
  kable_styling(full_width = FALSE)
cat(latex_table)

# ---- Special-model tables (Section 5.1) ------------------------------------
res_spe[is.na(res_spe)] = ""
res_spe = res_spe[, -c(1,9:14)]  # keep oracle/special and the refined block only

# Two-column wide layout
res_spe1 = cbind(res_spe[1:20, ], res_spe[right_rows, ])
latex_table <- kable(res_spe1, format = "latex", booktabs = TRUE) %>%
  kable_styling(full_width = FALSE)
cat(latex_table)

# RMSE-only view of res_spe (coverage values in parentheses are stripped out).
remove_parentheses <- function(string) {
  gsub("\\([^\\)]+\\)", "", string)
}

df <- as.data.frame(lapply(res_spe, remove_parentheses))
latex_table <- kable(df, format = "latex", booktabs = TRUE) %>%
  kable_styling(full_width = FALSE)
cat(latex_table)

# ---- GE / SPE / MIS comparison table (Section 5.3) -------------------------
res_mis[is.na(res_mis)] = ""
res_mis1 = cbind(res_mis[1:20, -1], res_mis[right_rows, -1])
latex_table <- kable(res_mis1, format = "latex", booktabs = TRUE) %>%
  kable_styling(full_width = FALSE)
cat(latex_table)