################################################################################
# RE:THINC: Residual Estimation: Tracing Hints to Infer Novoentropic Causality
################################################################################

library(ranger)
library(NHANES)
library(dplyr)
library(ggplot2)
library(doParallel)
library(foreach)
library(parallel)
library(FITSio)
library(nhanesA)
library(scales)

################################################################################
# TEST 1: Medical Data (NHANES) - FULL RE:THINC PIPELINE
################################################################################

#Loading, cleaning, and scaling data
df_med <- NHANES %>%
  select(BMI, BPSysAve, Age) %>%
  filter(complete.cases(.)) %>%
  mutate(across(everything(), ~ as.numeric(scale(.))))

N_samples <- nrow(df_med)
print(N_samples)
set.seed(1)

# Analytical ROPE (H0-Margin for Effect Size)
pooled_sd <- function(x, y) {
  nx <- length(x); ny <- length(y)
  sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
}

rf_full_A <- ranger(BPSysAve ~ BMI, data = df_med)
rf_full_B <- ranger(BMI ~ BPSysAve, data = df_med)

res_full_A <- df_med$BPSysAve - rf_full_A$predictions
res_full_B <- df_med$BMI - rf_full_B$predictions

NN_full_A <- residuals(lm(res_full_A ~ res_full_B))
NN_full_B <- residuals(lm(res_full_B ~ res_full_A))

rope_margin <- 0.1 * pooled_sd(NN_full_A, NN_full_B)
rope_lower <- -rope_margin
rope_upper <- rope_margin



# RE:THINC: Bivariate Base Causal Direction
n_iter <- 1000
sample_size <- floor(0.8 * nrow(df_med))
n_cores <- max(1, detectCores() - 1)

cat("Start Fractional Resampling Layer 1 (", n_iter, " Iterations on ", n_cores, " cores)...\n", sep="")

boot_diffs <- mclapply(1:n_iter, function(i) {
  df_sub <- df_med[sample(nrow(df_med), size = sample_size, replace = FALSE), ]
  
  rf_A <- ranger(BPSysAve ~ BMI, data = df_sub, num.threads = 1)
  rf_B <- ranger(BMI ~ BPSysAve, data = df_sub, num.threads = 1)
  
  res_A <- df_sub$BPSysAve - rf_A$predictions
  res_B <- df_sub$BMI - rf_B$predictions
  
  df_res <- data.frame(res_A = res_A, res_B = res_B)
  
  NN_A <- res_A - ranger(res_A ~ res_B, data = df_res, num.threads = 1)$predictions
  NN_B <- res_B - ranger(res_B ~ res_A, data = df_res, num.threads = 1)$predictions
  
  mean(abs(NN_A)) - mean(abs(NN_B))
}, mc.cores = n_cores) |> unlist()

boot_lower <- quantile(boot_diffs, 0.025)
boot_upper <- quantile(boot_diffs, 0.975)

#Decision
cat("\n======================================================\n")
cat("95% Analytical ROPE (H0): [", round(rope_lower, 5), ",", round(rope_upper, 5), "]\n")
cat("95% Bootstrap Signal:     [", round(boot_lower, 5), ",", round(boot_upper, 5), "]\n\n")

if(boot_upper < 0 | boot_lower > 0) {
  if(boot_upper < rope_lower | boot_lower > rope_upper) {
    if(boot_upper < 0) print("Result: BMI -> Blood Pressure (MEANINGFUL True Causal Direction!)")
    if(boot_lower > 0) print("Result: Blood Pressure -> BMI (MEANINGFUL True Causal Direction!)")
  } else {
    if(boot_upper < 0) print("Result: BMI -> Blood Pressure (Direction identified, but NOT practically meaningful / overlaps ROPE)")
    if(boot_lower > 0) print("Result: Blood Pressure -> BMI (Direction identified, but NOT practically meaningful / overlaps ROPE)")
  }
} else {
  print("Result: Undecided (Confidence interval contains 0 / Spurious Correlation)")
}
cat("======================================================\n\n")


# Layer 2 & 3: Confounder entropy reduction and causal cyclic analysis

n_boot <- 1000 
cat("    Generating topological distributions Layer 2 & 3 (", n_boot, " Iterations on ", n_cores, " cores)...\n", sep="")

results_list <- mclapply(1:n_boot, function(i) {
  idx <- sample(N_samples, size = sample_size, replace = FALSE)
  df_sub <- df_med[idx, ]
  
  #Base Extraction
  rf_A_t <- ranger(BPSysAve ~ BMI, data = df_sub, num.threads = 1)
  rf_B_t <- ranger(BMI ~ BPSysAve, data = df_sub, num.threads = 1)
  
  res_A_t <- df_sub$BPSysAve - rf_A_t$predictions
  res_B_t <- df_sub$BMI - rf_B_t$predictions
  
  df_res_t <- data.frame(res_A = res_A_t, res_B = res_B_t)
  
  #Projection (NN-Isolation)
  NN_A_t <- res_A_t - ranger(res_A ~ res_B, data = df_res_t, num.threads = 1)$predictions
  NN_B_t <- res_B_t - ranger(res_B ~ res_A, data = df_res_t, num.threads = 1)$predictions
  
  df_model <- cbind(df_res_t, df_sub)
  
  df_model_perm <- df_model
  df_model_perm$res_A <- sample(df_model_perm$res_A)
  df_model_perm$res_B <- sample(df_model_perm$res_B)
  
  NN_A_perm <- res_A_t - ranger(res_A ~ res_B, data = df_model_perm, num.threads = 1)$predictions
  NN_B_perm <- res_B_t - ranger(res_B ~ res_A, data = df_model_perm, num.threads = 1)$predictions
  
  # Topological test
  NNN_A_true <- NN_A_t - ranger(NN_A_t ~ BPSysAve, data = df_sub, num.threads = 1)$predictions
  NNN_B_true <- NN_B_t - ranger(NN_B_t ~ BMI, data = df_sub, num.threads = 1)$predictions
  
  # Permutation Null
  df_model_null <- df_sub
  df_model_null$BPSysAve <- sample(df_model_null$BPSysAve)
  df_model_null$BMI <- sample(df_model_null$BMI)
  
  NNN_A_null <- NN_A_t - ranger(NN_A_t ~ BPSysAve, data = df_model_null, num.threads = 1)$predictions
  NNN_B_null <- NN_B_t - ranger(NN_B_t ~ BMI, data = df_model_null, num.threads = 1)$predictions
  
  return(c(
    NN_perm  = mean(abs(NN_A_perm)) - mean(abs(NN_B_perm)),
    MAD_L1_A = mean(abs(NN_A_perm)),
    MAD_L2_A = mean(abs(NN_A_t)),
    MAD_L1_B = mean(abs(NN_B_perm)),
    MAD_L2_B = mean(abs(NN_B_t)),
    NNN_true = mean(abs(NNN_A_true)) - mean(abs(NNN_B_true)),
    NNN_null = mean(abs(NNN_A_null)) - mean(abs(NNN_B_null))
  ))
}, mc.cores = n_cores)

#Integrating results into a vector
results_matrix <- do.call(rbind, results_list)

NN_perm  <- results_matrix[, "NN_perm"]
MAD_L1_A <- results_matrix[, "MAD_L1_A"]
MAD_L2_A <- results_matrix[, "MAD_L2_A"]
MAD_L1_B <- results_matrix[, "MAD_L1_B"]
MAD_L2_B <- results_matrix[, "MAD_L2_B"]
NNN_true <- results_matrix[, "NNN_true"]
NNN_null <- results_matrix[, "NNN_null"]

cat("\n... done!\n")

#Confounder Entropy Reduction (CER)
drop_A_vec <- (1 - (MAD_L2_A / MAD_L1_A)) * 100
drop_B_vec <- (1 - (MAD_L2_B / MAD_L1_B)) * 100

ci_drop_A <- quantile(drop_A_vec, probs = c(0.025, 0.5, 0.975))
ci_drop_B <- quantile(drop_B_vec, probs = c(0.025, 0.5, 0.975))

cat("\n    [CER-Check] Confounder-Proportion in BPSysAve: ", round(ci_drop_A[2], 2), 
    "%  (95% CI: [", round(ci_drop_A[1], 2), " , ", round(ci_drop_A[3], 2), "])")

cat("\n    [CER-Check] Confounder-Proportion in BMI:      ", round(ci_drop_B[2], 2), 
    "%  (95% CI: [", round(ci_drop_B[1], 2), " , ", round(ci_drop_B[3], 2), "])\n")

#Statistical Inference Criterion for a) CER and b) cyclic causality
n_resamples <- 10000
relative_ratio_vec <- numeric(n_resamples)

for (i in 1:n_resamples) {
  idx <- sample(1:n_boot, size = n_boot, replace = TRUE)
  mad_true_total <- mean(abs(NNN_true[idx])) 
  mad_null_total <- mean(abs(NNN_null[idx]))
  relative_ratio_vec[i] <- mad_true_total / mad_null_total
}

ci_relative <- quantile(relative_ratio_vec, probs = c(0.025, 0.975))

cat("\n    95% CI of the relative proportion (True / Null): [", round(ci_relative[1], 4), ",", round(ci_relative[2], 4), "]\n\n")

if (ci_relative[2] < 1) {
  cat("Layer 3 conclusion: vicious cycle / bidirectional feedback loop\n")
} else if (ci_relative[1] > 1) {
  cat("Layer 3 conclusion: unidirectional + adversarial noise / heteroscedasticity\n")
} else {
  cat("Layer 3 conclusion: unidirection + symmetric noise\n")
}

# Paper plot
plot_data <- data.frame(Value = boot_diffs)

p_med <- ggplot(plot_data, aes(x = Value)) +
  annotate("rect", xmin = rope_lower, xmax = rope_upper, ymin = 0, ymax = Inf, 
           fill = "darkgray", alpha = 0.3) +
  geom_density(fill = "#3498db", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
  geom_vline(xintercept = c(rope_lower, rope_upper), linetype = "dashed", color = "#2c3e50") +
  theme_minimal() +
  labs(title = "RE:THINC: Causal Direction Validation (NHANES)", 
       subtitle = "Blue: Subsampled Causal Signal | Gray Zone: Analytical ROPE (0.1 SD)",
       x = "Difference in Noise Distortion (Delta NN)", y = "Density")

print(p_med)

################################################################################
# TEST 2: Environmental Chemistry (Temperature vs. Ozone)
################################################################################

# Loading, cleaning, scaling, and sampling data
set.seed(1)
df_chem <- airquality %>% 
  select(Temp, Ozone) %>% 
  filter(complete.cases(.)) %>% 
  mutate(across(everything(), ~ as.numeric(scale(.))))

N_samples_chem <- nrow(df_chem)
cat("\nStarting TEST 2: Environmental Chemistry (N =", N_samples_chem, ")\n")

# Analytical ROPE (H0-Margin for Effect Size)
rf_full_A_chem <- ranger(Ozone ~ Temp, data = df_chem)
rf_full_B_chem <- ranger(Temp ~ Ozone, data = df_chem)

res_full_A_chem <- df_chem$Ozone - rf_full_A_chem$predictions
res_full_B_chem <- df_chem$Temp - rf_full_B_chem$predictions

NN_full_A_chem <- residuals(lm(res_full_A_chem ~ res_full_B_chem))
NN_full_B_chem <- residuals(lm(res_full_B_chem ~ res_full_A_chem))

rope_margin_chem <- 0.1 * pooled_sd(NN_full_A_chem, NN_full_B_chem)
rope_lower_chem <- -rope_margin_chem
rope_upper_chem <- rope_margin_chem

# RE:THINC L1: Bivariate Base Causal Direction
cat("Start Fractional Resampling Layer 1 (", n_iter, " Iterations on ", n_cores, " cores)...\n", sep="")

boot_diffs_chem <- mclapply(1:n_iter, function(i) {
  df_sub <- df_chem[sample(nrow(df_chem), size = floor(0.8 * nrow(df_chem)), replace = FALSE), ]
  
  rf_A <- ranger(Ozone ~ Temp, data = df_sub, num.threads = 1)
  rf_B <- ranger(Temp ~ Ozone, data = df_sub, num.threads = 1)
  
  res_A <- df_sub$Ozone - rf_A$predictions
  res_B <- df_sub$Temp - rf_B$predictions
  
  df_res <- data.frame(res_A = res_A, res_B = res_B)
  
  NN_A <- res_A - ranger(res_A ~ res_B, data = df_res, num.threads = 1)$predictions
  NN_B <- res_B - ranger(res_B ~ res_A, data = df_res, num.threads = 1)$predictions
  
  mean(abs(NN_A)) - mean(abs(NN_B))
}, mc.cores = n_cores) |> unlist()

boot_lower_chem <- quantile(boot_diffs_chem, 0.025)
boot_upper_chem <- quantile(boot_diffs_chem, 0.975)

# Decision L1
cat("\n======================================================\n")
cat("95% Analytical ROPE (H0): [", round(rope_lower_chem, 5), ",", round(rope_upper_chem, 5), "]\n")
cat("95% Bootstrap Signal:     [", round(boot_lower_chem, 5), ",", round(boot_upper_chem, 5), "]\n\n")

if(boot_upper_chem < 0 | boot_lower_chem > 0) {
  if(boot_upper_chem < rope_lower_chem | boot_lower_chem > rope_upper_chem) {
    if(boot_upper_chem < 0) print("Result: Temp -> Ozone (MEANINGFUL True Causal Direction!)")
    if(boot_lower_chem > 0) print("Result: Ozone -> Temp (MEANINGFUL True Causal Direction!)")
  } else {
    if(boot_upper_chem < 0) print("Result: Temp -> Ozone (Direction identified, but NOT practically meaningful / overlaps ROPE)")
    if(boot_lower_chem > 0) print("Result: Ozone -> Temp (Direction identified, but NOT practically meaningful / overlaps ROPE)")
  }
} else {
  print("Result: Undecided (Confidence interval contains 0 / Spurious Correlation)")
}
cat("======================================================\n\n")

# Layer 2 & 3
cat("    Generating topological distributions Layer 2 & 3 (", n_boot, " Iterations on ", n_cores, " cores)...\n", sep="")

results_list_chem <- mclapply(1:n_boot, function(i) {
  df_sub <- df_chem[sample(nrow(df_chem), size = floor(0.8 * nrow(df_chem)), replace = FALSE), ]
  
  rf_A_t <- ranger(Ozone ~ Temp, data = df_sub, num.threads = 1)
  rf_B_t <- ranger(Temp ~ Ozone, data = df_sub, num.threads = 1)
  
  res_A_t <- df_sub$Ozone - rf_A_t$predictions
  res_B_t <- df_sub$Temp - rf_B_t$predictions
  
  df_res_t <- data.frame(res_A = res_A_t, res_B = res_B_t)
  
  NN_A_t <- res_A_t - ranger(res_A ~ res_B, data = df_res_t, num.threads = 1)$predictions
  NN_B_t <- res_B_t - ranger(res_B ~ res_A, data = df_res_t, num.threads = 1)$predictions
  
  df_model_perm <- cbind(df_res_t, df_sub)
  df_model_perm$res_A <- sample(df_model_perm$res_A)
  df_model_perm$res_B <- sample(df_model_perm$res_B)
  
  NN_A_perm <- res_A_t - ranger(res_A ~ res_B, data = df_model_perm, num.threads = 1)$predictions
  NN_B_perm <- res_B_t - ranger(res_B ~ res_A, data = df_model_perm, num.threads = 1)$predictions
  
  NNN_A_true <- NN_A_t - ranger(NN_A_t ~ Ozone, data = df_sub, num.threads = 1)$predictions
  NNN_B_true <- NN_B_t - ranger(NN_B_t ~ Temp, data = df_sub, num.threads = 1)$predictions
  
  df_model_null <- df_sub
  df_model_null$Ozone <- sample(df_model_null$Ozone)
  df_model_null$Temp <- sample(df_model_null$Temp)
  
  NNN_A_null <- NN_A_t - ranger(NN_A_t ~ Ozone, data = df_model_null, num.threads = 1)$predictions
  NNN_B_null <- NN_B_t - ranger(NN_B_t ~ Temp, data = df_model_null, num.threads = 1)$predictions
  
  return(c(
    MAD_L1_A = mean(abs(NN_A_perm)), MAD_L2_A = mean(abs(NN_A_t)),
    MAD_L1_B = mean(abs(NN_B_perm)), MAD_L2_B = mean(abs(NN_B_t)),
    NNN_true = mean(abs(NNN_A_true)) - mean(abs(NNN_B_true)),
    NNN_null = mean(abs(NNN_A_null)) - mean(abs(NNN_B_null))
  ))
}, mc.cores = n_cores)

res_mat_chem <- do.call(rbind, results_list_chem)
cat("\n... done!\n")

ci_drop_A_chem <- quantile((1 - (res_mat_chem[, "MAD_L2_A"] / res_mat_chem[, "MAD_L1_A"])) * 100, probs = c(0.025, 0.5, 0.975))
ci_drop_B_chem <- quantile((1 - (res_mat_chem[, "MAD_L2_B"] / res_mat_chem[, "MAD_L1_B"])) * 100, probs = c(0.025, 0.5, 0.975))

cat("\n    [CER-Check] Confounder-Proportion in Ozone: ", round(ci_drop_A_chem[2], 2), "%  (95% CI: [", round(ci_drop_A_chem[1], 2), " , ", round(ci_drop_A_chem[3], 2), "])\n")
cat("    [CER-Check] Confounder-Proportion in Temp: ", round(ci_drop_B_chem[2], 2), "%  (95% CI: [", round(ci_drop_B_chem[1], 2), " , ", round(ci_drop_B_chem[3], 2), "])\n")

ratio_vec_chem <- numeric(10000)
for (i in 1:10000) {
  idx <- sample(1:n_boot, replace = TRUE)
  ratio_vec_chem[i] <- mean(abs(res_mat_chem[idx, "NNN_true"])) / mean(abs(res_mat_chem[idx, "NNN_null"]))
}
ci_rel_chem <- quantile(ratio_vec_chem, probs = c(0.025, 0.975))

cat("\n    95% CI of the relative proportion (True / Null): [", round(ci_rel_chem[1], 4), ",", round(ci_rel_chem[2], 4), "]\n")
if (ci_rel_chem[2] < 1) cat("Layer 3 conclusion: vicious cycle / bidirectional feedback loop\n") else if (ci_rel_chem[1] > 1) cat("Layer 3 conclusion: unidirectional + adversarial noise / heteroscedasticity\n") else cat("Layer 3 conclusion: unidirection + symmetric noise\n")

# Paper plot
plot_data_chem <- data.frame(Value = boot_diffs_chem)
p_chem <- ggplot(plot_data_chem, aes(x = Value)) +
  annotate("rect", xmin = rope_lower_chem, xmax = rope_upper_chem, ymin = 0, ymax = Inf, fill = "darkgray", alpha = 0.3) +
  geom_density(fill = "#27ae60", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
  geom_vline(xintercept = c(rope_lower_chem, rope_upper_chem), linetype = "dashed", color = "#2c3e50") +
  theme_minimal() +
  labs(title = "Environmental Chemistry: Temperature vs. Ozone", subtitle = "Green: Subsampled Causal Signal | Gray Zone: Analytical ROPE (0.1 SD)", x = "Difference in Noise Distortion (Delta NN)", y = "Density")
print(p_chem)


################################################################################
# TEST 3: Geophysics (Eruptions vs. Waiting time)
################################################################################

# Loading, cleaning, scaling, and sampling data
set.seed(1)
df_geo <- faithful %>% 
  select(eruptions, waiting) %>% 
  filter(complete.cases(.)) %>% 
  mutate(across(everything(), ~ as.numeric(scale(.))))

N_samples_geo <- nrow(df_geo)
cat("\nStarting TEST 3: Geophysics (N =", N_samples_geo, ")\n")

# Analytical ROPE (H0-Margin for Effect Size)
rf_full_A_geo <- ranger(waiting ~ eruptions, data = df_geo)
rf_full_B_geo <- ranger(eruptions ~ waiting, data = df_geo)

res_full_A_geo <- df_geo$waiting - rf_full_A_geo$predictions
res_full_B_geo <- df_geo$eruptions - rf_full_B_geo$predictions

NN_full_A_geo <- residuals(lm(res_full_A_geo ~ res_full_B_geo))
NN_full_B_geo <- residuals(lm(res_full_B_geo ~ res_full_A_geo))

rope_margin_geo <- 0.1 * pooled_sd(NN_full_A_geo, NN_full_B_geo)
rope_lower_geo <- -rope_margin_geo
rope_upper_geo <- rope_margin_geo

# RE:THINC L1: Bivariate Base Causal Direction
cat("Start Fractional Resampling Layer 1 (", n_iter, " Iterations on ", n_cores, " cores)...\n", sep="")

boot_diffs_geo <- mclapply(1:n_iter, function(i) {
  df_sub <- df_geo[sample(nrow(df_geo), size = floor(0.8 * nrow(df_geo)), replace = FALSE), ]
  
  rf_A <- ranger(waiting ~ eruptions, data = df_sub, num.threads = 1)
  rf_B <- ranger(eruptions ~ waiting, data = df_sub, num.threads = 1)
  
  res_A <- df_sub$waiting - rf_A$predictions
  res_B <- df_sub$eruptions - rf_B$predictions
  
  df_res <- data.frame(res_A = res_A, res_B = res_B)
  
  NN_A <- res_A - ranger(res_A ~ res_B, data = df_res, num.threads = 1)$predictions
  NN_B <- res_B - ranger(res_B ~ res_A, data = df_res, num.threads = 1)$predictions
  
  mean(abs(NN_A)) - mean(abs(NN_B))
}, mc.cores = n_cores) |> unlist()

boot_lower_geo <- quantile(boot_diffs_geo, 0.025)
boot_upper_geo <- quantile(boot_diffs_geo, 0.975)

# Decision L1
cat("\n======================================================\n")
cat("95% Analytical ROPE (H0): [", round(rope_lower_geo, 5), ",", round(rope_upper_geo, 5), "]\n")
cat("95% Bootstrap Signal:     [", round(boot_lower_geo, 5), ",", round(boot_upper_geo, 5), "]\n\n")

if(boot_upper_geo < 0 | boot_lower_geo > 0) {
  if(boot_upper_geo < rope_lower_geo | boot_lower_geo > rope_upper_geo) {
    if(boot_upper_geo < 0) print("Result: Eruptions -> Waiting Time (MEANINGFUL True Causal Direction!)")
    if(boot_lower_geo > 0) print("Result: Waiting Time -> Eruptions (MEANINGFUL True Causal Direction!)")
  } else {
    if(boot_upper_geo < 0) print("Result: Eruptions -> Waiting Time (Direction identified, but NOT practically meaningful / overlaps ROPE)")
    if(boot_lower_geo > 0) print("Result: Waiting Time -> Eruptions (Direction identified, but NOT practically meaningful / overlaps ROPE)")
  }
} else {
  print("Result: Undecided (Confidence interval contains 0 / Spurious Correlation)")
}
cat("======================================================\n\n")

# Layer 2 & 3
cat("    Generating topological distributions Layer 2 & 3 (", n_boot, " Iterations on ", n_cores, " cores)...\n", sep="")

results_list_geo <- mclapply(1:n_boot, function(i) {
  df_sub <- df_geo[sample(nrow(df_geo), size = floor(0.8 * nrow(df_geo)), replace = FALSE), ]
  
  rf_A_t <- ranger(waiting ~ eruptions, data = df_sub, num.threads = 1)
  rf_B_t <- ranger(eruptions ~ waiting, data = df_sub, num.threads = 1)
  
  res_A_t <- df_sub$waiting - rf_A_t$predictions
  res_B_t <- df_sub$eruptions - rf_B_t$predictions
  
  df_res_t <- data.frame(res_A = res_A_t, res_B = res_B_t)
  
  NN_A_t <- res_A_t - ranger(res_A ~ res_B, data = df_res_t, num.threads = 1)$predictions
  NN_B_t <- res_B_t - ranger(res_B ~ res_A, data = df_res_t, num.threads = 1)$predictions
  
  df_model_perm <- cbind(df_res_t, df_sub)
  df_model_perm$res_A <- sample(df_model_perm$res_A)
  df_model_perm$res_B <- sample(df_model_perm$res_B)
  
  NN_A_perm <- res_A_t - ranger(res_A ~ res_B, data = df_model_perm, num.threads = 1)$predictions
  NN_B_perm <- res_B_t - ranger(res_B ~ res_A, data = df_model_perm, num.threads = 1)$predictions
  
  NNN_A_true <- NN_A_t - ranger(NN_A_t ~ waiting, data = df_sub, num.threads = 1)$predictions
  NNN_B_true <- NN_B_t - ranger(NN_B_t ~ eruptions, data = df_sub, num.threads = 1)$predictions
  
  df_model_null <- df_sub
  df_model_null$waiting <- sample(df_model_null$waiting)
  df_model_null$eruptions <- sample(df_model_null$eruptions)
  
  NNN_A_null <- NN_A_t - ranger(NN_A_t ~ waiting, data = df_model_null, num.threads = 1)$predictions
  NNN_B_null <- NN_B_t - ranger(NN_B_t ~ eruptions, data = df_model_null, num.threads = 1)$predictions
  
  return(c(
    MAD_L1_A = mean(abs(NN_A_perm)), MAD_L2_A = mean(abs(NN_A_t)),
    MAD_L1_B = mean(abs(NN_B_perm)), MAD_L2_B = mean(abs(NN_B_t)),
    NNN_true = mean(abs(NNN_A_true)) - mean(abs(NNN_B_true)),
    NNN_null = mean(abs(NNN_A_null)) - mean(abs(NNN_B_null))
  ))
}, mc.cores = n_cores)

res_mat_geo <- do.call(rbind, results_list_geo)
cat("\n... done!\n")

ci_drop_A_geo <- quantile((1 - (res_mat_geo[, "MAD_L2_A"] / res_mat_geo[, "MAD_L1_A"])) * 100, probs = c(0.025, 0.5, 0.975))
ci_drop_B_geo <- quantile((1 - (res_mat_geo[, "MAD_L2_B"] / res_mat_geo[, "MAD_L1_B"])) * 100, probs = c(0.025, 0.5, 0.975))

cat("\n    [CER-Check] Confounder-Proportion in Waiting Time: ", round(ci_drop_A_geo[2], 2), "%  (95% CI: [", round(ci_drop_A_geo[1], 2), " , ", round(ci_drop_A_geo[3], 2), "])\n")
cat("    [CER-Check] Confounder-Proportion in Eruptions: ", round(ci_drop_B_geo[2], 2), "%  (95% CI: [", round(ci_drop_B_geo[1], 2), " , ", round(ci_drop_B_geo[3], 2), "])\n")

ratio_vec_geo <- numeric(10000)
for (i in 1:10000) {
  idx <- sample(1:n_boot, replace = TRUE)
  ratio_vec_geo[i] <- mean(abs(res_mat_geo[idx, "NNN_true"])) / mean(abs(res_mat_geo[idx, "NNN_null"]))
}
ci_rel_geo <- quantile(ratio_vec_geo, probs = c(0.025, 0.975))

cat("\n    95% CI of the relative proportion (True / Null): [", round(ci_rel_geo[1], 4), ",", round(ci_rel_geo[2], 4), "]\n")
if (ci_rel_geo[2] < 1) cat("Layer 3 conclusion: vicious cycle / bidirectional feedback loop\n") else if (ci_rel_geo[1] > 1) cat("Layer 3 conclusion: unidirectional + adversarial noise / heteroscedasticity\n") else cat("Layer 3 conclusion: unidirection + symmetric noise\n")

# Paper plot
plot_data_geo <- data.frame(Value = boot_diffs_geo)
p_geo <- ggplot(plot_data_geo, aes(x = Value)) +
  annotate("rect", xmin = rope_lower_geo, xmax = rope_upper_geo, ymin = 0, ymax = Inf, fill = "darkgray", alpha = 0.3) +
  geom_density(fill = "#e67e22", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
  geom_vline(xintercept = c(rope_lower_geo, rope_upper_geo), linetype = "dashed", color = "#2c3e50") +
  theme_minimal() +
  labs(title = "Geophysics: Eruption Duration vs. Waiting Time", subtitle = "Orange: Subsampled Causal Signal | Gray Zone: Analytical ROPE (0.1 SD)", x = "Difference in Noise Distortion (Delta NN)", y = "Density")
print(p_geo)


################################################################################
# TEST 4: Socioeconomics (Carat vs. Price)
################################################################################

#Loading, cleaning, scaling, and sampling data
set.seed(1)
df_econ <- diamonds %>% 
  select(carat, price) %>% 
  filter(complete.cases(.)) %>%
  mutate(across(everything(), ~ as.numeric(scale(.))))

N_samples_econ <- nrow(df_econ)
cat("\nStarting TEST 4: Socioeconomics (N =", N_samples_econ, ")\n")

# Analytical ROPE (H0-Margin for Effect Size)
rf_full_A_econ <- ranger(price ~ carat, data = df_econ)
rf_full_B_econ <- ranger(carat ~ price, data = df_econ)

res_full_A_econ <- df_econ$price - rf_full_A_econ$predictions
res_full_B_econ <- df_econ$carat - rf_full_B_econ$predictions

NN_full_A_econ <- residuals(lm(res_full_A_econ ~ res_full_B_econ))
NN_full_B_econ <- residuals(lm(res_full_B_econ ~ res_full_A_econ))

rope_margin_econ <- 0.1 * pooled_sd(NN_full_A_econ, NN_full_B_econ)
rope_lower_econ <- -rope_margin_econ
rope_upper_econ <- rope_margin_econ

# RE:THINC L1: Bivariate Base Causal Direction
cat("Start Fractional Resampling Layer 1 (", n_iter, " Iterations on ", n_cores, " cores)...\n", sep="")

boot_diffs_econ <- mclapply(1:n_iter, function(i) {
  df_sub <- df_econ[sample(nrow(df_econ), size = floor(0.8 * nrow(df_econ)), replace = FALSE), ]
  
  rf_A <- ranger(price ~ carat, data = df_sub, num.threads = 1)
  rf_B <- ranger(carat ~ price, data = df_sub, num.threads = 1)
  
  res_A <- df_sub$price - rf_A$predictions
  res_B <- df_sub$carat - rf_B$predictions
  
  df_res <- data.frame(res_A = res_A, res_B = res_B)
  
  NN_A <- res_A - ranger(res_A ~ res_B, data = df_res, num.threads = 1)$predictions
  NN_B <- res_B - ranger(res_B ~ res_A, data = df_res, num.threads = 1)$predictions
  
  mean(abs(NN_A)) - mean(abs(NN_B))
}, mc.cores = n_cores) |> unlist()

boot_lower_econ <- quantile(boot_diffs_econ, 0.025)
boot_upper_econ <- quantile(boot_diffs_econ, 0.975)

#Decision L1
cat("\n======================================================\n")
cat("95% Analytical ROPE (H0): [", round(rope_lower_econ, 5), ",", round(rope_upper_econ, 5), "]\n")
cat("95% Bootstrap Signal:     [", round(boot_lower_econ, 5), ",", round(boot_upper_econ, 5), "]\n\n")

if(boot_upper_econ < 0 | boot_lower_econ > 0) {
  if(boot_upper_econ < rope_lower_econ | boot_lower_econ > rope_upper_econ) {
    if(boot_upper_econ < 0) print("Result: Carat -> Price (MEANINGFUL True Causal Direction!)")
    if(boot_lower_econ > 0) print("Result: Price -> Carat (MEANINGFUL True Causal Direction!)")
  } else {
    if(boot_upper_econ < 0) print("Result: Carat -> Price (Direction identified, but NOT practically meaningful / overlaps ROPE)")
    if(boot_lower_econ > 0) print("Result: Price -> Carat (Direction identified, but NOT practically meaningful / overlaps ROPE)")
  }
} else {
  print("Result: Undecided (Confidence interval contains 0 / Spurious Correlation)")
}
cat("======================================================\n\n")

# Layer 2 & 3
cat("    Generating topological distributions Layer 2 & 3 (", n_boot, " Iterations on ", n_cores, " cores)...\n", sep="")

results_list_econ <- mclapply(1:n_boot, function(i) {
  df_sub <- df_econ[sample(nrow(df_econ), size = floor(0.8 * nrow(df_econ)), replace = FALSE), ]
  
  rf_A_t <- ranger(price ~ carat, data = df_sub, num.threads = 1)
  rf_B_t <- ranger(carat ~ price, data = df_sub, num.threads = 1)
  
  res_A_t <- df_sub$price - rf_A_t$predictions
  res_B_t <- df_sub$carat - rf_B_t$predictions
  
  df_res_t <- data.frame(res_A = res_A_t, res_B = res_B_t)
  
  NN_A_t <- res_A_t - ranger(res_A ~ res_B, data = df_res_t, num.threads = 1)$predictions
  NN_B_t <- res_B_t - ranger(res_B ~ res_A, data = df_res_t, num.threads = 1)$predictions
  
  df_model_perm <- cbind(df_res_t, df_sub)
  df_model_perm$res_A <- sample(df_model_perm$res_A)
  df_model_perm$res_B <- sample(df_model_perm$res_B)
  
  NN_A_perm <- res_A_t - ranger(res_A ~ res_B, data = df_model_perm, num.threads = 1)$predictions
  NN_B_perm <- res_B_t - ranger(res_B ~ res_A, data = df_model_perm, num.threads = 1)$predictions
  
  NNN_A_true <- NN_A_t - ranger(NN_A_t ~ price, data = df_sub, num.threads = 1)$predictions
  NNN_B_true <- NN_B_t - ranger(NN_B_t ~ carat, data = df_sub, num.threads = 1)$predictions
  
  df_model_null <- df_sub
  df_model_null$price <- sample(df_model_null$price)
  df_model_null$carat <- sample(df_model_null$carat)
  
  NNN_A_null <- NN_A_t - ranger(NN_A_t ~ price, data = df_model_null, num.threads = 1)$predictions
  NNN_B_null <- NN_B_t - ranger(NN_B_t ~ carat, data = df_model_null, num.threads = 1)$predictions
  
  return(c(
    MAD_L1_A = mean(abs(NN_A_perm)), MAD_L2_A = mean(abs(NN_A_t)),
    MAD_L1_B = mean(abs(NN_B_perm)), MAD_L2_B = mean(abs(NN_B_t)),
    NNN_true = mean(abs(NNN_A_true)) - mean(abs(NNN_B_true)),
    NNN_null = mean(abs(NNN_A_null)) - mean(abs(NNN_B_null))
  ))
}, mc.cores = n_cores)

res_mat_econ <- do.call(rbind, results_list_econ)
cat("\n... done!\n")

ci_drop_A_econ <- quantile((1 - (res_mat_econ[, "MAD_L2_A"] / res_mat_econ[, "MAD_L1_A"])) * 100, probs = c(0.025, 0.5, 0.975))
ci_drop_B_econ <- quantile((1 - (res_mat_econ[, "MAD_L2_B"] / res_mat_econ[, "MAD_L1_B"])) * 100, probs = c(0.025, 0.5, 0.975))

cat("\n    [CER-Check] Confounder-Proportion in Price: ", round(ci_drop_A_econ[2], 2), "%  (95% CI: [", round(ci_drop_A_econ[1], 2), " , ", round(ci_drop_A_econ[3], 2), "])\n")
cat("    [CER-Check] Confounder-Proportion in Carat: ", round(ci_drop_B_econ[2], 2), "%  (95% CI: [", round(ci_drop_B_econ[1], 2), " , ", round(ci_drop_B_econ[3], 2), "])\n")

ratio_vec_econ <- numeric(10000)
for (i in 1:10000) {
  idx <- sample(1:n_boot, replace = TRUE)
  ratio_vec_econ[i] <- mean(abs(res_mat_econ[idx, "NNN_true"])) / mean(abs(res_mat_econ[idx, "NNN_null"]))
}
ci_rel_econ <- quantile(ratio_vec_econ, probs = c(0.025, 0.975))

cat("\n    95% CI of the relative proportion (True / Null): [", round(ci_rel_econ[1], 4), ",", round(ci_rel_econ[2], 4), "]\n")
if (ci_rel_econ[2] < 1) cat("Layer 3 conclusion: vicious cycle / bidirectional feedback loop\n") else if (ci_rel_econ[1] > 1) cat("Layer 3 conclusion: unidirectional + adversarial noise / heteroscedasticity\n") else cat("Layer 3 conclusion: unidirection + symmetric noise\n")

# Paper plot
plot_data_econ <- data.frame(Value = boot_diffs_econ)
p_econ <- ggplot(plot_data_econ, aes(x = Value)) +
  annotate("rect", xmin = rope_lower_econ, xmax = rope_upper_econ, ymin = 0, ymax = Inf, fill = "darkgray", alpha = 0.3) +
  geom_density(fill = "#2980b9", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
  geom_vline(xintercept = c(rope_lower_econ, rope_upper_econ), linetype = "dashed", color = "#2c3e50") +
  theme_minimal() +
  labs(title = "Socioeconomics: Physical Weight vs. Market Value", subtitle = "Blue: Subsampled Causal Signal | Gray Zone: Analytical ROPE (0.1 SD)", x = "Difference in Noise Distortion (Delta NN)", y = "Density")
print(p_econ)


################################################################################
# TEST 5: Physical Chemistry (Alcohol Concentration vs. Density | Control: Sugar)
################################################################################

url_wine <- "https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-white.csv"
df_chem <- read.csv(url_wine, sep=";")
df_chem <- df_chem %>%
  dplyr::select("alcohol", "density", "residual.sugar") %>% 
  dplyr::rename(sugar = residual.sugar) %>% 
  filter(complete.cases(.)) %>%
  mutate(across(everything(), ~ as.numeric(scale(.))))

set.seed(1)
N_samples_chem <- nrow(df_chem)
cat("\nStarting TEST 5: Physical Chemistry (Alcohol vs Density, N =", N_samples_chem, ")\n")

# Analytical ROPE (H0-Margin for Effect Size) - WITH SUGAR CONTROL
rf_full_A_chem <- ranger(density ~ alcohol + sugar, data = df_chem)
rf_full_B_chem <- ranger(alcohol ~ density + sugar, data = df_chem)

res_full_A_chem <- df_chem$density - rf_full_A_chem$predictions
res_full_B_chem <- df_chem$alcohol - rf_full_B_chem$predictions

NN_full_A_chem <- residuals(lm(res_full_A_chem ~ res_full_B_chem + df_chem$sugar))
NN_full_B_chem <- residuals(lm(res_full_B_chem ~ res_full_A_chem + df_chem$sugar))

rope_margin_chem <- 0.1 * pooled_sd(NN_full_A_chem, NN_full_B_chem)
rope_lower_chem <- -rope_margin_chem
rope_upper_chem <- rope_margin_chem

# RE:THINC L1
cat("Start Fractional Resampling Layer 1...\n")
boot_diffs_chem <- mclapply(1:n_iter, function(i) {
  df_sub <- df_chem[sample(nrow(df_chem), size = floor(0.8 * nrow(df_chem)), replace = FALSE), ]
  
  rf_A <- ranger(density ~ alcohol + sugar, data = df_sub, num.threads = 1)
  rf_B <- ranger(alcohol ~ density + sugar, data = df_sub, num.threads = 1)
  
  res_A <- df_sub$density - rf_A$predictions
  res_B <- df_sub$alcohol - rf_B$predictions
  
  df_res <- data.frame(res_A = res_A, res_B = res_B, sugar = df_sub$sugar)
  
  NN_A <- res_A - ranger(res_A ~ res_B + sugar, data = df_res, num.threads = 1)$predictions
  NN_B <- res_B - ranger(res_B ~ res_A + sugar, data = df_res, num.threads = 1)$predictions
  
  mean(abs(NN_A)) - mean(abs(NN_B))
}, mc.cores = n_cores) |> unlist()

boot_lower_chem <- quantile(boot_diffs_chem, 0.025)
boot_upper_chem <- quantile(boot_diffs_chem, 0.975)

#Decision L1
cat("\n======================================================\n")
cat("95% Analytical ROPE (H0): [", round(rope_lower_chem, 5), ",", round(rope_upper_chem, 5), "]\n")
cat("95% Bootstrap Signal:     [", round(boot_lower_chem, 5), ",", round(boot_upper_chem, 5), "]\n\n")

if(boot_upper_chem < 0 | boot_lower_chem > 0) {
  if(boot_upper_chem < rope_lower_chem | boot_lower_chem > rope_upper_chem) {
    if(boot_upper_chem < 0) print("Result: Alcohol -> Density (MEANINGFUL True Causal Direction!)")
    if(boot_lower_chem > 0) print("Result: Density -> Alcohol (MEANINGFUL True Causal Direction!)")
  } else {
    if(boot_upper_chem < 0) print("Result: Alcohol -> Density (Direction identified, but NOT practically meaningful)")
    if(boot_lower_chem > 0) print("Result: Density -> Alcohol (Direction identified, but NOT practically meaningful)")
  }
} else {
  print("Result: Undecided (Confidence interval contains 0 / Spurious Correlation)")
}
cat("======================================================\n\n")

# Layer 2 & 3
cat("    Generating topological distributions Layer 2 & 3...\n")
results_list_chem <- mclapply(1:n_boot, function(i) {
  df_sub <- df_chem[sample(nrow(df_chem), size = floor(0.8 * nrow(df_chem)), replace = FALSE), ]
  
  rf_A_t <- ranger(density ~ alcohol + sugar, data = df_sub, num.threads = 1)
  rf_B_t <- ranger(alcohol ~ density + sugar, data = df_sub, num.threads = 1)
  
  res_A_t <- df_sub$density - rf_A_t$predictions
  res_B_t <- df_sub$alcohol - rf_B_t$predictions
  
  df_res_t <- data.frame(res_A = res_A_t, res_B = res_B_t, sugar = df_sub$sugar)
  
  NN_A_t <- res_A_t - ranger(res_A ~ res_B + sugar, data = df_res_t, num.threads = 1)$predictions
  NN_B_t <- res_B_t - ranger(res_B ~ res_A + sugar, data = df_res_t, num.threads = 1)$predictions
  
  df_model_perm <- cbind(df_res_t, df_sub)
  df_model_perm$res_A <- sample(df_model_perm$res_A)
  df_model_perm$res_B <- sample(df_model_perm$res_B)
  # Sugar is NOT permuted to keep the marginal structural integrity!
  
  NN_A_perm <- res_A_t - ranger(res_A ~ res_B + sugar, data = df_model_perm, num.threads = 1)$predictions
  NN_B_perm <- res_B_t - ranger(res_B ~ res_A + sugar, data = df_model_perm, num.threads = 1)$predictions
  
  NNN_A_true <- NN_A_t - ranger(NN_A_t ~ density + sugar, data = df_sub, num.threads = 1)$predictions
  NNN_B_true <- NN_B_t - ranger(NN_B_t ~ alcohol + sugar, data = df_sub, num.threads = 1)$predictions
  
  df_model_null <- df_sub
  df_model_null$density <- sample(df_model_null$density)
  df_model_null$alcohol <- sample(df_model_null$alcohol)
  
  NNN_A_null <- NN_A_t - ranger(NN_A_t ~ density + sugar, data = df_model_null, num.threads = 1)$predictions
  NNN_B_null <- NN_B_t - ranger(NN_B_t ~ alcohol + sugar, data = df_model_null, num.threads = 1)$predictions
  
  return(c(
    MAD_L1_A = mean(abs(NN_A_perm)), MAD_L2_A = mean(abs(NN_A_t)),
    MAD_L1_B = mean(abs(NN_B_perm)), MAD_L2_B = mean(abs(NN_B_t)),
    NNN_true = mean(abs(NNN_A_true)) - mean(abs(NNN_B_true)),
    NNN_null = mean(abs(NNN_A_null)) - mean(abs(NNN_B_null))
  ))
}, mc.cores = n_cores)

res_mat_chem <- do.call(rbind, results_list_chem)
cat("\n... done!\n")

ci_drop_A_chem <- quantile((1 - (res_mat_chem[, "MAD_L2_A"] / res_mat_chem[, "MAD_L1_A"])) * 100, probs = c(0.025, 0.5, 0.975))
ci_drop_B_chem <- quantile((1 - (res_mat_chem[, "MAD_L2_B"] / res_mat_chem[, "MAD_L1_B"])) * 100, probs = c(0.025, 0.5, 0.975))

cat("\n    [CER-Check] Confounder-Proportion in Density: ", round(ci_drop_A_chem[2], 2), "%  (95% CI: [", round(ci_drop_A_chem[1], 2), " , ", round(ci_drop_A_chem[3], 2), "])\n")
cat("    [CER-Check] Confounder-Proportion in Alcohol: ", round(ci_drop_B_chem[2], 2), "%  (95% CI: [", round(ci_drop_B_chem[1], 2), " , ", round(ci_drop_B_chem[3], 2), "])\n")

ratio_vec_chem <- numeric(10000)
for (i in 1:10000) {
  idx <- sample(1:n_boot, replace = TRUE)
  ratio_vec_chem[i] <- mean(abs(res_mat_chem[idx, "NNN_true"])) / mean(abs(res_mat_chem[idx, "NNN_null"]))
}
ci_rel_chem <- quantile(ratio_vec_chem, probs = c(0.025, 0.975))

cat("\n    95% CI of the relative proportion (True / Null): [", round(ci_rel_chem[1], 4), ",", round(ci_rel_chem[2], 4), "]\n")
if (ci_rel_chem[2] < 1) cat("Layer 3 conclusion: vicious cycle / bidirectional feedback loop\n") else if (ci_rel_chem[1] > 1) cat("Layer 3 conclusion: unidirectional + adversarial noise / heteroscedasticity\n") else cat("Layer 3 conclusion: unidirection + symmetric noise\n")

# Paper plot
plot_data_chem <- data.frame(Value = boot_diffs_chem)
p_chem <- ggplot(plot_data_chem, aes(x = Value)) +
  annotate("rect", xmin = rope_lower_chem, xmax = rope_upper_chem, ymin = 0, ymax = Inf, fill = "darkgray", alpha = 0.3) +
  geom_density(fill = "#16a085", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
  geom_vline(xintercept = c(rope_lower_chem, rope_upper_chem), linetype = "dashed", color = "#2c3e50") +
  theme_minimal() +
  labs(title = "Physical Chemistry: Alcohol Concentration vs. Density", subtitle = "Teal: Subsampled Causal Signal | Gray Zone: Analytical ROPE (0.1 SD)", x = "Difference in Noise Distortion (Delta NN)", y = "Density")
print(p_chem)


################################################################################
# TEST 6: Biochemistry (Sugar vs. Alcohol)
################################################################################

df_6 <- df_chem %>% select(alcohol, sugar)
N_samples_6 <- nrow(df_6)
cat("\nStarting TEST 6: Biochemistry (Sugar vs Alcohol, N =", N_samples_6, ")\n")

# Analytical ROPE (H0-Margin for Effect Size)
rf_full_A_6 <- ranger(alcohol ~ sugar, data = df_6)
rf_full_B_6 <- ranger(sugar ~ alcohol, data = df_6)

res_full_A_6 <- df_6$alcohol - rf_full_A_6$predictions
res_full_B_6 <- df_6$sugar - rf_full_B_6$predictions

NN_full_A_6 <- residuals(lm(res_full_A_6 ~ res_full_B_6))
NN_full_B_6 <- residuals(lm(res_full_B_6 ~ res_full_A_6))

rope_margin_6 <- 0.1 * pooled_sd(NN_full_A_6, NN_full_B_6)
rope_lower_6 <- -rope_margin_6
rope_upper_6 <- rope_margin_6

# RE:THINC L1
cat("Start Fractional Resampling Layer 1...\n")
boot_diffs_6 <- mclapply(1:n_iter, function(i) {
  df_sub <- df_6[sample(nrow(df_6), size = floor(0.8 * nrow(df_6)), replace = FALSE), ]
  
  rf_A <- ranger(alcohol ~ sugar, data = df_sub, num.threads = 1)
  rf_B <- ranger(sugar ~ alcohol, data = df_sub, num.threads = 1)
  
  res_A <- df_sub$alcohol - rf_A$predictions
  res_B <- df_sub$sugar - rf_B$predictions
  
  df_res <- data.frame(res_A = res_A, res_B = res_B)
  
  NN_A <- res_A - ranger(res_A ~ res_B, data = df_res, num.threads = 1)$predictions
  NN_B <- res_B - ranger(res_B ~ res_A, data = df_res, num.threads = 1)$predictions
  
  mean(abs(NN_A)) - mean(abs(NN_B))
}, mc.cores = n_cores) |> unlist()

boot_lower_6 <- quantile(boot_diffs_6, 0.025)
boot_upper_6 <- quantile(boot_diffs_6, 0.975)

#Decision L1
cat("\n======================================================\n")
cat("95% Analytical ROPE (H0): [", round(rope_lower_6, 5), ",", round(rope_upper_6, 5), "]\n")
cat("95% Bootstrap Signal:     [", round(boot_lower_6, 5), ",", round(boot_upper_6, 5), "]\n\n")

if(boot_upper_6 < 0 | boot_lower_6 > 0) {
  if(boot_upper_6 < rope_lower_6 | boot_lower_6 > rope_upper_6) {
    if(boot_upper_6 < 0) print("Result: Sugar -> Alcohol (MEANINGFUL True Causal Direction!)")
    if(boot_lower_6 > 0) print("Result: Alcohol -> Sugar (MEANINGFUL True Causal Direction!)")
  } else {
    if(boot_upper_6 < 0) print("Result: Sugar -> Alcohol (Direction identified, but NOT practically meaningful)")
    if(boot_lower_6 > 0) print("Result: Alcohol -> Sugar (Direction identified, but NOT practically meaningful)")
  }
} else {
  print("Result: Undecided (Confidence interval contains 0 / Spurious Correlation)")
}
cat("======================================================\n\n")

# Layer 2 & 3
cat("    Generating topological distributions Layer 2 & 3...\n")
results_list_6 <- mclapply(1:n_boot, function(i) {
  df_sub <- df_6[sample(nrow(df_6), size = floor(0.8 * nrow(df_6)), replace = FALSE), ]
  
  rf_A_t <- ranger(alcohol ~ sugar, data = df_sub, num.threads = 1)
  rf_B_t <- ranger(sugar ~ alcohol, data = df_sub, num.threads = 1)
  
  res_A_t <- df_sub$alcohol - rf_A_t$predictions
  res_B_t <- df_sub$sugar - rf_B_t$predictions
  
  df_res_t <- data.frame(res_A = res_A_t, res_B = res_B_t)
  
  NN_A_t <- res_A_t - ranger(res_A ~ res_B, data = df_res_t, num.threads = 1)$predictions
  NN_B_t <- res_B_t - ranger(res_B ~ res_A, data = df_res_t, num.threads = 1)$predictions
  
  df_model_perm <- cbind(df_res_t, df_sub)
  df_model_perm$res_A <- sample(df_model_perm$res_A)
  df_model_perm$res_B <- sample(df_model_perm$res_B)
  
  NN_A_perm <- res_A_t - ranger(res_A ~ res_B, data = df_model_perm, num.threads = 1)$predictions
  NN_B_perm <- res_B_t - ranger(res_B ~ res_A, data = df_model_perm, num.threads = 1)$predictions
  
  NNN_A_true <- NN_A_t - ranger(NN_A_t ~ alcohol, data = df_sub, num.threads = 1)$predictions
  NNN_B_true <- NN_B_t - ranger(NN_B_t ~ sugar, data = df_sub, num.threads = 1)$predictions
  
  df_model_null <- df_sub
  df_model_null$alcohol <- sample(df_model_null$alcohol)
  df_model_null$sugar <- sample(df_model_null$sugar)
  
  NNN_A_null <- NN_A_t - ranger(NN_A_t ~ alcohol, data = df_model_null, num.threads = 1)$predictions
  NNN_B_null <- NN_B_t - ranger(NN_B_t ~ sugar, data = df_model_null, num.threads = 1)$predictions
  
  return(c(
    MAD_L1_A = mean(abs(NN_A_perm)), MAD_L2_A = mean(abs(NN_A_t)),
    MAD_L1_B = mean(abs(NN_B_perm)), MAD_L2_B = mean(abs(NN_B_t)),
    NNN_true = mean(abs(NNN_A_true)) - mean(abs(NNN_B_true)),
    NNN_null = mean(abs(NNN_A_null)) - mean(abs(NNN_B_null))
  ))
}, mc.cores = n_cores)

res_mat_6 <- do.call(rbind, results_list_6)
cat("\n... done!\n")

ci_drop_A_6 <- quantile((1 - (res_mat_6[, "MAD_L2_A"] / res_mat_6[, "MAD_L1_A"])) * 100, probs = c(0.025, 0.5, 0.975))
ci_drop_B_6 <- quantile((1 - (res_mat_6[, "MAD_L2_B"] / res_mat_6[, "MAD_L1_B"])) * 100, probs = c(0.025, 0.5, 0.975))

cat("\n    [CER-Check] Confounder-Proportion in Alcohol: ", round(ci_drop_A_6[2], 2), "%  (95% CI: [", round(ci_drop_A_6[1], 2), " , ", round(ci_drop_A_6[3], 2), "])\n")
cat("    [CER-Check] Confounder-Proportion in Sugar: ", round(ci_drop_B_6[2], 2), "%  (95% CI: [", round(ci_drop_B_6[1], 2), " , ", round(ci_drop_B_6[3], 2), "])\n")

ratio_vec_6 <- numeric(10000)
for (i in 1:10000) {
  idx <- sample(1:n_boot, replace = TRUE)
  ratio_vec_6[i] <- mean(abs(res_mat_6[idx, "NNN_true"])) / mean(abs(res_mat_6[idx, "NNN_null"]))
}
ci_rel_6 <- quantile(ratio_vec_6, probs = c(0.025, 0.975))

cat("\n    95% CI of the relative proportion (True / Null): [", round(ci_rel_6[1], 4), ",", round(ci_rel_6[2], 4), "]\n")
if (ci_rel_6[2] < 1) cat("Layer 3 conclusion: vicious cycle / bidirectional feedback loop\n") else if (ci_rel_6[1] > 1) cat("Layer 3 conclusion: unidirectional + adversarial noise / heteroscedasticity\n") else cat("Layer 3 conclusion: unidirection + symmetric noise\n")

# Paper plot
plot_data_6 <- data.frame(Value = boot_diffs_6)
p_6 <- ggplot(plot_data_6, aes(x = Value)) +
  annotate("rect", xmin = rope_lower_6, xmax = rope_upper_6, ymin = 0, ymax = Inf, fill = "darkgray", alpha = 0.3) +
  geom_density(fill = "#8e44ad", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
  geom_vline(xintercept = c(rope_lower_6, rope_upper_6), linetype = "dashed", color = "#2c3e50") +
  theme_minimal() +
  labs(title = "Biochemistry: Fermentation (Sugar vs. Alcohol)", subtitle = "Purple: Subsampled Causal Signal | Gray Zone: Analytical ROPE (0.1 SD)", x = "Difference in Noise Distortion (Delta NN)", y = "Density")
print(p_6)

################################################################################
# TEST 7: Solid State Physics (Spatial Volume vs. Mass)
################################################################################

# 1. Data Prep (building a deterministic closed system)
set.seed(1)
df_phys <- diamonds %>% 
  # Filtering error measures
  filter(x > 0 & y > 0 & z > 0) %>% 
  # Volume estimation
  mutate(volume = x * y * z) %>% 
  select(carat, volume) %>% 
  filter(complete.cases(.)) %>%
  mutate(across(everything(), ~ as.numeric(scale(.))))

df_7 <- df_phys
N_samples_7 <- nrow(df_7)
cat("\nStarting TEST 7: Solid State Physics (Spatial Volume vs. Mass, N =", N_samples_7, ")\n")

# Analytical ROPE (H0-Margin for Effect Size)
rf_full_A_7 <- ranger(carat ~ volume, data = df_7)
rf_full_B_7 <- ranger(volume ~ carat, data = df_7)

res_full_A_7 <- df_7$carat - rf_full_A_7$predictions
res_full_B_7 <- df_7$volume - rf_full_B_7$predictions

NN_full_A_7 <- residuals(lm(res_full_A_7 ~ res_full_B_7))
NN_full_B_7 <- residuals(lm(res_full_B_7 ~ res_full_A_7))

# Assuming pooled_sd is loaded in your environment (as in 6)
rope_margin_7 <- 0.1 * pooled_sd(NN_full_A_7, NN_full_B_7)
rope_lower_7 <- -rope_margin_7
rope_upper_7 <- rope_margin_7

# RE:THINC L1
cat("Start Fractional Resampling Layer 1...\n")
boot_diffs_7_raw <- mclapply(1:n_iter, function(i) {
  df_sub <- df_7[sample(nrow(df_7), size = floor(0.8 * nrow(df_7)), replace = FALSE), ]
  
  rf_A <- ranger(carat ~ volume, data = df_sub, num.threads = 1)
  rf_B <- ranger(volume ~ carat, data = df_sub, num.threads = 1)
  
  res_A <- df_sub$carat - rf_A$predictions
  res_B <- df_sub$volume - rf_B$predictions
  
  df_res <- data.frame(res_A = res_A, res_B = res_B)
  
  NN_A <- res_A - ranger(res_A ~ res_B, data = df_res, num.threads = 1)$predictions
  NN_B <- res_B - ranger(res_B ~ res_A, data = df_res, num.threads = 1)$predictions
  
  #Single MAD output for cover art
  c(MAD_A = mean(abs(NN_A)), MAD_B = mean(abs(NN_B)))
}, mc.cores = n_cores)

res_L1 <- do.call(rbind, boot_diffs_7_raw)
boot_diffs_7 <- res_L1[, "MAD_A"] - res_L1[, "MAD_B"]

boot_lower_7 <- quantile(boot_diffs_7, 0.025)
boot_upper_7 <- quantile(boot_diffs_7, 0.975)

# Decision L1
cat("\n======================================================\n")
cat("95% Analytical ROPE (H0): [", round(rope_lower_7, 5), ",", round(rope_upper_7, 5), "]\n")
cat("95% Bootstrap Signal:     [", round(boot_lower_7, 5), ",", round(boot_upper_7, 5), "]\n\n")

if(boot_upper_7 < 0 | boot_lower_7 > 0) {
  if(boot_upper_7 < rope_lower_7 | boot_lower_7 > rope_upper_7) {
    if(boot_upper_7 < 0) print("Result: Volume -> Mass (Carat) (MEANINGFUL True Causal Direction!)")
    if(boot_lower_7 > 0) print("Result: Mass (Carat) -> Volume (MEANINGFUL True Causal Direction!)")
  } else {
    if(boot_upper_7 < 0) print("Result: Volume -> Mass (Carat) (Direction identified, but NOT practically meaningful)")
    if(boot_lower_7 > 0) print("Result: Mass (Carat) -> Volume (Direction identified, but NOT practically meaningful)")
  }
} else {
  print("Result: Undecided (Confidence interval contains 0 / Spurious Correlation)")
}
cat("======================================================\n\n")

# Layer 2 & 3
cat("    Generating topological distributions Layer 2 & 3...\n")
results_list_7 <- mclapply(1:n_boot, function(i) {
  df_sub <- df_7[sample(nrow(df_7), size = floor(0.8 * nrow(df_7)), replace = FALSE), ]
  
  rf_A_t <- ranger(carat ~ volume, data = df_sub, num.threads = 1)
  rf_B_t <- ranger(volume ~ carat, data = df_sub, num.threads = 1)
  
  res_A_t <- df_sub$carat - rf_A_t$predictions
  res_B_t <- df_sub$volume - rf_B_t$predictions
  
  df_res_t <- data.frame(res_A = res_A_t, res_B = res_B_t)
  
  NN_A_t <- res_A_t - ranger(res_A ~ res_B, data = df_res_t, num.threads = 1)$predictions
  NN_B_t <- res_B_t - ranger(res_B ~ res_A, data = df_res_t, num.threads = 1)$predictions
  
  df_model_perm <- cbind(df_res_t, df_sub)
  df_model_perm$res_A <- sample(df_model_perm$res_A)
  df_model_perm$res_B <- sample(df_model_perm$res_B)
  
  NN_A_perm <- res_A_t - ranger(res_A ~ res_B, data = df_model_perm, num.threads = 1)$predictions
  NN_B_perm <- res_B_t - ranger(res_B ~ res_A, data = df_model_perm, num.threads = 1)$predictions
  
  NNN_A_true <- NN_A_t - ranger(NN_A_t ~ carat, data = df_sub, num.threads = 1)$predictions
  NNN_B_true <- NN_B_t - ranger(NN_B_t ~ volume, data = df_sub, num.threads = 1)$predictions
  
  df_model_null <- df_sub
  df_model_null$carat <- sample(df_model_null$carat)
  df_model_null$volume <- sample(df_model_null$volume)
  
  NNN_A_null <- NN_A_t - ranger(NN_A_t ~ carat, data = df_model_null, num.threads = 1)$predictions
  NNN_B_null <- NN_B_t - ranger(NN_B_t ~ volume, data = df_model_null, num.threads = 1)$predictions
  
  return(c(
    MAD_L1_A = mean(abs(NN_A_perm)), MAD_L2_A = mean(abs(NN_A_t)),
    MAD_L1_B = mean(abs(NN_B_perm)), MAD_L2_B = mean(abs(NN_B_t)),
    NNN_true = mean(abs(NNN_A_true)) - mean(abs(NNN_B_true)),
    NNN_null = mean(abs(NNN_A_null)) - mean(abs(NNN_B_null)),
    MAD_NNN_true_A = mean(abs(NNN_A_true)),
    MAD_NNN_null_A = mean(abs(NNN_A_null))
  ))
}, mc.cores = n_cores)

res_mat_7 <- do.call(rbind, results_list_7)
cat("\n... done!\n")

ci_drop_A_7 <- quantile((1 - (res_mat_7[, "MAD_L2_A"] / res_mat_7[, "MAD_L1_A"])) * 100, probs = c(0.025, 0.5, 0.975))
ci_drop_B_7 <- quantile((1 - (res_mat_7[, "MAD_L2_B"] / res_mat_7[, "MAD_L1_B"])) * 100, probs = c(0.025, 0.5, 0.975))

cat("\n    [CER-Check] Confounder-Proportion in Mass (Carat): ", round(ci_drop_A_7[2], 2), "%  (95% CI: [", round(ci_drop_A_7[1], 2), " , ", round(ci_drop_A_7[3], 2), "])\n")
cat("    [CER-Check] Confounder-Proportion in Volume: ", round(ci_drop_B_7[2], 2), "%  (95% CI: [", round(ci_drop_B_7[1], 2), " , ", round(ci_drop_B_7[3], 2), "])\n")

ratio_vec_7 <- numeric(10000)
for (i in 1:10000) {
  # Sampling with replacement using the exact number of rows in res_mat_7
  idx <- sample(1:nrow(res_mat_7), replace = TRUE)
  ratio_vec_7[i] <- mean(abs(res_mat_7[idx, "NNN_true"])) / mean(abs(res_mat_7[idx, "NNN_null"]))
}
ci_rel_7 <- quantile(ratio_vec_7, probs = c(0.025, 0.975))

cat("\n    95% CI of the relative proportion (True / Null): [", round(ci_rel_7[1], 4), ",", round(ci_rel_7[2], 4), "]\n")
if (ci_rel_7[2] < 1) cat("Layer 3 conclusion: vicious cycle / bidirectional feedback loop\n") else if (ci_rel_7[1] > 1) cat("Layer 3 conclusion: unidirectional + adversarial noise / heteroscedasticity\n") else cat("Layer 3 conclusion: unidirection + symmetric noise\n")

################################################################################
# NOVOENTROPIC SYNTAX KEY GRAPHIC GENERATOR
################################################################################

# PLOT 1: Acyclic Asymmetry (Pink vs Gray on Black)
plot_cover_1 <- data.frame(
  MAD = c(res_L1[, "MAD_A"], res_L1[, "MAD_B"]),
  Direction = factor(rep(c("True Direction", "Anticausal"), each = nrow(res_L1)), 
                     levels = c("Anticausal", "True Direction"))
)

p_cover_1 <- ggplot(plot_cover_1, aes(x = MAD, fill = Direction, color = Direction)) +
  geom_density(alpha = 0.6, linewidth = 0.8) +
  scale_fill_manual(values = c("Anticausal" = "#888888", "True Direction" = "#F69")) +
  scale_color_manual(values = c("Anticausal" = "#888888", "True Direction" = "#F69")) +
  theme_void() + 
  theme(
    plot.background = element_rect(fill = "#0a0a0a", color = NA),
    panel.background = element_rect(fill = "#0a0a0a", color = NA),
    legend.position = "bottom",
    legend.text = element_text(color = "white", size = 12, face = "bold"),
    legend.title = element_blank(),
    plot.margin = margin(30, 30, 30, 30)
  )

# PLOT 2: Cyclic Feedback Loop
plot_cover_2 <- data.frame(
  MAD = c(res_mat_7[, "MAD_NNN_true_A"], res_mat_7[, "MAD_NNN_null_A"]),
  Topology = factor(rep(c("True Topology", "Permuted Topology"), each = nrow(res_mat_7)),
                    levels = c("Permuted Topology", "True Topology"))
)

p_cover_2 <- ggplot(plot_cover_2, aes(x = MAD, fill = Topology, color = Topology)) +
  geom_density(alpha = 0.6, linewidth = 0.8) +
  scale_fill_manual(values = c("Permuted Topology" = "#888888", "True Topology" = "#89cff0")) +
  scale_color_manual(values = c("Permuted Topology" = "#888888", "True Topology" = "#89cff0")) +
  theme_void() + 
  theme(
    plot.background = element_rect(fill = "#0a0a0a", color = NA),
    panel.background = element_rect(fill = "#0a0a0a", color = NA),
    legend.position = "bottom",
    legend.text = element_text(color = "white", size = 12, face = "bold"),
    legend.title = element_blank(),
    plot.margin = margin(30, 30, 30, 30)
  )

# Output Plots in R
print(p_cover_1)
print(p_cover_2)

################################################################################
# 8. Power Analysis (Monte Carlo Simulation for RE:THINC)
################################################################################

# Defining simulation parameters
n_simulations <- 1000   
n_iter <- 100          

sample_sizes <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200,
                  1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200,
                  2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200,
                  3300, 3400, 3500)
effect_sizes <- c(0.20, 0.50, 0.80) # Small, Medium, Large Asymmetry

# Dataframe for final results
power_results <- data.frame()

# Parallelisation setup
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

#Plasmode setup
rf_base_model <- ranger(BPSysAve ~ BMI, data = df_med)

#true biological curve
true_signal_std <- as.numeric(scale(rf_base_model$predictions)) 

#true noise
true_oob_residuals <- as.numeric(scale(df_med$BPSysAve - rf_base_model$predictions))

#real x
real_X_std <- as.numeric(scale(df_med$BMI))

for(es in effect_sizes) {
  for(n in sample_sizes) {
    
    cat("Running Simulation for N =", n, "| Effect Size =", es, "...\n")
    
    sim_results <- foreach(sim = 1:n_simulations, .combine = rbind, .packages = c("ranger", "dplyr")) %dopar% {
      
    #Semi-synthetic data generation
      
      #Draw random row indices
      sim_indices <- sample(1:length(real_X_std), size = n, replace = TRUE)
      
      #Draw X and the corresponding true signal and noise
      X_sim <- real_X_std[sim_indices]
      Signal_sim <- true_signal_std[sim_indices]
      Noise_sim <- true_oob_residuals[sim_indices]
      
      #Causal injection
      Y_sim <- (es * Signal_sim) + sqrt(1 - es^2) * Noise_sim
      
      #Final dataset
      df <- data.frame(X = X_sim, Y = Y_sim) %>% 
        mutate(across(everything(), ~ as.numeric(scale(.))))
      
      # Fractional Resampling
      boot_diffs <- numeric(n_iter)
      sub_size <- floor(0.8 * nrow(df)) 
      
      for(i in 1:n_iter) {
        df_sub <- df[sample(nrow(df), size = sub_size, replace = FALSE), ]
        
        rf_boot_A <- ranger(Y ~ X, data = df_sub, num.trees = 150, num.threads = 1)
        rf_boot_B <- ranger(X ~ Y, data = df_sub, num.trees = 150, num.threads = 1)
        
        res_A <- df_sub$Y - rf_boot_A$predictions
        res_B <- df_sub$X - rf_boot_B$predictions
        
        df_res <- data.frame(res_A = res_A, res_B = res_B)
        
        NN_A <- res_A - ranger(res_A ~ res_B, data = df_res, num.threads = 1)$predictions
        NN_B <- res_B - ranger(res_B ~ res_A, data = df_res, num.threads = 1)$predictions
        
        mad_A <- mean(abs(NN_A))
        mad_B <- mean(abs(NN_B))
        
        boot_diffs[i] <- mad_A - mad_B
      }
      
      boot_lower <- quantile(boot_diffs, 0.025)
      boot_upper <- quantile(boot_diffs, 0.975)
      
      # Decision Rule
      if(boot_upper < 0) {
        decision <- "X->Y"
      } else if(boot_lower > 0) {
        decision <- "Y->X"
      } else {
        decision <- "Undecided"
      }
      
      # Return simple vector for speed and memory
      c(Correct = ifelse(decision == "X->Y", 1, 0),
        Wrong = ifelse(decision == "Y->X", 1, 0),
        Undecided = ifelse(decision == "Undecided", 1, 0))
    }
    
    # 3. Aggregate Results for this grid point
    power <- sum(sim_results[, "Correct"]) / n_simulations
    false_rate <- sum(sim_results[, "Wrong"]) / n_simulations
    undecided_rate <- sum(sim_results[, "Undecided"]) / n_simulations
    
    power_results <- rbind(power_results, data.frame(
      N = n,
      EffectSize = es,
      Power = power,
      FalseRate = false_rate,
      UndecidedRate = undecided_rate
    ))
  }
}

# Closing the clusters
stopCluster(cl)
print(power_results)

# 4. Plotting the Power Curves for the Paper
power_results$EffectLabel <- factor(power_results$EffectSize, 
                                    levels = c(0.2, 0.5, 0.8), 
                                    labels = c("Small Effect (d ~ 0.20)", 
                                               "Medium Effect (d ~ 0.50)", 
                                               "Large Effect (d ~ 0.80)"))

p_power <- ggplot(power_results, aes(x = N, y = Power, color = EffectLabel, group = EffectLabel)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "black", linewidth = 0.8) +
  annotate("text", x = max(sample_sizes) * 0.9, y = 0.83, label = "80% Power Threshold", fontface = "italic") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_x_continuous(breaks = sample_sizes) +
  scale_color_manual(values = c("#e74c3c", "#f39c12", "#27ae60")) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Statistical Power of the RE:THINC Algorithm",
    subtitle = "Semi-Synthetic Plasmode Simulation (True Causal Recovery: BMI \u2192 Blood Pressure)",
    x = "Sample Size (N)",
    y = "Statistical Power (Correct Identification Rate)"
  )

print(p_power)

################################################################################
#Separate Power Plots
################################################################################
df <- read.csv("")
df$N <- as.numeric(as.character(df$N))

#Creating new labels
df$EffectLabel <- factor(df$EffectSize, 
                         levels = c(0.2, 0.5, 0.8), 
                         labels = c("20% Signal Strength (w = 0.2)", 
                                    "50% Signal Strength (w = 0.5)", 
                                    "80% Signal Strength (w = 0.8)"))

#Plot with three separate windows
p_power <- ggplot(df, aes(x = N, y = Power, color = EffectLabel)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "black", linewidth = 0.8) +
  
  #Creates three plots
  facet_wrap(~ EffectLabel, ncol = 3) + 
  
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 3500, by = 1000)) +
  scale_color_manual(values = c("#e74c3c", "#f39c12", "#27ae60")) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 13, face = "bold"),
    panel.spacing = unit(2, "lines")
  ) +
  labs(
    title = "Statistical Power of the RE:THINC Algorithm",
    subtitle = "Plasmode Simulation: Attenuation of the Empirical BMI \u2192 Blood Pressure Signal",
    x = "Sample Size (N)",
    y = "Statistical Power (Correct Identification Rate)"
  )

#Saves the plots
print(p_power)
ggsave("plot_power_results.pdf", plot = p_power, width = 15, height = 6)

################################################################################
# TEST 9: Host Galaxies and Super Massive Black Holes
################################################################################

# 1. Data Prep (building a deterministic closed system)
set.seed(1)
df_raw <- readFrameFromFITS("")

df_astro <- df_raw %>%
  select(LOGMBH, OIII5007.3) %>%
  dplyr::rename(Mass_BH = LOGMBH, Mass_Galaxy = OIII5007.3) %>%
  mutate(
    Mass_BH = as.numeric(Mass_BH),
    Mass_Galaxy = as.numeric(Mass_Galaxy)
  ) %>%
  filter(is.finite(Mass_BH) & is.finite(Mass_Galaxy)) %>%
  filter(Mass_BH > 0 & Mass_Galaxy > 0) %>%
  dplyr::slice_sample(n = 10000) %>%
  mutate(across(everything(), ~ as.numeric(scale(.)))) %>%
  filter(complete.cases(.))

rm(df_raw); gc()

df_9 <- df_astro
N_samples_9 <- nrow(df_9)
cat("\nStarting TEST 9: Host Galaxies and Super Massive Black Holes (N =", N_samples_9, ")\n")

# Analytical ROPE (H0-Margin for Effect Size)
rf_full_A_9 <- ranger(Mass_Galaxy ~ Mass_BH, data = df_9)
rf_full_B_9 <- ranger(Mass_BH ~ Mass_Galaxy, data = df_9)

res_full_A_9 <- df_9$Mass_Galaxy - rf_full_A_9$predictions
res_full_B_9 <- df_9$Mass_BH - rf_full_B_9$predictions

NN_full_A_9 <- residuals(lm(res_full_A_9 ~ res_full_B_9))
NN_full_B_9 <- residuals(lm(res_full_B_9 ~ res_full_A_9))

# Assuming pooled_sd is loaded in your environment (as in 6/7)
rope_margin_9 <- 0.1 * pooled_sd(NN_full_A_9, NN_full_B_9)
rope_lower_9 <- -rope_margin_9
rope_upper_9 <- rope_margin_9

# RE:THINC L1
cat("Start Fractional Resampling Layer 1...\n")
boot_diffs_9 <- mclapply(1:n_iter, function(i) {
  df_sub <- df_9[sample(nrow(df_9), size = floor(0.8 * nrow(df_9)), replace = FALSE), ]
  
  rf_A <- ranger(Mass_Galaxy ~ Mass_BH, data = df_sub, num.threads = 1)
  rf_B <- ranger(Mass_BH ~ Mass_Galaxy, data = df_sub, num.threads = 1)
  
  res_A <- df_sub$Mass_Galaxy - rf_A$predictions
  res_B <- df_sub$Mass_BH - rf_B$predictions
  
  df_res <- data.frame(res_A = res_A, res_B = res_B)
  
  NN_A <- res_A - ranger(res_A ~ res_B, data = df_res, num.threads = 1)$predictions
  NN_B <- res_B - ranger(res_B ~ res_A, data = df_res, num.threads = 1)$predictions
  
  mean(abs(NN_A)) - mean(abs(NN_B))
}, mc.cores = n_cores) |> unlist()

boot_lower_9 <- quantile(boot_diffs_9, 0.025)
boot_upper_9 <- quantile(boot_diffs_9, 0.975)

# Decision L1
cat("\n======================================================\n")
cat("95% Analytical ROPE (H0): [", round(rope_lower_9, 5), ",", round(rope_upper_9, 5), "]\n")
cat("95% Bootstrap Signal:     [", round(boot_lower_9, 5), ",", round(boot_upper_9, 5), "]\n\n")

if(boot_upper_9 < 0 | boot_lower_9 > 0) {
  if(boot_upper_9 < rope_lower_9 | boot_lower_9 > rope_upper_9) {
    if(boot_upper_9 < 0) print("Result: Mass_BH -> Mass_Galaxy (MEANINGFUL True Causal Direction!)")
    if(boot_lower_9 > 0) print("Result: Mass_Galaxy -> Mass_BH (MEANINGFUL True Causal Direction!)")
  } else {
    if(boot_upper_9 < 0) print("Result: Mass_BH -> Mass_Galaxy (Direction identified, but NOT practically meaningful)")
    if(boot_lower_9 > 0) print("Result: Mass_Galaxy -> Mass_BH (Direction identified, but NOT practically meaningful)")
  }
} else {
  print("Result: Undecided (Confidence interval contains 0 / Spurious Correlation)")
}
cat("======================================================\n\n")

# Layer 2 & 3
cat("    Generating topological distributions Layer 2 & 3...\n")
results_list_9 <- mclapply(1:n_boot, function(i) {
  df_sub <- df_9[sample(nrow(df_9), size = floor(0.8 * nrow(df_9)), replace = FALSE), ]
  
  rf_A_t <- ranger(Mass_Galaxy ~ Mass_BH, data = df_sub, num.threads = 1)
  rf_B_t <- ranger(Mass_BH ~ Mass_Galaxy, data = df_sub, num.threads = 1)
  
  res_A_t <- df_sub$Mass_Galaxy - rf_A_t$predictions
  res_B_t <- df_sub$Mass_BH - rf_B_t$predictions
  
  df_res_t <- data.frame(res_A = res_A_t, res_B = res_B_t)
  
  NN_A_t <- res_A_t - ranger(res_A ~ res_B, data = df_res_t, num.threads = 1)$predictions
  NN_B_t <- res_B_t - ranger(res_B ~ res_A, data = df_res_t, num.threads = 1)$predictions
  
  df_model_perm <- cbind(df_res_t, df_sub)
  df_model_perm$res_A <- sample(df_model_perm$res_A)
  df_model_perm$res_B <- sample(df_model_perm$res_B)
  
  NN_A_perm <- res_A_t - ranger(res_A ~ res_B, data = df_model_perm, num.threads = 1)$predictions
  NN_B_perm <- res_B_t - ranger(res_B ~ res_A, data = df_model_perm, num.threads = 1)$predictions
  
  NNN_A_true <- NN_A_t - ranger(NN_A_t ~ Mass_Galaxy, data = df_sub, num.threads = 1)$predictions
  NNN_B_true <- NN_B_t - ranger(NN_B_t ~ Mass_BH, data = df_sub, num.threads = 1)$predictions
  
  df_model_null <- df_sub
  df_model_null$Mass_Galaxy <- sample(df_model_null$Mass_Galaxy)
  df_model_null$Mass_BH <- sample(df_model_null$Mass_BH)
  
  NNN_A_null <- NN_A_t - ranger(NN_A_t ~ Mass_Galaxy, data = df_model_null, num.threads = 1)$predictions
  NNN_B_null <- NN_B_t - ranger(NN_B_t ~ Mass_BH, data = df_model_null, num.threads = 1)$predictions
  
  return(c(
    MAD_L1_A = mean(abs(NN_A_perm)), MAD_L2_A = mean(abs(NN_A_t)),
    MAD_L1_B = mean(abs(NN_B_perm)), MAD_L2_B = mean(abs(NN_B_t)),
    NNN_true = mean(abs(NNN_A_true)) - mean(abs(NNN_B_true)),
    NNN_null = mean(abs(NNN_A_null)) - mean(abs(NNN_B_null))
  ))
}, mc.cores = n_cores)

res_mat_9 <- do.call(rbind, results_list_9)
cat("\n... done!\n")

ci_drop_A_9 <- quantile((1 - (res_mat_9[, "MAD_L2_A"] / res_mat_9[, "MAD_L1_A"])) * 100, probs = c(0.025, 0.5, 0.975))
ci_drop_B_9 <- quantile((1 - (res_mat_9[, "MAD_L2_B"] / res_mat_9[, "MAD_L1_B"])) * 100, probs = c(0.025, 0.5, 0.975))

cat("\n    [CER-Check] Confounder-Proportion in Mass Galaxy: ", round(ci_drop_A_9[2], 2), "%  (95% CI: [", round(ci_drop_A_9[1], 2), " , ", round(ci_drop_A_9[3], 2), "])\n")
cat("    [CER-Check] Confounder-Proportion in Mass BH: ", round(ci_drop_B_9[2], 2), "%  (95% CI: [", round(ci_drop_B_9[1], 2), " , ", round(ci_drop_B_9[3], 2), "])\n")

ratio_vec_9 <- numeric(10000)
for (i in 1:10000) {
  # Sampling with replacement using the exact number of rows in res_mat_9
  idx <- sample(1:nrow(res_mat_9), replace = TRUE)
  ratio_vec_9[i] <- mean(abs(res_mat_9[idx, "NNN_true"])) / mean(abs(res_mat_9[idx, "NNN_null"]))
}
ci_rel_9 <- quantile(ratio_vec_9, probs = c(0.025, 0.975))

cat("\n    95% CI of the relative proportion (True / Null): [", round(ci_rel_9[1], 4), ",", round(ci_rel_9[2], 4), "]\n")
if (ci_rel_9[2] < 1) cat("Layer 3 conclusion: vicious cycle / bidirectional feedback loop\n") else if (ci_rel_9[1] > 1) cat("Layer 3 conclusion: unidirectional + adversarial noise / heteroscedasticity\n") else cat("Layer 3 conclusion: unidirection + symmetric noise\n")

# Paper plot
plot_data_9 <- data.frame(Value = boot_diffs_9)
p_9 <- ggplot(plot_data_9, aes(x = Value)) +
  annotate("rect", xmin = rope_lower_9, xmax = rope_upper_9, ymin = 0, ymax = Inf, fill = "darkgray", alpha = 0.3) +
  geom_density(fill = "#c0392b", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
  geom_vline(xintercept = c(rope_lower_9, rope_upper_9), linetype = "dashed", color = "#2c3e50") +
  theme_minimal() +
  labs(title = "Astrophysics: Host Galaxies vs. Super Massive Black Holes", 
       subtitle = "Red: Subsampled Causal Signal | Gray Zone: Analytical ROPE (0.1 SD)", 
       x = "Difference in Noise Distortion (Delta NN)", 
       y = "Density")
print(p_9)

################################################################################
# RE:THINC - NHANES Causal Discovery Pipeline
# Target: Depression (PHQ-9) vs. Phyiological Biomarker
################################################################################

#Loads depression data
df_dep_raw <- nhanes('DPQ_J', translated = FALSE)

df_dep <- df_dep_raw %>%
  mutate(across(DPQ010:DPQ090, ~ ifelse(. >= 7, NA, as.numeric(.)))) %>%
  mutate(PHQ9_Score = rowSums(select(., DPQ010:DPQ090), na.rm = FALSE)) %>%
  select(SEQN, PHQ9_Score)

rm(df_dep_raw); gc()

################################################################################
#Configuration of test variables
################################################################################
biomarker_tasks <- list(
 # list(table = "HSCRP_J", col = "LBXHSCRP", name = "Inflammation (CRP)"),
#  list(table = "VID_J",   col = "LBXVIDMS", name = "Vitamin D"),
#  list(table = "BMX_J",   col = "BMXBMI",   name = "Obesity (BMI)"),
#  list(table = "GHB_J",   col = "LBXGH",    name = "Blood sugar (HbA1c)"),
  list(table = "COT_J",   col = "LBXCOT",   name = "Cotinine")
#  list(table = "PBCD_J",  col = "LBXBPB",   name = "Heavy Metal (Blood Lead)")
)

################################################################################
# L2-Orthogonal Projection
################################################################################
set.seed(1)
run_rethinc <- function(df_dep, task, n_iter = 100) {
  
  cat("\n=======================================================================\n")
  cat("Testing Hypothesis:", task$name, "vs. Depression (PHQ-9)\n")
  cat("Loading table:", task$table, "...\n")
  
  # Loading data
  df_raw <- suppressWarnings(nhanes(task$table, translated = FALSE))
  
  df_bio <- df_raw %>%
    mutate(Bio_Level = as.numeric(!!sym(task$col))) %>%
    select(SEQN, Bio_Level)
  
  # 2. Merging & Scaling
  df_merged <- inner_join(df_bio, df_dep, by = "SEQN") %>%
    filter(is.finite(Bio_Level) & is.finite(PHQ9_Score)) %>%
    mutate(
      Bio_Level = as.numeric(scale(Bio_Level)),
      PHQ9_Score = as.numeric(scale(PHQ9_Score))
    )
  
  N_samples <- nrow(df_merged)
  cat("Data (N):", N_samples, "\n")
  
  if(N_samples < 200) {
    cat("Warning: not enough data.\n")
    return(NULL)
  }
  
  # 3. Fractional Resampling (Layer 1 & 2 Quick-Scan)
  sample_size <- floor(0.8 * N_samples)
  n_cores <- max(1, detectCores() - 1)
  
  cat("Starting L2 Projection (", n_iter, " Iterations )...\n", sep="")
  
  boot_diffs <- mclapply(1:n_iter, function(i) {
    df_sub <- df_merged[sample(N_samples, size = sample_size, replace = FALSE), ]
    
    rf_A <- ranger(PHQ9_Score ~ Bio_Level, data = df_sub, num.threads = 1)
    rf_B <- ranger(Bio_Level ~ PHQ9_Score, data = df_sub, num.threads = 1)
    
    res_A <- df_sub$PHQ9_Score - rf_A$predictions
    res_B <- df_sub$Bio_Level - rf_B$predictions
    
    df_res <- data.frame(res_A = res_A, res_B = res_B)
    
    NN_A <- res_A - ranger(res_A ~ res_B, data = df_res, num.threads = 1)$predictions
    NN_B <- res_B - ranger(res_B ~ res_A, data = df_res, num.threads = 1)$predictions
    
    return(mean(abs(NN_A)) - mean(abs(NN_B)))
    
  }, mc.cores = n_cores) |> unlist()
  
  boot_lower <- quantile(boot_diffs, 0.025)
  boot_upper <- quantile(boot_diffs, 0.975)
  
  # ROPE Definition: 0.1 * Standard deviation of Bootstrap-differences
  rope_val <- 0.1 * sd(boot_diffs)
  rope_lower <- -rope_val
  rope_upper <- rope_val
  
  # 4. Decision Logic (Primary Driver)
  cat("\n95% Bootstrap Signal:    [", round(boot_lower, 5), ",", round(boot_upper, 5), "]\n")
  cat("95% ROPE (0.1 * SD):     [", round(rope_lower, 5), ",", round(rope_upper, 5), "]\n")
  
  if(boot_upper < rope_lower) {
    cat("Result:      ", task$name, " -> Depression (PHQ-9)\n")
    cat("CONCLUSION: DRIVER!\n")
  } else if(boot_lower > rope_upper) {
    cat("Result:      Depression (PHQ-9) ->", task$name, "\n")
    cat("CONCLUSION: PASSENGER!\n")
  } else {
    cat("Result:      Undecided (Bidirectional Loop or confounding)\n")
  }
  
  # ==========================================================================
  # DEEP DIVE: Structural Variance Ratio Test & CER-Check
  # ==========================================================================
  cat("\n--> Starting Deep Dive (CER-Check & Permutations-Null distribution)...\n")
  
  n_boot <- 1000 
  NNN_true <- numeric(n_boot)
  NNN_null <- numeric(n_boot) 
  NN_perm <- numeric(n_boot)
  
  #Initialising vectors for CER check
  MAD_L1_A <- numeric(n_boot)
  MAD_L2_A <- numeric(n_boot)
  MAD_L1_B <- numeric(n_boot)
  MAD_L2_B <- numeric(n_boot)
  
  cat("Generating topological distribution (", n_boot, " Iterations)...\n", sep="")
  
  for(i in 1:n_boot) {
    idx <- sample(N_samples, size = sample_size, replace = FALSE)
    df_sub <- df_merged[idx, ]
    
    # Base extraction
    rf_A_t <- ranger(PHQ9_Score ~ Bio_Level, data = df_sub, num.threads = 1)
    rf_B_t <- ranger(Bio_Level ~ PHQ9_Score, data = df_sub, num.threads = 1)
    
    res_A_t <- df_sub$PHQ9_Score - rf_A_t$predictions
    res_B_t <- df_sub$Bio_Level - rf_B_t$predictions
    
    #Data frame for layer 2
    df_res_t <- data.frame(res_A = res_A_t, res_B = res_B_t)
    
    # Cross prediction
    NN_A_t <- res_A_t - ranger(res_A ~ res_B, data = df_res_t, num.threads = 1)$predictions
    NN_B_t <- res_B_t - ranger(res_B ~ res_A, data = df_res_t, num.threads = 1)$predictions
    
    # Data frame for layer 3
    df_model <- cbind(df_res_t, df_sub)
    
    df_model_perm <- df_model
    df_model_perm$res_A <- sample(df_model_perm$res_A)
    df_model_perm$res_B <- sample(df_model_perm$res_B)
    
    NN_A_perm <- res_A_t - ranger(res_A ~ res_B, data = df_model_perm, num.threads = 1)$predictions
    NN_B_perm <- res_B_t - ranger(res_B ~ res_A, data = df_model_perm, num.threads = 1)$predictions
    
    NN_perm[i] <- mean(abs(NN_A_perm)) - mean(abs(NN_B_perm))
    
    #Asymmetry layer 1 vs 2
    MAD_L1_A[i] <- mean(abs(NN_A_perm))
    MAD_L2_A[i] <- mean(abs(NN_A_t))
    MAD_L1_B[i] <- mean(abs(NN_B_perm))
    MAD_L2_B[i] <- mean(abs(NN_B_t))
    
    #Feedback loop check
    NNN_A_true <- NN_A_t - ranger(NN_A_t ~ PHQ9_Score, data = df_model, num.threads = 1)$predictions
    NNN_B_true <- NN_B_t - ranger(NN_B_t ~ Bio_Level, data = df_model, num.threads = 1)$predictions
    
    NNN_true[i] <- mean(abs(NNN_A_true)) - mean(abs(NNN_B_true))
    
    #Permutation layer 3
    df_model_null <- df_model
    df_model_null$Bio_Level <- sample(df_model_null$Bio_Level)
    df_model_null$PHQ9_Score <- sample(df_model_null$PHQ9_Score)
    
    NNN_A_null <- NN_A_t - ranger(NN_A_t ~ PHQ9_Score, data = df_model_null, num.threads = 1)$predictions
    NNN_B_null <- NN_B_t - ranger(NN_B_t ~ Bio_Level, data = df_model_null, num.threads = 1)$predictions
    
    NNN_null[i] <- mean(abs(NNN_A_null)) - mean(abs(NNN_B_null))
  }
  
  #Confounder entropy reduction
  drop_A_vec <- (1 - (MAD_L2_A / MAD_L1_A)) * 100
  drop_B_vec <- (1 - (MAD_L2_B / MAD_L1_B)) * 100
  
  #Median and confidence interval
  ci_drop_A <- quantile(drop_A_vec, probs = c(0.025, 0.5, 0.975))
  ci_drop_B <- quantile(drop_B_vec, probs = c(0.025, 0.5, 0.975))
  
  cat("\n    [CER-Check] Noise reduction in depression: ", round(ci_drop_A[2], 2), 
      "%  (95% CI: [", round(ci_drop_A[1], 2), " , ", round(ci_drop_A[3], 2), "])")
  
  cat("\n    [CER-Check] Noise reduction in biomarker:  ", round(ci_drop_B[2], 2), 
      "%  (95% CI: [", round(ci_drop_B[1], 2), " , ", round(ci_drop_B[3], 2), "])\n")
  n_resamples <- 10000
  relative_ratio_vec <- numeric(n_resamples)
  
  for (i in 1:n_resamples) {
    idx <- sample(1:n_boot, size = n_boot, replace = TRUE)
    
    #MAD
    mad_true_total <- mean(abs(NNN_true[idx])) 
    mad_null_total <- mean(abs(NNN_null[idx]))
    
    # Ratio
    relative_ratio_vec[i] <- mad_true_total / mad_null_total
  }
  
  #Inference criteria
  ci_relative <- quantile(relative_ratio_vec, probs = c(0.025, 0.975))
  
  cat("\n    95% CI of the ratio (True / Null): [", round(ci_relative[1], 4), ",", round(ci_relative[2], 4), "]\n\n")
  
  if (ci_relative[2] < 1) {
    cat("LAYER 3 Conclusion: Vicious cycle!\n")
    
  } else if (ci_relative[1] > 1) {
    cat("LAYER 3 Conclusion: Unidirectional + ADVERSARIAL NOISE!\n")
    
  } else {
    cat("LAYER 3 Conclusion: Unidirectional + Symmetric Noise!\n")
  }
  
  rm(df_raw, df_bio, df_merged); gc()
}

################################################################################
#Running the main loop
################################################################################

for (task in biomarker_tasks) {
  run_rethinc(df_dep = df_dep, task = task, n_iter = 1000)
}

cat("\n=======================================================================\n")
cat("Pipeline Completed\n")