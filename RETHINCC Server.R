################################################################################
# RE:THINC - COMPLETE CAUSAL DISCOVERY PIPELINE (SERVER VERSION)
################################################################################

################################################################################
#Libraries and global settings
################################################################################
library(parallel)
library(doParallel)
library(foreach)
library(ranger)
library(dplyr)
library(NHANES)
library(ggplot2)

set.seed(1)
n_iter_global <- 1000
n_cores <- max(1, detectCores() - 1)
cat("Server initialized. Running on", n_cores, "cores with set.seed(1) and", n_iter_global, "iterations.\n")

pooled_sd <- function(x, y) {
  nx <- length(x); ny <- length(y)
  sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
}


################################################################################
# PART 3: NHANES BMI vs. Blood Pressure (With and Without Age)
################################################################################

df_med <- NHANES %>%
  select(BMI, BPSysAve, Age) %>%
  filter(complete.cases(.)) %>%
  mutate(across(everything(), ~ as.numeric(scale(.))))

cat("Loaded df_med (N =", nrow(df_med), ")\n")
### This part is unimportant for the power analysis
run_bmi_bp_test <- function(df, control_age = FALSE) {
  if (control_age) {
    cat("\n>>> RUNNING VERSION 2: BMI vs BPSysAve (CONTROLLED FOR AGE) <<<\n")
    form_A <- BPSysAve ~ BMI + Age
    form_B <- BMI ~ BPSysAve + Age
  } else {
    cat("\n>>> RUNNING VERSION 1: BMI vs BPSysAve (UNADJUSTED) <<<\n")
    form_A <- BPSysAve ~ BMI
    form_B <- BMI ~ BPSysAve
  }
  
  sample_size <- floor(0.8 * nrow(df))
  
  # Layer 1
  cat("Start Layer 1...\n")
  boot_diffs <- mclapply(1:n_iter_global, function(i) {
    df_sub <- df[sample(nrow(df), size = sample_size, replace = FALSE), ]
    rf_A <- ranger(form_A, data = df_sub, num.threads = 1)
    rf_B <- ranger(form_B, data = df_sub, num.threads = 1)
    res_A <- df_sub$BPSysAve - rf_A$predictions
    res_B <- df_sub$BMI - rf_B$predictions
    df_res <- data.frame(res_A = res_A, res_B = res_B)
    NN_A <- res_A - ranger(res_A ~ res_B, data = df_res, num.threads = 1)$predictions
    NN_B <- res_B - ranger(res_B ~ res_A, data = df_res, num.threads = 1)$predictions
    mean(abs(NN_A)) - mean(abs(NN_B))
  }, mc.cores = n_cores) |> unlist()
  
  cat("L1 95% CI: [", quantile(boot_diffs, 0.025), ",", quantile(boot_diffs, 0.975), "]\n")
  
  # Layer 2 & 3
  cat("Start Layer 2 & 3...\n")
  res_l23 <- mclapply(1:100, function(i) {
    df_sub <- df[sample(nrow(df), size = sample_size, replace = FALSE), ]
    rf_A_t <- ranger(form_A, data = df_sub, num.threads = 1)
    rf_B_t <- ranger(form_B, data = df_sub, num.threads = 1)
    res_A_t <- df_sub$BPSysAve - rf_A_t$predictions
    res_B_t <- df_sub$BMI - rf_B_t$predictions
    df_res_t <- data.frame(res_A = res_A_t, res_B = res_B_t)
    NN_A_t <- res_A_t - ranger(res_A ~ res_B, data = df_res_t, num.threads = 1)$predictions
    NN_B_t <- res_B_t - ranger(res_B ~ res_A, data = df_res_t, num.threads = 1)$predictions
    
    # Layer 2 (Permute residuals)
    df_model_perm <- cbind(df_res_t, df_sub)
    df_model_perm$res_A <- sample(df_model_perm$res_A)
    df_model_perm$res_B <- sample(df_model_perm$res_B)
    NN_A_perm <- res_A_t - ranger(res_A ~ res_B, data = df_model_perm, num.threads = 1)$predictions
    NN_B_perm <- res_B_t - ranger(res_B ~ res_A, data = df_model_perm, num.threads = 1)$predictions
    
    # Layer 3 (True loop)
    form_NN_A <- if(control_age) {NN_A_t ~ BPSysAve + Age} else {NN_A_t ~ BPSysAve}
    form_NN_B <- if(control_age) {NN_B_t ~ BMI + Age} else {NN_B_t ~ BMI}
    
    NNN_A_true <- NN_A_t - ranger(form_NN_A, data = cbind(NN_A_t = NN_A_t, df_sub), num.threads = 1)$predictions
    NNN_B_true <- NN_B_t - ranger(form_NN_B, data = cbind(NN_B_t = NN_B_t, df_sub), num.threads = 1)$predictions
    
    # Layer 3 (Null loop) - SHUFFLE ONLY TARGET VARS, KEEP AGE INTACT
    df_model_null <- df_sub
    df_model_null$BMI <- sample(df_model_null$BMI)
    df_model_null$BPSysAve <- sample(df_model_null$BPSysAve)
    
    NNN_A_null <- NN_A_t - ranger(form_NN_A, data = cbind(NN_A_t = NN_A_t, df_model_null), num.threads = 1)$predictions
    NNN_B_null <- NN_B_t - ranger(form_NN_B, data = cbind(NN_B_t = NN_B_t, df_model_null), num.threads = 1)$predictions
    
    c(MAD_L1_A = mean(abs(NN_A_perm)), MAD_L2_A = mean(abs(NN_A_t)),
      MAD_L1_B = mean(abs(NN_B_perm)), MAD_L2_B = mean(abs(NN_B_t)),
      NNN_true = mean(abs(NNN_A_true)) - mean(abs(NNN_B_true)),
      NNN_null = mean(abs(NNN_A_null)) - mean(abs(NNN_B_null)))
  }, mc.cores = n_cores)
  
  res_mat <- do.call(rbind, res_l23)
  
  ci_drop_A <- quantile((1 - (res_mat[, "MAD_L2_A"] / res_mat[, "MAD_L1_A"])) * 100, probs = c(0.025, 0.5, 0.975))
  ci_drop_B <- quantile((1 - (res_mat[, "MAD_L2_B"] / res_mat[, "MAD_L1_B"])) * 100, probs = c(0.025, 0.5, 0.975))
  cat("L2 CER BPSysAve: ", round(ci_drop_A[2], 2), "%\n")
  cat("L2 CER BMI:      ", round(ci_drop_B[2], 2), "%\n")
  
  ratio_vec <- numeric(10000)
  for(i in 1:10000) {
    idx <- sample(1:nrow(res_mat), replace = TRUE)
    ratio_vec[i] <- mean(abs(res_mat[idx, "NNN_true"])) / mean(abs(res_mat[idx, "NNN_null"]))
  }
  ci_rel <- quantile(ratio_vec, probs = c(0.025, 0.975))
  cat("L3 Ratio CI: [", round(ci_rel[1], 4), ",", round(ci_rel[2], 4), "]\n\n")
}

run_bmi_bp_test(df_med, control_age = FALSE)
run_bmi_bp_test(df_med, control_age = TRUE)


################################################################################
# PART 4: Power Analysis (Monte Carlo Simulation)
################################################################################

n_simulations <- 1000   
n_iter_power <- 100 
sample_sizes <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200,
                  1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200,
                  2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200,
                  3300, 3400, 3500)
effect_sizes <- c(0.20, 0.50, 0.80) 

power_results <- data.frame()

# Parallelisation setup for foreach
cl <- makeCluster(n_cores)
registerDoParallel(cl)

rf_base_model <- ranger(BPSysAve ~ BMI, data = df_med)
true_signal_std <- as.numeric(scale(rf_base_model$predictions)) 
true_oob_residuals <- as.numeric(scale(df_med$BPSysAve - rf_base_model$predictions))
real_X_std <- as.numeric(scale(df_med$BMI))

for(es in effect_sizes) {
  for(n in sample_sizes) {
    cat("Running Simulation for N =", n, "| Effect Size =", es, "...\n")
    
    sim_results <- foreach(sim = 1:n_simulations, .combine = rbind, .packages = c("ranger", "dplyr")) %dopar% {
      sim_indices <- sample(1:length(real_X_std), size = n, replace = TRUE)
      
      X_sim <- real_X_std[sim_indices]
      Signal_sim <- true_signal_std[sim_indices]
      Noise_sim <- true_oob_residuals[sim_indices]
      
      Y_sim <- (es * Signal_sim) + sqrt(1 - es^2) * Noise_sim
      
      df_sim <- data.frame(X = X_sim, Y = Y_sim) %>% 
        mutate(across(everything(), ~ as.numeric(scale(.))))
      
      boot_diffs_sim <- numeric(n_iter_power)
      sub_size <- floor(0.8 * nrow(df_sim)) 
      
      for(i in 1:n_iter_power) {
        df_sub <- df_sim[sample(nrow(df_sim), size = sub_size, replace = FALSE), ]
        rf_boot_A <- ranger(Y ~ X, data = df_sub, num.trees = 150, num.threads = 1)
        rf_boot_B <- ranger(X ~ Y, data = df_sub, num.trees = 150, num.threads = 1)
        res_A <- df_sub$Y - rf_boot_A$predictions
        res_B <- df_sub$X - rf_boot_B$predictions
        df_res <- data.frame(res_A = res_A, res_B = res_B)
        NN_A <- res_A - ranger(res_A ~ res_B, data = df_res, num.threads = 1)$predictions
        NN_B <- res_B - ranger(res_B ~ res_A, data = df_res, num.threads = 1)$predictions
        boot_diffs_sim[i] <- mean(abs(NN_A)) - mean(abs(NN_B))
      }
      
      boot_lower <- quantile(boot_diffs_sim, 0.025)
      boot_upper <- quantile(boot_diffs_sim, 0.975)
      
      if(boot_upper < 0) {
        decision <- "X->Y"
      } else if(boot_lower > 0) {
        decision <- "Y->X"
      } else {
        decision <- "Undecided"
      }
      
      c(Correct = ifelse(decision == "X->Y", 1, 0),
        Wrong = ifelse(decision == "Y->X", 1, 0),
        Undecided = ifelse(decision == "Undecided", 1, 0))
    }
    
    power <- sum(sim_results[, "Correct"]) / n_simulations
    false_rate <- sum(sim_results[, "Wrong"]) / n_simulations
    undecided_rate <- sum(sim_results[, "Undecided"]) / n_simulations
    
    power_results <- rbind(power_results, data.frame(
      N = n, EffectSize = es, Power = power,
      FalseRate = false_rate, UndecidedRate = undecided_rate
    ))
  }
}

stopCluster(cl)

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
  theme(legend.position = "bottom", legend.title = element_blank(), panel.grid.minor = element_blank()) +
  labs(title = "Statistical Power of the RE:THINCC Algorithm",
       subtitle = "Semi-Synthetic Plasmode Simulation (True Causal Recovery: BMI \u2192 Blood Pressure)",
       x = "Sample Size (N)", y = "Statistical Power (Correct Identification Rate)")

ggsave("plot_power_results.pdf", plot = p_power, width = 10, height = 7)
write.csv(power_results, "power_results.csv", row.names = FALSE)
saveRDS(power_results, "power_results.rds")