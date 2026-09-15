###############################################################################
## 04_simulation_poc.R
## Proof-of-concept simulation: verify the MCMC sampler recovers known
## parameters from synthetic data (small scale)
###############################################################################

set.seed(20260328)

base_dir <- getwd()
source(file.path(base_dir, "R", "01_lib_core.R"))
source(file.path(base_dir, "R", "simulation", "lib_simulation.R"))

cat("=== Proof-of-Concept Simulation ===\n")
cat("Small-scale test: n=10, p=100, T=3, K_true=3\n\n")

# ------------------------------------------------------------------
# 1. Generate synthetic data under M1 (tracking)
# ------------------------------------------------------------------
sim <- generate_sim_data(
  n = 10, p = 100, T_max = 3, K_true = 3,
  s_k = 10, kappa = 20, alpha_dp = 1,
  w_N_mean = 0.1, model = "M1", seed = 123
)

cat("True kappa:", sim$kappa_true, "\n")
cat("True pi:", round(sim$pi_true, 3), "\n")
cat("True wN (first 5):", round(sim$wN_true[1:5], 3), "\n")
cat("True K+:", sim$K_true, "\n\n")

# ------------------------------------------------------------------
# 2. Fit M1
# ------------------------------------------------------------------
cat("--- Fitting M1 ---\n")
fit1 <- fit_dp_concordance(
  Y_g = sim$Y_g, Y_cf = sim$Y_cf,
  patient = sim$patient, time = sim$time,
  K = 10, n_iter = 2000, n_burn = 1000,
  model = "M1", seed = 42, thin = 1, verbose = TRUE,
  kappa_method = "grid"
)

# ------------------------------------------------------------------
# 3. Fit M0
# ------------------------------------------------------------------
cat("\n--- Fitting M0 ---\n")
fit0 <- fit_dp_concordance(
  Y_g = sim$Y_g, Y_cf = sim$Y_cf,
  patient = sim$patient, time = sim$time,
  K = 10, n_iter = 2000, n_burn = 1000,
  model = "M0", seed = 43, thin = 1, verbose = TRUE
)

# ------------------------------------------------------------------
# 4. Check parameter recovery
# ------------------------------------------------------------------
cat("\n=== Parameter Recovery ===\n\n")

# kappa
cat(sprintf("kappa: true=%.1f, post_mean=%.1f, post_median=%.1f, 95%% CI=(%.1f, %.1f)\n",
            sim$kappa_true,
            mean(fit1$samples$kappa),
            median(fit1$samples$kappa),
            quantile(fit1$samples$kappa, 0.025),
            quantile(fit1$samples$kappa, 0.975)))

# alpha
cat(sprintf("alpha: true=%.1f, post_mean=%.2f, 95%% CI=(%.2f, %.2f)\n",
            sim$alpha_true,
            mean(fit1$samples$alpha_dp),
            quantile(fit1$samples$alpha_dp, 0.025),
            quantile(fit1$samples$alpha_dp, 0.975)))

# K+
cat(sprintf("K+: true=%d, post_mean=%.1f, post_median=%d\n",
            sim$K_true,
            mean(fit1$samples$K_plus),
            median(fit1$samples$K_plus)))

# Contamination
# Per-OBSERVATION contamination against the per-observation truth. wN_true is
# the per-patient mean (length n) and cannot be differenced against a length
# N_obs posterior summary.
wN_post <- colMeans(fit1$samples$omega_0)
w0_true <- as.numeric(sim$omega_0_true)
stopifnot(length(wN_post) == length(w0_true))
cat(sprintf("omega_0 bias: %.4f, RMSE: %.4f, cor: %.3f\n",
            mean(wN_post - w0_true),
            sqrt(mean((wN_post - w0_true)^2)),
            suppressWarnings(stats::cor(wN_post, w0_true))))

# Tracking correlation
alpha_m <- mean(fit1$samples$alpha_dp)
kappa_m <- mean(fit1$samples$kappa)
corr_post <- (kappa_m + 1) / (alpha_m + kappa_m + 1)
corr_true <- (sim$kappa_true + 1) / (sim$alpha_true + sim$kappa_true + 1)
cat(sprintf("Tracking correlation: true=%.3f, estimated=%.3f\n",
            corr_true, corr_post))

# Model comparison via WAIC
cat("\nModel comparison: use WAIC/LOO via R/08_model_comparison.R for formal comparison.\n")

cat("\n=== POC Complete ===\n")
cat("Sampler is functional. Proceed with full analysis.\n")

# Save POC results
poc_results <- list(
  sim = sim,
  fit1 = fit1, fit0 = fit0,
  metrics = compute_sim_metrics(fit1, sim)
)
saveRDS(poc_results, file.path(base_dir, "results", "poc_results.rds"))
