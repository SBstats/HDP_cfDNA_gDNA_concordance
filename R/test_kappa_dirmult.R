#!/usr/bin/env Rscript
###############################################################################
## test_kappa_dirmult.R
## Local validation: Dirichlet-Multinomial kappa update recovers true kappa
##
## Generates synthetic data from the model with known kappa_true, runs a
## short MCMC chain, and checks kappa posterior coverage.
###############################################################################

cat("=== Validation: Dirichlet-Multinomial Kappa Update ===\n\n")

base_dir <- getwd()
source(file.path(base_dir, "R", "01_lib_core.R"))

set.seed(42)

# --- Test parameters ---
K_true <- 3
p <- 100
n_subjects <- 5
T_max <- 2
kappa_true <- 5.0
alpha_dp <- 1.0

# Generate patient-time structure
patient <- rep(1:n_subjects, each = T_max)
time_idx <- rep(1:T_max, times = n_subjects)
N_obs <- length(patient)

cat(sprintf("Setup: K=%d, p=%d, n=%d, T=%d, N_obs=%d, kappa_true=%.1f\n",
            K_true, p, n_subjects, T_max, N_obs, kappa_true))

# --- Generate true parameters ---
# Stick-breaking weights
v_true <- rbeta(K_true, 1, alpha_dp)
v_true[K_true] <- 1
pi_true <- stick_break(v_true)

cat(sprintf("True pi: %s\n", paste(round(pi_true, 3), collapse=", ")))

# Subclonal signatures: sparse (only s_k active genes per component)
theta_true <- array(0, dim = c(p, K_true, T_max))
for (k in 1:K_true) {
  s_k <- max(3, round(p * 0.1))
  active_genes <- sample(p, s_k)
  for (tt in 1:T_max) {
    theta_true[active_genes, k, tt] <- rnorm(s_k, 0, 2)
  }
}

sigma_k2_true <- rep(0.5, K_true)
sigma_N2_true <- 2.0
wN_true <- rep(0.01, n_subjects)

# --- Generate data ---
Y_g <- matrix(0, p, N_obs)
Y_cf <- matrix(0, p, N_obs)

for (obs in 1:N_obs) {
  # Source-specific weights from Dir(kappa * pi)
  omega_g_obs <- rdirichlet_one(kappa_true * pi_true)
  omega_cf_obs <- rdirichlet_one(kappa_true * pi_true)

  for (j in 1:p) {
    # gDNA: mixture
    k_g <- sample(K_true, 1, prob = omega_g_obs)
    tt <- time_idx[obs]
    Y_g[j, obs] <- rnorm(1, theta_true[j, k_g, tt], sqrt(sigma_k2_true[k_g]))

    # cfDNA: contaminated mixture
    if (runif(1) < wN_true[patient[obs]]) {
      Y_cf[j, obs] <- rnorm(1, 0, sqrt(sigma_N2_true))
    } else {
      k_cf <- sample(K_true, 1, prob = omega_cf_obs)
      Y_cf[j, obs] <- rnorm(1, theta_true[j, k_cf, tt], sqrt(sigma_k2_true[k_cf]))
    }
  }
}

cat(sprintf("Data generated: Y_g [%d x %d], Y_cf [%d x %d]\n\n",
            nrow(Y_g), ncol(Y_g), nrow(Y_cf), ncol(Y_cf)))

# --- Run sampler with new DirMult kappa update ---
cat("--- Running sampler (DirMult kappa, K=5, 3000 iter, 1500 burn-in) ---\n")
fit <- fit_dp_concordance(
  Y_g, Y_cf, patient, time_idx,
  K = 5, n_iter = 3000, n_burn = 1500, thin = 1,
  model = "M1", kappa_method = "grid",
  seed = 123, verbose = TRUE,
  kappa_init = 10.0
)

kappa_samples <- fit$samples$kappa
kappa_mean <- mean(kappa_samples)
kappa_median <- median(kappa_samples)
kappa_ci <- quantile(kappa_samples, c(0.025, 0.975))
K_plus_samples <- fit$samples$K_plus

cat("\n=== RESULTS ===\n")
cat(sprintf("True kappa:      %.1f\n", kappa_true))
cat(sprintf("Posterior mean:   %.2f\n", kappa_mean))
cat(sprintf("Posterior median: %.2f\n", kappa_median))
cat(sprintf("95%% CI:          [%.2f, %.2f]\n", kappa_ci[1], kappa_ci[2]))
cat(sprintf("CI covers truth:  %s\n", ifelse(kappa_ci[1] <= kappa_true && kappa_true <= kappa_ci[2], "YES", "NO")))
cat(sprintf("K+ mode:         %d (true K=%d)\n", as.integer(names(which.max(table(K_plus_samples)))), K_true))

# Kappa grid values visited
kappa_tab <- table(kappa_samples)
cat("\nKappa grid usage:\n")
print(kappa_tab)

# --- Summary ---
cat("\n=== VALIDATION SUMMARY ===\n")
relative_error <- abs(kappa_median - kappa_true) / kappa_true
cat(sprintf("Relative error (median): %.1f%%\n", 100 * relative_error))

if (kappa_ci[1] <= kappa_true && kappa_true <= kappa_ci[2]) {
  cat("PASS: 95%% CI covers the true kappa\n")
} else {
  cat("WARN: 95%% CI does NOT cover the true kappa\n")
}

if (relative_error < 0.5) {
  cat("PASS: Relative error < 50%%\n")
} else {
  cat("WARN: Relative error >= 50%%\n")
}

cat("\nDone.\n")
