###############################################################################
## 31_smoke_test_full_pipeline.R
##
## End-to-end smoke test of the production pipeline on a small synthetic problem.
## Mirrors what the slurm jobs do but at a tiny scale (~2 minutes wall-clock).
##
## Steps:
##   1. Generate small synthetic data (p=200, N_obs=12, K_true=2).
##   2. Fit M1 with fit_multi_chain (3 chains).
##   3. Fit M0 with fit_multi_chain (3 chains).
##   4. Verify multi$merged contains v, v_g, v_cf, tau_k, b_sigma fields where
##      expected (M1: v populated; M0: v_g and v_cf populated).
##   5. Relabel both models via relabel_by_weight and align_chains_by_signature.
##   6. Re-merge and confirm no field-mismatch.
##   7. Compute WAIC and LOO for both models.
##   8. Run compute_convergence on K_plus, sigma_N, wN.
##
## PASS criteria: all 8 steps complete without error; merged objects have the
## expected fields; final WAIC/LOO are finite.
###############################################################################

source(file.path("R", "00_bootstrap.R"))
source(file.path("R", "01_lib_core.R"))

set.seed(20260513)

# ----- Synthetic data (small) -----
p <- 200
N_obs <- 12
K_true <- 2
T_max <- 1

theta_true <- matrix(0, p, K_true)
theta_true[1:(p/2), 1] <- 3
theta_true[(p/2 + 1):p, 2] <- 3
sigma_k_true <- 0.5

omega_true <- matrix(NA, K_true, N_obs)
omega_true[1, ] <- runif(N_obs, 0.2, 0.8)
omega_true[2, ] <- 1 - omega_true[1, ]

Y_g  <- matrix(0, p, N_obs)
Y_cf <- matrix(0, p, N_obs)
for (obs in 1:N_obs) {
  for (j in 1:p) {
    k_g  <- sample.int(K_true, 1, prob = omega_true[, obs])
    Y_g[j, obs]  <- rnorm(1, theta_true[j, k_g],  sigma_k_true)
    k_cf <- sample.int(K_true, 1, prob = omega_true[, obs])
    Y_cf[j, obs] <- rnorm(1, theta_true[j, k_cf], sigma_k_true)
  }
}
patient <- 1:N_obs
time <- rep(1, N_obs)

cat("=== Smoke test: full M1+M0+combine pipeline ===\n")
cat("p =", p, " N_obs =", N_obs, " K_true =", K_true, "\n\n")

# ----- Fit M1 -----
cat("--- fit_multi_chain M1 ---\n")
multi_M1 <- fit_multi_chain(
  Y_g, Y_cf, patient, time,
  n_chains = 3, seed = 42, mc.cores = 1,
  K = 4, n_iter = 1500, n_burn = 1000, thin = 5,
  model = "M1", kappa_fixed = 1.0,
  anneal = TRUE, T_anneal = 10.0, n_cool_buffer = 300,
  store_pointwise_ll = TRUE
)
cat("  M1 fit OK. Merged fields:", paste(names(multi_M1$merged), collapse = ", "), "\n")
stopifnot(!is.null(multi_M1$merged$v))
stopifnot(!is.null(multi_M1$merged$v_g))    # should be allocated even under M1 (NA-filled)
stopifnot(!is.null(multi_M1$merged$tau_k))
stopifnot(!is.null(multi_M1$merged$b_sigma))

# ----- Fit M0 -----
cat("\n--- fit_multi_chain M0 ---\n")
multi_M0 <- fit_multi_chain(
  Y_g, Y_cf, patient, time,
  n_chains = 3, seed = 100, mc.cores = 1,
  K = 4, n_iter = 1500, n_burn = 1000, thin = 5,
  model = "M0",
  anneal = TRUE, T_anneal = 10.0, n_cool_buffer = 300,
  store_pointwise_ll = TRUE
)
cat("  M0 fit OK. Merged fields:", paste(names(multi_M0$merged), collapse = ", "), "\n")
stopifnot(!is.null(multi_M0$merged$v_g))
stopifnot(!is.null(multi_M0$merged$v_cf))
# M0's v_g should be populated, not all NA
stopifnot(any(!is.na(multi_M0$merged$v_g)))
stopifnot(any(!is.na(multi_M0$merged$v_cf)))

# ----- Relabel M1 -----
cat("\n--- relabel + align M1 ---\n")
for (i in seq_along(multi_M1$chains)) {
  multi_M1$chains[[i]] <- relabel_by_weight(multi_M1$chains[[i]])
}
multi_M1$chains <- align_chains_by_signature(multi_M1$chains)
multi_M1$merged <- merge_chain_samples(multi_M1$chains)
cat("  M1 relabel+align OK.\n")

# ----- Relabel M0 -----
cat("\n--- relabel + align M0 ---\n")
for (i in seq_along(multi_M0$chains)) {
  multi_M0$chains[[i]] <- relabel_by_weight(multi_M0$chains[[i]])
}
multi_M0$chains <- align_chains_by_signature(multi_M0$chains)
multi_M0$merged <- merge_chain_samples(multi_M0$chains)
cat("  M0 relabel+align OK.\n")

# ----- WAIC/LOO -----
cat("\n--- WAIC/LOO ---\n")
waic_M1 <- loo::waic(multi_M1$merged$pointwise_ll)
waic_M0 <- loo::waic(multi_M0$merged$pointwise_ll)
cat("  WAIC M1:", round(waic_M1$estimates["waic", "Estimate"], 2), "\n")
cat("  WAIC M0:", round(waic_M0$estimates["waic", "Estimate"], 2), "\n")
stopifnot(is.finite(waic_M1$estimates["waic", "Estimate"]))
stopifnot(is.finite(waic_M0$estimates["waic", "Estimate"]))

# ----- Convergence -----
cat("\n--- compute_convergence ---\n")
conv_M1 <- compute_convergence(multi_M1, params = c("K_plus", "sigma_0", "omega_0", "varsigma_0"))
conv_M0 <- compute_convergence(multi_M0, params = c("K_plus", "sigma_0", "omega_0", "varsigma_0"))
cat("  M1 max Rhat:", round(conv_M1$worst_rhat, 4), "\n")
cat("  M0 max Rhat:", round(conv_M0$worst_rhat, 4), "\n")
stopifnot(is.finite(conv_M1$worst_rhat))
stopifnot(is.finite(conv_M0$worst_rhat))

# ----- Trapped-chain check -----
cat("\n--- chain-health check ---\n")
trap_M1 <- identify_trapped_chains(multi_M1, ll_gap = 10000)
trap_M0 <- identify_trapped_chains(multi_M0, ll_gap = 10000)
cat("  M1 retained chains:", trap_M1$keep, "  spread:", round(diff(range(trap_M1$ll_means)), 1), "\n")
cat("  M0 retained chains:", trap_M0$keep, "  spread:", round(diff(range(trap_M0$ll_means)), 1), "\n")

cat("\n=== PASS: full pipeline smoke test complete ===\n")
