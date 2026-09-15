###############################################################################
## 28_collapsed_validation.R
##
## Small DGP-based validation of the collapsed (Rao-Blackwellized) sampler.
## Generates synthetic data with known truth K+_true = 3 and fits multiple
## chains, then reports parameter recovery and chain-spread diagnostics.
##
## Designed to run in ~5-15 minutes on a laptop.
###############################################################################

source(file.path("R", "01_lib_core.R"))
source(file.path("R", "simulation", "lib_simulation.R"))

set.seed(20260425)

cat("=== Generating synthetic data (K_true = 3, denser signal) ===\n")
sim <- generate_sim_data(n = 10, p = 200, T_max = 3, K_true = 3,
                         s_k = 100, kappa = 1, alpha_dp = 1,
                         w_N_mean = 0.05, sigma_k = 0.5, sigma_N = sqrt(2),
                         model = "M1", seed = 1)
cat("  n = 10, p = 200, T = 3, K_true = 3, N_obs =", ncol(sim$Y_g), "\n")
cat("  true pi:", round(sim$pi_true, 3), "\n")
cat("  true w_N range:", round(range(sim$wN_true), 4), "\n\n")

cat("=== Fit: collapsed sampler, 2 chains, K=4 ===\n")
fits <- list()
for (i in 1:2) {
  cat("  Chain", i, "...\n")
  fits[[i]] <- fit_dp_concordance(
    Y_g = sim$Y_g, Y_cf = sim$Y_cf,
    patient = sim$patient, time = sim$time,
    K = 4, n_iter = 3000, n_burn = 2000, thin = 5,
    model = "M1", seed = 100 + i, verbose = FALSE,
    anneal = TRUE, T_anneal = 10.0, n_cool_buffer = 500,
    store_pointwise_ll = FALSE
  )
}

ll_means <- sapply(fits, function(f) mean(f$samples$loglik))
cat("\n  Mean post-burn-in log-likelihood per chain:", round(ll_means, 1), "\n")
cat("  Chain spread (max - min):", round(diff(range(ll_means)), 1), "log-units\n")

for (i in 1:2) {
  cat(sprintf("\n  Chain %d:\n", i))
  cat(sprintf("    K+ posterior:  %s\n",
              paste(names(table(fits[[i]]$samples$K_plus)),
                    "=", table(fits[[i]]$samples$K_plus), collapse = "  ")))
  cat(sprintf("    sigma_N mean:   %.3f  (true %.3f)\n",
              mean(fits[[i]]$samples$sigma_0), sqrt(2)))
  cat(sprintf("    sigma_k mean per k: %s\n",
              paste(round(colMeans(fits[[i]]$samples$sigma_k), 2), collapse = " ")))
  cat(sprintf("    wN mean per patient: %s\n",
              paste(round(colMeans(fits[[i]]$samples$omega_0), 4), collapse = " ")))
}

# Check that K+ recovers truth
K_plus_modes <- sapply(fits, function(f) {
  tbl <- table(f$samples$K_plus)
  as.integer(names(tbl)[which.max(tbl)])
})
cat("\n  Posterior modal K+ per chain:", K_plus_modes,
    "  (true K+ = ", 3, ")\n")

cat("\n=== VALIDATION DONE ===\n")
