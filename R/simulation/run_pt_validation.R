###############################################################################
## 27_pt_validation.R
##
## Small DGP-based validation of the parallel-tempering sampler. Uses the
## simulation infrastructure in R/03_simulation_dgp.R to generate synthetic
## data with known truth, then fits both the standard sampler and the PT
## sampler and compares parameter recovery and chain mixing.
##
## Designed to run in ~5-15 minutes on a laptop. Not for production use.
###############################################################################

source(file.path("R", "01_lib_core.R"))
source(file.path("R", "simulation", "lib_simulation.R"))

set.seed(20260424)

# ---------------------------------------------------------------------------
# Generate synthetic data with known K+_true = 3
# ---------------------------------------------------------------------------
cat("=== Generating synthetic data (K_true = 3) ===\n")
sim <- generate_sim_data(n = 10, p = 200, T_max = 3, K_true = 3,
                         s_k = 20, kappa = 1, alpha_dp = 1,
                         w_N_mean = 0.05, sigma_k = 1, sigma_N = sqrt(2),
                         model = "M1", seed = 1)
cat("  n =", 10, ", p =", 200, ", T =", 3, ", K_true =", 3, "\n")
cat("  N_obs =", ncol(sim$Y_g), "\n")
cat("  true pi:", round(sim$pi_true, 3), "\n\n")

# ---------------------------------------------------------------------------
# Fit 1: standard sampler (current default with symmetric annealing)
# ---------------------------------------------------------------------------
cat("=== Fit 1: standard sampler, 2 chains, K=4 ===\n")
fits_std <- list()
for (i in 1:2) {
  cat("  Chain", i, "...\n")
  fits_std[[i]] <- fit_dp_concordance(
    Y_g = sim$Y_g, Y_cf = sim$Y_cf,
    patient = sim$patient, time = sim$time,
    K = 4, n_iter = 2000, n_burn = 1500, thin = 5,
    model = "M1", seed = 100 + i, verbose = FALSE,
    anneal = TRUE, T_anneal = 10.0, n_cool_buffer = 200,
    store_pointwise_ll = FALSE
  )
}
ll_std <- sapply(fits_std, function(f) mean(f$samples$loglik))
cat("  Mean post-burn log-lik per chain:", round(ll_std, 1), "\n")
cat("  Chain spread:", round(diff(range(ll_std)), 1), "log-units\n")

cat("  K+ posterior (chain 1):", paste(table(fits_std[[1]]$samples$K_plus), collapse = "/"), "\n")
cat("  K+ posterior (chain 2):", paste(table(fits_std[[2]]$samples$K_plus), collapse = "/"), "\n\n")

# ---------------------------------------------------------------------------
# Fit 2: parallel-tempered sampler (3 rungs)
# ---------------------------------------------------------------------------
cat("=== Fit 2: parallel-tempered sampler, 2 chains, K=4, 3 rungs ===\n")
fits_pt <- list()
for (i in 1:2) {
  cat("  Chain", i, "...\n")
  fits_pt[[i]] <- fit_parallel_tempered(
    Y_g = sim$Y_g, Y_cf = sim$Y_cf,
    patient = sim$patient, time = sim$time,
    temperatures = c(1, 1.5, 2.25),
    swap_every = 10,
    K = 4, n_iter = 2000, n_burn = 1500, thin = 10,
    model = "M1", seed = 100 + i, verbose = FALSE,
    store_pointwise_ll = FALSE
  )
}
ll_pt <- sapply(fits_pt, function(f) mean(f$samples$loglik))
cat("  Mean post-burn log-lik (T=1 rung) per chain:", round(ll_pt, 1), "\n")
cat("  Chain spread:", round(diff(range(ll_pt)), 1), "log-units\n")
for (i in 1:2) {
  cat("  Chain", i, "swap acceptance rates:",
      round(fits_pt[[i]]$pt_diagnostics$acceptance_rate, 3), "\n")
}
cat("  K+ posterior (chain 1):", paste(table(fits_pt[[1]]$samples$K_plus), collapse = "/"), "\n")
cat("  K+ posterior (chain 2):", paste(table(fits_pt[[2]]$samples$K_plus), collapse = "/"), "\n\n")

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
cat("=== Summary ===\n")
cat("Standard sampler chain spread:", round(diff(range(ll_std)), 1), "log-units\n")
cat("PT sampler       chain spread:", round(diff(range(ll_pt)),  1), "log-units\n")
if (diff(range(ll_std)) > 0) {
  cat("Reduction factor:", round(diff(range(ll_std)) / max(diff(range(ll_pt)), 1), 2),
      "x (smaller = PT helped)\n")
}
cat("\nVALIDATION DONE\n")
