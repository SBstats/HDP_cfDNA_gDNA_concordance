###############################################################################
## 30_m0_prior_validation.R
##
## Validate the new M0 specification: independent source-specific DPs.
##   omega^g_{it}  ~ Dir(kappa * pi^g),  pi^g  ~ GEM(alpha)
##   omega^cf_{it} ~ Dir(kappa * pi^cf), pi^cf ~ GEM(alpha)
## with pi^g and pi^cf independent.
##
## Goals:
##   (1) Sanity-check that the new M0 sampler runs cleanly end-to-end on a
##       small synthetic problem (no crashes, all chains complete, finite LL).
##   (2) Check between-chain log-likelihood spread — the previous Dir(1,...,1)
##       M0 had spreads in the thousands; with source-specific pi the
##       Polya-urn predictive is anchored across observations within a source,
##       so spreads should be substantially smaller.
##   (3) Confirm the model recovers a clear two-subclone truth under M0:
##       K+ posterior near 2, sigma_k near the true within-cluster SD.
##   (4) Inspect that pi_g and pi_cf are independently estimated (they may
##       legitimately differ across chains by label permutation, which the
##       relabel-by-weight machinery handles).
##
## Setup: small p=300, N_obs=16, K_true=2, well-separated signatures so a
## working sampler should converge in a few thousand iterations.
###############################################################################

source(file.path("R", "01_lib_core.R"))

set.seed(20260513)

# Synthetic two-subclone problem
p <- 300
N_obs <- 16
K_true <- 2
T_max <- 1

theta_true <- matrix(0, p, K_true)
theta_true[1:(p/2), 1] <- 3
theta_true[(p/2 + 1):p, 2] <- 3
sigma_k_true <- 0.6

# Per-observation mixtures that genuinely differ between sources (the
# source-DP M0 should accommodate this, but the data also won't strongly
# discriminate against M1 because the means happen to be the same)
true_omega_g  <- matrix(NA, K_true, N_obs)
true_omega_cf <- matrix(NA, K_true, N_obs)
true_omega_g[1, ]  <- runif(N_obs, 0.2, 0.8)
true_omega_g[2, ]  <- 1 - true_omega_g[1, ]
true_omega_cf[1, ] <- runif(N_obs, 0.2, 0.8)
true_omega_cf[2, ] <- 1 - true_omega_cf[1, ]

Y_g  <- matrix(0, p, N_obs)
Y_cf <- matrix(0, p, N_obs)
for (obs in 1:N_obs) {
  for (j in 1:p) {
    k_g  <- sample.int(K_true, 1, prob = true_omega_g[, obs])
    Y_g[j, obs]  <- rnorm(1, theta_true[j, k_g], sigma_k_true)
    k_cf <- sample.int(K_true, 1, prob = true_omega_cf[, obs])
    Y_cf[j, obs] <- rnorm(1, theta_true[j, k_cf], sigma_k_true)
  }
}
patient <- 1:N_obs
time <- rep(1, N_obs)

cat("=== Source-DP M0 validation ===\n")
cat("p =", p, " N_obs =", N_obs, " K_true =", K_true, "\n")
cat("Y_g sample sd:", round(sd(Y_g), 2), "\n\n")

n_chains <- 3
fits <- list()
for (i in 1:n_chains) {
  cat(sprintf("--- M0 chain %d ---\n", i))
  fits[[i]] <- fit_dp_concordance(
    Y_g = Y_g, Y_cf = Y_cf,
    patient = patient, time = time,
    K = 4, n_iter = 2500, n_burn = 1500, thin = 5,
    model = "M0", seed = 300 + i, verbose = FALSE,
    anneal = TRUE, T_anneal = 10.0, n_cool_buffer = 400
  )
}

ll_means <- sapply(fits, function(f) mean(f$samples$loglik))
ll_spread <- diff(range(ll_means))
k_plus_means <- sapply(fits, function(f) mean(f$samples$K_plus))

cat("\n--- Summary across chains ---\n")
cat("  Mean log-likelihood per chain:", paste(round(ll_means, 1), collapse = "  "), "\n")
cat("  LL spread (max - min):", round(ll_spread, 1), "log-units\n")
cat("  K+ posterior mean per chain:", paste(round(k_plus_means, 2), collapse = "  "), "\n")

cat("\n--- Per-source stick-breaking (chain 1) ---\n")
v_g_means  <- colMeans(fits[[1]]$samples$v_g)
v_cf_means <- colMeans(fits[[1]]$samples$v_cf)
cat("  v_g posterior mean: ", paste(round(v_g_means, 3),  collapse = "  "), "\n")
cat("  v_cf posterior mean:", paste(round(v_cf_means, 3), collapse = "  "), "\n")

# Reconstruct pi_g, pi_cf at the posterior mean of v
pi_from_v <- function(v_mean) {
  v_full <- c(v_mean, 1)
  stick_break(v_full)
}
pi_g_post  <- pi_from_v(v_g_means)
pi_cf_post <- pi_from_v(v_cf_means)
cat("  pi_g posterior mean:", paste(round(pi_g_post, 3), collapse = "  "), "\n")
cat("  pi_cf posterior mean:", paste(round(pi_cf_post, 3), collapse = "  "), "\n")

# Apply relabel_by_weight to chain 1 and verify the trace works
cat("\n--- Relabel-by-weight (chain 1) ---\n")
fits[[1]] <- relabel_by_weight(fits[[1]])
cat("  v_g post-relabel mean:",  paste(round(colMeans(fits[[1]]$samples$v_g),  3), collapse = "  "), "\n")
cat("  v_cf post-relabel mean:", paste(round(colMeans(fits[[1]]$samples$v_cf), 3), collapse = "  "), "\n")

cat("\n=== Pass criteria ===\n")
cat("  (a) All chains finite LL:                 ",
    if (all(is.finite(ll_means))) "PASS" else "FAIL", "\n")
cat("  (b) LL spread across chains < 500 units:  ",
    if (ll_spread < 500) "PASS" else paste("FAIL (spread =", round(ll_spread, 0), ")"), "\n")
cat("  (c) Posterior K+ in [1.8, 2.5] (truth=2): ",
    if (all(k_plus_means >= 1.8 & k_plus_means <= 2.5)) "PASS" else "FAIL", "\n")
cat("  (d) Top-2 pi_g sum > 0.95 (sparsity):    ",
    if (sum(sort(pi_g_post,  decreasing = TRUE)[1:2]) > 0.95) "PASS" else "FAIL", "\n")
cat("  (e) Top-2 pi_cf sum > 0.95 (sparsity):   ",
    if (sum(sort(pi_cf_post, decreasing = TRUE)[1:2]) > 0.95) "PASS" else "FAIL", "\n")

cat("\n=== DONE ===\n")
