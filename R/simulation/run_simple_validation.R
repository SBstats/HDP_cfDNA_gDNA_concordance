###############################################################################
## 29_simple_validation.R
##
## Very simple validation: two well-separated subclones, larger p, lower noise.
## Just to confirm the collapsed sampler can recover a clear truth.
###############################################################################

source(file.path("R", "01_lib_core.R"))

set.seed(20260425)

# Two truly separated clusters: half the genes are "active" in subclone 1
# (theta = 3), the other half in subclone 2 (theta = -3). sigma_k = 0.5 (low
# within-cluster noise). 12 observations, each with one of the two
# compositions: 6 are mostly subclone 1, 6 are mostly subclone 2.
p <- 200
n_obs <- 12
N_obs <- n_obs
K_true <- 2
T_max <- 1

theta_true <- matrix(0, p, K_true)
theta_true[1:(p/2), 1] <- 3
theta_true[(p/2 + 1):p, 2] <- 3
sigma_k_true <- 0.5

# Each observation is a 90/10 or 10/90 mixture of subclones 1 and 2.
true_omega <- matrix(c(rep(c(0.9, 0.1), length.out = N_obs),
                       rep(c(0.1, 0.9), length.out = N_obs)),
                     nrow = K_true, byrow = TRUE)
true_omega[, seq(1, N_obs, 2)] <- c(0.9, 0.1)
true_omega[, seq(2, N_obs, 2)] <- c(0.1, 0.9)

Y_g  <- matrix(0, p, N_obs)
Y_cf <- matrix(0, p, N_obs)
for (obs in 1:N_obs) {
  # Sample z then Y for both sources.
  for (j in 1:p) {
    k_g  <- sample.int(K_true, 1, prob = true_omega[, obs])
    Y_g[j, obs]  <- rnorm(1, theta_true[j, k_g], sigma_k_true)
    # cfDNA: no contamination in this simple test.
    k_cf <- sample.int(K_true, 1, prob = true_omega[, obs])
    Y_cf[j, obs] <- rnorm(1, theta_true[j, k_cf], sigma_k_true)
  }
}
patient <- 1:N_obs
time <- rep(1, N_obs)

cat("=== Simple test: K_true=2, well-separated, p =", p, "===\n")
cat("Y_g sample sd:", round(sd(Y_g), 2), "\n")

cat("\n=== Fit: collapsed sampler, K=4, 2 chains, T_anneal=10 ===\n")
fits <- list()
for (i in 1:2) {
  cat("  Chain", i, "...\n")
  fits[[i]] <- fit_dp_concordance(
    Y_g = Y_g, Y_cf = Y_cf,
    patient = patient, time = time,
    K = 4, n_iter = 3000, n_burn = 2000, thin = 5,
    model = "M1", seed = 100 + i, verbose = FALSE,
    anneal = TRUE, T_anneal = 10.0, n_cool_buffer = 500
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
  cat(sprintf("    sigma_k mean per k: %s  (true %.2f)\n",
              paste(round(colMeans(fits[[i]]$samples$sigma_k), 2), collapse = " "),
              sigma_k_true))
  # omega mean per component
  K <- 4
  omega_means <- numeric(K)
  for (k in 1:K) {
    cols <- ((k - 1) * N_obs + 1):(k * N_obs)
    omega_means[k] <- mean(fits[[i]]$samples$omega_g_trace[, cols])
  }
  cat(sprintf("    omega_g mean per k: %s\n", paste(round(omega_means, 3), collapse = " ")))
}

cat("\n=== DONE ===\n")
