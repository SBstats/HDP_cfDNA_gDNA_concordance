###############################################################################
## 01_mcmc_sampler.R
## MCMC sampler for the Bayesian Hierarchical DP Mixture Model
##
## Implements the 8-step blocked Gibbs sampler from Algorithm 1 of the
## manuscript, with the truncated stick-breaking representation (K components).
##
## Key design choices for Biometrics reproducibility:
##   - Pure R (no C++/Rcpp) for portability
##   - Vectorized operations for speed
##   - Explicit seed control
##   - All hyperparameters as function arguments
###############################################################################

suppressPackageStartupMessages({
  library(mvtnorm)
  library(coda)
})

# ===========================================================================
# MAIN SAMPLER FUNCTION
# ===========================================================================

#' Fit the DP Mixture Concordance Model (M1 or M0)
#'
#' @param Y_g     Matrix [p x N_obs]: normalized gDNA signals
#' @param Y_cf    Matrix [p x N_obs]: normalized cfDNA signals
#' @param patient Integer vector [N_obs]: patient index for each observation
#' @param time    Integer vector [N_obs]: time index for each observation
#' @param K       Integer: truncation level (default 30)
#' @param n_iter  Integer: total MCMC iterations
#' @param n_burn  Integer: burn-in iterations
#' @param model   Character: "M1" (tracking) or "M0" (non-tracking)
#' @param seed    Integer: random seed
#' @param hyperparams List of hyperparameters (see defaults)
#' @param thin    Integer: thinning interval
#' @param verbose Logical: print progress
#'
#' @return List with posterior samples and diagnostics
fit_dp_concordance <- function(Y_g, Y_cf, patient, time,
                                K = 15,  # smaller default for speed
                                n_iter = 5000,
                                n_burn = 2500,
                                model = "M1",
                                seed = 42,
                                hyperparams = list(),
                                thin = 1,
                                verbose = TRUE) {

  set.seed(seed)

  # --- Dimensions ---
  p <- nrow(Y_g)
  N_obs <- ncol(Y_g)
  n <- max(patient)
  T_max <- max(time)

  # --- Hyperparameters (with defaults) ---
  hp <- list(
    a_alpha = 1, b_alpha = 1,       # DP concentration prior
    a_kappa = 2, b_kappa = 0.1,     # tracking precision prior: E[kappa]=20
    a_N = 1, b_N = 9,               # contamination prior
    a_sigma = 2, b_sigma = 1,       # component variance prior
    a_sigmaN = 2, b_sigmaN = 1,     # background variance prior
    sigma_kappa_mh = 0.3            # MH proposal SD for log(kappa)
  )
  hp[names(hyperparams)] <- hyperparams

  # --- Storage for posterior samples ---
  n_save <- floor((n_iter - n_burn) / thin)
  samples <- list(
    kappa    = numeric(n_save),
    alpha_dp = numeric(n_save),
    K_plus   = integer(n_save),
    wN       = matrix(NA, n_save, n),
    sigma_k  = matrix(NA, n_save, K),
    sigma_N  = numeric(n_save),
    # Store v (stick-breaking) rather than full pi for memory
    v        = matrix(NA, n_save, K - 1),
    # Track log-likelihood for model comparison
    loglik   = numeric(n_save)
  )

  # --- Initialize parameters ---

  # Stick-breaking variables (uniform initialization)
  v <- rep(0.5, K)
  v[K] <- 1
  pi_vec <- stick_break(v)

  # Source-specific weights: uniform initialization [K x N_obs]
  omega_g  <- matrix(1/K, K, N_obs)
  omega_cf <- matrix(1/K, K, N_obs)

  # Component variances
  sigma_k2 <- rep(1.0, K)
  sigma_N2 <- 2.0

  # Contamination fractions
  wN <- rep(0.1, n)

  # Subclonal signatures: initialize from data with noise
  # For speed, use a random subset of genes for initialization
  theta <- matrix(0, p, K)  # [p x K] -- shared across time for now
  # Initialize first few components from data clusters
  if (N_obs >= K) {
    km <- tryCatch({
      # Use gene means across samples, cluster into K groups
      gene_means <- rowMeans(Y_g[, sample(N_obs, min(N_obs, 20))])
      kmeans(gene_means, centers = min(K, 5), nstart = 1)
    }, error = function(e) NULL)
  }

  # Horseshoe parameters
  tau_k2   <- rep(0.1, K)          # global shrinkage
  lambda2  <- matrix(1, p, K)      # local shrinkage [p x K]
  nu_aux   <- matrix(1, p, K)      # local auxiliary
  xi_aux   <- rep(1, K)            # global auxiliary

  # Concentration parameters
  alpha_dp <- 1.0
  kappa    <- 10.0

  # Allocation variables [N_obs x p] -- stored as integer matrices
  # For memory: only store counts, not full allocations
  # n_g[k, obs] = number of genes allocated to component k for gDNA obs
  # We'll compute these on the fly

  # MH acceptance tracking
  kappa_accept <- 0
  kappa_total  <- 0

  # --- Helper: compute pi from v ---
  # (defined below as standalone function)

  # =====================================================================
  # MCMC LOOP
  # =====================================================================
  if (verbose) cat("Starting MCMC:", n_iter, "iterations,",
                   K, "components,", p, "genes\n")

  for (iter in 1:n_iter) {

    # -------------------------------------------------------------------
    # STEP 1: Update latent allocations z_g, z_cf
    # -------------------------------------------------------------------
    # For each observation (i,t) and gene j, allocate to component k
    # For gDNA: P(z=k) ∝ omega_g[k,obs] * dnorm(Y, theta[j,k], sigma_k)
    # For cfDNA: P(z=0) ∝ wN[i] * dnorm(Y, 0, sigma_N)
    #            P(z=k) ∝ (1-wN[i]) * omega_cf[k,obs] * dnorm(Y, theta[j,k], sigma_k)

    # Compute allocation counts efficiently
    alloc <- compute_allocations(Y_g, Y_cf, theta, sigma_k2, sigma_N2,
                                 omega_g, omega_cf, wN, patient, K, p)
    n_g_counts  <- alloc$n_g   # [K x N_obs]
    n_cf_counts <- alloc$n_cf  # [K x N_obs]
    N0_vec      <- alloc$N0    # [n] normal component counts per patient
    Nplus_vec   <- alloc$Nplus # [n] tumor component counts per patient

    # -------------------------------------------------------------------
    # STEP 2: Update subclonal signatures theta[j,k]
    # -------------------------------------------------------------------
    theta <- update_theta(Y_g, Y_cf, alloc$z_g, alloc$z_cf,
                          sigma_k2, tau_k2, lambda2, K, p, N_obs)

    # -------------------------------------------------------------------
    # STEP 3: Update stick-breaking variables v_k (M1 only)
    # -------------------------------------------------------------------
    if (model == "M1") {
      m_k <- rowSums(n_g_counts) + rowSums(n_cf_counts)  # total alloc per component
      for (k in 1:(K-1)) {
        a_post <- 1 + m_k[k]
        b_post <- alpha_dp + sum(m_k[(k+1):K])
        v[k] <- rbeta(1, a_post, b_post)
        # Clamp to avoid exact 0 or 1
        v[k] <- max(min(v[k], 1 - 1e-10), 1e-10)
      }
      v[K] <- 1
      pi_vec <- stick_break(v)
    }

    # -------------------------------------------------------------------
    # STEP 4: Update source-specific weights omega_g, omega_cf
    # -------------------------------------------------------------------
    for (obs in 1:N_obs) {
      if (model == "M1") {
        alpha_g_post <- kappa * pi_vec + n_g_counts[, obs]
        alpha_cf_post <- kappa * pi_vec + n_cf_counts[, obs]
      } else {
        alpha_g_post  <- rep(alpha_dp/K, K) + n_g_counts[, obs]
        alpha_cf_post <- rep(alpha_dp/K, K) + n_cf_counts[, obs]
      }
      omega_g[, obs]  <- rdirichlet_one(alpha_g_post)
      omega_cf[, obs] <- rdirichlet_one(alpha_cf_post)
    }

    # -------------------------------------------------------------------
    # STEP 5: Update variance and nuisance parameters
    # -------------------------------------------------------------------
    # 5a. Contamination fractions wN[i]
    for (i in 1:n) {
      wN[i] <- rbeta(1, hp$a_N + N0_vec[i], hp$b_N + Nplus_vec[i])
      wN[i] <- max(min(wN[i], 1 - 1e-10), 1e-10)
    }

    # 5b. Background variance sigma_N^2
    ss_N <- alloc$ss_normal  # sum of Y_cf^2 for normal-allocated obs
    N_total_0 <- sum(N0_vec)
    sigma_N2 <- 1 / rgamma(1, hp$a_sigmaN + N_total_0/2,
                            hp$b_sigmaN + ss_N/2)

    # 5c. Component variances sigma_k^2
    for (k in 1:K) {
      ss_k <- alloc$ss_component[k]  # sum of (y - theta)^2 for component k
      N_k  <- sum(n_g_counts[k, ]) + sum(n_cf_counts[k, ])
      sigma_k2[k] <- 1 / rgamma(1, hp$a_sigma + N_k/2,
                                  hp$b_sigma + ss_k/2)
    }

    # -------------------------------------------------------------------
    # STEP 6: Update horseshoe shrinkage parameters
    # -------------------------------------------------------------------
    for (k in 1:K) {
      # Local shrinkage lambda_{jk}^2
      # Since theta is shared across time (simplified), T_eff = 1
      theta_k_sq <- theta[, k]^2
      for (j in 1:p) {
        rate_lam <- 1/nu_aux[j, k] + theta_k_sq[j] / (2 * tau_k2[k])
        lambda2[j, k] <- 1 / rgamma(1, 1, rate_lam)
      }

      # Local auxiliary nu_{jk}
      for (j in 1:p) {
        nu_aux[j, k] <- 1 / rgamma(1, 1, 1 + 1/lambda2[j, k])
      }

      # Global shrinkage tau_k^2
      sum_theta_lam <- sum(theta_k_sq / lambda2[, k])
      tau_k2[k] <- 1 / rgamma(1, (p + 1)/2,
                                1/xi_aux[k] + sum_theta_lam/2)

      # Global auxiliary xi_k
      xi_aux[k] <- 1 / rgamma(1, 1, 1 + 1/tau_k2[k])
    }

    # -------------------------------------------------------------------
    # STEP 7: Update kappa (M1 only, Metropolis-Hastings)
    # -------------------------------------------------------------------
    if (model == "M1") {
      kappa_total <- kappa_total + 1
      log_kappa_prop <- rnorm(1, log(kappa), hp$sigma_kappa_mh)
      kappa_prop <- exp(log_kappa_prop)

      # Compute full Dirichlet log-density ratio directly for numerical stability
      # log Dir(x | alpha) = lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum((alpha-1)*log(x))
      pi_safe <- pmax(pi_vec, 1e-10)
      pi_safe <- pi_safe / sum(pi_safe)

      log_ratio_prior <- (hp$a_kappa - 1) * (log(kappa_prop) - log(kappa)) -
                         hp$b_kappa * (kappa_prop - kappa)

      log_ratio_dir <- 0
      for (obs in 1:N_obs) {
        og <- pmax(omega_g[, obs], 1e-300)
        oc <- pmax(omega_cf[, obs], 1e-300)

        # For gDNA: log Dir(og | kappa_new*pi) - log Dir(og | kappa_old*pi)
        alpha_new <- kappa_prop * pi_safe
        alpha_old <- kappa * pi_safe

        ld_g_new <- lgamma(sum(alpha_new)) - sum(lgamma(alpha_new)) +
                    sum((alpha_new - 1) * log(og))
        ld_g_old <- lgamma(sum(alpha_old)) - sum(lgamma(alpha_old)) +
                    sum((alpha_old - 1) * log(og))

        # For cfDNA
        ld_cf_new <- lgamma(sum(alpha_new)) - sum(lgamma(alpha_new)) +
                     sum((alpha_new - 1) * log(oc))
        ld_cf_old <- lgamma(sum(alpha_old)) - sum(lgamma(alpha_old)) +
                     sum((alpha_old - 1) * log(oc))

        log_ratio_dir <- log_ratio_dir +
          (ld_g_new - ld_g_old) + (ld_cf_new - ld_cf_old)
      }

      # Jacobian for log-normal proposal
      log_alpha <- log_ratio_prior + log_ratio_dir +
                   log(kappa_prop) - log(kappa)

      if (!is.finite(log_alpha)) log_alpha <- -Inf

      if (log(runif(1)) < log_alpha) {
        kappa <- kappa_prop
        kappa_accept <- kappa_accept + 1
      }
    }

    # -------------------------------------------------------------------
    # STEP 8: Update alpha_dp (conjugate Gamma)
    # -------------------------------------------------------------------
    if (model == "M1") {
      sum_log_1mv <- sum(log(1 - v[1:(K-1)]))
      alpha_dp <- rgamma(1, hp$a_alpha + K - 1,
                          hp$b_alpha - sum_log_1mv)
    }

    # -------------------------------------------------------------------
    # Save samples (after burn-in, with thinning)
    # -------------------------------------------------------------------
    if (iter > n_burn && (iter - n_burn) %% thin == 0) {
      idx <- (iter - n_burn) %/% thin
      samples$kappa[idx]    <- kappa
      samples$alpha_dp[idx] <- alpha_dp
      samples$K_plus[idx]   <- sum(colSums(n_g_counts + n_cf_counts) > 0 |
                                    pi_vec > 0.01)
      samples$wN[idx, ]     <- wN
      samples$sigma_k[idx, ] <- sqrt(sigma_k2)
      samples$sigma_N[idx]  <- sqrt(sigma_N2)
      samples$v[idx, ]      <- v[1:(K-1)]

      # Compute log-likelihood for bridge sampling
      samples$loglik[idx] <- compute_loglik(Y_g, Y_cf, theta, sigma_k2,
                                             sigma_N2, omega_g, omega_cf,
                                             wN, patient, K, p)
    }

    # Progress
    if (verbose && iter %% 500 == 0) {
      K_active <- sum(pi_vec > 0.01)
      cat(sprintf("  iter %d/%d | kappa=%.1f | alpha=%.2f | K+=%.0f | wN_mean=%.3f\n",
                  iter, n_iter, kappa, alpha_dp, K_active, mean(wN)))
    }
  }

  if (verbose && model == "M1") {
    cat(sprintf("  kappa MH acceptance rate: %.1f%%\n",
                100 * kappa_accept / kappa_total))
  }

  # --- Return ---
  list(
    samples = samples,
    n_save  = n_save,
    model   = model,
    K       = K,
    n_iter  = n_iter,
    n_burn  = n_burn,
    kappa_accept_rate = if (model == "M1") kappa_accept/kappa_total else NA,
    final_theta = theta,
    final_pi    = pi_vec,
    hyperparams = hp
  )
}

# ===========================================================================
# HELPER FUNCTIONS
# ===========================================================================

#' Stick-breaking transformation: v -> pi
stick_break <- function(v) {
  K <- length(v)
  pi_vec <- numeric(K)
  cum_prod <- 1
  for (k in 1:K) {
    pi_vec[k] <- v[k] * cum_prod
    cum_prod <- cum_prod * (1 - v[k])
  }
  # Normalize to handle numerical issues
  pi_vec <- pi_vec / sum(pi_vec)
  pi_vec
}

#' Single Dirichlet draw (avoids overhead of full package)
rdirichlet_one <- function(alpha) {
  x <- rgamma(length(alpha), alpha, 1)
  s <- sum(x)
  if (s < 1e-300) return(rep(1/length(alpha), length(alpha)))
  x / s
}

#' Compute allocations for all observations
#'
#' Returns allocation counts and sufficient statistics.
#' Uses vectorized operations over genes for speed.
compute_allocations <- function(Y_g, Y_cf, theta, sigma_k2, sigma_N2,
                                 omega_g, omega_cf, wN, patient, K, p) {
  N_obs <- ncol(Y_g)
  n <- max(patient)

  # Allocation matrices: z_g[j, obs] in {1,...,K}, z_cf[j, obs] in {0,...,K}
  z_g  <- matrix(1L, p, N_obs)
  z_cf <- matrix(0L, p, N_obs)

  # Counts
  n_g_counts  <- matrix(0, K, N_obs)
  n_cf_counts <- matrix(0, K, N_obs)
  N0_vec      <- numeric(n)
  Nplus_vec   <- numeric(n)
  ss_normal   <- 0        # sum of Y_cf^2 for normal-allocated
  ss_component <- numeric(K)  # sum of (y - theta_k)^2 for each component

  # For each observation (patient-timepoint pair)
  for (obs in 1:N_obs) {
    i <- patient[obs]

    # --- gDNA allocations ---
    # Log-probabilities for each component [K x p]
    log_probs_g <- matrix(NA, K, p)
    for (k in 1:K) {
      log_probs_g[k, ] <- log(omega_g[k, obs]) +
        dnorm(Y_g[, obs], theta[, k], sqrt(sigma_k2[k]), log = TRUE)
    }
    # Normalize and sample
    for (j in 1:p) {
      lp <- log_probs_g[, j]
      lp <- lp - max(lp)  # log-sum-exp trick
      probs <- exp(lp)
      probs <- probs / sum(probs)
      z_g[j, obs] <- sample.int(K, 1, prob = probs)
    }
    # Count allocations
    tab_g <- tabulate(z_g[, obs], nbins = K)
    n_g_counts[, obs] <- tab_g

    # --- cfDNA allocations ---
    # Normal component probability
    log_p0 <- log(wN[i]) + dnorm(Y_cf[, obs], 0, sqrt(sigma_N2), log = TRUE)

    # Tumor component probabilities
    log_probs_cf <- matrix(NA, K, p)
    for (k in 1:K) {
      log_probs_cf[k, ] <- log(1 - wN[i]) + log(omega_cf[k, obs]) +
        dnorm(Y_cf[, obs], theta[, k], sqrt(sigma_k2[k]), log = TRUE)
    }

    for (j in 1:p) {
      lp <- c(log_p0[j], log_probs_cf[, j])
      lp <- lp - max(lp)
      probs <- exp(lp)
      probs <- probs / sum(probs)
      z_cf[j, obs] <- sample.int(K + 1, 1, prob = probs) - 1L  # 0-indexed
    }

    # Count cfDNA allocations
    cf_alloc <- z_cf[, obs]
    n_normal <- sum(cf_alloc == 0)
    N0_vec[i] <- N0_vec[i] + n_normal
    Nplus_vec[i] <- Nplus_vec[i] + (p - n_normal)
    for (k in 1:K) {
      n_cf_counts[k, obs] <- sum(cf_alloc == k)
    }

    # Sufficient statistics for normal component
    ss_normal <- ss_normal + sum(Y_cf[cf_alloc == 0, obs]^2)
  }

  # Compute ss_component: sum of (y - theta_k)^2 over all allocations
  for (k in 1:K) {
    ss_k <- 0
    for (obs in 1:N_obs) {
      # gDNA allocated to k
      idx_g <- which(z_g[, obs] == k)
      if (length(idx_g) > 0) {
        ss_k <- ss_k + sum((Y_g[idx_g, obs] - theta[idx_g, k])^2)
      }
      # cfDNA allocated to k
      idx_cf <- which(z_cf[, obs] == k)
      if (length(idx_cf) > 0) {
        ss_k <- ss_k + sum((Y_cf[idx_cf, obs] - theta[idx_cf, k])^2)
      }
    }
    ss_component[k] <- ss_k
  }

  list(n_g = n_g_counts, n_cf = n_cf_counts,
       N0 = N0_vec, Nplus = Nplus_vec,
       ss_normal = ss_normal, ss_component = ss_component,
       z_g = z_g, z_cf = z_cf)
}

#' Update theta[j,k] (conjugate Normal)
update_theta <- function(Y_g, Y_cf, z_g, z_cf, sigma_k2,
                          tau_k2, lambda2, K, p, N_obs) {
  theta <- matrix(0, p, K)
  for (k in 1:K) {
    for (j in 1:p) {
      # Collect all observations allocated to component k at gene j
      y_vals <- c()
      for (obs in 1:N_obs) {
        if (z_g[j, obs] == k)  y_vals <- c(y_vals, Y_g[j, obs])
        if (z_cf[j, obs] == k) y_vals <- c(y_vals, Y_cf[j, obs])
      }
      n_jk <- length(y_vals)

      # Posterior precision and mean
      prior_prec <- 1 / (tau_k2[k] * lambda2[j, k])
      data_prec  <- n_jk / sigma_k2[k]
      post_prec  <- data_prec + prior_prec
      post_var   <- 1 / post_prec

      if (n_jk > 0) {
        post_mean <- post_var * (sum(y_vals) / sigma_k2[k])
      } else {
        post_mean <- 0
      }

      theta[j, k] <- rnorm(1, post_mean, sqrt(post_var))
    }
  }
  theta
}

#' Compute observed-data log-likelihood (for bridge sampling)
compute_loglik <- function(Y_g, Y_cf, theta, sigma_k2, sigma_N2,
                            omega_g, omega_cf, wN, patient, K, p) {
  N_obs <- ncol(Y_g)
  ll <- 0

  for (obs in 1:N_obs) {
    i <- patient[obs]
    # Sum over a subset of genes for computational tractability
    gene_subset <- seq(1, p, by = max(1, p %/% 500))
    for (j in gene_subset) {
      # gDNA
      log_mix_g <- -Inf
      for (k in 1:K) {
        log_mix_g <- log_sum_exp(log_mix_g,
          log(omega_g[k, obs]) + dnorm(Y_g[j, obs], theta[j, k],
                                        sqrt(sigma_k2[k]), log = TRUE))
      }
      ll <- ll + log_mix_g

      # cfDNA
      log_normal <- log(wN[i]) + dnorm(Y_cf[j, obs], 0, sqrt(sigma_N2), log = TRUE)
      log_tumor <- -Inf
      for (k in 1:K) {
        log_tumor <- log_sum_exp(log_tumor,
          log(1 - wN[i]) + log(omega_cf[k, obs]) +
          dnorm(Y_cf[j, obs], theta[j, k], sqrt(sigma_k2[k]), log = TRUE))
      }
      ll <- ll + log_sum_exp(log_normal, log_tumor)
    }
  }

  # Scale up to full p
  ll * (p / length(gene_subset))
}

#' Log-sum-exp for two values
log_sum_exp <- function(a, b) {
  m <- max(a, b)
  if (is.infinite(m)) return(-Inf)
  m + log(exp(a - m) + exp(b - m))
}
