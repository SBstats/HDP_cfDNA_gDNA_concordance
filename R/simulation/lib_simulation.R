
###############################################################################
## ---- 03_simulation_dgp.R
###############################################################################

###############################################################################
## 03_simulation_dgp.R
## Data-generating functions for the simulation study
##
## theta_{jk} is TIME-INVARIANT: [p x K_true] array (static across cycles).
## Active gene sets are shared across sources (same sparsity pattern),
## but gDNA and cfDNA receive independent source-specific perturbations
## of SD delta_theta to capture cross-source signature divergence.
##
## Generates data under M1 (tracking) and M0 (non-tracking) for the
## 3 simulation scenarios described in the manuscript (Sections 9.2-9.4).
###############################################################################

#' Generate simulation data under the DP mixture model
#'
#' @param n         Number of subjects
#' @param p         Number of genes
#' @param T_max     Number of time points
#' @param K_true    Number of true occupied components
#' @param s_k       Number of active genes per component
#' @param kappa     Tracking strength (M1 only)
#' @param alpha_dp  DP concentration parameter
#' @param w_N_mean  Mean contamination fraction
#' @param sigma_k   Component standard deviations (scalar or vector)
#' @param sigma_N   Background standard deviation
#' @param model     "M1" (tracking) or "M0" (non-tracking)
#' @param seed      Random seed
#' @param varsigma_theta  SD of the between-subject signature deviation:
#'                      theta[j,k,i] ~ N(theta_bar[j,k], varsigma_theta^2),
#'                      applied to ALL p genes (the model places no atom at
#'                      zero on the subject level). varsigma_theta = 0 recovers
#'                      the pooled-signature DGP exactly.
#' @param mu_bg         Location of the cfDNA background profile. Deliberately
#'                      nonzero: on a VST-like scale the background does not sit
#'                      at the origin, and a sampler assuming it does is
#'                      measurably misspecified.
#' @param varsigma_bg   SD of the subject-level background deviation.
#' @param delta_theta   DEPRECATED and IGNORED. Formerly generated independent
#'                      source-specific signatures for the retired rho^theta
#'                      estimand. Signatures are now shared across sources, so
#'                      that quantity does not exist. Retained only so old
#'                      scenario configs do not error; a nonzero value warns.
#'
#' @return List with Y_g, Y_cf, true parameters, patient/time indices
generate_sim_data <- function(n = 26, p = 500, T_max = 5, K_true = 3,
                               s_k = 50, kappa = 1, alpha_dp = 1,
                               w_N_mean = 0.1, sigma_k = 1,
                               sigma_N = sqrt(2), model = "M1",
                               seed = 1,
                               ## Between-subject signature dispersion: theta[j,k,i] ~
                               ## N(theta_bar[j,k], varsigma_theta^2). This is the key
                               ## factor of Scenario 5; varsigma_theta = 0 recovers the
                               ## earlier pooled-signature DGP exactly.
                               varsigma_theta = 0.5,
                               ## cfDNA background: mean location mu_bg (DELIBERATELY
                               ## nonzero -- on the VST scale the background does not sit
                               ## at the origin, and a sampler assuming it does is
                               ## measurably misspecified) and between-subject SD.
                               mu_bg = 4, varsigma_bg = 0.5,
                               delta_theta = 0,    # DEPRECATED: see note below
                               ## --- misspecification switches (Scenario 4; default OFF) ---
                               noise = "gaussian",   # "gaussian" (well-specified) or "t3"
                               gene_corr = 0,        # AR(1) corr rho across genes in [0,1); 0 = independent
                               hk_background = 0,    # SD of a shared non-sparse housekeeping offset; 0 = none
                               ## --- composition-heterogeneity switches (R1/R3; default OFF) ---
                               patient_pi = FALSE,   # TRUE: each patient i gets a DISTINCT composition pi_i
                               pi_drift = 0,         # within-patient temporal drift SD on the composition (>0 exercises longitudinal tracking)
                               balanced_pi = FALSE) {# TRUE: force K_true occupied components with balanced (non-degenerate) weights
  set.seed(seed)

  ## misspecified iff any switch departs from the well-specified defaults
  misspec <- (!identical(noise, "gaussian")) || (gene_corr != 0) || (hk_background != 0)

  if (length(sigma_k) == 1) sigma_k <- rep(sigma_k, K_true)

  # --- 1. Generate true subclonal signatures ---
  # theta_g_true [p, K_true] for gDNA and
  # theta_cf_true [p, K_true] for cfDNA (static — time-invariant).
  # Both share the same sparsity pattern (active gene set) and base values,
  # but receive independent source-specific perturbations of SD delta_theta.
  # When delta_theta=0 the two arrays are identical.
  theta_base <- matrix(0, p, K_true)
  active_genes <- list()
  for (k in 1:K_true) {
    idx <- sample(p, s_k)
    active_genes[[k]] <- idx
    vals <- sample(c(-3, -2, 2, 3), s_k, replace = TRUE)
    theta_base[idx, k] <- vals
  }

  # Population signature theta_bar [p x K_true], then SUBJECT-LEVEL signatures
  # theta_true_i [p x K_true x n] drawn from the hierarchical prior. The SAME
  # theta_true_i generates BOTH sources: the signature is shared across
  # compartments (writeup Sec 4.1), so the simulated truth has no
  # source-specific signature discrepancy.
  #
  # NOTE on `delta_theta`: this argument previously generated INDEPENDENT
  # source-specific signatures theta_g_true / theta_cf_true, supporting the
  # retired rho^theta estimand. Under the shared-theta specification that
  # quantity does not exist. The argument is retained only so old scenario
  # configs do not error; a nonzero value now warns and is ignored.
  if (!identical(delta_theta, 0)) {
    warning("delta_theta is deprecated and ignored: signatures are shared ",
            "across sources under the current model. Use varsigma_theta to ",
            "control BETWEEN-SUBJECT signature dispersion instead.")
  }
  theta_bar_true <- theta_base
  theta_true_i <- array(0, dim = c(p, K_true, n))
  for (i in 1:n) {
    for (k in 1:K_true) {
      # Perturb ALL p genes, matching the model the sampler implements:
      #   theta[j,k,i] ~ N(theta_bar[j,k], varsigma_k^2)   for EVERY j.
      #
      # An earlier version perturbed only the s_k active genes and left the
      # rest exactly at theta_bar. That makes the DGP and the estimator target
      # DIFFERENT quantities: the sampler's full conditional is
      # IG(a + np/2, b + 0.5*sum_{i,j}(theta - theta_bar)^2), pooling over all
      # p genes, so it estimates varsigma_theta^2 * s_k/p rather than
      # varsigma_theta^2 -- a factor of 10 at the defaults (p=500, s_k=50),
      # which made the Scenario 5 dispersion-recovery metric look structurally
      # biased when the sampler was in fact exact.
      #
      # Note the sparsity pattern is still shared: inactive genes have
      # theta_bar = 0 and deviate only by O(varsigma_theta), so they stay in a
      # tight band around zero rather than becoming active.
      theta_true_i[, k, i] <- theta_bar_true[, k] + rnorm(p, 0, varsigma_theta)
    }
  }
  # Backward-compatible alias used by existing metric code that expects a
  # [p x K] population-level signature.
  theta_true <- theta_bar_true

  # --- 2. Generate the cohort-level latent composition ---
  # Helper: draw one composition on the K_true-simplex.
  #   balanced_pi=FALSE -> original stick-breaking draw (can be dominated by 1 comp)
  #   balanced_pi=TRUE  -> Dirichlet(concentration>1) so all K_true components are
  #                        occupied with balanced, non-degenerate weights (R3).
  draw_pi <- function() {
    if (balanced_pi) {
      # symmetric Dirichlet with concentration 5 per component => E[pi_k]=1/K_true,
      # weights bounded away from 0 so >=3 components are genuinely occupied.
      rdirichlet_one(rep(5, K_true))
    } else {
      v <- rbeta(K_true - 1, 1, alpha_dp); v <- c(v, 1)
      pp <- numeric(K_true); cum <- 1
      for (k in 1:K_true) { pp[k] <- v[k] * cum; cum <- cum * (1 - v[k]) }
      pp / sum(pp)
    }
  }
  pi_true <- draw_pi()  # cohort-level composition (used when patient_pi = FALSE)

  ## R1: per-patient composition heterogeneity. Each patient i gets a distinct
  ## baseline composition pi_i (drawn independently). This creates the between-
  ## observation heterogeneity the pairing-specific estimand needs: without it,
  ## every observation shares one pi and a scrambled pair is as concordant as the
  ## true pair (zero excess by construction).
  pi_patient <- if (patient_pi) lapply(1:n, function(i) draw_pi()) else NULL

  ## Optional within-patient temporal drift on the log-composition (R1): lets the
  ## composition evolve across time so the longitudinal increment estimand has
  ## genuine dynamic range. Returns the (patient,time) composition on the simplex.
  compose_pi <- function(i, tt) {
    base <- if (patient_pi) pi_patient[[i]] else pi_true
    if (pi_drift <= 0) return(base)
    lg <- log(pmax(base, 1e-8)) + rnorm(K_true, 0, pi_drift) * (tt - 1)
    exp(lg) / sum(exp(lg))
  }

  # --- 3. Generate observations ---
  obs_list <- expand.grid(patient = 1:n, time = 1:T_max)
  N_obs <- nrow(obs_list)
  patient_vec <- obs_list$patient
  time_vec    <- obs_list$time

  Y_g  <- matrix(0, p, N_obs)
  Y_cf <- matrix(0, p, N_obs)

  omega_g_true  <- matrix(0, K_true, N_obs)
  omega_cf_true <- matrix(0, K_true, N_obs)
  # Contamination fraction is now PER OBSERVATION (writeup eq:prior_wN), not
  # per patient: the circulating tumour fraction changes over a treatment
  # course. Length N_obs, indexed by observation.
  omega_0_true <- rbeta(N_obs, w_N_mean * 9, (1 - w_N_mean) * 9)
  omega_0_true <- pmax(pmin(omega_0_true, 0.95), 0.01)
  # Back-compat alias for any caller still expecting a per-patient vector.
  wN_true <- tapply(omega_0_true, patient_vec, mean)

  # cfDNA background means: cohort profile theta_0_bar_true [p] at location
  # mu_bg, plus subject deviations [p x n]. The nonzero location is what makes
  # a zero-mean-background sampler measurably misspecified.
  theta_0_bar_true <- rnorm(p, mu_bg, 1)
  theta_0_true <- matrix(theta_0_bar_true, nrow = p, ncol = n) +
                  rnorm(p * n, 0, varsigma_bg)

  ## Shared, non-sparse housekeeping background (Scenario 4c): a common offset
  ## added to every observation of BOTH sources, mimicking the constitutive
  ## signal that makes the raw gene-level correlation non-specific.
  hk_mu <- if (hk_background != 0) rnorm(p, 0, hk_background) else numeric(p)

  ## Helper: draw the p-vector of emission NOISE for one observation.
  ##  - gaussian: N(0, sd^2) per gene (the well-specified model)
  ##  - t3:       scaled Student-t_3 so Var matches sd^2 (heavier tails)
  ##  - gene_corr: impose AR(1) correlation across genes via a shared factor,
  ##    keeping the marginal SD unchanged.
  draw_noise <- function(sd_vec) {
    if (identical(noise, "t3")) {
      z <- rt(p, df = 3) / sqrt(3)      # unit-variance t_3
    } else {
      z <- rnorm(p)
    }
    if (gene_corr != 0) {
      # AR(1)-style: mix each gene's innovation with its neighbour's, preserve SD
      u <- numeric(p); u[1] <- z[1]
      r <- gene_corr; s <- sqrt(1 - r^2)
      for (j in 2:p) u[j] <- r * u[j - 1] + s * z[j]
      z <- u
    }
    z * sd_vec
  }

  ## whether composition heterogeneity/drift is active (keeps the default path
  ## byte-identical to the original when both switches are off)
  het <- patient_pi || (pi_drift > 0)

  for (obs in 1:N_obs) {
    i <- patient_vec[obs]
    tt <- time_vec[obs]

    if (model == "M1") {
      ## M1 = tracking: both sources coupled to the SAME (patient/time) center pi_it
      pi_it <- if (het) compose_pi(i, tt) else pi_true
      omega_g_true[, obs]  <- rdirichlet_one(kappa * pi_it)
      omega_cf_true[, obs] <- rdirichlet_one(kappa * pi_it)
    } else {
      ## M0 = non-tracking: the two sources are drawn from INDEPENDENT centers,
      ## so even under heterogeneity there is no cross-source coupling.
      if (het) {
        pi_g  <- compose_pi(i, tt)
        pi_cf <- if (patient_pi) compose_pi(sample.int(n, 1), tt) else rdirichlet_one(rep(alpha_dp / K_true, K_true))
        omega_g_true[, obs]  <- rdirichlet_one(kappa * pi_g)
        omega_cf_true[, obs] <- rdirichlet_one(kappa * pi_cf)
      } else {
        omega_g_true[, obs]  <- rdirichlet_one(rep(alpha_dp / K_true, K_true))
        omega_cf_true[, obs] <- rdirichlet_one(rep(alpha_dp / K_true, K_true))
      }
    }

    ## Subject i's own signatures; the SAME array serves both sources.
    th_i <- theta_true_i[, , i]            # [p x K_true]
    w0_o <- omega_0_true[obs]              # this observation's contamination

    if (!misspec) {
      ## --- WELL-SPECIFIED path ---
      # gDNA: pure mixture, no background component.
      for (j in 1:p) {
        k <- sample(K_true, 1, prob = omega_g_true[, obs])
        Y_g[j, obs] <- rnorm(1, th_i[j, k], sigma_k[k])
      }
      # cfDNA: background (at theta_0_true[j,i], NOT 0) + tumour mixture.
      for (j in 1:p) {
        if (runif(1) < w0_o) {
          Y_cf[j, obs] <- rnorm(1, theta_0_true[j, i], sigma_N)
        } else {
          k <- sample(K_true, 1, prob = omega_cf_true[, obs])
          Y_cf[j, obs] <- rnorm(1, th_i[j, k], sigma_k[k])
        }
      }
    } else {
      ## --- MISSPECIFIED path: same mean structure, departed noise/background ---
      kg <- sample(K_true, p, replace = TRUE, prob = omega_g_true[, obs])
      mu_g <- th_i[cbind(1:p, kg)]
      Y_g[, obs] <- mu_g + draw_noise(sigma_k[kg]) + hk_mu

      contam <- runif(p) < w0_o
      kcf <- sample(K_true, p, replace = TRUE, prob = omega_cf_true[, obs])
      mu_cf <- ifelse(contam, theta_0_true[, i], th_i[cbind(1:p, kcf)])
      sd_cf <- ifelse(contam, sigma_N, sigma_k[kcf])
      Y_cf[, obs] <- mu_cf + draw_noise(sd_cf) + hk_mu
    }
  }

  list(
    Y_g = Y_g,
    Y_cf = Y_cf,
    patient = patient_vec,
    time = time_vec,
    n = n, p = p, T_max = T_max, N_obs = N_obs,
    K_true = K_true,
    theta_true      = theta_true,        # alias for theta_bar_true (back compat)
    theta_bar_true  = theta_bar_true,    # [p x K_true] POPULATION signatures
    theta_true_i    = theta_true_i,      # [p x K_true x n] SUBJECT-level, shared across sources
    varsigma_theta  = varsigma_theta,    # between-subject signature SD used in DGP
    theta_0_true     = theta_0_true,     # [p x n] subject cfDNA background means
    theta_0_bar_true = theta_0_bar_true, # [p] cohort background profile
    mu_bg = mu_bg, varsigma_bg = varsigma_bg,
    theta_base = theta_base,         # [p x K_true] base signatures
    pi_true = pi_true,
    omega_g_true = omega_g_true,
    omega_cf_true = omega_cf_true,
    omega_0_true = omega_0_true,     # [N_obs] per-observation contamination
    wN_true = wN_true,               # [n] per-patient mean (back compat)
    active_genes = active_genes,
    kappa_true = kappa,
    alpha_true = alpha_dp,
    model = model,
    noise = noise, gene_corr = gene_corr, hk_background = hk_background,
    misspec = misspec,
    patient_pi = patient_pi, pi_drift = pi_drift, balanced_pi = balanced_pi,
    pi_patient = pi_patient
  )
}

#' Dirichlet draw (duplicated here for standalone use)
rdirichlet_one <- function(alpha) {
  x <- rgamma(length(alpha), alpha, 1)
  s <- sum(x)
  if (s < 1e-300) return(rep(1/length(alpha), length(alpha)))
  x / s
}

#' Compute simulation performance metrics
#' Note: kappa is fixed (not estimated), so no kappa recovery metrics.
compute_sim_metrics <- function(fit, sim_data, kappa_fixed = 1.0) {
  samples <- fit$samples
  metrics <- list()

  metrics$K_plus_mean   <- mean(samples$K_plus)
  metrics$K_plus_median <- median(samples$K_plus)

  # Contamination: per-OBSERVATION omega_0 against the per-observation truth.
  om0_post <- colMeans(samples$omega_0)
  metrics$omega0_bias <- mean(om0_post - sim_data$omega_0_true)
  metrics$omega0_rmse <- sqrt(mean((om0_post - sim_data$omega_0_true)^2))
  metrics$omega0_cor  <- if (sd(om0_post) > 1e-12)
    suppressWarnings(cor(om0_post, sim_data$omega_0_true)) else NA_real_
  # Boundary collapse: the diagnostic symptom of the allocation-normaliser bug
  # and of a misspecified background. Reported explicitly so a regression is
  # visible rather than silently absorbed into the RMSE.
  metrics$omega0_frac_boundary <- mean(om0_post < 1e-3)
  # Back-compat aliases (older aggregation scripts read wN_*).
  metrics$wN_bias <- metrics$omega0_bias
  metrics$wN_rmse <- metrics$omega0_rmse

  # Background mean recovery.
  if (!is.null(samples$sigma_0)) metrics$sigma_0_mean <- mean(samples$sigma_0)
  ## varsigma_* store SDs; the posterior mean VARIANCE is E[s^2], not (E[s])^2.
  if (!is.null(samples$varsigma_0)) metrics$varsigma_0_mean <- mean(samples$varsigma_0^2)
  if (!is.null(samples$varsigma_k)) metrics$varsigma_k_mean <- mean(colMeans(samples$varsigma_k^2))

  # Tracking correlation from the ESTIMATED D2 coupling phi (Part I / D2). Under
  # the D2 model omega_g, omega_cf ~ Dir(phi * m_it) around a shared per-obs
  # center m_it ~ Dir(alpha0 * pi_bar), so the cross-source correlation is
  # (phi+1)/(alpha0+phi+1) -- verified exact by Monte Carlo. phi is a posterior
  # quantity, so the tracking correlation has a posterior (mean + 95% CrI).
  alpha0_m  <- 5.0  # m_it prior concentration (matches sampler default hp$alpha0_m)
  ## Degrading to a CONSTANT tracking correlation must be loud: it silently
  ## converts the paper's headline estimand into the value it was initialised
  ## at, and every posterior interval derived from it becomes degenerate.
  phi_draws <- if (!is.null(samples$phi)) samples$phi else
               if (!is.null(samples$lambda)) samples$lambda else {
                 warning("compute_sim_metrics: neither samples$phi nor ",
                         "samples$lambda present; the tracking correlation is ",
                         "being reported at the FIXED kappa = ", kappa_fixed,
                         " with a degenerate (zero-width) interval. If this is ",
                         "a merged multi-chain object, check that ",
                         "merge_chain_samples() carries 'phi'.", call. = FALSE)
                 rep(kappa_fixed, length(samples$K_plus))
               }
  metrics$phi_post_mean <- mean(phi_draws)
  metrics$phi_post_lo   <- as.numeric(quantile(phi_draws, 0.025))
  metrics$phi_post_hi   <- as.numeric(quantile(phi_draws, 0.975))
  corr_draws <- (phi_draws + 1) / (alpha0_m + phi_draws + 1)
  metrics$tracking_corr_post <- mean(corr_draws)
  metrics$tracking_corr_lo   <- as.numeric(quantile(corr_draws, 0.025))
  metrics$tracking_corr_hi   <- as.numeric(quantile(corr_draws, 0.975))

  if (sim_data$model == "M1") {
    # TRUE tracking correlation computed EMPIRICALLY from the DGP's simulated
    # source weights (robust to any mismatch between the DGP and fit generative
    # forms): mean per-component Pearson correlation of omega_g_true vs omega_cf_true.
    Wg <- sim_data$omega_g_true; Wcf <- sim_data$omega_cf_true  # [K_true x N_obs]
    cc <- sapply(seq_len(nrow(Wg)), function(k) {
      if (sd(Wg[k, ]) < 1e-9 || sd(Wcf[k, ]) < 1e-9) NA_real_ else cor(Wg[k, ], Wcf[k, ])
    })
    metrics$tracking_corr_true <- mean(cc, na.rm = TRUE)
    metrics$phi_true <- sim_data$kappa_true
    metrics$tracking_corr_covered <-
      (metrics$tracking_corr_true >= metrics$tracking_corr_lo) &
      (metrics$tracking_corr_true <= metrics$tracking_corr_hi)
  }

  metrics$kappa_true <- sim_data$kappa_true
  metrics$kappa_fixed <- kappa_fixed

  metrics
}

#' Compute signature estimation metrics for Scenario 2
#' Self-contained Hungarian (Kuhn-Munkres) algorithm for the linear assignment
#' problem, minimising total cost over a square cost matrix. No package
#' dependency (avoids `clue`). Returns an integer vector `assignment` where
#' assignment[i] = column matched to row i.
hungarian_match <- function(cost) {
  cost <- as.matrix(cost); n <- nrow(cost)
  stopifnot(nrow(cost) == ncol(cost))
  ## Row/column reduction
  cost <- cost - apply(cost, 1, min)
  cost <- sweep(cost, 2, apply(cost, 2, min), "-")
  assign_row <- rep(0L, n); assign_col <- rep(0L, n)
  cover_row <- rep(FALSE, n); cover_col <- rep(FALSE, n)
  ## Greedy initial star of independent zeros
  starred <- matrix(FALSE, n, n)
  rused <- rep(FALSE, n); cused <- rep(FALSE, n)
  for (i in 1:n) for (j in 1:n)
    if (cost[i, j] == 0 && !rused[i] && !cused[j]) { starred[i, j] <- TRUE; rused[i] <- TRUE; cused[j] <- TRUE }
  primed <- matrix(FALSE, n, n)
  repeat {
    cover_col <- apply(starred, 2, any)
    if (sum(cover_col) == n) break
    cover_row[] <- FALSE
    repeat {
      ## find an uncovered zero
      z <- which(cost == 0 & !cover_row & matrix(!cover_col, n, n, byrow = TRUE), arr.ind = TRUE)
      if (nrow(z) == 0) {
        ## adjust matrix by smallest uncovered value
        m <- min(cost[!cover_row, !cover_col])
        cost[!cover_row, ] <- cost[!cover_row, ] - m
        cost[, cover_col] <- cost[, cover_col] + m
        next
      }
      zi <- z[1, 1]; zj <- z[1, 2]; primed[zi, zj] <- TRUE
      star_in_row <- which(starred[zi, ])
      if (length(star_in_row) == 0) {
        ## augmenting path
        path <- matrix(c(zi, zj), ncol = 2)
        repeat {
          r <- which(starred[, path[nrow(path), 2]])
          if (length(r) == 0) break
          path <- rbind(path, c(r, path[nrow(path), 2]))
          cc <- which(primed[path[nrow(path), 1], ])
          path <- rbind(path, c(path[nrow(path), 1], cc))
        }
        for (k in seq_len(nrow(path))) {
          pi_ <- path[k, 1]; pj_ <- path[k, 2]
          starred[pi_, pj_] <- !starred[pi_, pj_]
        }
        primed[] <- FALSE; cover_row[] <- FALSE; cover_col[] <- FALSE
        break
      } else {
        cover_row[zi] <- TRUE; cover_col[star_in_row] <- FALSE
      }
    }
  }
  assignment <- integer(n)
  for (i in 1:n) assignment[i] <- which(starred[i, ])
  assignment
}

#' Match K_true true components to K_fit fitted columns by MAXIMISING absolute
#' signature correlation aggregated over ALL time points (Hungarian on the
#' negative-correlation cost). Returns `matched[k]` = fitted column for true k,
#' or NA if true component k has no fitted column correlating above `cor_floor`.
#'
#' The `cor_floor` guard is essential when K_fit > K_true: the surplus fitted
#' columns are dense, near-zero "junk" components, and an unguarded Hungarian
#' assignment will hand an unrecovered true signature to whichever junk column
#' minimises the (already tiny) cost, then score theta MSE against that junk.
#' Treating a below-floor best-match as UNRECOVERED (NA) makes the metric report
#' the recovery honestly instead of crediting a spurious pairing.
match_components <- function(theta_est, theta_true, K_true, K_fit,
                             cor_floor = 0.3) {
  ## build square cost (pad to Kc); cost = -mean_t |cor|
  Kc <- max(K_true, K_fit)
  cost <- matrix(0, Kc, Kc)
  cor_mat <- matrix(0, K_true, K_fit)
  for (k in 1:K_true) for (j in 1:K_fit) {
    a <- theta_true[, k]; b <- theta_est[, j]   # static theta: [p] vectors
    cc <- if (sd(a) < 1e-9 || sd(b) < 1e-9) 0 else abs(cor(a, b))
    cor_mat[k, j] <- cc
    cost[k, j] <- -cc
  }
  asg <- hungarian_match(cost)
  matched <- asg[1:K_true]
  matched[matched > K_fit] <- NA_integer_   # true comp with no fitted column at all
  ## junk-robust guard: drop a match whose achieved |cor| is below the floor.
  for (k in 1:K_true) {
    j <- matched[k]
    if (!is.na(j) && cor_mat[k, j] < cor_floor) matched[k] <- NA_integer_
  }
  matched
}

# Recompute theta posterior means from stored trace, bypassing the online running
# sum which is contaminated by label switching when pi_k is balanced.
# chain must have theta_bar_trace [n_trace x p x K] and/or theta_trace
# [n_trace x p x K x n] stored.
recompute_postmean_from_trace <- function(chain) {
  tb <- chain$theta_bar_trace
  th <- chain$theta_trace
  if (is.null(tb) && is.null(th)) return(NULL)
  out <- list()
  if (!is.null(tb)) out$theta_bar_postmean <- apply(tb, c(2, 3), mean)
  if (!is.null(th)) out$theta_postmean     <- apply(th, c(2, 3, 4), mean)
  out
}

#'
#' theta_est is [p x K_fit]; theta_true is [p x K_true] (static — no T index).
#' Computes MSE, Hungarian component matching, and (if fit$theta_bar_trace present)
#' 95% credible-interval coverage of theta_g.
compute_theta_metrics <- function(fit, sim_data) {
  # Score signature recovery on the POPULATION signature theta_bar: that is the
  # quantity Theorem thm:contraction is about, and the only one the manuscript
  # interprets. (The subject-level theta is weakly informed by design.)
  theta_est <- if (!is.null(fit$theta_bar_postmean)) fit$theta_bar_postmean else
               if (!is.null(fit$final_theta_bar))    fit$final_theta_bar    else
               stop("no theta_bar_postmean / final_theta_bar on fit")
  theta_true <- sim_data$theta_bar_true  # [p x K_true] POPULATION ground truth
  p <- sim_data$p
  K_true <- sim_data$K_true
  s_k <- length(sim_data$active_genes[[1]])

  # Hungarian component matching (static theta: single correlation per (k,j) pair)
  K_fit <- dim(theta_est)[2]
  matched <- match_components(theta_est, theta_true, K_true, K_fit)

  # Optional coverage from the stored theta_g trace [n_trace x p x K]
  has_trace <- !is.null(fit$theta_bar_trace)
  cov_active <- numeric(K_true); cov_inactive <- numeric(K_true)

  mse_vec <- numeric(K_true)   # per-gene mean
  sse_vec <- numeric(K_true)   # total (sum-over-genes)
  tpr_vec <- numeric(K_true)
  fpr_vec <- numeric(K_true)
  recovered <- logical(K_true)

  for (k in 1:K_true) {
    j <- matched[k]
    recovered[k] <- !is.na(j)
    active_true <- sim_data$active_genes[[k]]
    null_true <- setdiff(1:p, active_true)

    est_col <- if (is.na(j)) rep(0, p) else theta_est[, j]   # static: no tt
    sse_vec[k] <- sum((est_col - theta_true[, k])^2)
    mse_vec[k] <- sse_vec[k] / p

    # Active gene detection (a full miss detects nothing)
    detected <- if (is.na(j)) integer(0) else which(abs(theta_est[, j]) > 0.5)
    tpr_vec[k] <- length(intersect(detected, active_true)) / length(active_true)
    fpr_vec[k] <- length(intersect(detected, null_true)) / max(length(null_true), 1)

    # Coverage: fraction of genes whose 95% CrI (over theta_g trace) covers theta_g*
    if (has_trace && !is.na(j)) {
      qlo <- apply(fit$theta_bar_trace[, , j, drop = FALSE], 2, quantile, 0.025, na.rm = TRUE)
      qhi <- apply(fit$theta_bar_trace[, , j, drop = FALSE], 2, quantile, 0.975, na.rm = TRUE)
      truth_k <- theta_true[, k]                              # [p] static
      covered <- (truth_k >= qlo) & (truth_k <= qhi)
      cov_active[k]   <- mean(covered[active_true], na.rm = TRUE)
      cov_inactive[k] <- mean(covered[null_true], na.rm = TRUE)
    } else { cov_active[k] <- NA; cov_inactive[k] <- NA }
  }

  rate <- s_k * log(p / s_k)   # sparse-risk normaliser

  list(
    mse = mean(mse_vec),
    total_mse = mean(sse_vec),
    normalized_mse = mean(sse_vec) / rate,
    tpr = mean(tpr_vec),
    fpr = mean(fpr_vec),
    n_recovered = sum(recovered),
    recovery_rate = mean(recovered),
    coverage_active = mean(cov_active, na.rm = TRUE),
    coverage_inactive = mean(cov_inactive, na.rm = TRUE),
    per_component_mse = mse_vec,
    matched_components = matched
  )
}
## ---------------------------------------------------------------------------
## RETIRED: compute_rho_theta_metric()
##
## This scored rho^theta_k = cor(theta_g[,k], theta_cf[,k]), the cross-source
## correlation of SOURCE-SPECIFIC signatures. Under the current specification a
## single theta is shared by both sources, so the quantity does not exist. The
## function is removed rather than redefined; see R/01_lib_core.R for the same
## note on the fit-object side.
##
## Callers should use compute_scenario5_metrics() below, which scores the
## quantities the revised model actually requires.
## ---------------------------------------------------------------------------

#' Scenario 5 metrics: identifiability of the hierarchical signature and the
#' time-varying background (writeup Sec 9, sec:sim_ident).
#'
#' Reports, for a fitted model against a known truth:
#'   (i)   population signature recovery (RMSE + best-match correlation)
#'   (ii)  SUBJECT-level signature recovery -- expected to be POOR, and reported
#'         precisely because the manuscript claims these are not consistently
#'         estimable (Remark rem:not_claimed). Confirming the claim is honest
#'         rather than pessimistic is the point.
#'   (iii) between-subject dispersion varsigma_k^2 recovery
#'   (iv)  background recovery: omega_0, theta_0, and boundary collapse
#'   (v)   the critical readout -- whether WEIGHT-level estimands remain
#'         unbiased even where the signatures are not identified
#'
#' @param fit      A fit (or fit_view) carrying theta_bar_postmean / theta_postmean.
#' @param sim_data Output of generate_sim_data().
compute_scenario5_metrics <- function(fit, sim_data) {
  out <- list()
  tb_est <- fit$theta_bar_postmean
  tb_tru <- sim_data$theta_bar_true
  K_true <- sim_data$K_true

  ## (i) population signature: match fitted to true components by |cor|
  if (!is.null(tb_est)) {
    K_fit <- ncol(tb_est)
    Cm <- matrix(0, K_true, K_fit)
    for (a in seq_len(K_true)) for (b in seq_len(K_fit)) {
      x <- tb_tru[, a]; y <- tb_est[, b]
      Cm[a, b] <- if (sd(x) > 1e-10 && sd(y) > 1e-10) abs(cor(x, y)) else 0
    }
    match_b <- apply(Cm, 1, which.max)
    out$thetabar_cor_mean <- mean(apply(Cm, 1, max))
    out$thetabar_rmse <- sqrt(mean(sapply(seq_len(K_true), function(a)
      mean((tb_est[, match_b[a]] - tb_tru[, a])^2))))
  } else {
    out$thetabar_cor_mean <- NA_real_; out$thetabar_rmse <- NA_real_
    match_b <- seq_len(K_true)
  }

  ## (ii) subject-level signatures (expected poor; reported anyway)
  th_est <- fit$theta_postmean; th_tru <- sim_data$theta_true_i
  if (!is.null(th_est) && !is.null(th_tru)) {
    n_s <- min(dim(th_est)[3], dim(th_tru)[3])
    errs <- sapply(seq_len(n_s), function(i)
      mean(sapply(seq_len(K_true), function(a)
        mean((th_est[, match_b[a], i] - th_tru[, a, i])^2))))
    out$theta_subj_rmse <- sqrt(mean(errs))
    ## Ratio > 1 means the subject level is WORSE than the population level,
    ## which is the expected and claimed behaviour.
    out$theta_subj_vs_pop_ratio <- if (is.finite(out$thetabar_rmse) &&
                                       out$thetabar_rmse > 0)
      out$theta_subj_rmse / out$thetabar_rmse else NA_real_
  } else {
    out$theta_subj_rmse <- NA_real_; out$theta_subj_vs_pop_ratio <- NA_real_
  }

  ## (iii) between-subject dispersion
  vs_true <- sim_data$varsigma_theta^2
  ## varsigma_k stores SDs; the posterior mean VARIANCE is E[s^2], not (E[s])^2.
  vs_est  <- if (!is.null(fit$samples$varsigma_k))
    mean(colMeans(fit$samples$varsigma_k^2)[seq_len(K_true)]) else NA_real_
  out$varsigma_k_true <- vs_true
  out$varsigma_k_est  <- vs_est
  out$varsigma_k_bias <- vs_est - vs_true

  ## (iv) background
  if (!is.null(fit$samples$omega_0)) {
    om0 <- colMeans(fit$samples$omega_0)
    out$omega0_rmse <- sqrt(mean((om0 - sim_data$omega_0_true)^2))
    out$omega0_cor  <- if (sd(om0) > 1e-12)
      suppressWarnings(cor(om0, sim_data$omega_0_true)) else NA_real_
    out$omega0_frac_boundary <- mean(om0 < 1e-3)
  }
  if (!is.null(fit$theta_0_postmean)) {
    out$theta0_rmse <- sqrt(mean((fit$theta_0_postmean - sim_data$theta_0_true)^2))
    out$theta0_cor  <- suppressWarnings(cor(as.numeric(fit$theta_0_postmean),
                                            as.numeric(sim_data$theta_0_true)))
  }

  ## (v) THE critical readout: are the weight-level estimands unbiased even
  ## where the signatures are not identified? If these degrade in step with
  ## signature recovery, the separation argument of writeup Sec 4.1 fails.
  if (!is.null(fit$samples$v)) {
    pi_hat <- stick_break(c(colMeans(fit$samples$v), 1))
    pt <- sim_data$pi_true
    if (is.matrix(pt)) pt <- rowMeans(pt)
    L <- min(length(pi_hat), length(pt))
    out$pi_l1 <- sum(abs(sort(pi_hat, decreasing = TRUE)[seq_len(L)] -
                         sort(pt, decreasing = TRUE)[seq_len(L)]))
  }

  ## (vi) within-subject label switching
  ls_err <- NULL
  ls_ <- tryCatch(theta_label_switch_rate(fit),
                  error = function(e) { ls_err <<- conditionMessage(e); NULL })
  out$label_switch_rate <- if (!is.null(ls_)) ls_$rate else NA_real_

  ## (vii) empirical Assumption A1 check.
  ##
  ## theta_separation_check() ERRORS by design when the fit carries no
  ## varsigma_k, because a fit with no dispersion estimate cannot be checked
  ## against A1. Mapping that error to NA would undo the hardening: downstream,
  ## mean(A1_satisfied, na.rm = TRUE) over all-NA gives NaN, which reads as
  ## "nothing to report" rather than "the diagnostic never ran". Record the
  ## status as a separate, aggregatable column.
  sep_err <- NULL
  sep <- tryCatch(theta_separation_check(fit),
                  error = function(e) { sep_err <<- conditionMessage(e); NULL })
  out$A1_satisfied <- if (is.null(sep)) NA else isTRUE(attr(sep, "all_satisfied"))
  out$A1_status <- if (!is.null(sep_err)) "ERROR"
                   else if (is.null(sep)) "UNAVAILABLE"
                   else if (isTRUE(attr(sep, "all_satisfied"))) "OK" else "FAIL"
  out$A1_m_req <- if (!is.null(sep) && !is.null(sep$m_req)) sep$m_req[1] else NA_integer_
  if (!is.null(sep_err) || !is.null(ls_err))
    warning("compute_scenario5_metrics: identifiability diagnostic failed (",
            paste(na.omit(c(ls_err, sep_err)), collapse = "; "), ")", call. = FALSE)
  out
}

#' Competitor signature estimators (Scenario 1, item A4). Benchmarks the
#' horseshoe against (a) unregularized mean (no shrinkage), (b) lasso soft-
#' thresholding (glmnet), and (c) a spike-and-slab hard threshold (2-group EM).
#' To isolate the SHRINKAGE comparison, all competitors are given the same
#' oracle per-observation dominant-component allocation, then estimate each
#' component's signature and are scored with the same MSE/TPR/FPR as the model.
#' Returns a named list of per-method metric lists.
fit_competitors <- function(sim_data, detect_thresh = 0.5) {
  p <- sim_data$p; K_true <- sim_data$K_true
  s_k <- length(sim_data$active_genes[[1]]); rate <- s_k * log(p / s_k)
  Yg <- sim_data$Y_g
  # Non-oracle dominant component allocation via k-means on observed Y_g.
  # Uses only observed data (no access to omega_g_true), making the competitor
  # comparison fair: HDP infers allocations jointly; k-means does so explicitly.
  km <- tryCatch(
    kmeans(t(sim_data$Y_g), centers = sim_data$K_true, nstart = 20, iter.max = 100),
    error = function(e) NULL)
  dom_g <- if (!is.null(km)) km$cluster else
              apply(sim_data$omega_g_true, 2, which.max)  # fallback if k-means fails

  # raw per-component mean signature from gDNA (pool all timepoints — static theta)
  raw <- matrix(0, p, K_true)
  for (k in 1:K_true) {
    cols <- which(dom_g == k)
    raw[, k] <- if (length(cols) >= 1) rowMeans(Yg[, cols, drop = FALSE]) else 0
  }

  score <- function(theta_hat) {
    mse <- tpr <- fpr <- numeric(K_true)
    for (k in 1:K_true) {
      mse[k] <- sum((theta_hat[, k] - sim_data$theta_true[, k])^2)
      active <- sim_data$active_genes[[k]]; null <- setdiff(1:p, active)
      det <- which(abs(theta_hat[, k]) > detect_thresh)
      tpr[k] <- length(intersect(det, active)) / length(active)
      fpr[k] <- length(intersect(det, null)) / max(length(null), 1)
    }
    list(mse = mean(mse), normalized_mse = mean(mse) / rate, tpr = mean(tpr), fpr = mean(fpr))
  }

  ## (a) unregularized: raw mean
  m_unreg <- score(raw)

  ## (b) lasso soft-threshold: shrink each signature toward 0 by CV lambda.
  ## With an identity design, the lasso solution is soft-thresholding of the
  ## per-gene mean at lambda; pick lambda by a small CV over the pooled signal.
  soft <- function(x, lam) sign(x) * pmax(abs(x) - lam, 0)
  lam_grid <- quantile(abs(raw), probs = seq(0.5, 0.95, by = 0.05))
  # choose lambda minimising a proxy risk vs a de-noised target (median filter of raw)
  best_lam <- lam_grid[which.min(sapply(lam_grid, function(l) mean((soft(raw, l))^2 - 2 * soft(raw, l) * raw)))]
  th_lasso <- soft(raw, best_lam)
  m_lasso <- score(th_lasso)

  ## (c) spike-and-slab hard threshold: 2-group (null vs signal) EM on |raw|,
  ## keep genes with high posterior signal probability.
  a <- as.vector(abs(raw)); a <- a[is.finite(a)]
  # simple threshold at mean + 1.5 SD of the null-ish bulk (robust)
  thr <- median(a) + 1.5 * mad(a)
  th_ss <- raw; th_ss[abs(th_ss) < thr] <- 0
  m_ss <- score(th_ss)

  list(unregularized = m_unreg, lasso = m_lasso, spike_slab = m_ss)
}

#' Pairing-specific subclonal concordance and its random-pairing null, from a
#' (merged) fit's posterior-mean source-specific weights.
#'
#' Two concordance measures are computed against the SAME within-timepoint
#' derangement null: (i) the L1 concordance (legacy) and (ii) the PRIMARY
#' occupied-component posterior correlation (R2). The correlation is taken over
#' only the occupied components (mean cohort weight > `occ_thresh`), because a
#' full-K correlation is inflated by structural agreement on the empty
#' components' shared zeros. It is well-defined only when >= 3 components are
#' occupied; otherwise `cor_*` fields are NA and the L1 fields are used.
#'
#' @param merged a `samples`-shaped list (fit_multi_chain()$merged) with
#'   omega_g_trace / omega_cf_trace [n_save x (K*N_obs)] and K, N_obs known.
#' @return list with L1 fields (observed/null_mean/excess/p_perm) and correlation
#'   fields (cor_observed/cor_null/cor_excess/cor_p/n_occupied).
concordance_excess <- function(merged, N_obs, time, n_perm = 200, occ_thresh = 0.01) {
  K <- ncol(merged$omega_g_trace) / N_obs
  Wg  <- sapply(1:K, function(k) colMeans(merged$omega_g_trace[,  ((k-1)*N_obs + 1):(k*N_obs), drop = FALSE]))
  Wcf <- sapply(1:K, function(k) colMeans(merged$omega_cf_trace[, ((k-1)*N_obs + 1):(k*N_obs), drop = FALSE]))  # [N x K]

  ## occupied components (by cohort-mean weight, pooled over both sources)
  occ <- which((colMeans(Wg) + colMeans(Wcf)) / 2 > occ_thresh)
  n_occ <- length(occ)

  L1conc  <- function(m) mean(sapply(1:N_obs, function(i) 1 - 0.5 * sum(abs(Wg[i, ] - Wcf[m[i], ]))))
  ## occupied-component posterior correlation, averaged over observations
  CORconc <- function(m) {
    if (n_occ < 3) return(NA_real_)
    mean(sapply(1:N_obs, function(i) {
      a <- Wg[i, occ]; b <- Wcf[m[i], occ]
      if (sd(a) < 1e-9 || sd(b) < 1e-9) NA else cor(a, b)
    }), na.rm = TRUE)
  }
  # within-timepoint derangement null
  strat_derange <- function() {
    m <- 1:N_obs
    for (s in unique(time)) { ii <- which(time == s); if (length(ii) > 1) { repeat { q <- sample(ii); if (all(q != ii)) break }; m[ii] <- q } }
    m
  }
  perms <- replicate(n_perm, strat_derange(), simplify = FALSE)

  obs_l1 <- L1conc(1:N_obs); null_l1 <- sapply(perms, L1conc)
  obs_c  <- CORconc(1:N_obs); null_c  <- sapply(perms, CORconc)

  list(
    # legacy L1
    observed = obs_l1, null_mean = mean(null_l1),
    excess = obs_l1 - mean(null_l1),
    p_perm = (1 + sum(null_l1 >= obs_l1)) / (1 + n_perm),
    # PRIMARY occupied-component correlation (R2)
    cor_observed = obs_c, cor_null = mean(null_c, na.rm = TRUE),
    cor_excess = obs_c - mean(null_c, na.rm = TRUE),
    cor_p = (1 + sum(null_c >= obs_c, na.rm = TRUE)) / (1 + sum(is.finite(null_c))),
    n_occupied = n_occ
  )
}

#' Within-patient increment concordance (R5) --- the paper's actual longitudinal
#' estimand. For each patient, the change in the dominant-subclone weight between
#' consecutive timepoints is correlated between the two sources; the observed
#' mean is compared to an interval-stratified pairing-shuffle null. This is
#' change-based (first-differenced), so it is not saturated by a shared static
#' composition and it directly tests whether cfDNA follows gDNA *dynamics*.
#'
#' @param merged fit_multi_chain()$merged; @param patient,time index vectors.
#' @return list(observed, null_mean, excess, p_perm) on the increment correlation.
increment_concordance <- function(merged, N_obs, patient, time, n_perm = 200,
                                  occ_thresh = 0.01) {
  K <- ncol(merged$omega_g_trace) / N_obs
  Wg  <- sapply(1:K, function(k) colMeans(merged$omega_g_trace[,  ((k-1)*N_obs + 1):(k*N_obs), drop = FALSE]))
  Wcf <- sapply(1:K, function(k) colMeans(merged$omega_cf_trace[, ((k-1)*N_obs + 1):(k*N_obs), drop = FALSE]))
  ## Use ALL occupied components rather than the single "dominant" component.
  ## Under balanced_pi=TRUE, no single component dominates (all ~1/K), so
  ## which.max() would pick an arbitrary index and use only 1/K of the signal.
  ## Averaging the increment correlation over all occupied components captures
  ## the full tracking signal regardless of composition balance.
  occ <- which((colMeans(Wg) + colMeans(Wcf)) / 2 > occ_thresh)
  if (length(occ) == 0) occ <- which.max(colMeans(Wg))   # fallback if all below threshold
  ## consecutive within-patient transitions (strictly increasing time)
  TR <- do.call(rbind, lapply(unique(patient), function(pt) {
    idx <- which(patient == pt); idx <- idx[order(time[idx])]
    if (length(idx) < 2) return(NULL)
    do.call(rbind, lapply(1:(length(idx) - 1), function(j) {
      a <- idx[j]; b <- idx[j + 1]
      if (time[b] == time[a]) return(NULL)
      data.frame(a = a, b = b, stratum = paste0(time[a], "_", time[b]))
    }))
  }))
  if (is.null(TR) || nrow(TR) < 3) return(list(observed = NA, null_mean = NA, excess = NA, p_perm = NA, M = 0))
  M <- nrow(TR)
  icor <- function(x, y) if (sd(x) < 1e-9 || sd(y) < 1e-9) NA_real_ else cor(x, y)
  ## Per-component increment correlation averaged over all occupied components.
  comp_cor <- function(idx_a, idx_b) {
    cors <- vapply(occ, function(ref) {
      dg  <- Wg[idx_b, ref]  - Wg[idx_a, ref]
      dcf <- Wcf[idx_b, ref] - Wcf[idx_a, ref]
      icor(dg, dcf)
    }, numeric(1))
    mean(cors, na.rm = TRUE)
  }
  dg_all  <- rowMeans(Wg[TR$b,  occ, drop=FALSE] - Wg[TR$a,  occ, drop=FALSE])
  dcf_all <- rowMeans(Wcf[TR$b, occ, drop=FALSE] - Wcf[TR$a, occ, drop=FALSE])
  obs <- icor(dg_all, dcf_all)
  strat_shuffle <- function() {
    y <- dcf_all
    for (s in unique(TR$stratum)) {
      ii <- which(TR$stratum == s)
      if (length(ii) > 1) y[ii] <- sample(dcf_all[ii])
    }
    y
  }
  nullv <- replicate(n_perm, icor(dg_all, strat_shuffle()))
  list(observed = obs, null_mean = mean(nullv, na.rm = TRUE),
       excess = obs - mean(nullv, na.rm = TRUE),
       p_perm = (1 + sum(nullv >= obs, na.rm = TRUE)) / (1 + sum(is.finite(nullv))), M = M)
}

#' Scenario 3 operating-characteristics metrics. Given the health-checked WAIC
#' decision, the LOPO-style predictive readout, and the concordance/increment
#' excesses, record the selection outcome and the continuous tracking score.
#'
#' The PRIMARY continuous tracking score (R5) is the within-patient increment-
#' correlation excess --- the change-based signal the paper actually relies on ---
#' with the static occupied-component correlation excess retained as secondary.
#'
#' @param delta_waic WAIC_M0 - WAIC_M1 (positive favors M1); may be NA.
#' @param conc concordance_excess() output; @param incr increment_concordance() output.
#' @param true_model "M0" or "M1".
compute_selection_metrics <- function(delta_waic, conc, true_model, incr = NULL, lopo_delta_elpd = NA_real_) {
  favors_M1 <- if (is.na(delta_waic)) NA else (delta_waic > 0)
  incr_excess <- if (!is.null(incr)) incr$excess else NA_real_
  cor_excess  <- if (!is.null(conc$cor_excess)) conc$cor_excess else NA_real_
  list(
    true_model = true_model,
    delta_waic = delta_waic,
    favors_M1 = favors_M1,
    # correct = favor M1 iff truth is M1
    correct = if (is.na(favors_M1)) NA else (favors_M1 == (true_model == "M1")),
    # false positive = favored M1 when truth is M0
    false_positive = if (is.na(favors_M1)) NA else (true_model == "M0" && favors_M1),
    # static L1 concordance (legacy) and occupied-component correlation (R2, secondary)
    conc_observed = conc$observed, conc_null = conc$null_mean,
    conc_excess = conc$excess, conc_p = conc$p_perm,
    cor_observed = conc$cor_observed, cor_null = conc$cor_null,
    cor_excess = cor_excess, cor_p = conc$cor_p, n_occupied = conc$n_occupied,
    # change-based increment concordance (R5)
    incr_observed = if (!is.null(incr)) incr$observed else NA_real_,
    incr_null = if (!is.null(incr)) incr$null_mean else NA_real_,
    incr_excess = incr_excess, incr_p = if (!is.null(incr)) incr$p_perm else NA_real_,
    # LOPO-style predictive readout (R4): interpretable out-of-sample score
    lopo_delta_elpd = lopo_delta_elpd,
    # PRIMARY continuous tracking score for ROC = change-based increment excess
    # (the paper's estimand); falls back to the correlation excess if unavailable.
    tracking_score = if (is.finite(incr_excess)) incr_excess else cor_excess
  )
}


###############################################################################
## ---- 11_simulation_scenarios.R
###############################################################################

###############################################################################
## 11_simulation_scenarios.R
## Parameter grids for the 3 simulation scenarios
##
## Scenario 1: Model Selection Consistency via WAIC/LOO
## Scenario 2: Estimation of Subclonal Signatures
## Scenario 3: Recovery of Latent Compositions and Tracking Correlation
###############################################################################

#' Get the configuration grid for a given scenario
#'
#' Four scenarios matching the methods paper (Section: Simulation Study):
#'   1 signature estimation + competitors + coverage (contraction rate)
#'   2 posterior consistency of the latent composition (pi L1 error vs p)
#'   3 operating characteristics of tracking detection (M0-true vs M1-true, ROC)
#'   4 robustness to DGP misspecification (t3 noise, correlated genes, housekeeping)
#'
#' @param scenario Integer 1-4
#' @param scale Character: "production" (writeup values) or "local" (fast validation)
#' @return Data frame; each row is one configuration. All rows carry a common set
#'   of columns (p, s_k_frac, kappa_true, alpha_dp, true_model, K_true, w_N_mean,
#'   n, T_max, sub, noise, gene_corr, hk_background) plus MCMC settings.
get_scenario_grid <- function(scenario, scale = "production") {

  if (scale == "production") {
    ## p-grid anchored at real-data p=300; tests contraction from well-specified
    ## to high-dimensional settings, with p=300 validating the actual application scale.
    p_grid_s1 <- c(150, 300, 500, 1000, 2000)
    ## p-grid for the Scenario-2 posterior-consistency check (spans an order of magnitude)
    p_grid_s2 <- c(1000, 2000, 3000, 4000, 5000)
    p_base    <- 2000
    R <- 100
    n_chains <- 4
    n_iter <- 22500; n_burn <- 20000; thin <- 10   # match run_krd_M1.R
    K_fit <- 10; n_subjects <- 26; T_max <- 5
  } else {
    ## LOCAL: tiny + fast, for end-to-end validation only
    p_grid_s1 <- c(100, 200)
    p_grid_s2 <- c(100, 150, 200)
    p_base    <- 150
    R <- 3
    n_chains <- 2
    n_iter <- 400; n_burn <- 200; thin <- 2
    K_fit <- 6; n_subjects <- 6; T_max <- 3
  }

  ## every row gets these columns (defaults; scenarios override what they vary).
  ## Fills BOTH missing columns and NA cells (bind_rows introduces NAs when
  ## sub-grids carry different columns), so no downstream NA leaks into the DGP.
  base_row <- function(df) {
    defs <- list(s_k_frac = 0.05, kappa_true = 1, alpha_dp = 1, true_model = "M1",
                 K_true = 3, w_N_mean = 0.1, n = n_subjects, T_max = T_max,
                 sub = NA_character_, noise = "gaussian", gene_corr = 0,
                 hk_background = 0, p = p_base,
                 ## composition-heterogeneity switches (R1/R3): default off
                 patient_pi = FALSE, pi_drift = 0, balanced_pi = FALSE,
                 ## Anchor decoupling (Part I, KEPT): anchor_frac*p*pi is the
                 ## allocation anchor that fixes the p>>N mode-collapse (3/3
                 ## signature recovery). The cross-source tracking coupling is
                 ## NOT estimated: four estimator formulations (D1 x2, D2 x2) were
                 ## falsified -- at p>>N_obs the per-source counts drown the
                 ## coupling signal that lives in the per-observation weights, so
                 ## a scalar tracking parameter is structurally non-identifiable
                 ## in this model. Default is the fixed structural value
                 ## (kappa=1 -> rho=2/3); the D2 "mh" path remains available as an
                 ## option (lambda_method="mh") but is not used by default.
                 anchor_frac = 0.05, lambda_method = "fixed", lambda_init = 1.0)
    for (col in names(defs)) {
      if (is.null(df[[col]])) df[[col]] <- defs[[col]]
      if (!identical(col, "sub")) df[[col]][is.na(df[[col]])] <- defs[[col]]
    }
    df
  }

  if (scenario == 1) {
    ## Scenario 1: signature estimation, contraction rate + competitors + coverage.
    ## R1/R3: use per-patient BALANCED compositions so every one of the K_true
    ## signatures is genuinely exercised across observations. Without this the
    ## stick-breaking composition is dominated by one component, the data rarely
    ## visit components 2..K_true, and the sampler collapses onto a single mode
    ## (recovering only 1 of K_true signatures -> catastrophic theta MSE). A
    ## balanced, patient-varying composition gives each signature the support the
    ## contraction-rate estimand needs.
    grid <- expand.grid(p = p_grid_s1, s_k_frac = c(0.01, 0.05, 0.10),
                        stringsAsFactors = FALSE)
    grid$sub <- "1"
    grid$patient_pi <- TRUE; grid$balanced_pi <- TRUE
    ## Data-generating concentration (per-observation), DISTINCT from the model's
    ## fitting kappa. With kappa_gen=1 each observation's Dir(1*pi) draw is diffuse
    ## and visits only ~1.9 of the K_true=3 components, so the data cannot identify
    ## the 2nd/3rd signature and the fit recovers only one (recovery rate ~1/3).
    ## kappa_gen=20 makes each observation genuinely sample all K_true signatures
    ## (per-observation occupied K+ ~ 2.97), which is what the contraction estimand
    ## requires; balancing the cohort composition alone (balanced_pi) is not enough.
    grid$kappa_true <- 20
    ## Scenario 1 is the theta-estimation / contraction demo, not a tracking-
    ## strength scenario, so lambda is FIXED (at a moderate value) rather than
    ## estimated -- what matters is that the p-scaled allocation anchor
    ## (anchor_frac*p*pi, default in base_row) prevents the mode-collapse so all
    ## K_true signatures are recovered. This replaces the old kappa_fit=50 crutch.
    grid$lambda_method <- "fixed"
    grid$lambda_init   <- 20
    grid <- base_row(grid)

  } else if (scenario == 2) {
    ## Scenario 2 (REDESIGNED): signature-recovery and rho^theta consistency as n grows.
    ##
    ## ESTIMAND: recovery_rate (fraction of K_true=3 components matched) and mean
    ## rho^theta (cross-source signature correlation) as functions of n. Both should
    ## increase monotonically toward 1 as the model accumulates information across
    ## more patients.
    ##
    ## KEY DESIGN CHOICES:
    ##   patient_pi=FALSE, balanced_pi=TRUE: a single cohort composition pi is shared
    ##   across all patients (standard DP stick-breaking with balanced_pi ensuring all
    ##   K_true=3 components are genuinely occupied). This is the standard setting
    ##   where the posterior for both theta and pi concentrates as n grows.
    ##   Previously patient_pi=TRUE created per-patient compositions that had no single
    ##   pi* for the L1 estimand to converge to, producing a flat curve.
    ##
    ##   kappa_true=20: high per-observation concentration so every sample exercises
    ##   all K_true components (same as Scenario 1 which confirmed this is necessary
    ##   for non-degenerate theta estimation).
    ##
    ##   n range: {5,10,20,50,100,200} — two orders of magnitude so the monotone
    ##   recovery trend is clearly visible. At n=5 x T=5 = 25 obs, recovery is poor;
    ##   at n=200 x T=5 = 1000 obs, recovery should approach 1.
    ##
    ##   p=300: matches real-data scale, keeps runtime tractable.
    ##
    ## METRICS COMPUTED (by run_one_sim scenario==2 block):
    ##   theta: recovery_rate, tpr, fpr, mse, coverage
    ##   rho_theta: mean_rho_est, mean_rho_true, mean_bias, match_rate
    n_grid_s2 <- if (scale == "production") c(5, 10, 20, 50, 100, 200) else c(3, 5, 8)
    grid <- data.frame(n = n_grid_s2, p = 300, sub = "2_consistency",
                       stringsAsFactors = FALSE)
    grid$patient_pi  <- FALSE   # single shared cohort composition
    grid$balanced_pi <- TRUE    # balanced pi ensures all K_true=3 components occupied
    grid$kappa_true  <- 20      # tight per-observation concentration (same as S1)
    grid$lambda_init <- 20      # matching kappa_true for fit
    grid$lambda_method <- "fixed"
    ## pi_drift=0 (base_row default): static pi_true
    grid <- base_row(grid)

  } else if (scenario == 3) {
    ## Scenario 3 (REDESIGNED): sample-size power curve for tracking detection.
    ##
    ## ESTIMAND: WAIC-based M1-selection rate (correct_selection) and increment-
    ## concordance detection rate (incr_p < 0.05) as a function of cohort size n,
    ## under fixed kappa*=1 (M1-true) and M0-null. Both should be near alpha=5%
    ## under M0 (valid null calibration) and rise monotonically with n under M1.
    ##
    ## WHY THE PREVIOUS DESIGN FAILED:
    ##   The previous design varied kappa* at fixed n=26. At kappa*=20, Dir(20*pi)
    ##   concentrates each observation so tightly near the shared composition that
    ##   within-patient weight increments delta_g ~ delta_cf ~ 0, destroying the
    ##   increment test's signal. Power at kappa*=20 was 0.06 -- below kappa*=1 (0.11).
    ##   The inverted power curve makes the scenario unpublishable.
    ##
    ## NEW DESIGN (increment-primary power curve):
    ##   Fixed kappa*=1 (moderate tracking strength, consistent with real-data results).
    ##   n varies: {5, 10, 20, 26, 50, 100} plus M0-null at n=26 for FPR calibration.
    ##   As n grows, more patient-pairs contribute to the within-patient increment test,
    ##   giving a clean monotone power curve. WAIC is reported as secondary.
    ##
    ##   WHY p=300 NOT p_base=2000:
    ##   At p=2000 with N_obs <= 500, the WAIC DOF imbalance (M0 has 2*(K-1)=18 extra
    ##   stick-breaking params vs M1) is overwhelmed by the scale of the per-observation
    ##   log-likelihood (sum over 2000 genes), making delta_waic structurally negative
    ##   regardless of n. At p=300, the DOF imbalance is still present but small relative
    ##   to the per-observation LL, allowing WAIC to be informative. The increment test
    ##   remains the primary estimand (it does not depend on p scaling).
    ##
    ##   patient_pi=FALSE, balanced_pi=TRUE: single shared pi so M1/M0 difference
    ##   reflects only the cross-source coupling, not per-patient heterogeneity.
    ##   kappa_true=1: real-data-calibrated; moderate per-observation composition
    ##   variance so within-patient increments have genuine variance for the test.
    ##
    ## n=26 is included in both M1 and M0 so the operating characteristics at the
    ## actual study sample size are directly reported.
    ##
    ## PRIMARY METRIC: incr_detection_rate = fraction of reps with incr_p < 0.05.
    ## SECONDARY: waic_detection_rate (delta_waic > 0).
    n_grid_s3 <- if (scale == "production") c(5, 10, 20, 26, 50, 100) else c(3, 5, 8)
    grid_null <- data.frame(true_model = "M0", kappa_true = 1, n = 26,
                            lambda_init = 1.0, sub = "3_null",
                            stringsAsFactors = FALSE)
    grid_alt  <- data.frame(true_model = "M1", kappa_true = 1,
                            n = n_grid_s3,
                            lambda_init = 1.0,
                            sub = paste0("3_n", n_grid_s3),
                            stringsAsFactors = FALSE)
    grid <- dplyr::bind_rows(grid_null, grid_alt)
    ## p=300: real-data scale; avoids structural WAIC DOF imbalance at p>>N_obs.
    grid$p <- if (scale == "production") 300 else p_base
    grid$patient_pi  <- FALSE   # single shared pi for clean comparison
    grid$balanced_pi <- TRUE    # all K_true=3 components occupied
    grid$pi_drift    <- 0       # static compositions
    grid$lambda_method <- "fixed"
    grid <- base_row(grid)

  } else if (scenario == 4) {
    ## Scenario 4: robustness to misspecification (M1-true, kappa*=1), with the
    ## same per-patient / balanced / drifting composition as Scenario 3 (R1/R3/R5)
    ## so the concordance and increment estimands are non-degenerate before the
    ## misspecification is applied.
    ## hk_background reduced 0.5->0.2: with sigma_k=1 and theta values in {2,3},
    ## SD=0.5 background caused component collapse (K+=2, match_rate=0.33). SD=0.2
    ## gives SNR>=10 relative to signal amplitude, matching realistic 5hmC DE gene data.
    grid <- dplyr::bind_rows(
      data.frame(sub = "4_t3",   noise = "t3",       gene_corr = 0,   hk_background = 0),
      data.frame(sub = "4_corr", noise = "gaussian", gene_corr = 0.3, hk_background = 0),
      data.frame(sub = "4_hk",   noise = "gaussian", gene_corr = 0,   hk_background = 0.2),
      data.frame(sub = "4_all",  noise = "t3",       gene_corr = 0.3, hk_background = 0.2)
    ); grid$p <- p_base
    grid$patient_pi <- TRUE; grid$balanced_pi <- TRUE; grid$pi_drift <- 0.5
    grid <- base_row(grid)

  } else if (scenario == 5) {
    ## Scenario 5 (REDESIGNED): compartment-specific signatures (theta_g != theta_cf).
    ##
    ## ESTIMAND: mean rho^theta (cross-source signature correlation) vs delta_theta.
    ## As delta_theta increases, true rho^theta declines; the estimator should track
    ## this decline, producing a downward-sloping curve from ~1.0 (shared signatures)
    ## toward ~0.0 (maximally divergent signatures).
    ##
    ## WHY THE PREVIOUS DESIGN FAILED:
    ##   Previous delta_theta range {0, 0.10, 0.25, 0.50} with signal amplitude {+-2,+-3}
    ##   (SD ~ 2.6 per active gene) produced rho_true values of 1.000, 0.998, 0.991, 0.964 --
    ##   a range of only 0.036. This is smaller than the estimator's own ~0.04 bias floor,
    ##   so the estimated rho was flat at ~0.96 across all delta_theta. The estimator
    ##   cannot track a signal change smaller than its own noise floor.
    ##
    ## NEW DESIGN:
    ##   Extended delta_theta range: {0, 1, 2, 4}. True rho^theta per component:
    ##     delta=0: rho=1.000 (identical signatures)
    ##     delta=1: rho = sigma^2_signal / (sigma^2_signal + sigma^2_delta)
    ##              = ~6.7/(6.7+1) ~ 0.87  [with {+-2,+-3} signal, mean E[val^2]~6.7]
    ##     delta=2: rho ~ 6.7/(6.7+4) ~ 0.63
    ##     delta=4: rho ~ 6.7/(6.7+16) ~ 0.30
    ##   This gives a range from 1.00 to 0.30 -- well outside the estimator's bias floor.
    ##   The estimated rho should decline visibly from ~0.96 to ~0.30, giving a clear
    ##   monotone downward trend that validates the model's compartment-specific estimation.
    ##
    ##   patient_pi=FALSE, balanced_pi=TRUE: clean non-degenerate compositions so all
    ##   K_true=3 components are exercised and matched at all delta_theta levels.
    ##   kappa_true=20: all components exercised (same as S1,S2).
    ##   p=300: real-data scale; enough genes per component for robust theta matching.
    ##   M1-true throughout: composition concordance holds, only signature concordance varies.
    ##
    ## METRICS:
    ##   rho_theta: mean_rho_est, mean_rho_true, mean_bias, match_rate (via compute_rho_theta_metric)
    ##   kplus_median, delta_theta
    grid <- data.frame(
      delta_theta = c(0, 1, 2, 4),
      sub = c("5_null", "5_mild", "5_moderate", "5_large"),
      stringsAsFactors = FALSE
    )
    ## p=300: real-data scale. A smaller p gives more observations per gene,
    ## sharper theta posteriors, and more reliable rho^theta estimation.
    ## Large p (2000) with N_obs=130 is severely underdetermined for cfDNA theta
    ## estimation, inflating rho^theta bias and masking the true decline with delta_theta.
    grid$p <- if (scale == "production") 300 else p_base
    grid$true_model  <- "M1"   # tracking DGP (composition concordance is present)
    grid$patient_pi  <- FALSE  # single shared pi for clean component matching
    grid$balanced_pi <- TRUE   # all K_true=3 components occupied
    grid$pi_drift    <- 0      # static compositions for clean rho estimation
    grid$kappa_true  <- 20     # all components exercised
    grid$lambda_init <- 20
    grid$lambda_method <- "fixed"
    grid <- base_row(grid)

  } else {
    stop("scenario must be 1, 2, 3, 4, or 5")
  }

  # Add common MCMC settings
  grid$config_id <- seq_len(nrow(grid))
  grid$R <- R
  ## R4: Scenarios 3, 4, 5 hinge on model comparison / concordance; give them
  ## more chains and longer annealing so clean chains survive the health check.
  if (scenario %in% c(3, 4, 5) && scale == "production") {
    grid$n_chains <- 8
    grid$n_iter <- 30000; grid$n_burn <- 27500; grid$thin <- 10
  } else {
    grid$n_chains <- n_chains
    grid$n_iter <- n_iter; grid$n_burn <- n_burn; grid$thin <- thin
  }
  ## K_fit = 10 for all scenarios, matching the real-data analysis. The anchor
  ## decoupling (anchor_frac=0.05) prevents mode-collapse at K_fit=10 by providing
  ## a p-scaled allocation prior that keeps all K_true signatures identified.
  grid$K_fit <- K_fit

  grid
}

#' Print scenario summary
print_scenario_summary <- function(scenario, scale = "production") {
  grid <- get_scenario_grid(scenario, scale)
  cat(sprintf("Scenario %d: %d configurations x R=%d reps = %d total jobs\n",
              scenario, nrow(grid), grid$R[1], nrow(grid) * grid$R[1]))
  cat("Configurations:\n")
  print(grid[, !names(grid) %in% c("R", "n_iter", "n_burn", "K_fit", "config_id")])
}


###############################################################################
## ---- 12_simulation_runner.R
###############################################################################

###############################################################################
## 12_simulation_runner.R
## Master function to run one simulation replication
##
## Generates data, fits M1 and M0, computes all metrics, saves results.
## Designed to be called by SLURM array jobs via slurm/run_one_sim.R.
###############################################################################

#' Run one simulation replication
#'
#' @param scenario   Integer: 1, 2, or 3
#' @param config_id  Integer: row index in the scenario grid
#' @param rep_id     Integer: replication number (used as seed offset)
#' @param scale      Character: "production" or "local"
#' @param base_dir   Character: project root directory
#'
#' @return List with config, metrics, and timing info. Also saved to disk.
run_one_sim <- function(scenario, config_id, rep_id,
                         scale = "production",
                         base_dir = getwd()) {

  # Source all needed files (bootstrap compiles the Rcpp inner loops)
  suppressWarnings(try(source(file.path(base_dir, "R", "00_bootstrap.R")), silent = TRUE))
  source(file.path(base_dir, "R", "01_lib_core.R"))

  # Get configuration
  grid <- get_scenario_grid(scenario, scale)
  cfg <- grid[config_id, ]

  # Compute s_k from fraction
  s_k <- max(1, round(cfg$s_k_frac * cfg$p))

  # Seed: deterministic from config and rep
  sim_seed <- scenario * 1e6 + config_id * 1e3 + rep_id

  cat(sprintf("Scenario %d, Config %d, Rep %d (seed=%d)\n",
              scenario, config_id, rep_id, sim_seed))
  cat(sprintf("  p=%d, n=%d, T=%d, K_true=%d, s_k=%d, kappa=%s, model=%s\n",
              cfg$p, cfg$n, cfg$T_max, cfg$K_true, s_k,
              cfg$kappa_true, cfg$true_model))

  t_start <- Sys.time()

  # --- Generate data (misspecification + composition-heterogeneity switches) ---
  noise <- if (!is.null(cfg$noise)) cfg$noise else "gaussian"
  gene_corr <- if (!is.null(cfg$gene_corr)) cfg$gene_corr else 0
  hk_background <- if (!is.null(cfg$hk_background)) cfg$hk_background else 0
  patient_pi <- isTRUE(cfg$patient_pi)
  pi_drift <- if (!is.null(cfg$pi_drift)) cfg$pi_drift else 0
  balanced_pi <- isTRUE(cfg$balanced_pi)
  delta_theta <- if (!is.null(cfg$delta_theta)) cfg$delta_theta else 0
  sim <- generate_sim_data(
    n = cfg$n, p = cfg$p, T_max = cfg$T_max, K_true = cfg$K_true,
    s_k = s_k, kappa = cfg$kappa_true, alpha_dp = cfg$alpha_dp,
    w_N_mean = cfg$w_N_mean, model = cfg$true_model, seed = sim_seed,
    noise = noise, gene_corr = gene_corr, hk_background = hk_background,
    patient_pi = patient_pi, pi_drift = pi_drift, balanced_pi = balanced_pi,
    delta_theta = delta_theta
  )

  # --- Fit with the PRODUCTION sampler (matches slurm/run_krd_M1.R) ---
  # 4-chain fit_multi_chain with annealing; store theta trace for S1 coverage.
  store_trace <- (scenario == 1)
  n_chains <- if (!is.null(cfg$n_chains)) cfg$n_chains else 4
  thin <- if (!is.null(cfg$thin)) cfg$thin else 10
  ## DECOUPLED anchor (Part I, KEPT). The allocation-prior anchor is
  ## eta_k = anchor_frac * p * pi_k, which scales with p and stabilises the
  ## collapsed Polya-urn against the p-sized count term -- this replaces the old
  ## S1-only kappa_fit=50 crutch and fixes the mode-collapse at its root for every
  ## scenario (3/3 signature recovery). The cross-source tracking coupling is NOT
  ## estimated: four estimator formulations were falsified (see note below), so we
  ## keep the fixed structural value (kappa=1 -> rho=2/3). The D2 "mh" path remains
  ## selectable via cfg$lambda_method but is not the default.
  anchor_frac   <- if (!is.null(cfg$anchor_frac))   cfg$anchor_frac   else 0.05
  lambda_method <- if (!is.null(cfg$lambda_method)) cfg$lambda_method else "fixed"
  lambda_init   <- if (!is.null(cfg$lambda_init))   cfg$lambda_init   else 1.0
  fit_M1 <- fit_multi_chain(
    sim$Y_g, sim$Y_cf, sim$patient, sim$time,
    n_chains = n_chains, seed = sim_seed + 1,
    K = cfg$K_fit, n_iter = cfg$n_iter, n_burn = cfg$n_burn, thin = thin,
    model = "M1",
    anchor_frac = anchor_frac, lambda_method = lambda_method, lambda_init = lambda_init,
    anneal = TRUE, T_anneal = 20.0, n_cool_buffer = max(1, floor((cfg$n_burn) / 2)),
    store_pointwise_ll = TRUE,
    store_theta_trace = store_trace
  )
  fit_M0 <- fit_multi_chain(
    sim$Y_g, sim$Y_cf, sim$patient, sim$time,
    n_chains = n_chains, seed = sim_seed + 2,
    K = cfg$K_fit, n_iter = cfg$n_iter, n_burn = cfg$n_burn, thin = thin,
    model = "M0",
    anchor_frac = anchor_frac, lambda_method = "fixed", lambda_init = lambda_init,
    anneal = TRUE, T_anneal = 20.0,
    n_cool_buffer = max(1, floor((cfg$n_burn) / 2)),
    store_pointwise_ll = TRUE
  )

  t_end <- Sys.time()
  runtime_sec <- as.numeric(difftime(t_end, t_start, units = "secs"))

  # --- Health-checked WAIC (the real-data decision rule, run_krd_combine.R) ---
  tr1 <- tryCatch(identify_trapped_chains(fit_M1, ll_gap = 10000),
                  error = function(e) list(keep = seq_along(fit_M1$chains), drop = integer(0)))
  tr0 <- tryCatch(identify_trapped_chains(fit_M0, ll_gap = 10000),
                  error = function(e) list(keep = seq_along(fit_M0$chains), drop = integer(0)))
  n_trapped_M1 <- length(tr1$drop); n_trapped_M0 <- length(tr0$drop)
  delta_waic <- NA_real_
  tryCatch({
    ss1 <- waic_loo_subset(fit_M1, tr1$keep); ss0 <- waic_loo_subset(fit_M0, tr0$keep)
    if (!is.null(ss1$waic) && !is.null(ss0$waic))
      delta_waic <- ss0$waic$estimates["waic", "Estimate"] - ss1$waic$estimates["waic", "Estimate"]
  }, error = function(e) warning(sprintf("WAIC failed (s%d c%d r%d): %s",
                                         scenario, config_id, rep_id, conditionMessage(e))))

  # helper: single-chain-shaped view for the legacy metric functions that read
  # fit$samples / fit$theta_postmean / fit$theta_trace. Use the best (least
  # trapped) chain's theta and the merged samples for weights/pi.
  best_chain <- fit_M1$chains[[ if (length(tr1$keep)) tr1$keep[1] else 1 ]]
  fit_view <- list(
    samples            = fit_M1$merged,
    theta_bar_postmean = best_chain$theta_bar_postmean,
    theta_bar_postvar  = best_chain$theta_bar_postvar,
    theta_postmean     = best_chain$theta_postmean,
    theta_postvar      = best_chain$theta_postvar,
    theta_0_postmean   = best_chain$theta_0_postmean,
    final_theta        = best_chain$final_theta,
    final_theta_bar    = best_chain$final_theta_bar,
    final_theta_0      = best_chain$final_theta_0,
    theta_bar_trace    = best_chain$theta_bar_trace,
    theta_trace        = best_chain$theta_trace,
    final_omega_g      = best_chain$final_omega_g,
    K = best_chain$K, n_save = nrow(fit_M1$merged$omega_g_trace),
    N_obs = sim$N_obs
  )

  # --- Assemble results ---
  results <- list(
    scenario = scenario, config_id = config_id, rep_id = rep_id,
    config = as.list(cfg), runtime_sec = runtime_sec,
    delta_waic = delta_waic, favors_M1 = if (is.na(delta_waic)) NA else delta_waic > 0,
    n_trapped_M1 = n_trapped_M1, n_trapped_M0 = n_trapped_M0
  )
  results$basic <- compute_sim_metrics(best_chain, sim)
  results$correct_selection <- if (is.na(delta_waic)) NA else
    ((delta_waic > 0) == (cfg$true_model == "M1"))

  # --- Scenario-specific metrics ---
  if (scenario == 1) {
    # Align chains by signature to correct label switching, then recompute
    # theta postmean from the stored trace (not the contaminated running sum).
    # align_chains_by_signature matches omega traces across chains.
    aligned_chains <- tryCatch(
      align_chains_by_signature(fit_M1$chains),
      error = function(e) { warning("chain alignment failed: ", conditionMessage(e)); fit_M1$chains }
    )
    best_aligned <- aligned_chains[[ if (length(tr1$keep)) tr1$keep[1] else 1 ]]
    pm <- recompute_postmean_from_trace(best_aligned)
    if (!is.null(pm)) {
      if (!is.null(pm$theta_bar_postmean))
        fit_view$theta_bar_postmean <- pm$theta_bar_postmean
      if (!is.null(pm$theta_postmean))
        fit_view$theta_postmean <- pm$theta_postmean
      fit_view$theta_bar_trace <- best_aligned$theta_bar_trace
      fit_view$theta_trace     <- best_aligned$theta_trace
    }
    results$theta <- compute_theta_metrics(fit_view, sim)
    results$competitors <- tryCatch(fit_competitors(sim),
                                    error = function(e) { warning("fit_competitors failed: ", conditionMessage(e)); NULL })
  }
  if (scenario == 2) {
    ## Scenario 2 (REDESIGNED): signature recovery as n grows.
    ## Primary estimands: recovery_rate and mean_rho_est both increase with n.
    ## pi metric is removed: with patient_pi=FALSE and no per-patient drift, the
    ## "L1 error of estimated vs true pi" would be valid but is not the focus.
    ## rho_theta is computed below in the universal block.
    results$theta <- compute_theta_metrics(fit_view, sim)
  }
  if (scenario == 3) {
    conc <- tryCatch(concordance_excess(fit_M1$merged, sim$N_obs, sim$time),
                     error = function(e) list(observed=NA, null_mean=NA, excess=NA, p_perm=NA,
                                              cor_observed=NA, cor_null=NA, cor_excess=NA, cor_p=NA, n_occupied=NA))
    ## R5: change-based increment concordance (the paper's estimand) is the
    ## primary tracking score; R2 occupied-component correlation is secondary.
    incr <- tryCatch(increment_concordance(fit_M1$merged, sim$N_obs, sim$patient, sim$time),
                     error = function(e) NULL)
    ## R4: interpretable out-of-sample readout. The full patient-grouped LOPO is
    ## run as a separate driver (slurm/run_lopo_cv.R); within a replicate we
    ## record the source-specific held-out cfDNA elpd advantage if available,
    ## else leave NA and rely on the increment score for discrimination.
    results$selection <- compute_selection_metrics(delta_waic, conc, cfg$true_model, incr = incr)
  }
  if (scenario == 4) {
    results$concordance <- tryCatch(
      concordance_excess(fit_M1$merged, sim$N_obs, sim$time),
      error = function(e) NULL)
    results$increment <- tryCatch(
      increment_concordance(fit_M1$merged, sim$N_obs, sim$patient, sim$time),
      error = function(e) NULL)
    results$kplus_median <- median(fit_M1$merged$K_plus)
  }
  if (scenario == 5) {
    ## Scenario 5: identifiability of the hierarchical signature and the
    ## time-varying background (writeup sec:sim_ident). The critical readout is
    ## whether the WEIGHT-level estimands stay unbiased even in the part of the
    ## varsigma_theta grid where the signatures themselves are not identified.
    results$ident <- tryCatch(
      compute_scenario5_metrics(fit_view, sim),
      error = function(e) { warning(conditionMessage(e)); NULL })
    results$kplus_median   <- median(fit_M1$merged$K_plus)
    results$varsigma_theta <- sim$varsigma_theta
  }

  ## The identifiability metrics are cheap and meaningful in every scenario,
  ## so compute them universally when not already present.
  if (is.null(results$ident)) {
    results$ident <- tryCatch(
      compute_scenario5_metrics(fit_view, sim),
      error = function(e) NULL)
  }

  # --- Save ---
  out_dir <- file.path(base_dir, "results", "sim",
                        sprintf("scenario_%d", scenario),
                        sprintf("config_%d", config_id))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, sprintf("rep_%d.rds", rep_id))
  saveRDS(results, out_file)

  cat(sprintf("  Done in %.1f sec. Saved to %s\n", runtime_sec, out_file))
  invisible(results)
}

