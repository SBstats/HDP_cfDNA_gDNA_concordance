
###############################################################################
## ---- 01_mcmc_sampler.R
###############################################################################

###############################################################################
## 01_mcmc_sampler.R
## Canonical MCMC sampler for the Bayesian Hierarchical DP Mixture Model
##
## Implements the 6-step blocked Gibbs sampler from Algorithm 1 of the
## manuscript, with the truncated stick-breaking representation (K components).
## Both alpha and kappa are FIXED (alpha=1, kappa=1 by default).
##
## theta_{jk} is TIME-INVARIANT: [p x K] array — same subclonal signature
## across all treatment cycles. Longitudinal dynamics are absorbed entirely
## by the time-varying mixture weights omega_{ikt}. The horseshoe local
## shrinkage lambda_{jk} is shared across both sources.
##
## Hyperparameters match manuscript exactly:
##   K=8, a_kappa=2, b_kappa=0.2 (E[kappa]=10), a_N=1, b_N=9
##
## Features:
##   - kappa fixed at 1 (default); adjustable via kappa_fixed parameter
##   - Vectorized horseshoe updates
##   - Full log-likelihood computation (no gene subsampling)
##   - Pointwise log-likelihood storage for WAIC
##   - Clamping to prevent numerical underflow/overflow
##   - Rcpp acceleration for inner loops (use_rcpp=TRUE, default)
##
## When use_rcpp=TRUE, the three bottleneck functions (compute_allocations,
## update_theta, compute_loglik) are replaced with C++ implementations
## from R/cpp/sampler_core.cpp, giving ~50-100x speedup.
## Pure R fallbacks are retained with _R suffix for debugging.
###############################################################################

# Auto-install required packages
.required_pkgs <- c("coda", "Rcpp", "loo", "posterior")
for (.pkg in .required_pkgs) {
  if (!requireNamespace(.pkg, quietly = TRUE)) {
    message("Installing missing package: ", .pkg)
    install.packages(.pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
}
rm(.required_pkgs, .pkg)

suppressPackageStartupMessages({
  library(coda)
  library(Rcpp)
})

# Source Rcpp functions (compiled on first load).
#
# IMPORTANT: On some platforms (notably macOS), calling sourceCpp() from
# inside a source() frame triggers "C stack usage too close to the limit"
# errors. The recommended pattern is to bootstrap Rcpp at the TOP LEVEL of
# a fresh R session BEFORE sourcing this library:
#
#     source("R/00_bootstrap.R")
#     source("R/01_lib_core.R")
#
# If 00_bootstrap.R already set .rcpp_available, we skip the in-source
# attempt. Otherwise we try here and fall back to pure R on failure.
if (!exists(".rcpp_available", envir = .GlobalEnv)) {
  .rcpp_available <- FALSE
  .cpp_file <- file.path(getwd(), "R", "cpp", "sampler_core.cpp")
  if (file.exists(.cpp_file)) {
    tryCatch({
      sourceCpp(.cpp_file)
      .rcpp_available <- TRUE
    }, error = function(e) {
      message("Rcpp compilation failed: ", conditionMessage(e),
              "\nFalling back to pure R. ",
              "Run source('R/00_bootstrap.R') BEFORE this library to avoid ",
              "the macOS C-stack issue with sourceCpp inside source().")
    })
  }
  rm(.cpp_file)
}

# ===========================================================================
# MAIN SAMPLER FUNCTION
# ===========================================================================

#' Fit the DP Mixture Concordance Model (M1 or M0)
#'
#' @param Y_g     Matrix [p x N_obs]: normalized gDNA signals
#' @param Y_cf    Matrix [p x N_obs]: normalized cfDNA signals
#' @param patient Integer vector [N_obs]: patient index for each observation
#' @param time    Integer vector [N_obs]: time index (1,...,T_max) for each obs
#' @param K       Integer: truncation level (default 30, matching manuscript)
#' @param n_iter  Integer: total MCMC iterations
#' @param n_burn  Integer: burn-in iterations
#' @param model   Character: "M1" (tracking) or "M0" (non-tracking)
#' @param seed    Integer: random seed
#' @param hyperparams List of hyperparameters (see defaults)
#' @param thin    Integer: thinning interval
#' @param verbose Logical: print progress
#' @param kappa_method Character: "grid" (recommended) or "adaptive_mh"
#' @param store_pointwise_ll Logical: store pointwise log-lik for WAIC (memory!)
#'
#' @return List with posterior samples and diagnostics
fit_dp_concordance <- function(Y_g, Y_cf, patient, time,
                                K = 10,
                                n_iter = 15000,
                                n_burn = 10000,
                                model = "M1",
                                seed = 42,
                                hyperparams = list(),
                                thin = 20,
                                verbose = TRUE,
                                kappa_fixed = 1.0,
                                ## --- decoupled anchor / tracking (Part I) ---
                                ## anchor_frac: allocation-prior pseudocount is
                                ##   eta_k = anchor_frac * p * pi_k  (role: stabilise
                                ##   the collapsed Polya-urn against the p-sized count
                                ##   term, INDEPENDENT of tracking strength). When NULL,
                                ##   the legacy anchor eta_k = kappa_fixed * pi_k is used
                                ##   (byte-identical back-compat path).
                                anchor_frac = NULL,
                                ## lambda: the tracking/dispersion scalar in the ω
                                ##   reconstruction Dir(lambda*pi + n) and the tracking
                                ##   correlation (lambda+1)/(alpha+lambda+1). "fixed"
                                ##   holds it at lambda_init; "mh" estimates it (M1 only)
                                ##   by log-normal random-walk Metropolis-Hastings.
                                lambda_method = c("fixed", "mh"),
                                lambda_init = NULL,
                                store_pointwise_ll = FALSE,
                                use_rcpp = TRUE,
                                anneal = TRUE,
                                T_anneal = 20.0,
                                n_cool_buffer = 10000,
                                fixed_temp = NULL,
                                store_theta_trace = FALSE,
                                theta_trace_thin = 5,
                                ## Re-anchor components that hold zero allocations, during
                                ## burn-in only. Prevents the prior-only random walk described
                                ## at the re-anchor block below. Default ON: the drift is a
                                ## defect rather than a modelling choice. Set FALSE to
                                ## reproduce pre-fix behaviour exactly.
                                empty_reanchor = TRUE,
                                ## Annealing shape. "geometric" is linear in log T (i.e. beta
                                ## grows exponentially) and spends far more of the burn-in in
                                ## the critical T ~ 1-2 window where the mixture commits to a
                                ## K+. "linear" is the historical schedule.
                                anneal_schedule = c("linear", "geometric"),
                                init = NULL) {

  set.seed(seed)

  # --- Rcpp toggle ---
  if (use_rcpp && !exists(".rcpp_available")) use_rcpp <- FALSE
  if (use_rcpp && !.rcpp_available) {
    if (verbose) cat("  Rcpp not available, falling back to pure R\n")
    use_rcpp <- FALSE
  }

  # --- Input validation ---
  stopifnot(is.matrix(Y_g), is.matrix(Y_cf))
  stopifnot(nrow(Y_g) == nrow(Y_cf))
  stopifnot(ncol(Y_g) == ncol(Y_cf))
  stopifnot(length(patient) == ncol(Y_g))
  stopifnot(length(time) == ncol(Y_g))
  stopifnot(model %in% c("M1", "M0"))
  stopifnot(K >= 2)
  stopifnot(n_iter > n_burn)

  # --- Dimensions ---
  p <- nrow(Y_g)
  N_obs <- ncol(Y_g)
  n <- max(patient)
  T_max <- max(time)

  # --- Hyperparameters (manuscript defaults) ---
  hp <- list(
    a_alpha = 1, b_alpha = 1,         # DP concentration: E[alpha]=1
    a_kappa = 2, b_kappa = 0.2,       # Tracking precision: E[kappa]=10
    a_0 = 1, b_0 = 9,                 # Contamination: E[omega_0]=0.1  (eq:prior_wN)
    a_sigma = 2,                      # Component variance (IG shape)
    a_b_sigma = 2, b_b_sigma = 1,     # Gamma hyperprior on component-variance rate
                                      # (pools sigma_k^2 across components to prevent locking)
    a_sigma0 = 2, b_sigma0 = 1,       # Background variance sigma_0^2
    # Between-subject signature dispersion varsigma_k^2 ~ IG(a_vs, b_vs).
    # Weakly informative: E[varsigma_k^2] = 2 (eq:prior_varsigma).
    a_vs = 2, b_vs = 2,
    # Background dispersion varsigma_0^2 ~ IG(a_vs0, b_vs0). DELIBERATELY
    # informative (E = 0.25): shrinking theta_0 toward a common cohort profile
    # is what stops the background absorbing a genuine subclone. See
    # writeup Sec 5.3 "Background means".
    a_vs0 = 3, b_vs0 = 0.5,
    # Prior on the cohort background mean theta_0_bar[j] ~ N(m_0, s_0^2).
    # m_0 = NULL means "use the cohort mean of the observed cfDNA signal",
    # which anchors the background where non-tumour signal actually lies
    # rather than at an arbitrary origin.
    m_0 = NULL, s_0_2 = 10,
    sigma_kappa_mh = 0.3,             # Initial MH proposal SD for log(kappa)
    kappa_target_accept = 0.25,       # Target acceptance for adaptive MH
    kappa_grid = c(0.1, 0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 4, 5, 7, 10, 12, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200)
  )
  hp[names(hyperparams)] <- hyperparams
  # Back-compat: accept the old contamination/background hyperparameter names.
  if (!is.null(hyperparams$a_N))      hp$a_0      <- hyperparams$a_N
  if (!is.null(hyperparams$b_N))      hp$b_0      <- hyperparams$b_N
  if (!is.null(hyperparams$a_sigmaN)) hp$a_sigma0 <- hyperparams$a_sigmaN
  if (!is.null(hyperparams$b_sigmaN)) hp$b_sigma0 <- hyperparams$b_sigmaN
  if (is.null(hp$m_0)) hp$m_0 <- mean(Y_cf, na.rm = TRUE)

  # --- Storage for posterior samples ---
  n_save <- floor((n_iter - n_burn) / thin)
  K_trace <- K  # store omega traces for ALL K components (needed for posterior concordance metrics)

  samples <- list(
    kappa    = numeric(n_save),   # legacy name; now stores the tracking scalar (phi)
    lambda   = numeric(n_save),   # == kappa slot (back-compat)
    phi      = numeric(n_save),   # D2 cross-source coupling (tracking estimand)
    alpha_dp = numeric(n_save),
    K_plus   = integer(n_save),
    ## Occupancy count (any component with >= 1 allocated gene). Diagnostic
    ## only -- K_plus above is the reported quantity and uses a weight
    ## threshold, identically for M1 and M0.
    K_plus_occupied = integer(n_save),
    # Contamination fraction, now PER OBSERVATION (omega_0[i,t]); columns follow
    # the observation ordering of Y_g/Y_cf. Use model_data$paired_info to map a
    # column back to (patient, timepoint).
    omega_0  = matrix(NA, n_save, N_obs),
    sigma_k  = matrix(NA, n_save, K),
    sigma_0  = numeric(n_save),
    # Between-subject dispersions: varsigma_k^2 (signatures), varsigma_0^2 (background)
    varsigma_k = matrix(NA, n_save, K),
    varsigma_0 = numeric(n_save),
    b_sigma  = numeric(n_save),          # shared hyperparameter (pools sigma_k^2 across k)
    tau_k    = matrix(NA, n_save, K),    # horseshoe global scale per component
    # Under M1, v is a single stick-breaking vector for the shared pi.
    # Under M0 (source-DP formulation), the shared v is not used; instead
    # v_g and v_cf parameterize the two independent source-specific DPs.
    v        = matrix(NA, n_save, K - 1),
    v_g      = matrix(NA, n_save, K - 1), # M0 only: gDNA stick-breaking
    v_cf     = matrix(NA, n_save, K - 1), # M0 only: cfDNA stick-breaking
    loglik   = numeric(n_save),
    # Mixture weight traces: [n_save x N_obs x K_trace] stored as [n_save x (N_obs*K_trace)]
    omega_g_trace  = matrix(NA, n_save, N_obs * K_trace),
    omega_cf_trace = matrix(NA, n_save, N_obs * K_trace)
  )

  # Running posterior mean and sum-of-squares for the signature hierarchy,
  # accumulated over saved iterations in canonical component order.
  #   theta     [p x K x n]  subject-level signatures (shared across sources)
  #   theta_bar [p x K]      population signatures
  #   theta_0   [p x n]      subject-level cfDNA background means
  # Online accumulation: posterior mean = sum/n, variance = (sumsq - sum^2/n)/(n-1).
  theta_sum       <- array(0, dim = c(p, K, n))
  theta_sumsq     <- array(0, dim = c(p, K, n))
  theta_bar_sum   <- array(0, dim = c(p, K))
  theta_bar_sumsq <- array(0, dim = c(p, K))
  theta_0_sum     <- array(0, dim = c(p, n))
  theta_0_sumsq   <- array(0, dim = c(p, n))
  theta_nsave     <- 0L

  # OPTIONAL thinned traces for between-chain split-Rhat audit, in canonical
  # component order. Off by default.
  #
  # MEMORY NOTE: the subject-level trace is [n_tr x p x K x n]. At this
  # analysis's dimensions (p=318, K=10, n=26, 50 draws) that is ~0.03 GB per
  # chain -- fine. At genome-wide scale (p~11,837) the same structure is
  # ~123 GB per chain. Do NOT port this setting to a full-gene pipeline
  # without restricting to occupied components or a gene subset.
  store_theta_trace <- isTRUE(store_theta_trace)
  n_theta_trace <- if (store_theta_trace) floor(n_save / theta_trace_thin) else 0L
  theta_trace     <- if (store_theta_trace && n_theta_trace > 0)
                       array(NA_real_, dim = c(n_theta_trace, p, K, n)) else NULL
  theta_bar_trace <- if (store_theta_trace && n_theta_trace > 0)
                       array(NA_real_, dim = c(n_theta_trace, p, K)) else NULL
  theta_trace_i  <- 0L

  if (store_pointwise_ll) {
    samples$pointwise_ll <- matrix(NA, n_save, 2 * N_obs)
  }

  # --- Initialize parameters ---
  if (!is.null(init)) {
    # Data-driven initialization (e.g., from K-means)
    if (verbose) cat("  Using provided initialization\n")
    v <- init$v
    v[K] <- 1
    pi_vec <- stick_break(v)
    omega_g  <- init$omega_g
    omega_cf <- init$omega_cf
    sigma_k2 <- init$sigma_k2
    sigma_0_2 <- if (!is.null(init$sigma_0_2)) init$sigma_0_2 else
                 if (!is.null(init$sigma_N2)) init$sigma_N2 else 2.0
    # Population signature: from init if given (K-means centres), else zero.
    theta_bar <- if (!is.null(init$theta_bar)) init$theta_bar else
                 if (!is.null(init$theta_g))   init$theta_g   else array(0, dim = c(p, K))  # accepts an old-format init
    # Subject-level signatures start AT the population value; the sampler
    # separates them as each subject's data warrant (writeup Remark rem:anneal).
    theta <- if (!is.null(init$theta)) init$theta else
             array(rep(as.numeric(theta_bar), n), dim = c(p, K, n))
  } else {
    v <- rep(0.5, K)
    v[K] <- 1
    pi_vec <- stick_break(v)
    omega_g  <- matrix(1/K, K, N_obs)
    omega_cf <- matrix(1/K, K, N_obs)
    sigma_k2 <- rep(1.0, K)
    sigma_0_2 <- 2.0
    theta_bar <- array(0, dim = c(p, K))
    theta     <- array(0, dim = c(p, K, n))
  }

  # Source-specific stick-breaking state (used only under M0; allocated for
  # both models so warm-starts from M0 fits roundtrip cleanly).
  if (!is.null(init$v_g)) {
    v_g <- init$v_g
  } else {
    v_g <- v  # warm-start v_g from the shared v
  }
  v_g[K] <- 1
  pi_g_vec <- stick_break(v_g)
  if (!is.null(init$v_cf)) {
    v_cf <- init$v_cf
  } else {
    v_cf <- v
  }
  v_cf[K] <- 1
  pi_cf_vec <- stick_break(v_cf)

  # omega_0, background, b_sigma, horseshoe state: pull from init if provided
  # (used by parallel tempering warm-starts), else initialize from prior means.
  #
  # omega_0 is length N_obs (per observation), NOT length n. A warm start from
  # an older per-patient wN is expanded by patient.
  omega_0 <- if (!is.null(init$omega_0)) init$omega_0
             else if (!is.null(init$wN) && length(init$wN) == n) init$wN[patient]
             else rep(hp$a_0 / (hp$a_0 + hp$b_0), N_obs)
  stopifnot(length(omega_0) == N_obs)

  # Background means. theta_0_bar[j] is the cohort background profile;
  # theta_0[j,i] is subject i's deviation from it. Initialized at the observed
  # per-gene cfDNA mean so the background starts where non-tumour signal
  # actually lies (a zero start is badly misspecified on the VST scale).
  theta_0_bar <- if (!is.null(init$theta_0_bar)) init$theta_0_bar
                 else rowMeans(Y_cf, na.rm = TRUE)
  theta_0 <- if (!is.null(init$theta_0)) init$theta_0
             else matrix(rep(theta_0_bar, n), nrow = p, ncol = n)

  # Between-subject dispersions.
  varsigma_k2 <- if (!is.null(init$varsigma_k2)) init$varsigma_k2
                 else rep(hp$b_vs / (hp$a_vs - 1), K)      # prior mean
  varsigma_0_2 <- if (!is.null(init$varsigma_0_2)) init$varsigma_0_2
                  else hp$b_vs0 / (hp$a_vs0 - 1)

  # Shared hyperparameter on the rate of sigma_k^2's Inv-Gamma prior:
  #   sigma_k^2 | b_sigma ~ Inv-Gamma(a_sigma, b_sigma)
  #   b_sigma            ~ Gamma(a_b_sigma, b_b_sigma)
  # This lets the K components borrow strength on their variance scale,
  # preventing different chains from locking onto different per-component
  # variances for sparsely-allocated components.
  b_sigma <- if (!is.null(init$b_sigma)) init$b_sigma else hp$a_b_sigma / hp$b_b_sigma

  # Horseshoe parameters
  # lambda_{jk} is shared across time (manuscript Sec 4.1)
  tau_k2   <- if (!is.null(init$tau_k2))  init$tau_k2  else rep(0.1, K)
  lambda2  <- if (!is.null(init$lambda2)) init$lambda2 else matrix(1, p, K)
  nu_aux   <- if (!is.null(init$nu_aux))  init$nu_aux  else matrix(1, p, K)
  xi_aux   <- if (!is.null(init$xi_aux))  init$xi_aux  else rep(1, K)

  alpha_dp <- 1.0
  kappa    <- kappa_fixed  # legacy scalar (anchor+dispersion fused); see below

  # -------------------------------------------------------------------
  # Decoupled anchor / tracking (Part I).
  #   anchor_mass_k = anchor_frac * p * pi_k   -> allocation Polya-urn prior
  #   lambda_track  -> dispersion of the ω reconstruction Dir(lambda*pi + n)
  #                    and the tracking correlation (lambda+1)/(alpha+lambda+1)
  # Back-compat: anchor_frac=NULL reproduces the legacy fused anchor kappa*pi.
  # -------------------------------------------------------------------
  lambda_method <- match.arg(lambda_method)
  anneal_schedule <- match.arg(anneal_schedule)
  ## Per-gene data centre, used by the empty-component re-anchor below.
  ## Loop-invariant, so computed once.
  y_ctr_reanchor <- rowMeans(cbind(Y_g, Y_cf), na.rm = TRUE)   # [p]
  lambda_track  <- if (!is.null(lambda_init)) lambda_init else kappa_fixed
  use_p_anchor  <- !is.null(anchor_frac)
  # anchor multiplier applied to pi to form eta (either c*p, or legacy kappa)
  anchor_mult   <- if (use_p_anchor) anchor_frac * p else kappa
  # MH state for lambda (M1 only)
  lambda_accept <- 0L; lambda_total <- 0L
  sigma_lambda_mh <- if (!is.null(hp$sigma_kappa_mh)) hp$sigma_kappa_mh else 0.3

  # -------------------------------------------------------------------
  # D2 cross-source coupling (Part I, phi at the pi level).
  #   Each observation (i,t) has a latent SHARED center m_it on the simplex;
  #   both sources' weights concentrate around it: omega_s_it ~ Dir(phi * m_it).
  #   phi (tracking) is identified from how tightly the PAIRED (omega_g, omega_cf)
  #   agree through their common m_it across the N_obs observations -- a genuine
  #   cross-source signal (unlike the D1 per-source dispersion, which pinned to a
  #   p/K-driven constant). m_it prior center is the cohort pi (alpha0 * pi_bar).
  #   Active only under M1 with lambda_method=="mh"; M0 is the phi->0 / independent
  #   limit and keeps the source-specific pi_g/pi_cf path unchanged.
  # -------------------------------------------------------------------
  use_d2 <- (model == "M1" && lambda_method == "mh")
  phi_track   <- lambda_track          # reuse lambda_track as the phi state/init
  phi_accept  <- 0L; phi_total <- 0L
  alpha0_m    <- if (!is.null(hp$alpha0_m)) hp$alpha0_m else 5.0  # m_it prior conc
  sigma_phi_mh <- sigma_lambda_mh
  # NOTE: under the marginalized-m D2 update, m_it is NOT a persistent state
  # variable (it is integrated out of the phi likelihood). For downstream omega
  # reconstruction we draw a per-observation shared center from its posterior
  # given the pooled counts, then draw each source's omega around it.
  a_lambda <- if (!is.null(hp$a_kappa)) hp$a_kappa else 2.0
  b_lambda <- if (!is.null(hp$b_kappa)) hp$b_kappa else 0.2

  # =====================================================================
  # MCMC LOOP
  # =====================================================================
  if (verbose) cat("Starting MCMC:", n_iter, "iterations,",
                   K, "components,", p, "genes,", T_max, "timepoints, model:", model, "\n")
  if (verbose && model == "M1") cat("  kappa:", kappa_fixed, "(fixed)\n")
  n_cool_end <- max(1, n_burn - n_cool_buffer)
  if (verbose && anneal) cat("  Annealing: T_anneal =", T_anneal, ", cooling to T=1 at iter", n_cool_end,
                              "(buffer:", n_cool_buffer, "iters at T=1 before saving)\n")

  for (iter in 1:n_iter) {

    # Temperature: either fixed (parallel tempering rung) or annealed (default).
    if (!is.null(fixed_temp)) {
      temp <- fixed_temp
    } else {
      # Annealed: cool to T=1 with buffer before burn-in ends.
      n_cool <- max(1, n_burn - n_cool_buffer)
      temp <- if (anneal && iter <= n_cool) {
        ## The LINEAR-in-T schedule races through the critical region: at
        ## T_anneal = 20 with n_cool = 10,000 it spends only ~2.5% of burn-in
        ## at T <= 2, which is where the mixture commits to a number of
        ## occupied components. The geometric schedule is linear in log T and
        ## spends ~23% of burn-in there at identical cost.
        if (identical(anneal_schedule, "geometric")) {
          T_anneal^(1 - (iter - 1) / n_cool)
        } else {
          T_anneal * (1 - (iter - 1) / n_cool) + 1
        }
      } else {
        1.0
      }
    }

    # -------------------------------------------------------------------
    # STEP 1: Update latent allocations z_g, z_cf via the COLLAPSED /
    #   Rao-Blackwellized Polya-urn predictive (Section 5 of writeup).
    #   Source-specific weights omega_g, omega_cf are integrated out;
    #   each gene's allocation is sampled given the running counts of
    #   the other genes' allocations within the same observation.
    #
    #   Under M1 (tracking, hierarchical DP coupling): both sources share
    #     the global stick-breaking weights pi, so
    #       eta_g_k = eta_cf_k = kappa * pi_k.
    #
    #   Under M0 (independent source-DP non-tracking): the two sources have
    #     independent stick-breaking weights pi_g, pi_cf, so
    #       eta_g_k = kappa * pi_g_k,  eta_cf_k = kappa * pi_cf_k.
    #     This is the "two independent DPs" formulation; everything else
    #     (kappa, alpha, K, horseshoe, sigma_k^2, ...) is identical to M1.
    #
    #   Tempering (writeup Algorithm 1, Step 1): the Gaussian factor is raised
    #   to 1/T AND the Polya-urn running count is divided by T, i.e. the prior
    #   factor is (eta + n^{(-j)}/T). The anchor eta is not tempered. Tempering
    #   the count makes the collapsed urn the marginal of a tempered
    #   Dirichlet-multinomial, so both factors see the same likelihood power.
    # -------------------------------------------------------------------
    ## Allocation-prior pseudocount (anchor role). With anchor_frac set this is
    ## anchor_frac*p*pi (scales with p so it is not swamped by the p-sized count
    ## term); with anchor_frac=NULL it is the legacy kappa*pi. This is DECOUPLED
    ## from the tracking scalar lambda_track, which enters only the ω reconstruction.
    if (model == "M1") {
      eta_g  <- anchor_mult * pi_vec
      eta_cf <- anchor_mult * pi_vec
    } else {
      eta_g  <- anchor_mult * pi_g_vec
      eta_cf <- anchor_mult * pi_cf_vec
    }
    if (use_rcpp) {
      alloc <- compute_allocations_cpp(Y_g, Y_cf,
                                        as.numeric(theta), as.numeric(theta_0),
                                        sigma_k2, sigma_0_2,
                                        eta_g, eta_cf, omega_0,
                                        patient, time, K, p, temp)
    } else {
      alloc <- compute_allocations(Y_g, Y_cf, theta, theta_0, sigma_k2, sigma_0_2,
                                   eta_g, eta_cf, omega_0, patient, time, K, p, temp)
    }
    n_g_counts  <- alloc$n_g
    n_cf_counts <- alloc$n_cf
    N0_vec      <- alloc$N0
    Nplus_vec   <- alloc$Nplus

    # -------------------------------------------------------------------
    # SYMMETRIC TEMPERING (Steps 2-4)
    #   Scale data-driven sufficient statistics by 1/temp so that Steps 2-4
    #   see a tempered likelihood of the same power as Step 1. The Step 1
    #   Polya-urn prior factor is intentionally NOT tempered; only its
    #   Gaussian-density factor is.
    #   Horseshoe (Step 5) is a prior-side update and is NOT tempered.
    # -------------------------------------------------------------------
    inv_temp <- 1 / temp
    n_g_counts_T  <- n_g_counts  * inv_temp
    n_cf_counts_T <- n_cf_counts * inv_temp
    N0_vec_T      <- N0_vec      * inv_temp
    Nplus_vec_T   <- Nplus_vec   * inv_temp
    ss_component_T <- alloc$ss_component * inv_temp
    ss_normal_T   <- alloc$ss_normal * inv_temp

    # -------------------------------------------------------------------
    # STEP 2: Update the signature hierarchy (writeup Algorithm Step 2a-2c).
    #   2a  theta[j,k,i]  -- subject-level, pooling BOTH sources' allocated
    #       observations for that subject (the signature is shared across
    #       sources). Empty component reverts to theta_bar, not to zero.
    #   2b  theta_bar[j,k] -- population level, horseshoe prior. Not tempered.
    #   2c  theta_0[j,i]   -- subject-level cfDNA background mean.
    #   Data sufficient statistics are scaled by 1/temp inside the updater.
    # -------------------------------------------------------------------
    if (use_rcpp) {
      theta_list <- update_theta_cpp(Y_g, Y_cf, alloc$z_g, alloc$z_cf,
                                      sigma_k2, varsigma_k2, as.numeric(theta_bar),
                                      tau_k2, lambda2,
                                      sigma_0_2, varsigma_0_2, theta_0_bar,
                                      patient, K, p, N_obs, n, inv_temp)
    } else {
      theta_list <- update_theta(Y_g, Y_cf, alloc$z_g, alloc$z_cf,
                                  sigma_k2, varsigma_k2, theta_bar,
                                  tau_k2, lambda2,
                                  sigma_0_2, varsigma_0_2, theta_0_bar,
                                  patient, K, p, N_obs, n, inv_temp)
    }
    theta     <- theta_list$theta
    theta_bar <- theta_list$theta_bar
    theta_0   <- theta_list$theta_0

    # --- Hyperparameters of the signature hierarchy ---
    # varsigma_k^2 | . ~ IG(a_vs + np/2, b_vs + 0.5*sum_{i,j}(theta - theta_bar)^2)
    # Shape grows as np, so this is sharply determined even at modest n --
    # which is what lets the hierarchical level regularize effectively even
    # though any individual theta[j,k,i] is weakly informed.
    for (k in 1:K) {
      dev_k <- theta[, k, ] - matrix(theta_bar[, k], nrow = p, ncol = n)
      varsigma_k2[k] <- 1 / rgamma(1, hp$a_vs + n * p / 2,
                                    hp$b_vs + 0.5 * sum(dev_k^2))
      varsigma_k2[k] <- max(varsigma_k2[k], 1e-8)
    }

    # ---- EMPTY-COMPONENT RE-ANCHOR (burn-in only) ----------------------
    #
    # PROBLEM. A component k holding zero allocations has n_ijk = 0 for every
    # (i, j). Step 2a therefore draws theta[, k, ] from its prior
    # N(theta_bar[, k], varsigma_k^2), and Step 2b draws theta_bar[, k] back
    # from those same thetas under the horseshoe. That two-block Gibbs cycle
    # touches NO likelihood term: it is a random walk on a hierarchical normal
    # with no data, restrained only by a half-Cauchy tail. The horseshoe adds
    # positive feedback, since a larger theta_bar inflates lambda^2, which
    # lowers the prior precision, which permits a still larger theta_bar.
    #
    # On the production fit this drove |theta_bar| to 153 while the data lie in
    # [0.01, 11.49]. Two consequences, the second being the damaging one:
    #   (a) spurious between-chain variance in theta_bar, varsigma_k and tau_k;
    #   (b) the component becomes a USELESS BIRTH CANDIDATE. Occupying it costs
    #       roughly 3.2 nats per gene and needs ~1000 genes to move at once, so
    #       single-site Gibbs cannot do it. Chains then cannot agree on the
    #       number of occupied components, which is exactly what was observed
    #       (K+ modes 3, 2, 2, 2 with loglik split-Rhat 1.79).
    #
    # FIX. While a component is empty, re-anchor it on the per-gene data centre
    # rather than letting the prior-only cycle carry it away. It stays a
    # plausible birth candidate, so K+ can move.
    #
    # VALIDITY. Gated on iter <= n_burn, so it cannot touch any saved draw:
    # sampling of the reported posterior begins at iter > n_burn (see the
    # save block below). This is the same logical status as the existing
    # `anneal` block -- both alter the burn-in trajectory only, and neither
    # changes the post-burn-in transition kernel. The move fires only for
    # components with zero allocations, whose coordinates enter the joint
    # density through the prior alone.
    if (isTRUE(empty_reanchor) && iter <= n_burn) {
      N_k_all <- rowSums(n_g_counts) + rowSums(n_cf_counts)
      k_empty <- which(N_k_all == 0)
      if (length(k_empty)) {
        vs_prior_mean <- hp$b_vs / max(hp$a_vs - 1, 1e-8)
        for (k in k_empty) {
          sd_k <- sqrt(varsigma_k2[k])
          theta_bar[, k]  <- y_ctr_reanchor + rnorm(p, 0, sd_k)
          theta[, k, ]    <- rnorm(p * n, rep(theta_bar[, k], n), sd_k)
          varsigma_k2[k]  <- vs_prior_mean
          tau_k2[k]       <- 0.1
          lambda2[, k]    <- 1
          nu_aux[, k]     <- 1
          xi_aux[k]       <- 1
        }
      }
    }
    # --------------------------------------------------------------------

    # --- Background hierarchy: theta_0_bar[j] and varsigma_0^2 ---
    # theta_0_bar[j] | . ~ N( s2*(sum_i theta_0[j,i]/vs0 + m_0/s_0_2), s2 )
    prec_0bar   <- n / varsigma_0_2 + 1 / hp$s_0_2
    var_0bar    <- 1 / prec_0bar
    mean_0bar   <- var_0bar * (rowSums(theta_0) / varsigma_0_2 + hp$m_0 / hp$s_0_2)
    theta_0_bar <- rnorm(p, mean_0bar, sqrt(var_0bar))

    dev_0 <- theta_0 - matrix(theta_0_bar, nrow = p, ncol = n)
    varsigma_0_2 <- 1 / rgamma(1, hp$a_vs0 + n * p / 2,
                                hp$b_vs0 + 0.5 * sum(dev_0^2))
    varsigma_0_2 <- max(varsigma_0_2, 1e-8)

    # -------------------------------------------------------------------
    # STEP 3: Update stick-breaking variables.
    #   Under M1, one shared v parameterizes the global pi, with counts
    #     pooled across both sources: m_k = sum over (i,t) of (n_g + n_cf).
    #   Under M0, two independent stick-breaking vectors v_g, v_cf
    #     parameterize the source-specific pi_g, pi_cf, with counts taken
    #     from only the corresponding source: m_g_k uses n_g, m_cf_k uses n_cf.
    #   Counts are tempered: at T>1 the sampler sees effectively fewer
    #   allocations, so stick-breaking is closer to its prior.
    # -------------------------------------------------------------------
    if (model == "M1") {
      m_k <- rowSums(n_g_counts_T) + rowSums(n_cf_counts_T)
      for (k in 1:(K-1)) {
        a_post <- 1 + m_k[k]
        b_post <- alpha_dp + sum(m_k[(k+1):K])
        v[k] <- rbeta(1, a_post, b_post)
        v[k] <- max(min(v[k], 1 - 1e-10), 1e-10)
      }
      v[K] <- 1
      pi_vec <- stick_break(v)
    } else {
      # M0: two independent DPs (one per source).
      m_g_k  <- rowSums(n_g_counts_T)
      m_cf_k <- rowSums(n_cf_counts_T)
      for (k in 1:(K-1)) {
        a_g  <- 1 + m_g_k[k]
        b_g  <- alpha_dp + sum(m_g_k[(k+1):K])
        v_g[k] <- rbeta(1, a_g, b_g)
        v_g[k] <- max(min(v_g[k], 1 - 1e-10), 1e-10)
        a_cf <- 1 + m_cf_k[k]
        b_cf <- alpha_dp + sum(m_cf_k[(k+1):K])
        v_cf[k] <- rbeta(1, a_cf, b_cf)
        v_cf[k] <- max(min(v_cf[k], 1 - 1e-10), 1e-10)
      }
      v_g[K]  <- 1
      v_cf[K] <- 1
      pi_g_vec  <- stick_break(v_g)
      pi_cf_vec <- stick_break(v_cf)
    }

    # -------------------------------------------------------------------
    # STEP 3b (D2): Update the cross-source coupling phi (M1 only) with the
    #   per-observation shared center m_it MARGINALIZED OUT. This REPLACES both
    #   the D1 per-source dispersion (non-identified) AND the first D2 attempt
    #   that ESTIMATED m_it per observation -- that estimated m_it overfit the
    #   p-sized counts and absorbed the coupling, so phi pinned to a p/K constant
    #   (~40) irrespective of the truth. Integrating m_it over its prior instead
    #   restores identifiability: on realistic p=2000 counts the marginalized
    #   estimator discriminates kappa_true=2 (phi_hat~1.0) from 20 (~1.7), while
    #   the estimated-m version gave ~40 for both.
    #
    #   Model: m_it ~ Dir(alpha0 * pi);  omega_s_it ~ Dir(phi * m_it), s in {g,cf};
    #          n_s_it ~ Multinomial(omega_s_it)  =>  n_s_it | m_it ~ DirMult(phi*m_it).
    #   phi likelihood integrates m_it out by Monte Carlo over its prior:
    #     L(phi) = prod_it E_{m~Dir(alpha0*pi)}[ DM(n_g_it;phi*m) * DM(n_cf_it;phi*m) ].
    #   The PAIR shares the same m draw, so phi is identified by cross-source
    #   AGREEMENT (both count vectors point the same way only when phi is large).
    #   Common random numbers (shared m draws) are used for the current and
    #   proposed phi so the MH ratio is low-variance. Not tempered.
    # -------------------------------------------------------------------
    if (use_d2) {
      pi_c <- pmax(pi_vec, 1e-10); pi_c <- pi_c / sum(pi_c)
      n_mc_m <- if (!is.null(hp$n_mc_m)) hp$n_mc_m else 32L
      # pre-draw shared m samples (common random numbers across the MH ratio)
      m_draws <- replicate(n_mc_m, rdirichlet_one(alpha0_m * pi_c))  # K x n_mc_m
      lse <- function(v) { mx <- max(v); mx + log(mean(exp(v - mx))) }
      # collapsed paired DirMult log-lik at concentration a = phi*m, summed over obs
      phi_marg_ll <- function(phi) {
        tot <- 0
        for (obs in 1:N_obs) {
          ng_o <- n_g_counts[, obs]; ncf_o <- n_cf_counts[, obs]
          Ng <- sum(ng_o); Ncf <- sum(ncf_o)
          if (Ng <= 0 && Ncf <= 0) next
          ll_mc <- numeric(n_mc_m)
          for (s in 1:n_mc_m) {
            a <- phi * pmax(m_draws[, s], 1e-12); sa <- sum(a); lg_a <- sum(lgamma(a))
            ll <- 0
            if (Ng  > 0) ll <- ll + lgamma(sa) - lgamma(Ng  + sa) + sum(lgamma(ng_o  + a)) - lg_a
            if (Ncf > 0) ll <- ll + lgamma(sa) - lgamma(Ncf + sa) + sum(lgamma(ncf_o + a)) - lg_a
            ll_mc[s] <- ll
          }
          tot <- tot + lse(ll_mc)   # log E_m over the shared draws
        }
        tot
      }
      phi_total <- phi_total + 1L
      log_phi_prop <- rnorm(1, log(phi_track), sigma_phi_mh)
      phi_prop <- exp(log_phi_prop)
      log_ratio_prior <- (a_lambda - 1) * (log(phi_prop) - log(phi_track)) -
                         b_lambda * (phi_prop - phi_track)
      log_alpha <- log_ratio_prior + phi_marg_ll(phi_prop) - phi_marg_ll(phi_track) +
                   log(phi_prop) - log(phi_track)  # log-normal Jacobian
      if (!is.finite(log_alpha)) log_alpha <- -Inf
      if (log(runif(1)) < log_alpha) {
        phi_track <- phi_prop
        phi_accept <- phi_accept + 1L
      }
      lambda_track <- phi_track  # keep the reported tracking scalar in sync
    }

    # -------------------------------------------------------------------
    # (The omega Gibbs update has been marginalized into Step 1.)
    # omega_g, omega_cf are no longer state variables. We reconstruct them
    # from a single Dirichlet draw at each SAVED iteration only, so that
    # downstream concordance metrics still have an omega trace.
    # -------------------------------------------------------------------

    # -------------------------------------------------------------------
    # STEP 4: Update variance and nuisance parameters (tempered)
    # -------------------------------------------------------------------
    # Contamination fractions, now one per OBSERVATION (eq:app_fc_wN):
    #   omega_0[i,t] | z^cf_{it} ~ Beta(a_0 + N0_{it}, b_0 + Nplus_{it})
    # Counts are specific to (i,t) rather than pooled over a subject's
    # timepoints, which is what makes the fraction time-resolved. Each is
    # informed by p allocations from a single observation, so the posterior is
    # wider than the old per-patient version -- the honest cost of resolving
    # contamination in time.
    omega_0 <- rbeta(N_obs, hp$a_0 + N0_vec_T, hp$b_0 + Nplus_vec_T)
    omega_0 <- pmax(pmin(omega_0, 1 - 1e-10), 1e-10)

    # Background variance (eq:app_fc_sigmaN). ss_normal is a sum of squared
    # residuals about theta_0[j,i], NOT a raw second moment: under the old
    # mean-zero background this was sum(Y_cf^2), which conflated the
    # background's location with its dispersion.
    #
    # sigma_0^2 shares the rate hyperparameter b_sigma with the sigma_k^2, so
    # the background and tumour dispersions are on a common scale. This keeps
    # sigma_0 from drifting independently when the background holds few
    # allocations; it is a mild pooling choice, not load-bearing for
    # identification (the allocation normaliser above is what makes the
    # background identifiable).
    N_total_0_T <- sum(N0_vec_T)
    sigma_0_2 <- 1 / rgamma(1, hp$a_sigma0 + N_total_0_T/2,
                             b_sigma + ss_normal_T/2)

    for (k in 1:K) {
      ss_k <- ss_component_T[k]
      N_k  <- (sum(n_g_counts[k, ]) + sum(n_cf_counts[k, ])) * inv_temp
      sigma_k2[k] <- 1 / rgamma(1, hp$a_sigma + N_k/2,
                                  b_sigma + ss_k/2)
      ## FLOOR sigma_k^2, as varsigma_k^2 already is (see the varsigma update).
      ## Without it a component whose allocated genes happen to be near-identical
      ## can draw sigma_k^2 ~ 0; dnorm(y, theta, 0, log = TRUE) is then -Inf for
      ## EVERY component, log_sum_exp2 propagates -Inf, and a single such draw
      ## poisons loglik, WAIC, PSIS-LOO and every downstream Rhat.
      ##
      ## NOTE the mechanism: this is NOT floating-point underflow of the density.
      ## dnorm(..., log = TRUE) is exact even at |y - theta| = 1e4 (it returns
      ## -5e7, not -Inf). The -Inf requires a degenerate sigma_k, so flooring the
      ## variance is the correct and sufficient guard.
      sigma_k2[k] <- max(sigma_k2[k], 1e-8)
    }

    # Update shared rate b_sigma for the variance priors (hierarchical pooling).
    # sigma_0^2 now draws from the SAME rate as the sigma_k^2, so it contributes
    # one additional shape count and one additional 1/sigma^2 term:
    #   b_sigma | . ~ Gamma((K+1)*a_sigma + a_b_sigma,
    #                       sum_k 1/sigma_k^2 + 1/sigma_0^2 + b_b_sigma)
    # b_sigma is a prior-side hyperparameter and is NOT tempered.
    b_sigma <- rgamma(1,
                      (K + 1) * hp$a_sigma + hp$a_b_sigma,
                      sum(1 / sigma_k2) + 1 / sigma_0_2 + hp$b_b_sigma)

    # -------------------------------------------------------------------
    # STEP 5: Update horseshoe shrinkage parameters
    #   The horseshoe now acts on the single POPULATION signature:
    #     theta_bar[j,k] ~ N(0, tau_k^2 * lambda_{jk}^2)
    #   so each (j,k) contributes ONE squared term instead of two.
    #   lambda_{jk}^2 ~ IG(1,        1/nu + theta_bar_{jk}^2/(2 tau_k^2))
    #   tau_k^2       ~ IG((p+1)/2,  1/xi + sum_j theta_bar_{jk}^2/lambda_{jk}^2 / 2)
    #   Shapes were (3/2) and ((2p+1)/2) under the old source-specific,
    #   time-indexed specification; see writeup eq:app_fc_lambda / eq:app_fc_tau.
    #   Cross-source and cross-subject borrowing has NOT disappeared -- it now
    #   lives in the shared theta and in varsigma_k^2, rather than being routed
    #   through the shrinkage hyperparameter.
    #   This is a prior-side update and is NOT tempered.
    # -------------------------------------------------------------------
    for (k in 1:K) {
      theta_k_sq <- theta_bar[, k]^2   # [p] vector, one term per (j,k)

      # Local shrinkage lambda_{jk}^2 (vectorized over genes); shape 1
      rate_lam <- 1/nu_aux[, k] + theta_k_sq / (2 * tau_k2[k])
      rate_lam <- pmax(rate_lam, 1e-300)
      lambda2[, k] <- 1 / rgamma(p, 1, rate_lam)

      # Local auxiliary nu_{jk}
      rate_nu <- 1 + 1/lambda2[, k]
      rate_nu <- pmax(rate_nu, 1e-300)
      nu_aux[, k] <- 1 / rgamma(p, 1, rate_nu)

      # Global shrinkage tau_k^2; shape (p+1)/2
      sum_theta_lam <- sum(theta_k_sq / lambda2[, k])
      rate_tau <- 1/xi_aux[k] + sum_theta_lam/2
      rate_tau <- max(rate_tau, 1e-300)
      tau_k2[k] <- 1 / rgamma(1, (p + 1)/2, rate_tau)

      # Global auxiliary xi_k
      xi_aux[k] <- 1 / rgamma(1, 1, 1 + 1/tau_k2[k])
    }

    # -------------------------------------------------------------------
    # FIXED HYPERPARAMETERS: kappa and alpha_dp are both fixed at 1 and
    # are not updated. See manuscript Remark rem:kappa_grid.
    # -------------------------------------------------------------------

    # -------------------------------------------------------------------
    # Save samples (after burn-in, with thinning)
    # -------------------------------------------------------------------
    if (iter > n_burn && (iter - n_burn) %% thin == 0) {
      idx <- (iter - n_burn) %/% thin
      samples$kappa[idx]    <- lambda_track   # tracking scalar (legacy slot name)
      samples$lambda[idx]   <- lambda_track
      samples$phi[idx]      <- phi_track      # D2 cross-source coupling
      samples$alpha_dp[idx] <- alpha_dp
      samples$omega_0[idx, ] <- omega_0
      samples$sigma_0[idx]  <- sqrt(sigma_0_2)
      samples$varsigma_0[idx]   <- sqrt(varsigma_0_2)
      samples$b_sigma[idx]  <- b_sigma
      # Component-indexed quantities (sigma_k, varsigma_k, tau_k) are stored
      # below, AFTER ord_theta is computed, so that every component-indexed
      # slot refers to the same canonical component as theta and the omega
      # traces. Storing them in raw order here would silently desynchronise
      # them from theta.

      # Accumulate the signature hierarchy toward posterior mean and variance.
      # IMPORTANT: accumulate in a CANONICAL component order (decreasing
      # cohort weight), the same order applied to the omega traces below.
      # sigma_k / varsigma_k / lambda2 are indexed by the same canonical k.
      ord_theta <- if (model == "M1") {
        order(pi_vec, decreasing = TRUE)
      } else {
        order(rowSums(n_g_counts) + rowSums(n_cf_counts), decreasing = TRUE)
      }
      th_ord  <- theta[, ord_theta, , drop = FALSE]      # [p x K x n]
      tb_ord  <- theta_bar[, ord_theta, drop = FALSE]    # [p x K]
      theta_sum       <- theta_sum       + th_ord
      theta_sumsq     <- theta_sumsq     + th_ord^2
      theta_bar_sum   <- theta_bar_sum   + tb_ord
      theta_bar_sumsq <- theta_bar_sumsq + tb_ord^2
      # theta_0 has no component index, so it needs no reordering.
      theta_0_sum     <- theta_0_sum     + theta_0
      theta_0_sumsq   <- theta_0_sumsq   + theta_0^2
      theta_nsave <- theta_nsave + 1L

      # Also reorder the per-component variances so that every stored
      # component-indexed quantity refers to the same canonical component.
      samples$sigma_k[idx, ]    <- sqrt(sigma_k2[ord_theta])
      samples$varsigma_k[idx, ] <- sqrt(varsigma_k2[ord_theta])
      samples$tau_k[idx, ]      <- sqrt(tau_k2[ord_theta])

      # Store thinned theta traces in the same canonical order.
      if (store_theta_trace && theta_trace_i < n_theta_trace &&
          (idx %% theta_trace_thin == 0)) {
        theta_trace_i <- theta_trace_i + 1L
        theta_trace[theta_trace_i, , , ]  <- th_ord
        theta_bar_trace[theta_trace_i, , ] <- tb_ord
      }

      if (model == "M1") {
        samples$K_plus[idx] <- sum(pi_vec > 0.01)
        ## Store the stick-breaking variables in the SAME canonical order as
        ## theta / sigma_k / varsigma_k / tau_k. Previously v was stored raw
        ## while everything else was canonical, so relabel_by_weight() derived
        ## its permutation from the raw v and applied it a SECOND time to the
        ## already-canonical arrays -- desynchronising sigma_k[k] from
        ## theta[,k]. Storing canonically here makes stick_break(v) already
        ## sorted, so the downstream permutation is the identity.
        samples$v[idx, ] <- inv_stick_break(pi_vec[ord_theta])[1:(K-1)]
      } else {
        # K+ MUST use the same definition as M1 (weight above 0.01), otherwise
        # the two models are not comparable on this quantity. The previous
        # definition here counted any component with >= 1 allocated gene, which
        # is a much weaker criterion: on the production fit it reported K+ = 5
        # for M0 against 2-3 for M1, an apparent model difference that was
        # entirely an artifact of the two definitions. Under the common
        # weight-based rule M0 gives a consistent K+ = 3.
        #
        # The occupancy count is retained separately as K_plus_occupied, since
        # it is a useful diagnostic (it detects components holding a handful of
        # genes), but it is NOT the reported K+.
        m_g_total  <- rowSums(n_g_counts)
        m_cf_total <- rowSums(n_cf_counts)
        pi_g_s  <- stick_break(c(v_g[1:(K-1)],  1))
        pi_cf_s <- stick_break(c(v_cf[1:(K-1)], 1))
        samples$K_plus[idx] <- sum(pmax(pi_g_s, pi_cf_s) > 0.01)
        if (!is.null(samples$K_plus_occupied))
          samples$K_plus_occupied[idx] <- sum((m_g_total > 0) | (m_cf_total > 0))
        samples$v_g[idx, ]  <- v_g[1:(K-1)]
        samples$v_cf[idx, ] <- v_cf[1:(K-1)]
      }

      # POST-HOC reconstruction of source-specific weights from their
      # conditional Dirichlet posterior (used only for downstream concordance
      # metrics; does not feed back into the sampler).
      #   D2 (M1, coupling estimated): both sources concentrate around the SAME
      #     per-observation center m_it, so omega_s ~ Dir(phi*m_it + n_s). The
      #     shared m_it is what makes the reconstructed omega_g, omega_cf agree
      #     (the cross-source tracking the concordance metrics measure).
      #   M1 legacy (lambda fixed, no D2): omega_s ~ Dir(lambda*pi + n_s).
      #   M0 (source-DP): omega_s ~ Dir(lambda*pi_s + n_s), independent per source.
      if (use_d2) {
        pi_c <- pmax(pi_vec, 1e-10); pi_c <- pi_c / sum(pi_c)
        for (obs in 1:N_obs) {
          # shared per-obs center from its posterior given the pooled counts
          m_obs <- rdirichlet_one(pmax(alpha0_m * pi_c + n_g_counts[, obs] + n_cf_counts[, obs], 1e-10))
          am <- phi_track * m_obs
          omega_g[, obs]  <- rdirichlet_one(pmax(am + n_g_counts[, obs],  1e-10))
          omega_cf[, obs] <- rdirichlet_one(pmax(am + n_cf_counts[, obs], 1e-10))
        }
      } else {
        if (model == "M1") {
          prior_alpha_g  <- lambda_track * pi_vec
          prior_alpha_cf <- lambda_track * pi_vec
        } else {
          prior_alpha_g  <- lambda_track * pi_g_vec
          prior_alpha_cf <- lambda_track * pi_cf_vec
        }
        for (obs in 1:N_obs) {
          omega_g[, obs]  <- rdirichlet_one(pmax(prior_alpha_g  + n_g_counts[, obs],  1e-10))
          omega_cf[, obs] <- rdirichlet_one(pmax(prior_alpha_cf + n_cf_counts[, obs], 1e-10))
        }
      }

      # Store omega traces for the top K_trace components, in the SAME
      # canonical component order as theta / sigma_k / varsigma_k / tau_k / v.
      # Storing these raw while the scalars were canonical was the second half
      # of the desynchronisation bug: relabel_by_weight() then permuted the
      # traces (raw -> canonical, correct) and the scalars (canonical ->
      # canonical again, WRONG) in the same pass.
      for (kk in 1:K_trace) {
        col_start <- (kk - 1) * N_obs + 1
        col_end   <- kk * N_obs
        k_src <- ord_theta[kk]
        samples$omega_g_trace[idx, col_start:col_end]  <- omega_g[k_src, ]
        samples$omega_cf_trace[idx, col_start:col_end] <- omega_cf[k_src, ]
      }

      if (use_rcpp) {
        ll_result <- compute_loglik_cpp(Y_g, Y_cf,
                                         as.numeric(theta), as.numeric(theta_0),
                                         sigma_k2, sigma_0_2,
                                         omega_g, omega_cf, omega_0,
                                         patient, time, K, p,
                                         store_pointwise_ll)
      } else {
        ll_result <- compute_loglik(Y_g, Y_cf, theta, theta_0, sigma_k2,
                                     sigma_0_2, omega_g, omega_cf,
                                     omega_0, patient, time, K, p,
                                     return_pointwise = store_pointwise_ll)
      }
      samples$loglik[idx] <- ll_result$total

      if (store_pointwise_ll) {
        samples$pointwise_ll[idx, ] <- ll_result$pointwise
      }
    }

    # Progress
    if (verbose && iter %% 500 == 0) {
      if (model == "M1") {
        K_active <- sum(pi_vec > 0.01)
      } else {
        K_active <- sum((pi_g_vec > 0.01) | (pi_cf_vec > 0.01))
      }
      cat(sprintf("  iter %d/%d | T=%.2f | kappa=%.1f | K+=%.0f | om0_mean=%.3f | vs_k=%.2f\n",
                  iter, n_iter, temp, kappa, K_active, mean(omega_0),
                  sqrt(mean(varsigma_k2))))
    }
  }

  # --- Return ---
  # Reconstruct omega_g and omega_cf one last time from the post-iteration
  # allocation counts so that the returned `final_state` carries an
  # internally consistent set of state variables for downstream code (e.g.,
  # the parallel-tempering swap step that evaluates the observed-data
  # log-likelihood).
  if (exists("n_g_counts", inherits = FALSE) &&
      exists("n_cf_counts", inherits = FALSE)) {
    if (use_d2) {
      pi_c <- pmax(pi_vec, 1e-10); pi_c <- pi_c / sum(pi_c)
      for (obs in 1:N_obs) {
        m_obs <- rdirichlet_one(pmax(alpha0_m * pi_c + n_g_counts[, obs] + n_cf_counts[, obs], 1e-10))
        am <- phi_track * m_obs
        omega_g[, obs]  <- rdirichlet_one(pmax(am + n_g_counts[, obs],  1e-10))
        omega_cf[, obs] <- rdirichlet_one(pmax(am + n_cf_counts[, obs], 1e-10))
      }
    } else {
      if (model == "M1") {
        prior_alpha_g_final  <- lambda_track * pi_vec
        prior_alpha_cf_final <- lambda_track * pi_vec
      } else {
        prior_alpha_g_final  <- lambda_track * pi_g_vec
        prior_alpha_cf_final <- lambda_track * pi_cf_vec
      }
      for (obs in 1:N_obs) {
        omega_g[, obs]  <- rdirichlet_one(pmax(prior_alpha_g_final  + n_g_counts[, obs],  1e-10))
        omega_cf[, obs] <- rdirichlet_one(pmax(prior_alpha_cf_final + n_cf_counts[, obs], 1e-10))
      }
    }
  }

  # final_state contains everything needed to warm-start a subsequent run
  # via the `init` argument (used by parallel tempering driver).
  final_state <- list(
    theta       = theta,
    theta_bar   = theta_bar,
    theta_0     = theta_0,
    theta_0_bar = theta_0_bar,
    v         = v,
    pi_vec    = pi_vec,
    v_g       = v_g,
    v_cf      = v_cf,
    pi_g_vec  = pi_g_vec,
    pi_cf_vec = pi_cf_vec,
    omega_g   = omega_g,
    omega_cf  = omega_cf,
    sigma_k2  = sigma_k2,
    sigma_0_2 = sigma_0_2,
    varsigma_k2  = varsigma_k2,
    varsigma_0_2 = varsigma_0_2,
    omega_0   = omega_0,
    tau_k2    = tau_k2,
    lambda2   = lambda2,
    nu_aux    = nu_aux,
    xi_aux    = xi_aux,
    b_sigma   = b_sigma
  )
  # Compute posterior means and variances from accumulators
  pm <- function(sum_, fallback) if (theta_nsave > 0) sum_ / theta_nsave else fallback
  pv <- function(sumsq_, sum_, dims) if (theta_nsave > 1)
    (sumsq_ - sum_^2 / theta_nsave) / (theta_nsave - 1) else array(NA_real_, dim = dims)

  theta_postmean     <- pm(theta_sum,     theta)
  theta_bar_postmean <- pm(theta_bar_sum, theta_bar)
  theta_0_postmean   <- pm(theta_0_sum,   theta_0)
  theta_postvar      <- pv(theta_sumsq,     theta_sum,     c(p, K, n))
  theta_bar_postvar  <- pv(theta_bar_sumsq, theta_bar_sum, c(p, K))
  theta_0_postvar    <- pv(theta_0_sumsq,   theta_0_sum,   c(p, n))

  list(
    samples = samples,
    n_save  = n_save,
    model   = model,
    K       = K,
    n_iter  = n_iter,
    n_burn  = n_burn,
    thin    = thin,
    kappa_fixed = kappa_fixed,
    # Signature hierarchy (last iteration)
    final_theta     = theta,
    final_theta_bar = theta_bar,
    final_theta_0   = theta_0,
    # Posterior means and variances over saved iterations, in canonical
    # component order. theta is [p x K x n] (subject-level, shared across
    # sources); theta_bar is [p x K] (population); theta_0 is [p x n].
    # Interpret theta_bar substantively -- the subject-level signatures are
    # weakly informed by design and are not consistently estimable.
    theta_postmean     = theta_postmean,
    theta_bar_postmean = theta_bar_postmean,
    theta_0_postmean   = theta_0_postmean,
    theta_postvar      = theta_postvar,
    theta_bar_postvar  = theta_bar_postvar,
    theta_0_postvar    = theta_0_postvar,
    theta_nsave        = theta_nsave,
    # Thinned traces for the split-Rhat audit (NULL unless store_theta_trace=TRUE)
    theta_trace        = theta_trace,
    theta_bar_trace    = theta_bar_trace,
    theta_trace_thin   = if (store_theta_trace) theta_trace_thin else NA_integer_,
    final_pi    = pi_vec,
    final_pi_g  = pi_g_vec,
    final_pi_cf = pi_cf_vec,
    final_omega_g = omega_g,
    final_omega_cf = omega_cf,
    final_state = final_state,
    hyperparams = hp,
    K_trace = K_trace,
    p = p, N_obs = N_obs, n = n, T_max = T_max
  )
}

# ===========================================================================
# HELPER FUNCTIONS
# ===========================================================================

stick_break <- function(v) {
  K <- length(v)
  pi_vec <- numeric(K)
  cum_prod <- 1
  for (k in 1:K) {
    pi_vec[k] <- v[k] * cum_prod
    cum_prod <- cum_prod * (1 - v[k])
  }
  pi_vec <- pi_vec / sum(pi_vec)
  pi_vec
}

#' Invert the stick-breaking map: recover v from a weight vector.
#'
#' Right inverse of stick_break(): stick_break(inv_stick_break(pi)) == pi
#' up to the truncation of the final stick. Used to store the stick-breaking
#' variables in the same canonical (decreasing-weight) component order as
#' theta / sigma_k / varsigma_k, so that every component-indexed slot in
#' `samples` refers to the same component k.
#'
#' @param pi_vec Length-K weight vector (need not be sorted).
#' @return Length-K vector v with v[K] = 1.
inv_stick_break <- function(pi_vec) {
  K <- length(pi_vec)
  v <- numeric(K)
  cum <- 1
  if (K > 1) for (k in 1:(K - 1)) {
    v[k] <- if (cum > 0) pi_vec[k] / cum else 0.5
    v[k] <- max(min(v[k], 1 - 1e-10), 1e-10)
    cum  <- cum * (1 - v[k])
  }
  v[K] <- 1
  v
}

rdirichlet_one <- function(alpha) {
  x <- rgamma(length(alpha), alpha, 1)
  s <- sum(x)
  if (s < 1e-300) return(rep(1/length(alpha), length(alpha)))
  x / s
}

#' Compute allocations for all observations (COLLAPSED / Rao-Blackwellized)
#'
#' Pure-R fallback that mirrors the Rcpp `compute_allocations_cpp`. Each
#' gene's allocation is sampled from the Dirichlet-multinomial Polya-urn
#' predictive (\code{eta + n^{(-j)}/T}) * Gaussian likelihood^(1/T), with genes
#' scanned in a random order within each observation. BOTH the Gaussian factor
#' and the urn running count are tempered (the anchor \code{eta} is not); see
#' writeup Algorithm 1, Step 1.
#'
#' @param eta_g numeric vector of length K, the prior parameter for the
#'   collapsed Dirichlet on gDNA allocations (kappa*pi under M1, or
#'   kappa*pi_g under the source-DP M0).
#' @param eta_cf numeric vector of length K, the prior parameter for the
#'   collapsed Dirichlet on cfDNA allocations (kappa*pi under M1, or
#'   kappa*pi_cf under the source-DP M0).
compute_allocations <- function(Y_g, Y_cf, theta, theta_0, sigma_k2, sigma_0_2,
                                 eta_g, eta_cf, omega_0, patient, time,
                                 K, p, temperature = 1.0) {
  N_obs <- ncol(Y_g)
  n <- max(patient)

  z_g  <- matrix(1L, p, N_obs)
  z_cf <- matrix(0L, p, N_obs)

  n_g_counts  <- matrix(0, K, N_obs)
  n_cf_counts <- matrix(0, K, N_obs)
  # Background / tumour counts are PER OBSERVATION (omega_0 is per observation).
  N0_vec      <- numeric(N_obs)
  Nplus_vec   <- numeric(N_obs)
  ss_normal   <- 0
  ss_component <- numeric(K)

  log_sigma_k <- 0.5 * log(sigma_k2)
  inv_2sigma_k2 <- 1 / (2 * sigma_k2)
  log_sigma_0 <- 0.5 * log(sigma_0_2)
  inv_2sigma_0_2 <- 1 / (2 * sigma_0_2)

  inv_temp <- 1 / temperature

  for (obs in 1:N_obs) {
    i <- patient[obs]

    n_run_g  <- numeric(K)
    n_run_cf <- numeric(K)

    # --- gDNA allocations: pure mixture over subject i's own signatures ---
    ord_g <- sample.int(p)
    for (j in ord_g) {
      resid <- Y_g[j, obs] - theta[j, , i]   # [K] vector
      log_prior <- log(pmax(eta_g + inv_temp * n_run_g, 1e-300))
      log_gauss <- (-log_sigma_k - resid^2 * inv_2sigma_k2) * inv_temp
      lp <- log_prior + log_gauss
      lp <- lp - max(lp)
      probs <- exp(lp); probs <- probs / sum(probs)
      zk <- sample.int(K, 1, prob = probs)
      z_g[j, obs] <- zk
      n_run_g[zk] <- n_run_g[zk] + 1
    }
    n_g_counts[, obs] <- n_run_g

    # --- cfDNA allocations: background at slot 0 + K tumour components ---
    log_w0   <- log(max(omega_0[obs], 1e-300))
    log_1mw0 <- log(max(1 - omega_0[obs], 1e-300))
    ord_cf <- sample.int(p)
    for (j in ord_cf) {
      resid <- Y_cf[j, obs] - theta[j, , i]   # [K] vector
      # Normalise the tumour block: the background is a Bernoulli split
      # outside the Dirichlet, so it must be compared against the tumour
      # components' NORMALISED probabilities. See the extended note in
      # cpp/sampler_core.cpp -- omitting this gives every tumour component a
      # spurious +log(p) bonus and collapses omega_0 to the boundary.
      urn <- pmax(eta_cf + inv_temp * n_run_cf, 1e-300)
      log_prior <- log(urn) - log(sum(urn))
      log_gauss <- (-log_sigma_k - resid^2 * inv_2sigma_k2) * inv_temp
      lp_tumor <- log_1mw0 + log_prior + log_gauss
      # Background residual about the ESTIMATED mean theta_0[j,i], not 0.
      resid0 <- Y_cf[j, obs] - theta_0[j, i]
      lp_norm <- log_w0 + (-log_sigma_0 - resid0^2 * inv_2sigma_0_2) * inv_temp
      lp <- c(lp_norm, lp_tumor)
      lp <- lp - max(lp)
      probs <- exp(lp); probs <- probs / sum(probs)
      zk <- sample.int(K + 1, 1, prob = probs) - 1L  # 0 = background
      z_cf[j, obs] <- zk
      if (zk == 0) {
        N0_vec[obs] <- N0_vec[obs] + 1
        ss_normal <- ss_normal + resid0^2       # centred residual
      } else {
        Nplus_vec[obs] <- Nplus_vec[obs] + 1
        n_run_cf[zk] <- n_run_cf[zk] + 1
      }
    }
    n_cf_counts[, obs] <- n_run_cf
  }

  # Sufficient statistics for sigma_k^2. Both sources' residuals are taken
  # about the same subject-level signature theta[j,k,i].
  for (k in 1:K) {
    ss_k <- 0
    for (obs in 1:N_obs) {
      i <- patient[obs]
      idx_g <- which(z_g[, obs] == k)
      if (length(idx_g) > 0) {
        ss_k <- ss_k + sum((Y_g[idx_g, obs] - theta[idx_g, k, i])^2)
      }
      idx_cf <- which(z_cf[, obs] == k)
      if (length(idx_cf) > 0) {
        ss_k <- ss_k + sum((Y_cf[idx_cf, obs] - theta[idx_cf, k, i])^2)
      }
    }
    ss_component[k] <- ss_k
  }

  list(n_g = n_g_counts, n_cf = n_cf_counts,
       N0 = N0_vec, Nplus = Nplus_vec,
       ss_normal = ss_normal, ss_component = ss_component,
       z_g = z_g, z_cf = z_cf)
}

#' Update the signature hierarchy (conjugate Normal at three levels)
#'
#' Mirrors `update_theta_cpp`. Draws, in order:
#'   2a theta[j,k,i]   pooling BOTH sources' allocated obs for subject i.
#'                     Empty component reverts to theta_bar[j,k], NOT to zero.
#'   2b theta_bar[j,k] population level, horseshoe prior. Not tempered.
#'   2c theta_0[j,i]   subject-level cfDNA background mean.
#' Data sufficient statistics are multiplied by inv_temp for annealing.
update_theta <- function(Y_g, Y_cf, z_g, z_cf, sigma_k2,
                          varsigma_k2, theta_bar,
                          tau_k2, lambda2,
                          sigma_0_2, varsigma_0_2, theta_0_bar,
                          patient, K, p, N_obs, n,
                          inv_temp = 1.0) {
  theta         <- array(0, dim = c(p, K, n))
  theta_bar_new <- array(0, dim = c(p, K))
  theta_0       <- array(0, dim = c(p, n))

  # --- Accumulate (subject, gene, component) sufficient statistics, pooled
  #     over BOTH sources and all timepoints of that subject. ---
  sum_y <- array(0, dim = c(p, K, n))
  n_ijk <- array(0, dim = c(p, K, n))
  sum_y0 <- matrix(0, p, n)
  n_ij0  <- matrix(0, p, n)

  for (obs in 1:N_obs) {
    i <- patient[obs]
    for (k in 1:K) {
      mg <- (z_g[, obs] == k)
      if (any(mg)) {
        sum_y[mg, k, i] <- sum_y[mg, k, i] + Y_g[mg, obs]
        n_ijk[mg, k, i] <- n_ijk[mg, k, i] + 1
      }
      mc <- (z_cf[, obs] == k)
      if (any(mc)) {
        sum_y[mc, k, i] <- sum_y[mc, k, i] + Y_cf[mc, obs]
        n_ijk[mc, k, i] <- n_ijk[mc, k, i] + 1
      }
    }
    m0 <- (z_cf[, obs] == 0L)
    if (any(m0)) {
      sum_y0[m0, i] <- sum_y0[m0, i] + Y_cf[m0, obs]
      n_ij0[m0, i]  <- n_ij0[m0, i] + 1
    }
  }

  # --- 2a: subject-level signatures ---
  for (k in 1:K) {
    prior_prec_sub <- 1 / varsigma_k2[k]
    for (i in 1:n) {
      data_prec <- (n_ijk[, k, i] * inv_temp) / sigma_k2[k]
      post_prec <- data_prec + prior_prec_sub
      post_var  <- 1 / post_prec
      # n_ijk == 0 collapses to the prior mean theta_bar[, k].
      post_mean <- post_var * ((sum_y[, k, i] * inv_temp) / sigma_k2[k] +
                                 theta_bar[, k] * prior_prec_sub)
      theta[, k, i] <- rnorm(p, post_mean, sqrt(post_var))
    }
  }

  # --- 2b: population signatures (prior-side; not tempered) ---
  for (k in 1:K) {
    # Sum subject-level signatures for component k: [p x n] -> [p]
    sum_sub    <- rowSums(matrix(theta[, k, ], nrow = p, ncol = n))
    prior_prec <- 1 / (tau_k2[k] * lambda2[, k])
    post_prec  <- n / varsigma_k2[k] + prior_prec
    post_var   <- 1 / post_prec
    post_mean  <- post_var * (sum_sub / varsigma_k2[k])
    theta_bar_new[, k] <- rnorm(p, post_mean, sqrt(post_var))
  }

  # --- 2c: subject-level background means ---
  prior_prec_0 <- 1 / varsigma_0_2
  for (i in 1:n) {
    data_prec <- (n_ij0[, i] * inv_temp) / sigma_0_2
    post_prec <- data_prec + prior_prec_0
    post_var  <- 1 / post_prec
    post_mean <- post_var * ((sum_y0[, i] * inv_temp) / sigma_0_2 +
                               theta_0_bar * prior_prec_0)
    theta_0[, i] <- rnorm(p, post_mean, sqrt(post_var))
  }

  list(theta = theta, theta_bar = theta_bar_new, theta_0 = theta_0)
}

#' Compute observed-data log-likelihood (shared subject-indexed theta)
#'
#' gDNA is a pure mixture; cfDNA adds a background component with the
#' estimated mean theta_0[j,i]. Pointwise layout is [2 * N_obs]: entries
#' 1..N_obs are gDNA (summed over genes), (N_obs+1)..(2*N_obs) are cfDNA.
compute_loglik <- function(Y_g, Y_cf, theta, theta_0, sigma_k2, sigma_0_2,
                            omega_g, omega_cf, omega_0, patient, time,
                            K, p,
                            return_pointwise = FALSE) {
  N_obs <- ncol(Y_g)
  ll_total <- 0

  if (return_pointwise) {
    pw <- numeric(2 * N_obs)
  }

  for (obs in 1:N_obs) {
    i <- patient[obs]
    ll_g <- 0
    ll_cf <- 0

    log_w0   <- log(max(omega_0[obs], 1e-300))
    log_1mw0 <- log(max(1 - omega_0[obs], 1e-300))

    for (j in 1:p) {
      # gDNA: pure mixture over subject i's own signatures
      log_mix_g <- -Inf
      for (k in 1:K) {
        log_mix_g <- log_sum_exp(log_mix_g,
          log(max(omega_g[k, obs], 1e-300)) +
          dnorm(Y_g[j, obs], theta[j, k, i], sqrt(sigma_k2[k]), log = TRUE))
      }
      ll_g <- ll_g + log_mix_g

      # cfDNA: background (estimated mean theta_0[j,i]) + tumour mixture
      log_normal <- log_w0 +
        dnorm(Y_cf[j, obs], theta_0[j, i], sqrt(sigma_0_2), log = TRUE)
      log_tumor <- -Inf
      for (k in 1:K) {
        log_tumor <- log_sum_exp(log_tumor,
          log_1mw0 +
          log(max(omega_cf[k, obs], 1e-300)) +
          dnorm(Y_cf[j, obs], theta[j, k, i], sqrt(sigma_k2[k]), log = TRUE))
      }
      ll_cf <- ll_cf + log_sum_exp(log_normal, log_tumor)
    }

    ll_total <- ll_total + ll_g + ll_cf

    if (return_pointwise) {
      pw[obs] <- ll_g
      pw[N_obs + obs] <- ll_cf
    }
  }

  result <- list(total = ll_total)
  if (return_pointwise) result$pointwise <- pw
  result
}

log_sum_exp <- function(a, b) {
  m <- max(a, b)
  if (is.infinite(m)) return(-Inf)
  m + log(exp(a - m) + exp(b - m))
}


###############################################################################
## ---- 06_multi_chain.R
###############################################################################

###############################################################################
## 06_multi_chain.R
## Multi-chain MCMC infrastructure
##
## Runs multiple chains of fit_dp_concordance() in parallel using
## parallel::mclapply(). Chains are initialized from K-means clustering
## on the concatenated data and use simulated annealing during burn-in.
##
## Returns individual chain fits plus merged post-burn-in samples.
###############################################################################

#' Initialize sampler parameters from K-means clustering
#'
#' Runs K-means on the transposed concatenated data [Y_g | Y_cf]^T
#' and derives initial theta, v, omega, sigma_k2 from the clusters.
#'
#' @param Y_g   Matrix [p x N_obs]: gDNA signals
#' @param Y_cf  Matrix [p x N_obs]: cfDNA signals
#' @param K     Integer: number of mixture components (= number of clusters)
#' @param T_max Integer: number of time points
#' @param seed  Integer: random seed for K-means
#' @return List with theta, v, omega_g, omega_cf, sigma_k2
initialize_from_kmeans <- function(Y_g, Y_cf, K, T_max = 1, seed = 42) {
  set.seed(seed)
  p <- nrow(Y_g)
  N_obs <- ncol(Y_g)

  # Concatenate gDNA and cfDNA, transpose to [2*N_obs x p]
  Y_combined <- t(cbind(Y_g, Y_cf))

  # Run K-means
  km <- kmeans(Y_combined, centers = K, nstart = 25, iter.max = 100)

  # Cluster centers -> initial theta [p x K] (static, time-invariant)
  centers <- t(km$centers)  # [p x K]
  theta <- array(0, dim = c(p, K))
  theta <- centers

  # Cluster proportions -> initial pi -> initial v (inverse stick-breaking)
  pi_vec <- km$size / sum(km$size)
  # Sort by decreasing weight for consistency
  ord <- order(pi_vec, decreasing = TRUE)
  pi_vec <- pi_vec[ord]
  theta <- theta[, ord, drop = FALSE]

  v <- numeric(K)
  cum <- 1
  for (k in 1:(K-1)) {
    v[k] <- pi_vec[k] / cum
    v[k] <- max(min(v[k], 1 - 1e-10), 1e-10)
    cum <- cum * (1 - v[k])
  }
  v[K] <- 1

  # omega: use pi_vec for all observations (will be updated in first iteration)
  omega_g <- matrix(pi_vec, K, N_obs)
  omega_cf <- matrix(pi_vec, K, N_obs)

  # Within-cluster variance -> initial sigma_k2
  sigma_k2 <- numeric(K)
  for (k in 1:K) {
    members <- which(km$cluster == ord[k])  # account for reordering
    if (length(members) > 1) {
      sigma_k2[k] <- mean(apply(Y_combined[members, , drop = FALSE], 2, var))
    } else {
      sigma_k2[k] <- 1.0
    }
    sigma_k2[k] <- max(sigma_k2[k], 0.01)  # floor to avoid degeneracy
  }

  # K-means centres initialize the POPULATION signature; every subject starts
  # at the population value and the sampler separates them as the per-subject
  # data warrant (writeup Remark rem:anneal).
  list(theta_bar = theta,
       v = v[1:(K-1)],  # v has K-1 entries (v_K=1 set in sampler)
       omega_g = omega_g, omega_cf = omega_cf,
       sigma_k2 = sigma_k2)
}

#' Fit the DP mixture model with multiple parallel chains
#'
#' @param Y_g, Y_cf, patient, time  Data arguments (passed to fit_dp_concordance)
#' @param n_chains  Integer: number of parallel chains (default 4)
#' @param seed      Integer: base seed; chain c uses seed + (c-1)*1000
#' @param mc.cores  Integer: number of cores for parallel execution
#' @param ...       Additional arguments passed to fit_dp_concordance
#'
#' @return List with:
#'   $chains   — list of individual fit objects
#'   $merged   — merged posterior samples (all chains concatenated)
#'   $n_chains — number of chains
fit_multi_chain <- function(Y_g, Y_cf, patient, time,
                             n_chains = 4,
                             seed = 42,
                             mc.cores = min(n_chains, parallel::detectCores() - 1),
                             use_kmeans_init = TRUE,
                             ...) {

  cat(sprintf("Running %d chains on %d cores\n", n_chains, mc.cores))

  seeds <- seed + (seq_len(n_chains) - 1) * 1000

  # K-means initialization: per-chain with different seeds so chains are
  # overdispersed across starting basins. This makes split-Rhat meaningful
  # and reduces the chance that all chains land in the same local mode.
  dots <- list(...)
  K <- if (!is.null(dots$K)) dots$K else 10
  T_max_time <- max(time)

  km_inits <- NULL
  if (use_kmeans_init) {
    cat("  Computing per-chain K-means initializations (K =", K, ")...\n")
    km_inits <- vector("list", n_chains)
    for (i in seq_len(n_chains)) {
      km_inits[[i]] <- initialize_from_kmeans(
        Y_g, Y_cf, K, T_max = T_max_time, seed = seeds[i]
      )
      cat(sprintf("    Chain %d init cluster sizes: %s %%\n", i,
                  paste(round(stick_break(c(km_inits[[i]]$v, 1)) * 100, 1),
                        collapse = " ")))
    }
  }

  fit_one <- function(i) {
    fit_dp_concordance(
      Y_g = Y_g, Y_cf = Y_cf,
      patient = patient, time = time,
      seed = seeds[i],
      verbose = FALSE,
      init = if (!is.null(km_inits)) km_inits[[i]] else NULL,
      ...
    )
  }

  # Run chains in parallel
  chains <- parallel::mclapply(seq_len(n_chains), fit_one, mc.cores = mc.cores)

  # Check for errors
  for (i in seq_along(chains)) {
    if (inherits(chains[[i]], "try-error") || is.null(chains[[i]])) {
      stop(sprintf("Chain %d failed. Error: %s", i, as.character(chains[[i]])))
    }
  }

  # Merge samples across chains
  merged <- merge_chain_samples(chains)

  list(
    chains = chains,
    merged = merged,
    n_chains = n_chains,
    seeds = seeds
  )
}

#' Merge posterior samples from multiple chains
#'
#' Concatenates post-burn-in samples from all chains into a single
#' samples list with the same structure as a single-chain fit.
merge_chain_samples <- function(chains) {
  n_chains <- length(chains)
  ref <- chains[[1]]$samples

  merged <- list()

  # Scalar parameters: concatenate vectors
  ## "phi" and "lambda" are the D2 cross-source coupling draws -- the paper's
  ## tracking estimand. Omitting them made merged$phi silently NULL, and
  ## compute_sim_metrics() falls back through lambda to the FIXED kappa, which
  ## would report a constant tracking correlation with no warning.
  for (nm in c("kappa", "alpha_dp", "K_plus", "K_plus_occupied", "sigma_0",
               "varsigma_0", "b_sigma", "loglik", "phi", "lambda")) {
    if (!is.null(ref[[nm]])) {
      merged[[nm]] <- do.call(c, lapply(chains, function(ch) ch$samples[[nm]]))
    }
  }

  # Matrix parameters: rbind
  for (nm in c("omega_0", "sigma_k", "varsigma_k", "tau_k", "v", "v_g", "v_cf",
               "omega_g_trace", "omega_cf_trace")) {
    if (!is.null(ref[[nm]])) {
      merged[[nm]] <- do.call(rbind, lapply(chains, function(ch) ch$samples[[nm]]))
    }
  }

  # Pointwise log-likelihood (if stored)
  if (!is.null(ref$pointwise_ll)) {
    merged$pointwise_ll <- do.call(rbind,
      lapply(chains, function(ch) ch$samples$pointwise_ll))
  }

  merged
}

#' Extract mcmc.list object for convergence diagnostics (coda format)
#'
#' Creates a coda::mcmc.list from the multi-chain fit, suitable for
#' computing Rhat and ESS via the posterior package.
#'
#' @param multi_fit Output of fit_multi_chain()
#' @param params Character vector of parameter names to include
#' @return coda::mcmc.list object
extract_mcmc_list <- function(multi_fit, params = c("K_plus", "sigma_0")) {
  chains <- multi_fit$chains
  n_chains <- length(chains)

  mcmc_chains <- list()
  for (i in seq_len(n_chains)) {
    ch <- chains[[i]]
    n_save <- ch$n_save
    thin <- ch$thin

    # Build matrix of monitored parameters
    mat <- matrix(NA, n_save, 0)
    col_names <- c()

    for (nm in params) {
      val <- ch$samples[[nm]]
      if (is.matrix(val)) {
        # Matrix parameter (e.g., omega_0): include all columns
        mat <- cbind(mat, val)
        col_names <- c(col_names, paste0(nm, "[", 1:ncol(val), "]"))
      } else {
        mat <- cbind(mat, val)
        col_names <- c(col_names, nm)
      }
    }
    colnames(mat) <- col_names

    mcmc_chains[[i]] <- coda::mcmc(mat, start = ch$n_burn + 1,
                                     thin = thin)
  }

  coda::mcmc.list(mcmc_chains)
}


## ---------------------------------------------------------------------------
## RETIRED ESTIMAND
##
## `theta_concordance_summary()` and `theta_concordance_multi_chain()` computed
## rho^theta_k = cor(theta_g[,k], theta_cf[,k]) -- the cross-source correlation
## of SOURCE-SPECIFIC signatures. Under the present specification a single
## theta is shared by both sources (writeup Sec 4.1), so there is no cf-vs-g
## signature pair to correlate and the quantity does not exist. It has been
## removed rather than redefined: all cross-source information now lives in
## the omega weights, which is what makes weight concordance the estimand.
##
## Callers that previously used rho^theta should use the omega-based
## concordance metrics (R/04_concordance_metrics.R, R/21_concordance_hierarchy.R).
##
## The two functions below replace it with the diagnostics the revised
## specification actually requires.
## ---------------------------------------------------------------------------

#' Check the Assumption A1 separation condition empirically.
#'
#' Assumption A1 (writeup Sec 8) requires a quantitative MARGIN, not merely
#' that the population signatures differ. Specifically there must exist
#' delta > 0 such that for each occupied pair (k, k') the set
#'
#'   J_kk'(delta) = { j : |thetabar_jk - thetabar_jk'| >= delta }
#'
#' has at least 3K elements, and varsigma_k^2 < delta^2 / 4 for every
#' occupied k.
#'
#' This function inverts that: for each pair it finds the LARGEST delta whose
#' set still has >= 3K genes (i.e. the 3K-th largest absolute gap), and
#' compares varsigma_k^2 against delta^2/4. Reporting the achievable delta is
#' more informative than a bare pass/fail, because it says how much separation
#' the fit actually supports.
#'
#' An earlier version took an infimum over the set where the signatures merely
#' differ. That is generically 0 (continuous signatures produce arbitrarily
#' small gaps), so the check would essentially always fail regardless of the
#' fit -- an uninformative diagnostic.
#'
#' @param fit    A fit from fit_dp_concordance() (or one chain of fit_multi_chain()).
#' @param k_occ  Optional occupied components. Defaults to posterior-mean weight >= 0.01.
#' @param n_req  Cardinality required of J_kk'(delta). Defaults to 3*K_plus,
#'               where K_plus is the number of OCCUPIED components -- that is
#'               what the Kruskal-rank argument consumes (three groups, one
#'               discriminating gene per occupied pair per group). Using the
#'               truncation level K instead would demand discriminating genes
#'               for empty components, which are horseshoe-shrunk to zero and
#'               have no gaps, forcing delta down into the noise floor and
#'               making the check fail regardless of the fit.
#' @return data.frame, one row per occupied pair, with the achievable delta,
#'   the implied bound delta^2/4, max(varsigma^2) over the pair, and the verdict.
#'   Attribute "all_satisfied" summarises; NA if undeterminable.
theta_separation_check <- function(fit, k_occ = NULL, n_req = NULL) {
  tb <- fit$theta_bar_postmean
  if (is.null(tb)) stop("theta_bar_postmean not found; refit with the current sampler.")
  ## varsigma_k stores SDs; the posterior mean VARIANCE is E[s^2], not (E[s])^2.
  ## Absent field must ERROR, not fall through to NA: max(NA, na.rm=TRUE) is
  ## -Inf and -Inf < bound is TRUE, which would report the condition as
  ## SATISFIED on a fit carrying no dispersion estimate at all.
  if (is.null(fit$samples$varsigma_k))
    stop("theta_separation_check: fit has no samples$varsigma_k; ",
         "refit with the current sampler.")
  vs <- colMeans(fit$samples$varsigma_k^2, na.rm = TRUE)
  K  <- ncol(tb)

  if (is.null(k_occ)) {
    ## Occupancy MUST be judged in the same component order as theta_bar_postmean,
    ## i.e. the CANONICAL (decreasing-weight) order the sampler accumulates in.
    ##
    ## `final_omega_g` is the LAST RAW ITERATION in unrelabelled order, so using
    ## it selects the wrong components: on the production M1 fit it returned
    ## k = 4, 5, 10 while the occupied components in canonical order are 1, 2, 3.
    ## The check was then comparing theta_bar columns for horseshoe-shrunk empty
    ## components -- whose gaps are noise -- and reporting A1 as violated.
    ##
    ## Use the posterior-mean cohort weight from the omega trace, which is stored
    ## canonically (see the omega_g_trace block in the sampler).
    ogt <- fit$samples$omega_g_trace
    if (!is.null(ogt) && !is.null(fit$N_obs) && fit$N_obs > 0) {
      N_obs_f <- fit$N_obs
      K_tr    <- ncol(ogt) %/% N_obs_f
      wbar    <- vapply(seq_len(K_tr), function(k)
                   mean(ogt[, ((k - 1L) * N_obs_f + 1L):(k * N_obs_f)], na.rm = TRUE),
                   numeric(1))
      k_occ <- which(is.finite(wbar) & wbar >= 0.01)
    } else if (!is.null(fit$final_omega_g)) {
      ## Fallback for fits without a stored trace. Flag it: the ordering may not
      ## match theta_bar, so the result should be treated as unverified.
      warning("theta_separation_check: no omega_g_trace; falling back to ",
              "final_omega_g, whose component order may not match theta_bar. ",
              "Interpret the A1 verdict with caution.", call. = FALSE)
      k_occ <- which(rowMeans(fit$final_omega_g) >= 0.01)
    } else {
      k_occ <- seq_len(min(K, 3))
    }
  }
  k_occ <- k_occ[k_occ <= K & k_occ <= length(vs)]
  ## Cardinality requirement. Assumption A1 asks for |J_kk'(delta)| >= 3*m, where
  ## m is CALIBRATED to the cohort size and the target error, NOT set to K_plus:
  ##
  ##   m >= log(3 * n * choose(K_plus, 2) / eps) / log(1 / q),   q = Phi(-sqrt(2))
  ##
  ## The union bound in the identifiability proof runs over the three Kruskal
  ## groups x choose(K_plus,2) pairs x n SUBJECTS, because the label permutation
  ## must be common across subjects. Taking m = K_plus ignores the factor n and
  ## gives a vacuous bound: at n = 26, K_plus = 2 it is 3*26*q^2 = 0.48.
  ## Using 3*K_plus here would let a fit PASS this check while violating the
  ## assumption the theory actually uses (the 6th largest gap instead of the 12th).
  n_subj <- NULL
  if (!is.null(fit$samples$omega_0) && !is.null(fit$patient))
    n_subj <- length(unique(fit$patient))
  if (is.null(n_subj) && !is.null(fit$theta_0_postmean))
    n_subj <- ncol(fit$theta_0_postmean)   # theta_0 is [p x n]
  if (is.null(n_subj)) n_subj <- 26L       # cohort size of the application
  Kp  <- max(length(k_occ), 2L)
  eps <- 0.01
  q   <- stats::pnorm(-sqrt(2))            # 0.07865
  m_req <- ceiling(log(3 * n_subj * choose(Kp, 2) / eps) / log(1 / q))
  m_req <- max(m_req, 1L)
  if (is.null(n_req)) n_req <- 3L * as.integer(m_req)
  if (length(k_occ) < 2) {
    out <- data.frame(k = NA_integer_, k_prime = NA_integer_, delta = NA_real_,
                      sep_bound = NA_real_, varsigma2 = NA_real_, satisfied = NA)
    attr(out, "all_satisfied") <- NA
    return(out)
  }

  rows <- list()
  for (a in seq_along(k_occ)) for (b in seq_along(k_occ)) {
    if (a >= b) next
    ka <- k_occ[a]; kb <- k_occ[b]
    gaps <- sort(abs(tb[, ka] - tb[, kb]), decreasing = TRUE)
    if (length(gaps) < n_req) next
    ## Largest delta whose J_kk'(delta) still has >= n_req members.
    delta <- gaps[n_req]
    sep_bound <- delta^2 / 4
    vmax <- max(vs[c(ka, kb)])
    if (!is.finite(vmax) || !is.finite(sep_bound)) next
    rows[[length(rows) + 1L]] <- data.frame(
      k = ka, k_prime = kb,
      delta = delta,
      sep_bound = sep_bound,
      varsigma2 = vmax,
      satisfied = isTRUE(vmax < sep_bound)
    )
  }
  out <- if (length(rows)) do.call(rbind, rows) else
    data.frame(k = NA_integer_, k_prime = NA_integer_, delta = NA_real_,
               sep_bound = NA_real_, varsigma2 = NA_real_, satisfied = NA)
  out$n_req <- n_req
  out$m_req <- m_req
  out$n_subj <- n_subj
  attr(out, "all_satisfied") <- if (all(is.na(out$satisfied))) NA else
                                 all(out$satisfied, na.rm = TRUE)
  attr(out, "n_req") <- n_req
  attr(out, "m_req") <- m_req
  attr(out, "eps")   <- eps
  out
}

#' Within-subject label-switch rate.
#'
#' With subject-indexed signatures, a chain can settle on DIFFERENT component
#' orderings for different subjects (subject i's component 1 matching subject
#' i''s component 2). Every per-chain diagnostic would look healthy while the
#' component label loses its cross-subject meaning, invalidating cohort-level
#' statements about pi or thetabar. This is the diagnostic counterpart of the
#' Assumption A1 separation condition; see writeup Sec 10 and Scenario 5.
#'
#' For each subject we match that subject's components to the cohort population
#' signatures by maximum absolute correlation, giving a per-subject implied
#' permutation. We then report the fraction of subjects whose permutation
#' disagrees with the MODAL permutation across subjects.
#'
#' Comparing against the modal permutation rather than the identity is
#' deliberate. What invalidates a cohort statement is subjects disagreeing with
#' EACH OTHER; a single global relabel (every subject permuted the same way)
#' leaves the label perfectly consistent across subjects and is harmless, but
#' scores 1.0 against the identity. The modal comparison scores it 0, which is
#' the interpretable answer.
#'
#' Degenerate cases are reported as indeterminate (NA) rather than as
#' switching: if any subject has a constant (numerically empty) signature, or
#' if two cohort signatures are so similar that the matching is ambiguous, the
#' correlation argmax is arbitrary and a rate computed from it is meaningless.
#'
#' @param fit   A fit from fit_dp_concordance().
#' @param k_occ Occupied components to check (defaults as in theta_separation_check).
#' @param tie_tol Relative tolerance for declaring a correlation match ambiguous.
#' @return List with $rate (fraction of subjects disagreeing with the modal
#'   permutation; NA if indeterminate), $switched (integer subject indices),
#'   $n_subj, $n_indeterminate, and $modal_perm.
theta_label_switch_rate <- function(fit, k_occ = NULL, tie_tol = 0.02) {
  th <- fit$theta_postmean       # [p x K x n]
  tb <- fit$theta_bar_postmean   # [p x K]
  if (is.null(th) || is.null(tb))
    stop("theta_postmean / theta_bar_postmean not found; refit with the current sampler.")
  dims <- dim(th); K <- dims[2]; n <- dims[3]

  if (is.null(k_occ)) {
    om <- fit$final_omega_g
    k_occ <- if (!is.null(om)) which(rowMeans(om) >= 0.01) else seq_len(min(K, 3))
  }
  k_occ <- k_occ[k_occ <= K]
  if (length(k_occ) < 2) {
    return(list(rate = 0, switched = integer(0), n_subj = n,
                note = "fewer than 2 occupied components; switching undefined"))
  }

  ## Guard against ambiguous cohort signatures BEFORE matching: if two
  ## population signatures are nearly identical, which.max picks the first
  ## every time, forcing duplicates and a spurious rate of 1.0.
  L <- length(k_occ)
  for (a in seq_len(L)) for (b in seq_len(L)) {
    if (a >= b) next
    ya <- tb[, k_occ[a]]; yb <- tb[, k_occ[b]]
    if (sd(ya) <= 1e-10 || sd(yb) <= 1e-10 ||
        abs(abs(cor(ya, yb)) - 1) < tie_tol) {
      return(list(rate = NA_real_, switched = integer(0), n_subj = n,
                  n_indeterminate = n, modal_perm = NA,
                  note = paste0("cohort signatures ", k_occ[a], " and ", k_occ[b],
                                " are not distinguishable (|cor| ~ 1); the ",
                                "label-switch rate is undefined. This is itself ",
                                "an Assumption A1 failure.")))
    }
  }

  perms <- vector("list", n)
  indet <- integer(0)
  for (i in seq_len(n)) {
    # Match each of subject i's occupied signatures to each cohort signature by
    # maximum |correlation|; the argmax defines the subject's implied ordering.
    Cmat <- matrix(NA_real_, L, L)
    degenerate <- FALSE
    for (a in seq_len(L)) for (b in seq_len(L)) {
      x <- th[, k_occ[a], i]; y <- tb[, k_occ[b]]
      if (sd(x) <= 1e-10) degenerate <- TRUE
      Cmat[a, b] <- if (sd(x) > 1e-10 && sd(y) > 1e-10) abs(cor(x, y)) else NA_real_
    }
    if (degenerate || anyNA(Cmat)) { indet <- c(indet, i); next }
    perms[[i]] <- as.integer(apply(Cmat, 1, which.max))
  }

  ok <- which(!vapply(perms, is.null, logical(1)))
  if (!length(ok))
    return(list(rate = NA_real_, switched = integer(0), n_subj = n,
                n_indeterminate = length(indet), modal_perm = NA,
                note = "no subject had a well-defined signature match"))

  ## Modal permutation across subjects. A permutation shared by every subject
  ## is a global relabel and is harmless; only disagreement BETWEEN subjects
  ## breaks cohort-level statements about pi or thetabar.
  keys <- vapply(perms[ok], paste, character(1), collapse = "-")
  modal_key <- names(sort(table(keys), decreasing = TRUE))[1]
  switched <- ok[keys != modal_key]

  list(rate = length(switched) / length(ok),
       switched = switched,
       n_subj = n,
       n_indeterminate = length(indet),
       modal_perm = as.integer(strsplit(modal_key, "-", fixed = TRUE)[[1]]),
       note = if (length(indet))
         paste0(length(indet), " subject(s) indeterminate and excluded") else NULL)
}

###############################################################################
## ---- 07_convergence_diagnostics.R
###############################################################################

###############################################################################
## 07_convergence_diagnostics.R
## Convergence diagnostics for multi-chain MCMC output
##
## Computes rank-normalized split-Rhat (Vehtari et al. 2021),
## bulk ESS, and tail ESS using the posterior package.
##
## Returns pass/fail flag based on Rhat < 1.01 threshold.
###############################################################################

suppressPackageStartupMessages({
  library(posterior)
  library(coda)
})

#' Compute convergence diagnostics for a multi-chain fit
#'
#' @param multi_fit Output of fit_multi_chain()
#' @param params Character vector of scalar parameter names to monitor
#' @param rhat_threshold Numeric: maximum acceptable Rhat (default 1.01)
#'
#' @return List with:
#'   $summary     — data.frame with Rhat, bulk_ess, tail_ess per parameter
#'   $all_pass    — logical: TRUE if all Rhat < threshold
#'   $worst_rhat  — worst Rhat value across all parameters
#'   $diagnostics — full posterior::summarise_draws output
compute_convergence <- function(multi_fit,
                                 params = c("K_plus", "sigma_0"),
                                 rhat_threshold = 1.01) {

  chains <- multi_fit$chains
  n_chains <- length(chains)

  if (n_chains < 2) {
    warning("Need >= 2 chains for Rhat computation. Returning NA diagnostics.")
    return(list(
      summary = data.frame(parameter = params, rhat = NA, bulk_ess = NA, tail_ess = NA),
      all_pass = NA,
      worst_rhat = NA
    ))
  }

  # Build draws_array [iterations x chains x variables]
  n_save <- chains[[1]]$n_save

  # Collect parameter names and values
  var_names <- c()
  draws_list <- list()

  for (nm in params) {
    val <- chains[[1]]$samples[[nm]]
    if (is.matrix(val)) {
      ncols <- ncol(val)
      for (col in 1:ncols) {
        var_name <- paste0(nm, "[", col, "]")
        var_names <- c(var_names, var_name)
        arr <- matrix(NA, n_save, n_chains)
        for (ch in 1:n_chains) {
          arr[, ch] <- chains[[ch]]$samples[[nm]][, col]
        }
        draws_list[[var_name]] <- arr
      }
    } else {
      var_names <- c(var_names, nm)
      arr <- matrix(NA, n_save, n_chains)
      for (ch in 1:n_chains) {
        arr[, ch] <- chains[[ch]]$samples[[nm]]
      }
      draws_list[[nm]] <- arr
    }
  }

  # Build draws_array [iterations, chains, variables]
  n_vars <- length(var_names)
  draws_arr <- array(NA, dim = c(n_save, n_chains, n_vars),
                     dimnames = list(NULL, paste0("chain", 1:n_chains), var_names))
  for (v in 1:n_vars) {
    draws_arr[, , v] <- draws_list[[var_names[v]]]
  }

  # Convert to posterior draws object
  draws <- posterior::as_draws_array(draws_arr)

  # Compute diagnostics
  diag <- posterior::summarise_draws(draws,
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail
  )

  summary_df <- data.frame(
    parameter = diag$variable,
    rhat = diag$rhat,
    bulk_ess = diag$ess_bulk,
    tail_ess = diag$ess_tail,
    stringsAsFactors = FALSE
  )

  worst_rhat <- max(summary_df$rhat, na.rm = TRUE)
  all_pass <- all(summary_df$rhat < rhat_threshold, na.rm = TRUE)

  list(
    summary = summary_df,
    all_pass = all_pass,
    worst_rhat = worst_rhat,
    rhat_threshold = rhat_threshold,
    diagnostics = diag
  )
}

#' Print convergence diagnostic summary
print_convergence <- function(conv, verbose = TRUE) {
  if (verbose) {
    cat("\n=== Convergence Diagnostics ===\n")
    cat(sprintf("Rhat threshold: %.3f\n", conv$rhat_threshold))
    cat(sprintf("Worst Rhat: %.4f\n", conv$worst_rhat))
    cat(sprintf("All pass: %s\n\n", ifelse(conv$all_pass, "YES", "NO")))

    print(conv$summary, row.names = FALSE, digits = 4)
    cat("\n")

    # Flag any failures
    failures <- conv$summary[conv$summary$rhat >= conv$rhat_threshold, ]
    if (nrow(failures) > 0) {
      cat("WARNING: The following parameters have Rhat >= threshold:\n")
      print(failures, row.names = FALSE)
    }
  }
  invisible(conv)
}

#' Identify chains trapped in low-density posterior modes via log-likelihood spread
#'
#' Returns indices of chains whose mean post-burn-in log-likelihood is more than
#' `ll_gap` units below the best chain's mean. Useful for excluding trapped
#' chains from WAIC/LOO aggregation.
#'
#' @param multi_fit  Output of fit_multi_chain()
#' @param ll_gap     Threshold in log-units (default 10000)
#' @return List with $keep (indices to keep), $drop (indices to drop), $ll_means
identify_trapped_chains <- function(multi_fit, ll_gap = 10000) {
  chains <- multi_fit$chains
  n_chains <- length(chains)
  ll_means <- sapply(chains, function(ch) mean(ch$samples$loglik))
  best <- max(ll_means)
  drop <- which((best - ll_means) > ll_gap)
  keep <- setdiff(seq_len(n_chains), drop)
  list(keep = keep, drop = drop, ll_means = ll_means, best = best, ll_gap = ll_gap)
}

#' Print log-likelihood spread across chains and flag trapped chains
print_ll_spread <- function(trap, label = "") {
  pfx <- if (nchar(label) > 0) paste0(label, ": ") else ""
  cat(sprintf("\n=== %sChain Log-Likelihood Spread ===\n", pfx))
  for (i in seq_along(trap$ll_means)) {
    status <- if (i %in% trap$drop) "TRAPPED" else "OK"
    cat(sprintf("  Chain %d: mean LL = %.1f  (gap %.1f) [%s]\n",
                i, trap$ll_means[i], trap$best - trap$ll_means[i], status))
  }
  cat(sprintf("Threshold: gap > %.0f log-units flagged as trapped.\n", trap$ll_gap))
  cat(sprintf("Keeping chains: %s\n", paste(trap$keep, collapse = ",")))
  if (length(trap$drop) > 0) {
    cat(sprintf("Dropping chains: %s\n", paste(trap$drop, collapse = ",")))
  }
  invisible(trap)
}

#' Compute WAIC and LOO on a subset of chains (pointwise log-likelihood)
#'
#' Aggregates pointwise log-likelihoods from the selected chains only,
#' then runs loo::waic() and loo::loo() on the combined matrix.
#'
#' @param multi_fit       Output of fit_multi_chain() (must have pointwise_ll stored)
#' @param chain_indices   Integer vector: which chains to include
#' @return List with $waic, $loo (or NULL if computation fails)
waic_loo_subset <- function(multi_fit, chain_indices) {
  if (!requireNamespace("loo", quietly = TRUE)) {
    stop("loo package required")
  }
  pws <- lapply(chain_indices, function(i) multi_fit$chains[[i]]$samples$pointwise_ll)
  if (any(sapply(pws, is.null))) {
    stop("pointwise_ll not stored for all requested chains")
  }
  ll_mat <- do.call(rbind, pws)
  list(
    waic = tryCatch(loo::waic(ll_mat), error = function(e) {
      message("waic failed: ", conditionMessage(e)); NULL
    }),
    loo  = tryCatch(loo::loo(ll_mat),  error = function(e) {
      message("loo failed: ",  conditionMessage(e)); NULL
    }),
    n_chains = length(chain_indices),
    chain_indices = chain_indices
  )
}


###############################################################################
## ---- 08_model_comparison.R
###############################################################################

###############################################################################
## 08_model_comparison.R
## Model comparison tools: WAIC and PSIS-LOO
##
## WAIC and PSIS-LOO use the loo package.
## These are the primary model comparison criteria.
###############################################################################

suppressPackageStartupMessages({
  library(loo)
})

#' Compute WAIC for a fitted model
#'
#' Requires the model to have been fit with store_pointwise_ll = TRUE.
#' The pointwise log-likelihood matrix has dimensions [n_save x (2*N_obs)]:
#'   columns 1:N_obs = gDNA observations (sum over genes)
#'   columns (N_obs+1):(2*N_obs) = cfDNA observations (sum over genes)
#'
#' @param fit Output of fit_dp_concordance() with store_pointwise_ll=TRUE
#' @return loo::waic object
compute_waic <- function(fit) {
  if (is.null(fit$samples$pointwise_ll)) {
    stop("Pointwise log-likelihood not stored. Re-fit with store_pointwise_ll=TRUE.")
  }

  ll_matrix <- fit$samples$pointwise_ll
  loo::waic(ll_matrix)
}

#' Compute PSIS-LOO for a fitted model
#'
#' @param fit Output of fit_dp_concordance() with store_pointwise_ll=TRUE
#' @return loo::loo object
compute_psis_loo <- function(fit) {
  if (is.null(fit$samples$pointwise_ll)) {
    stop("Pointwise log-likelihood not stored. Re-fit with store_pointwise_ll=TRUE.")
  }

  ll_matrix <- fit$samples$pointwise_ll
  loo::loo(ll_matrix)
}

#' Compare two models via WAIC
#'
#' @param fit_M1 Fitted M1 model (with pointwise_ll)
#' @param fit_M0 Fitted M0 model (with pointwise_ll)
#' @return List with delta_waic, se, and interpretation
compare_models_waic <- function(fit_M1, fit_M0) {
  waic_M1 <- compute_waic(fit_M1)
  waic_M0 <- compute_waic(fit_M0)

  comp <- loo::loo_compare(list(M1 = waic_M1, M0 = waic_M0))

  delta_waic <- waic_M0$estimates["waic", "Estimate"] -
                waic_M1$estimates["waic", "Estimate"]

  list(
    waic_M1 = waic_M1,
    waic_M0 = waic_M0,
    delta_waic = delta_waic,  # positive favors M1
    comparison = comp,
    favors = if (delta_waic > 0) "M1 (tracking)" else "M0 (non-tracking)"
  )
}

#' Compare two models via PSIS-LOO
#'
#' @param fit_M1 Fitted M1 model (with pointwise_ll)
#' @param fit_M0 Fitted M0 model (with pointwise_ll)
#' @return List with delta_loo, se, and interpretation
compare_models_loo <- function(fit_M1, fit_M0) {
  loo_M1 <- compute_psis_loo(fit_M1)
  loo_M0 <- compute_psis_loo(fit_M0)

  comp <- loo::loo_compare(list(M1 = loo_M1, M0 = loo_M0))

  delta_loo <- loo_M0$estimates["looic", "Estimate"] -
               loo_M1$estimates["looic", "Estimate"]

  list(
    loo_M1 = loo_M1,
    loo_M0 = loo_M0,
    delta_loo = delta_loo,  # positive favors M1
    comparison = comp,
    favors = if (delta_loo > 0) "M1 (tracking)" else "M0 (non-tracking)"
  )
}

#' Full model comparison report
#'
#' @param fit_M1, fit_M0 Fitted models (with store_pointwise_ll=TRUE)
#' @param verbose Print results
#' @return List with all comparison metrics
compare_models <- function(fit_M1, fit_M0, verbose = TRUE) {
  results <- list()

  # WAIC and LOO (if pointwise ll available)
  has_pw <- !is.null(fit_M1$samples$pointwise_ll) &&
            !is.null(fit_M0$samples$pointwise_ll)

  if (has_pw) {
    results$waic <- compare_models_waic(fit_M1, fit_M0)
    results$loo  <- tryCatch(
      compare_models_loo(fit_M1, fit_M0),
      error = function(e) {
        if (verbose) cat("PSIS-LOO failed:", conditionMessage(e), "\n")
        NULL
      }
    )
  }

  if (verbose) {
    cat("\n=== Model Comparison ===\n")

    if (has_pw) {
      cat(sprintf("Delta WAIC (M0 - M1): %.2f (positive favors M1)\n",
                  results$waic$delta_waic))
      cat(sprintf("  WAIC favors: %s\n", results$waic$favors))

      if (!is.null(results$loo)) {
        cat(sprintf("Delta LOO (M0 - M1): %.2f (positive favors M1)\n",
                    results$loo$delta_loo))
        cat(sprintf("  LOO favors: %s\n", results$loo$favors))
      }
    } else {
      cat("(Pointwise log-lik not stored; WAIC/LOO unavailable)\n")
    }
  }

  invisible(results)
}


###############################################################################
## ---- 09_label_switching.R
###############################################################################

###############################################################################
## 09_label_switching.R
## Label switching resolution for mixture model posterior
##
## Uses a simple pivot-based relabeling: at each iteration, components
## are reordered by decreasing stick-breaking weight pi_k.
## This is appropriate when the primary estimands (kappa, omega_0, K+, WAIC)
## are permutation-invariant and only component-specific quantities
## (theta_k, pi_k) need relabeling.
##
## For more sophisticated approaches (Stephens 2000 KL), the
## label.switching package can be used as a post-processing step.
###############################################################################

#' Relabel components by decreasing stick-breaking weight
#'
#' At each saved iteration, reorders the K components so that
#' pi_1 >= pi_2 >= ... >= pi_K. This resolves the most common
#' form of label switching (weight permutation).
#'
#' Under M1 (shared global pi), gDNA and cfDNA share the same component
#' index space, so a single ordering is applied to both omega traces and to
#' sigma_k (which is also shared across sources by the model).
#'
#' Under M0 with independent source-DPs (pi_g, pi_cf), gDNA and cfDNA have
#' separate stick-breaking weights. Relabeling is performed independently per
#' source: each source's omega trace is reordered by its own pi (pi_g for
#' omega_g, pi_cf for omega_cf). sigma_k indices remain shared by the model
#' specification (same Inv-Gamma prior, same horseshoe), so we reorder
#' sigma_k by the gDNA pi as a convention; reviewers comparing sigma_k traces
#' across chains for M0 should be aware that sigma_k labels follow the gDNA
#' relabeling convention.
#'
#' @param fit Output of fit_dp_concordance()
#' @return Modified fit object with relabeled sigma_k, v (or v_g/v_cf), and
#'   omega samples.
relabel_by_weight <- function(fit) {
  n_save <- fit$n_save
  K <- fit$K
  K_trace <- if (!is.null(fit$K_trace)) fit$K_trace else min(K, 3)
  N_obs <- fit$N_obs
  is_M0 <- !is.null(fit$model) && fit$model == "M0" &&
           !is.null(fit$samples$v_g) && !is.null(fit$samples$v_cf)

  for (s in 1:n_save) {
    if (is_M0) {
      v_g_s  <- c(fit$samples$v_g[s, ],  1)
      v_cf_s <- c(fit$samples$v_cf[s, ], 1)
      pi_g_s  <- stick_break(v_g_s)
      pi_cf_s <- stick_break(v_cf_s)
      ord_g  <- order(pi_g_s,  decreasing = TRUE)
      ord_cf <- order(pi_cf_s, decreasing = TRUE)
    } else {
      v_s <- c(fit$samples$v[s, ], 1)
      pi_s <- stick_break(v_s)
      ord_g  <- order(pi_s, decreasing = TRUE)
      ord_cf <- ord_g
    }

    ## IDEMPOTENCE GUARD. The current sampler already stores every
    ## component-indexed slot -- theta, theta_bar, sigma_k, varsigma_k, tau_k,
    ## v, and the omega traces -- in canonical (decreasing-weight) order, so
    ## `ord_g` computed here is the identity and this loop is a no-op.
    ##
    ## It is retained only for fits produced BEFORE canonical-order storage was
    ## added. Applying it to an already-canonical fit permuted sigma_k /
    ## varsigma_k / tau_k / omega a SECOND time while leaving theta_bar alone
    ## (the fit-level block below could not detect it, because `v` had by then
    ## been overwritten with the sorted weights, making its own ordering the
    ## identity). The result was sigma_k[k] and theta[,k] referring to
    ## different components in every consumer of the standard
    ## relabel -> align -> merge pipeline.
    ##
    ## Skipping on the identity makes the function idempotent, which is the
    ## property the pipeline actually relies on.
    if (identical(as.integer(ord_g), seq_len(K)) &&
        identical(as.integer(ord_cf), seq_len(K))) next

    # Reorder ALL component-indexed scalars by the gDNA ordering
    # (M1: identical for both sources; M0: convention).
    # varsigma_k and tau_k were previously omitted, which left them
    # desynchronised from sigma_k after relabelling.
    fit$samples$sigma_k[s, ] <- fit$samples$sigma_k[s, ord_g]
    if (!is.null(fit$samples$varsigma_k))
      fit$samples$varsigma_k[s, ] <- fit$samples$varsigma_k[s, ord_g]
    if (!is.null(fit$samples$tau_k))
      fit$samples$tau_k[s, ] <- fit$samples$tau_k[s, ord_g]

    # Reorder omega traces by their own source's ordering.
    if (!is.null(fit$samples$omega_g_trace) && !is.null(N_obs)) {
      og_mat <- matrix(NA, K_trace, N_obs)
      oc_mat <- matrix(NA, K_trace, N_obs)
      for (kk in 1:K_trace) {
        col_start <- (kk - 1) * N_obs + 1
        col_end   <- kk * N_obs
        og_mat[kk, ] <- fit$samples$omega_g_trace[s, col_start:col_end]
        oc_mat[kk, ] <- fit$samples$omega_cf_trace[s, col_start:col_end]
      }
      ## When K_trace < K the permutation can point beyond the traced block.
      ## Leaving such a row at its UNPERMUTED value (as an earlier version did)
      ## makes the mapping non-bijective and silently desynchronises the omega
      ## traces from sigma_k / theta, which are permuted over the full K. Set
      ## the row to NA instead so the loss of information is visible rather
      ## than being papered over with a stale value.
      og_reord <- og_mat
      oc_reord <- oc_mat
      for (kk in 1:K_trace) {
        og_reord[kk, ] <- if (ord_g[kk]  <= K_trace) og_mat[ord_g[kk],  ] else NA_real_
        oc_reord[kk, ] <- if (ord_cf[kk] <= K_trace) oc_mat[ord_cf[kk], ] else NA_real_
      }
      for (kk in 1:K_trace) {
        col_start <- (kk - 1) * N_obs + 1
        col_end   <- kk * N_obs
        fit$samples$omega_g_trace[s, col_start:col_end]  <- og_reord[kk, ]
        fit$samples$omega_cf_trace[s, col_start:col_end] <- oc_reord[kk, ]
      }
    }

    # Reorder stick-breaking variables to match the sorted pi.
    inv_stickbreak <- function(pi_sorted) {
      v_new <- numeric(K)
      cum <- 1
      for (k in 1:(K-1)) {
        v_new[k] <- pi_sorted[k] / cum
        v_new[k] <- max(min(v_new[k], 1 - 1e-10), 1e-10)
        cum <- cum * (1 - v_new[k])
      }
      v_new[K] <- 1
      v_new
    }
    if (is_M0) {
      v_g_new  <- inv_stickbreak(pi_g_s[ord_g])
      v_cf_new <- inv_stickbreak(pi_cf_s[ord_cf])
      fit$samples$v_g[s, ]  <- v_g_new[1:(K-1)]
      fit$samples$v_cf[s, ] <- v_cf_new[1:(K-1)]
    } else {
      v_new <- inv_stickbreak(pi_s[ord_g])
      fit$samples$v[s, ] <- v_new[1:(K-1)]
    }
  }

  # ---------------------------------------------------------------------
  # Fit-level (not per-iteration) theta quantities.
  #
  # The sampler already accumulates these in canonical order, so under the
  # current sampler this is a no-op. It matters for (a) fits produced before
  # the canonical-order accumulation was added, and (b) keeping theta in step
  # with the per-iteration permutation applied above -- previously theta was
  # never permuted here at all, so theta[,k] and sigma_k[k] could refer to
  # different components (see writeup Sec 10 on label switching).
  #
  # We use the ordering implied by the posterior-mean stick-breaking weights,
  # which is the fit-level analogue of the per-iteration ord_g.
  # ---------------------------------------------------------------------
  if (!is.null(fit$theta_bar_postmean)) {
    v_bar <- if (is_M0 && !is.null(fit$samples$v_g))
               c(colMeans(fit$samples$v_g, na.rm = TRUE), 1)
             else c(colMeans(fit$samples$v, na.rm = TRUE), 1)
    ord_fit <- order(stick_break(v_bar), decreasing = TRUE)
    if (!identical(ord_fit, seq_len(K))) {
      perm_mat <- function(M) if (is.null(M)) NULL else M[, ord_fit, drop = FALSE]
      perm_arr <- function(A) if (is.null(A)) NULL else A[, ord_fit, , drop = FALSE]
      fit$theta_bar_postmean <- perm_mat(fit$theta_bar_postmean)
      fit$theta_bar_postvar  <- perm_mat(fit$theta_bar_postvar)
      fit$theta_postmean     <- perm_arr(fit$theta_postmean)
      fit$theta_postvar      <- perm_arr(fit$theta_postvar)
      if (!is.null(fit$theta_bar_trace))
        fit$theta_bar_trace <- fit$theta_bar_trace[, , ord_fit, drop = FALSE]
      if (!is.null(fit$theta_trace))
        fit$theta_trace <- fit$theta_trace[, , ord_fit, , drop = FALSE]
    }
  }

  fit
}

#' Align component labels across chains using posterior mean omega signatures
#'
#' After per-chain relabeling by weight, components across chains should
#' already be roughly aligned (both sorted by decreasing pi). This function
#' provides an additional alignment step by matching components across chains
#' based on their posterior mean omega_g signatures, using chain 1 as reference.
#'
#' @param chains List of fit objects (already relabeled by weight)
#' @return List of fit objects with cross-chain aligned labels
align_chains_by_signature <- function(chains) {
  n_chains <- length(chains)
  if (n_chains < 2) return(chains)

  ref <- chains[[1]]
  K_trace <- if (!is.null(ref$K_trace)) ref$K_trace else min(ref$K, 3)
  N_obs <- ref$N_obs
  n_save <- ref$n_save

  # Limit alignment to top K_align components to avoid K! permutation explosion
  K_align <- min(K_trace, 5)

  if (is.null(ref$samples$omega_g_trace) || K_align < 2) return(chains)

  # Compute reference chain's posterior mean omega_g for top K_align components
  ref_means <- matrix(NA, K_align, N_obs)
  for (kk in 1:K_align) {
    col_start <- (kk - 1) * N_obs + 1
    col_end   <- kk * N_obs
    ref_means[kk, ] <- colMeans(ref$samples$omega_g_trace[, col_start:col_end, drop = FALSE])
  }

  # Generate permutations of 1:K_align only (not 1:K_trace, which could be K!=10!)
  perms <- combinat_perms(K_align)

  for (ch_idx in 2:n_chains) {
    ch <- chains[[ch_idx]]

    # Compute this chain's posterior mean omega_g for top K_align components
    ch_means <- matrix(NA, K_align, N_obs)
    for (kk in 1:K_align) {
      col_start <- (kk - 1) * N_obs + 1
      col_end   <- kk * N_obs
      ch_means[kk, ] <- colMeans(ch$samples$omega_g_trace[, col_start:col_end, drop = FALSE])
    }

    # Find best permutation of top K_align components (using gDNA signatures).
    best_perm_align_g <- 1:K_align
    best_cost <- Inf
    for (p in perms) {
      cost <- sum((ch_means[p, ] - ref_means)^2)
      if (cost < best_cost) {
        best_cost <- cost
        best_perm_align_g <- p
      }
    }

    # Under M0 (source-DP), align cfDNA labels independently using cfDNA means.
    is_M0 <- !is.null(ch$model) && ch$model == "M0" &&
             !is.null(ch$samples$v_g) && !is.null(ch$samples$v_cf)
    best_perm_align_cf <- best_perm_align_g  # default for M1: shared labels
    if (is_M0) {
      ref_means_cf <- matrix(NA, K_align, N_obs)
      ch_means_cf <- matrix(NA, K_align, N_obs)
      for (kk in 1:K_align) {
        col_start <- (kk - 1) * N_obs + 1
        col_end   <- kk * N_obs
        ref_means_cf[kk, ] <- colMeans(ref$samples$omega_cf_trace[, col_start:col_end, drop = FALSE])
        ch_means_cf[kk, ]  <- colMeans(ch$samples$omega_cf_trace[, col_start:col_end, drop = FALSE])
      }
      best_cost_cf <- Inf
      for (p in perms) {
        cost <- sum((ch_means_cf[p, ] - ref_means_cf)^2)
        if (cost < best_cost_cf) {
          best_cost_cf <- cost
          best_perm_align_cf <- p
        }
      }
    }

    # Extend to full K_trace permutations (identity for positions > K_align).
    best_perm_g  <- 1:K_trace
    best_perm_cf <- 1:K_trace
    best_perm_g[1:K_align]  <- best_perm_align_g
    best_perm_cf[1:K_align] <- best_perm_align_cf

    apply_g  <- !identical(best_perm_g,  1:K_trace)
    apply_cf <- !identical(best_perm_cf, 1:K_trace)

    if (apply_g || apply_cf) {
      if (is_M0) {
        cat(sprintf("  Chain %d: relabeling gDNA top %d %s -> %s; cfDNA %s -> %s\n",
                    ch_idx, K_align,
                    paste(best_perm_g[1:K_align], collapse=","), paste(1:K_align, collapse=","),
                    paste(best_perm_cf[1:K_align], collapse=","), paste(1:K_align, collapse=",")))
      } else {
        cat(sprintf("  Chain %d: relabeling top %d components %s -> %s\n",
                    ch_idx, K_align, paste(best_perm_g[1:K_align], collapse=","),
                    paste(1:K_align, collapse=",")))
      }

      inv_stickbreak <- function(pi_sorted, K_total) {
        v_new <- numeric(K_total)
        cum <- 1
        for (k in 1:(K_total - 1)) {
          v_new[k] <- pi_sorted[k] / cum
          v_new[k] <- max(min(v_new[k], 1 - 1e-10), 1e-10)
          cum <- cum * (1 - v_new[k])
        }
        v_new[K_total] <- 1
        v_new
      }

      # Fit-level theta quantities must receive the SAME permutation as the
      # per-iteration quantities below. Omitting this (as earlier versions did)
      # left theta[,k] and sigma_k[k] referring to different components after
      # cross-chain alignment -- which silently corrupted any downstream code
      # that combined them, including the held-out predictive scorers.
      if (!identical(best_perm_g, 1:K_trace)) {
        pg <- best_perm_g
        if (!is.null(ch$theta_bar_postmean) && ncol(ch$theta_bar_postmean) >= K_trace) {
          ch$theta_bar_postmean[, 1:K_trace] <- ch$theta_bar_postmean[, pg, drop = FALSE]
          if (!is.null(ch$theta_bar_postvar))
            ch$theta_bar_postvar[, 1:K_trace] <- ch$theta_bar_postvar[, pg, drop = FALSE]
        }
        if (!is.null(ch$theta_postmean) && dim(ch$theta_postmean)[2] >= K_trace) {
          ch$theta_postmean[, 1:K_trace, ] <- ch$theta_postmean[, pg, , drop = FALSE]
          if (!is.null(ch$theta_postvar))
            ch$theta_postvar[, 1:K_trace, ] <- ch$theta_postvar[, pg, , drop = FALSE]
        }
        if (!is.null(ch$theta_bar_trace))
          ch$theta_bar_trace[, , 1:K_trace] <- ch$theta_bar_trace[, , pg, drop = FALSE]
        if (!is.null(ch$theta_trace))
          ch$theta_trace[, , 1:K_trace, ] <- ch$theta_trace[, , pg, , drop = FALSE]
      }

      for (s in 1:n_save) {
        # Reorder all component-indexed scalars by the gDNA permutation.
        ch$samples$sigma_k[s, 1:K_trace] <- ch$samples$sigma_k[s, best_perm_g]
        if (!is.null(ch$samples$varsigma_k))
          ch$samples$varsigma_k[s, 1:K_trace] <- ch$samples$varsigma_k[s, best_perm_g]
        if (!is.null(ch$samples$tau_k))
          ch$samples$tau_k[s, 1:K_trace] <- ch$samples$tau_k[s, best_perm_g]

        # Reorder omega traces per source.
        og_tmp <- matrix(NA, K_trace, N_obs)
        oc_tmp <- matrix(NA, K_trace, N_obs)
        for (kk in 1:K_trace) {
          col_start <- (kk - 1) * N_obs + 1
          col_end   <- kk * N_obs
          og_tmp[kk, ] <- ch$samples$omega_g_trace[s, col_start:col_end]
          oc_tmp[kk, ] <- ch$samples$omega_cf_trace[s, col_start:col_end]
        }
        og_perm <- og_tmp[best_perm_g,  , drop = FALSE]
        oc_perm <- oc_tmp[best_perm_cf, , drop = FALSE]
        for (kk in 1:K_trace) {
          col_start <- (kk - 1) * N_obs + 1
          col_end   <- kk * N_obs
          ch$samples$omega_g_trace[s, col_start:col_end]  <- og_perm[kk, ]
          ch$samples$omega_cf_trace[s, col_start:col_end] <- oc_perm[kk, ]
        }

        # Reorder stick-breaking variables: per source under M0, shared under M1.
        if (is_M0) {
          v_g_s  <- c(ch$samples$v_g[s, ],  1)
          v_cf_s <- c(ch$samples$v_cf[s, ], 1)
          pi_g_s  <- stick_break(v_g_s)
          pi_cf_s <- stick_break(v_cf_s)
          pi_g_s[1:K_trace]  <- pi_g_s[best_perm_g]
          pi_cf_s[1:K_trace] <- pi_cf_s[best_perm_cf]
          v_g_new  <- inv_stickbreak(pi_g_s,  ch$K)
          v_cf_new <- inv_stickbreak(pi_cf_s, ch$K)
          ch$samples$v_g[s, ]  <- v_g_new[1:(ch$K - 1)]
          ch$samples$v_cf[s, ] <- v_cf_new[1:(ch$K - 1)]
        } else {
          v_s <- c(ch$samples$v[s, ], 1)
          pi_s <- stick_break(v_s)
          pi_s[1:K_trace] <- pi_s[best_perm_g]
          v_new <- inv_stickbreak(pi_s, ch$K)
          ch$samples$v[s, ] <- v_new[1:(ch$K - 1)]
        }
      }

      chains[[ch_idx]] <- ch
    } else {
      cat(sprintf("  Chain %d: labels already aligned\n", ch_idx))
    }
  }

  chains
}

# Simple permutation generator (no external dependency)
# Returns list of all n! permutations of 1:n as integer vectors
combinat_perms <- function(n) {
  if (n == 1) return(list(1L))
  prev <- combinat_perms(n - 1)
  result <- list()
  for (p in prev) {
    for (pos in 1:n) {
      new_perm <- append(p, n, after = pos - 1)
      result <- c(result, list(as.integer(new_perm)))
    }
  }
  result
}

#' Check for label switching by examining bimodality of component weights
#'
#' If a component's weight trace has two distinct modes (e.g., alternating
#' between 0.3 and 0.7), label switching is occurring.
#'
#' @param fit Output of fit_dp_concordance()
#' @param K_check Number of top components to check (default: 5)
#' @return Logical: TRUE if label switching detected
detect_label_switching <- function(fit, K_check = 5) {
  n_save <- fit$n_save
  K <- fit$K
  is_M0 <- !is.null(fit$model) && fit$model == "M0" &&
           !is.null(fit$samples$v_g)

  pi_traces <- matrix(NA, n_save, K)
  for (s in 1:n_save) {
    if (is_M0) {
      # Use the gDNA stick-breaking by convention for detection under M0.
      v_s <- c(fit$samples$v_g[s, ], 1)
    } else {
      v_s <- c(fit$samples$v[s, ], 1)
    }
    pi_traces[s, ] <- stick_break(v_s)
  }

  # Check top K_check components for bimodality
  # Simple heuristic: if the range of pi_k is > 0.3, likely switching
  switching_detected <- FALSE
  for (k in 1:min(K_check, K)) {
    pi_k <- pi_traces[, k]
    iqr_k <- IQR(pi_k)
    range_k <- diff(range(pi_k))
    if (range_k > 0.4 && iqr_k > 0.15) {
      switching_detected <- TRUE
      break
    }
  }

  switching_detected
}


###############################################################################
## ---- 26_parallel_tempering.R
###############################################################################

###############################################################################
## 26_parallel_tempering.R
##
## Parallel-tempered MCMC driver (Priority 4 of the convergence remediation
## plan). Wraps the existing fit_dp_concordance() sampler with a temperature
## ladder and Metropolis-Hastings swap moves between adjacent rungs.
##
## Background:
## After Rounds 1-3 (symmetric annealing, per-chain K-means init, hierarchical
## sigma_k^2 prior, chain-health check) the M1 sampler still exhibits residual
## discrete-mode multimodality: in some runs one of the four chains lands in a
## populated-component subset that is ~10^4 log-units below the dominant mode
## and is correctly flagged for exclusion by the chain-health check. Parallel
## tempering is the standard escape from this kind of multimodality: each
## logical chain consists of R rungs at temperatures 1 = T_1 < T_2 < ... < T_R;
## the high-T rungs flatten the posterior and can move between modes, and
## periodic Metropolis-Hastings swaps push the state of the T=1 rung between
## modes as well. Only the T=1 rung's samples are reported.
##
## Implementation strategy:
## We run each rung as a short batch of MCMC iterations via the existing
## fit_dp_concordance() API with `fixed_temp = T_r` and warm-start it from
## the previous batch's final state. After each batch we propose a swap of
## states between adjacent rungs with MH acceptance ratio
##   alpha = min(1, exp((beta_r - beta_{r+1}) * (logL_{r+1} - logL_r)))
## where beta_r = 1/T_r and logL_r is the joint log-likelihood at rung r's
## current state. We swap *states* (theta, theta_bar, theta_0, v, omega,
## sigma_k2, sigma_0_2, varsigma, omega_0, b_sigma, horseshoe state), not
## labels, so the T=1 rung continues to
## target the posterior exactly.
##
## API:
##   fit_parallel_tempered(Y_g, Y_cf, patient, time,
##                         temperatures = c(1, 2, 4, 8),
##                         swap_every  = 10,
##                         n_iter, n_burn, thin, ...)
##
## Returns a list with the same structure as fit_dp_concordance() for the
## T=1 rung, plus a `pt_diagnostics` element with swap acceptance rates.
###############################################################################

# (Note: parallel-tempering used to source the standalone sampler file. It is
#  now part of this same library — no source() needed.)

#' Compute joint log-likelihood at a single sampler state.
#'
#' Wraps the existing compute_loglik / compute_loglik_cpp interface.
#' Used to evaluate the MH swap acceptance ratio between rungs.
pt_state_loglik <- function(state, Y_g, Y_cf, patient, time, K, p, T_max,
                            use_rcpp = TRUE) {
  N_obs <- ncol(Y_g)
  if (use_rcpp && exists(".rcpp_available") && .rcpp_available) {
    res <- compute_loglik_cpp(Y_g, Y_cf,
                              as.numeric(state$theta), as.numeric(state$theta_0),
                              state$sigma_k2, state$sigma_0_2,
                              state$omega_g, state$omega_cf,
                              state$omega_0, patient, time, K, p,
                              FALSE)
  } else {
    res <- compute_loglik(Y_g, Y_cf, state$theta, state$theta_0, state$sigma_k2,
                          state$sigma_0_2, state$omega_g, state$omega_cf,
                          state$omega_0, patient, time, K, p,
                          return_pointwise = FALSE)
  }
  res$total
}

#' Fit the DP mixture concordance model with parallel tempering.
#'
#' @param Y_g, Y_cf, patient, time  Data arguments (see fit_dp_concordance)
#' @param temperatures  Numeric vector of rung temperatures, ascending and
#'                      starting at 1. e.g., c(1, 2, 4, 8). Length R >= 2.
#' @param swap_every    Integer: number of iterations per batch before
#'                      proposing swaps between adjacent rungs. Default 10.
#' @param n_iter        Total iterations to run *at each rung*.
#' @param n_burn        Burn-in iterations (samples not saved). Swap proposals
#'                      still occur during burn-in.
#' @param thin          Thinning interval applied to T=1 rung saved samples.
#' @param K, model, seed, kappa_fixed, hyperparams, store_pointwise_ll,
#'        use_rcpp, init   Passed through to fit_dp_concordance.
#' @param verbose       Print swap acceptance rates and batch progress.
#' @return A list with the same structure as fit_dp_concordance() output for
#'         the T=1 rung, plus a `pt_diagnostics` element containing:
#'           $temperatures    the rung schedule
#'           $swap_attempts   matrix [n_batches x (R-1)]: 1 if a swap was
#'                            proposed between rungs (r, r+1) at this batch
#'           $swap_accepted   matrix [n_batches x (R-1)]: 1 if accepted
#'           $acceptance_rate per-pair acceptance rate
fit_parallel_tempered <- function(Y_g, Y_cf, patient, time,
                                   temperatures = c(1, 2, 4, 8),
                                   swap_every = 10,
                                   K = 10,
                                   n_iter = 22500,
                                   n_burn = 20000,
                                   thin = 10,
                                   model = "M1",
                                   seed = 42,
                                   hyperparams = list(),
                                   kappa_fixed = 1.0,
                                   store_pointwise_ll = FALSE,
                                   use_rcpp = TRUE,
                                   init = NULL,
                                   verbose = TRUE) {

  # ---------------------------------------------------------------------
  # Validation
  # ---------------------------------------------------------------------
  stopifnot(is.numeric(temperatures), length(temperatures) >= 2)
  stopifnot(abs(temperatures[1] - 1) < 1e-12)
  stopifnot(all(diff(temperatures) > 0))   # strictly ascending
  stopifnot(swap_every >= 1)
  stopifnot(n_iter %% swap_every == 0)
  stopifnot(n_burn  %% swap_every == 0)

  R_rungs   <- length(temperatures)
  betas     <- 1 / temperatures
  n_batches <- n_iter %/% swap_every
  burn_batches <- n_burn  %/% swap_every

  if (verbose) cat("Parallel tempering:", R_rungs, "rungs at T =",
                   paste(temperatures, collapse = ", "), "; swap every",
                   swap_every, "iters;", n_batches, "batches total\n")

  # ---------------------------------------------------------------------
  # Storage: only the T=1 rung's samples are kept as the official posterior.
  # We pre-allocate sample arrays of the right size.
  # ---------------------------------------------------------------------
  p     <- nrow(Y_g)
  N_obs <- ncol(Y_g)
  n     <- max(patient)
  T_max <- max(time)
  n_save_total <- floor((n_iter - n_burn) / thin)

  # Per-batch state: list of length R_rungs, each containing the warm-start
  # state for that rung. Initially all rungs warm-start from `init`.
  # NB: in R, `rung_state[[r]] <- NULL` would delete element r and shrink
  # the list, so we use single-bracket assignment with list() wrapping to
  # preserve NULL slots when `init` itself is NULL.
  rung_state <- rep(list(init), R_rungs)

  # Accumulators for the T=1 rung
  T1_samples_list <- list()  # filled in batch-by-batch after burn-in

  # Swap diagnostics
  swap_attempts <- matrix(0L, n_batches, R_rungs - 1)
  swap_accepted <- matrix(0L, n_batches, R_rungs - 1)

  # Per-rung running seed counters (so each rung gets distinct RNG streams)
  rung_seeds <- seed + (seq_len(R_rungs) - 1) * 100000

  # ---------------------------------------------------------------------
  # Main PT loop: run a batch at each rung, then propose swaps.
  # ---------------------------------------------------------------------
  for (b in seq_len(n_batches)) {

    # Decide whether to save samples from this batch (T=1 only, post burn-in).
    save_this_batch <- (b > burn_batches)

    # If we save, we need to align swap_every with thin to keep clean spacing.
    # Simplest: configure the sampler to run swap_every iters with n_burn=0
    # (since burn-in is handled at the outer level), and thin the saved
    # subset down to (swap_every / thin) samples per batch when on T=1.
    batch_thin <- if (save_this_batch) thin else swap_every  # avoid saving during burn

    # Run each rung for swap_every iterations.
    for (r in seq_len(R_rungs)) {
      # Force evaluation of init before the call (avoids any lazy-eval
      # interaction with the outer loop variable `r` and the inner Rcpp
      # bridge).
      init_r <- rung_state[[r]]
      seed_r <- rung_seeds[r] + b
      temp_r <- temperatures[r]
      store_ll_r <- save_this_batch && r == 1 && store_pointwise_ll

      rung_fit <- fit_dp_concordance(
        Y_g = Y_g, Y_cf = Y_cf,
        patient = patient, time = time,
        K = K,
        n_iter = swap_every,
        n_burn = 0,                  # burn handled at outer (batch) level
        model = model,
        seed = seed_r,
        hyperparams = hyperparams,
        thin = batch_thin,
        verbose = FALSE,
        kappa_fixed = kappa_fixed,
        store_pointwise_ll = store_ll_r,
        use_rcpp = use_rcpp,
        anneal = FALSE,
        fixed_temp = temp_r,
        init = init_r
      )
      rung_state[[r]] <- rung_fit$final_state

      # Stash T=1 samples for the official posterior
      if (save_this_batch && r == 1) {
        T1_samples_list[[length(T1_samples_list) + 1]] <- rung_fit$samples
      }
    }

    # Swap proposals between adjacent rungs.
    # We loop over R-1 adjacent pairs in random order to preserve detailed
    # balance.
    pair_order <- sample.int(R_rungs - 1)
    for (r in pair_order) {
      # Compute log-likelihood at each rung's current state.
      # We compute at temperature 1 (untempered) likelihood; the MH ratio
      # accounts for the tempering exponents.
      logL_r   <- pt_state_loglik(rung_state[[r]],   Y_g, Y_cf, patient, time,
                                  K, p, T_max, use_rcpp)
      logL_rp1 <- pt_state_loglik(rung_state[[r+1]], Y_g, Y_cf, patient, time,
                                  K, p, T_max, use_rcpp)

      # MH acceptance ratio for swap of states (r) and (r+1):
      #   alpha = min(1, exp((beta_r - beta_{r+1}) * (logL_{r+1} - logL_r)))
      # Note (beta_r - beta_{r+1}) > 0 since betas decrease with r.
      log_alpha <- (betas[r] - betas[r+1]) * (logL_rp1 - logL_r)
      swap_attempts[b, r] <- 1L
      if (log(runif(1)) < log_alpha) {
        # Swap states
        tmp <- rung_state[[r]]
        rung_state[[r]] <- rung_state[[r+1]]
        rung_state[[r+1]] <- tmp
        swap_accepted[b, r] <- 1L
      }
    }

    if (verbose && b %% max(1, n_batches %/% 20) == 0) {
      acc_rates <- colSums(swap_accepted[1:b, , drop = FALSE]) /
                   pmax(colSums(swap_attempts[1:b, , drop = FALSE]), 1)
      cat(sprintf("  batch %d/%d  swap acc rates: %s\n", b, n_batches,
                  paste(sprintf("%.2f", acc_rates), collapse = " ")))
    }
  }

  # ---------------------------------------------------------------------
  # Concatenate T=1 batch samples into a single fit object that looks like
  # fit_dp_concordance() output.
  # ---------------------------------------------------------------------
  if (length(T1_samples_list) == 0) {
    stop("No samples were saved at T=1; check n_iter > n_burn and thin settings.")
  }

  # Each batch produced (swap_every / thin) saved samples; concatenate.
  concat_field <- function(field) {
    pieces <- lapply(T1_samples_list, function(s) s[[field]])
    if (is.null(pieces[[1]])) return(NULL)
    if (is.matrix(pieces[[1]])) do.call(rbind, pieces)
    else                        do.call(c, pieces)
  }

  fields <- names(T1_samples_list[[1]])
  samples <- setNames(lapply(fields, concat_field), fields)

  # Final state = rung 1's last state
  final_state <- rung_state[[1]]

  # PT diagnostics
  pt_diagnostics <- list(
    temperatures = temperatures,
    betas        = betas,
    swap_every   = swap_every,
    n_batches    = n_batches,
    swap_attempts = swap_attempts,
    swap_accepted = swap_accepted,
    acceptance_rate = colSums(swap_accepted) / pmax(colSums(swap_attempts), 1)
  )

  list(
    samples = samples,
    n_save  = nrow(samples$omega_0),  # actual saved sample count
    model   = model,
    K       = K,
    n_iter  = n_iter,
    n_burn  = n_burn,
    thin    = thin,
    kappa_fixed = kappa_fixed,
    final_theta     = final_state$theta,
    final_theta_bar = final_state$theta_bar,
    final_theta_0   = final_state$theta_0,
    final_pi       = final_state$pi_vec,
    final_pi_g     = final_state$pi_g_vec,
    final_pi_cf    = final_state$pi_cf_vec,
    final_omega_g  = final_state$omega_g,
    final_omega_cf = final_state$omega_cf,
    final_state    = final_state,
    pt_diagnostics = pt_diagnostics,
    hyperparams    = list(),         # passed-through; downstream code expects this slot
    K_trace        = K,
    p = p, N_obs = N_obs, n = n, T_max = T_max
  )
}

