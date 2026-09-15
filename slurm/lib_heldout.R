###############################################################################
## slurm/lib_heldout.R -- shared held-out predictive scoring for K-fold and LOPO
##
## Previously `heldout_lpd`, `trained_params`, `score` and the convergence gate
## were duplicated byte-for-byte in run_kfold_aggregate.R and
## run_lopo_aggregate.R. They are factored out here so the two cannot drift.
##
## SCORING UNDER THE CURRENT MODEL (writeup Sec 7, eq:lopo_omega--eq:lopo_w0)
##
## Every subject-indexed parameter is unavailable for a patient who was never
## in the training set: the held-out patient has no theta[,,i], no theta_0[,i]
## and no omega_0[i,t]. The earlier implementation plugged in cohort means,
## which understates predictive uncertainty and implicitly scores a different
## model from the one fitted.
##
## We instead Monte-Carlo marginalise over the subject-level prior. For each of
## M draws we first sample a POSTERIOR index d from the training fit, then:
##
##   (sigma_k^2, sigma_0^2, varsigma_k^2, varsigma_0^2)  <- posterior draw d
##   omega^s_{i*t}  ~ Dir(kappa * pi_hat)            [M1] or Dir(kappa * pi_hat^s) [M0]
##   theta_{i*jk}   ~ N(theta_bar_hat[j,k], varsigma_k^2[d])
##   theta_0_{i*j}  ~ N(theta_0_bar_hat[j], varsigma_0^2[d])
##   omega_0_{i*t}  ~ Beta(a_0, b_0)
##
## then evaluate the observation model and average the densities on the log
## scale.
##
## WHAT IS AND IS NOT INTEGRATED OVER. The subject-level parameters and the
## population VARIANCE components are drawn, so both between-subject
## variability and posterior uncertainty in the dispersions enter the
## predictive. The population MEANS (pi_hat, theta_bar_hat, theta_0_bar_hat)
## remain plug-ins, because the sampler stores theta_bar only as a per-chain
## posterior mean and not as a full trace; enabling store_theta_trace and
## drawing from theta_bar_trace would remove this last approximation. The
## quantity computed is therefore p(y* | data) marginalised over the
## subject level and the variance components, and conditional on the
## population means -- stated explicitly so the reported SEs are not
## over-interpreted.
##
## This is intrinsically noisier than a full plug-in, which is appropriate
## rather than a defect: the earlier scoring reported a precision the design
## does not support.
###############################################################################

logsumexp <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }

#' Held-out log predictive density for one observation, one source.
#'
#' @param Y_obs   numeric [p]: the held-out observation for this source.
#' @param P       trained-parameter list from `trained_params()`.
#' @param source  "g" or "cf". cfDNA adds the background component.
#' @param M       Monte Carlo draws.
#' @param p       number of features (passed explicitly; do NOT rely on scope).
#' @return list(lpd, mcse) -- the log predictive density and its Monte Carlo SE.
heldout_lpd <- function(Y_obs, P, source = c("g", "cf"), M = 400, p = length(Y_obs)) {
  source <- match.arg(source)
  pi_vec <- if (source == "g") P$pig else P$picf
  K <- length(pi_vec)
  a <- P$kappa * pi_vec; a[a < 1e-6] <- 1e-6

  ## Population-level variance components: draw a posterior index per MC
  ## replicate so that population uncertainty is integrated over rather than
  ## conditioned away. Falls back to the point estimates if no draws were
  ## supplied (e.g. an older fit object).
  has_draws <- !is.null(P$draws) && nrow(as.matrix(P$draws$varsigma_k2)) > 0
  n_draw <- if (has_draws) nrow(as.matrix(P$draws$varsigma_k2)) else 0L

  lps <- numeric(M)
  for (m in seq_len(M)) {
    if (has_draws) {
      d      <- sample.int(n_draw, 1L)
      sd_th  <- sqrt(pmax(P$draws$varsigma_k2[d, ], 1e-12))
      sd_t0  <- sqrt(max(P$draws$varsigma_0_2[d], 1e-12))
      sig_k  <- sqrt(pmax(P$draws$sigma_k2[d, ], 1e-12))
      sig_0  <- sqrt(max(P$draws$sigma_0_2[d], 1e-12))
    } else {
      sd_th <- sqrt(pmax(P$varsigma_k2, 1e-12))
      sd_t0 <- sqrt(max(P$varsigma_0_2, 1e-12))
      sig_k <- P$sigk
      sig_0 <- P$sigma_0
    }

    ## (1) subject weights from the trained prior
    w <- rgamma(K, a, 1); w <- w / sum(w)

    ## (2) subject-level signatures from the trained hierarchical prior.
    ## This is the key change: theta is NOT plugged in at a cohort mean, it is
    ## drawn, so between-subject variability enters the predictive.
    th <- P$theta_bar + matrix(rnorm(p * K, 0, rep(sd_th, each = p)), p, K)

    logcomp <- matrix(0, p, K)
    for (k in seq_len(K)) {
      logcomp[, k] <- log(w[k] + 1e-300) +
        dnorm(Y_obs, th[, k], sig_k[k], log = TRUE)
    }
    lg <- apply(logcomp, 1, logsumexp)

    if (source == "cf") {
      ## (3) background mean and (4) contamination fraction, both drawn
      t0 <- P$theta_0_bar + rnorm(p, 0, sd_t0)
      w0 <- rbeta(1, P$a_0, P$b_0)
      lt <- lg + log(1 - w0)
      ln <- log(w0 + 1e-300) + dnorm(Y_obs, t0, sig_0, log = TRUE)
      mx <- pmax(lt, ln)
      lg <- mx + log(exp(lt - mx) + exp(ln - mx))
    }
    lps[m] <- sum(lg)
  }

  lpd <- logsumexp(lps) - log(M)

  ## Monte Carlo diagnostics.
  ##
  ## `mcse` is the delta-method SE of the log-mean-exp. The functional form is
  ## valid (lps is i.i.d. across m, since both the posterior index and the
  ## subject-level draws are resampled inside the loop), BUT it saturates in
  ## the regime this marginalisation creates: lps sums over p genes, so its
  ## across-replicate spread is O(sqrt(p)) log-units, and once sd(lps) is more
  ## than a few units a single replicate dominates exp(lps - max). Then
  ## sd(z)/mean(z) -> sqrt(M) and the reported MCSE approaches 1 regardless of
  ## M -- i.e. it stops being informative exactly when it matters.
  ##
  ## We therefore also return an importance-weight ESS and the observed
  ## sd(lps). A small `ess` (say < 0.1*M) means the estimate rests on a handful
  ## of draws and `mcse` is a LOWER BOUND on the true uncertainty; increase M
  ## or reduce the dimensionality being marginalised.
  z    <- exp(lps - max(lps))
  mcse <- if (M > 1) sd(z) / (sqrt(M) * mean(z)) else NA_real_
  wn   <- z / sum(z)
  ess  <- 1 / sum(wn^2)
  list(lpd = lpd, mcse = mcse, ess = ess, sd_lps = if (M > 1) sd(lps) else NA_real_)
}

#' Extract trained population-level parameters from a multi-chain fit.
#'
#' Applies the standard relabel -> cross-chain align -> merge sequence, then
#' pools the population-level posterior means across chains. Only
#' POPULATION-level quantities are extracted: subject-level theta / theta_0 /
#' omega_0 are deliberately NOT used, because a held-out patient has none.
trained_params <- function(fit, model) {
  for (i in seq_along(fit$chains)) fit$chains[[i]] <- relabel_by_weight(fit$chains[[i]])
  fit$chains <- align_chains_by_signature(fit$chains)
  m   <- merge_chain_samples(fit$chains)
  ch1 <- fit$chains[[1]]
  nc  <- length(fit$chains)

  if (is.null(ch1$theta_bar_postmean)) {
    stop("Fit has no theta_bar_postmean. This aggregator requires fits from the ",
         "current sampler (shared subject-indexed theta). Refit required.")
  }
  theta_bar   <- Reduce(`+`, lapply(fit$chains, function(c) c$theta_bar_postmean)) / nc
  theta_0_bar <- Reduce(`+`, lapply(fit$chains, function(c) rowMeans(c$theta_0_postmean))) / nc

  ## Per-draw population parameters, so the predictive can integrate over
  ## POPULATION uncertainty rather than conditioning on point estimates. With
  ## n ~ 25 training patients the posterior spread of varsigma_k^2 is not
  ## negligible, and plugging in the posterior mean of a VARIANCE is
  ## particularly bad (Jensen: it understates predictive spread). theta_bar is
  ## only available as a per-chain posterior mean (no full trace is stored),
  ## so it remains a plug-in; the variance components are drawn.
  draws <- list(
    sigma_k2    = m$sigma_k^2,      # [n_save x K]
    sigma_0_2   = m$sigma_0^2,      # [n_save]
    varsigma_k2 = m$varsigma_k^2,   # [n_save x K]
    varsigma_0_2= m$varsigma_0^2    # [n_save]
  )

  if (model == "M0") {
    pig  <- stick_break(c(colMeans(m$v_g),  1))
    picf <- stick_break(c(colMeans(m$v_cf), 1))
  } else {
    pib  <- stick_break(c(colMeans(m$v), 1)); pig <- pib; picf <- pib
  }
  stopifnot(all(is.finite(pig)), all(is.finite(picf)))

  hp <- ch1$hyperparams
  out <- list(
    theta_bar    = theta_bar,
    theta_0_bar  = theta_0_bar,
    ## NOTE ON SQUARING. samples$sigma_k / varsigma_k / varsigma_0 store
    ## STANDARD DEVIATIONS (the sampler writes sqrt(.) into them). The
    ## posterior mean of a variance is E[s^2], NOT (E[s])^2 -- by Jensen the
    ## latter is strictly smaller, by ~20% on an IG(2,2)-scale posterior. Using
    ## (E[s])^2 here would systematically narrow the held-out predictive, which
    ## is exactly the defect the prior-marginalisation redesign exists to avoid.
    ## So: average the SQUARES.
    sigk         = sqrt(colMeans(m$sigma_k^2)),
    sigma_0      = sqrt(mean(m$sigma_0^2)),
    varsigma_k2  = colMeans(m$varsigma_k^2),
    varsigma_0_2 = mean(m$varsigma_0^2),
    a_0          = if (!is.null(hp$a_0)) hp$a_0 else 1,
    b_0          = if (!is.null(hp$b_0)) hp$b_0 else 9,
    pig = pig, picf = picf,
    kappa = ch1$kappa_fixed,
    draws = draws
  )
  stopifnot(all(is.finite(out$sigk)), is.finite(out$sigma_0),
            all(is.finite(out$varsigma_k2)), is.finite(out$varsigma_0_2))
  out
}

#' Score one held-out observation under both sources.
score_obs <- function(P, Y_g, Y_cf, o, M = 400, p = nrow(Y_g)) {
  rg <- heldout_lpd(Y_g[, o],  P, "g",  M = M, p = p)
  rc <- heldout_lpd(Y_cf[, o], P, "cf", M = M, p = p)
  stopifnot(is.finite(rg$lpd), is.finite(rc$lpd))
  ## ess_* below 0.1*M indicates a degenerate importance-weight distribution,
  ## in which case mcse_* understates the true Monte Carlo error.
  c(lpd_g = rg$lpd, lpd_cf = rc$lpd, lpd = rg$lpd + rc$lpd,
    mcse_g = rg$mcse, mcse_cf = rc$mcse,
    ess_g = rg$ess, ess_cf = rc$ess,
    sdlps_g = rg$sd_lps, sdlps_cf = rc$sd_lps)
}

#' Pre-specified per-fold convergence gate.
#'
#' The within-fit chain-health check cannot see a fold in which ALL chains agree
#' at a bad mode, so this compares folds against each other. A fold is flagged
#' if EITHER (a) its excess training log-likelihood exceeds the spread of the
#' other folds (requires >= 2 others), or (b) its fraction-of-best falls below
#' 0.7. The OR rule was validated against six scenarios including the real
#' fold-3 pathology; a MAD z-score is degenerate at small F and a
#' ratio-to-best inverts when log-likelihoods are positive.
fold_gate <- function(ll_by_fold, frac_thresh = 0.7) {
  F_ <- length(ll_by_fold)
  flagged <- logical(F_)
  best <- max(ll_by_fold, na.rm = TRUE)
  for (f in seq_len(F_)) {
    others <- ll_by_fold[-f]
    others <- others[is.finite(others)]
    cond_spread <- FALSE
    if (length(others) >= 2) {
      spread <- diff(range(others))
      excess <- abs(ll_by_fold[f] - median(others))
      cond_spread <- is.finite(excess) && spread > 0 && (excess / spread) > 1
    }
    frac <- if (is.finite(ll_by_fold[f]) && best != 0) ll_by_fold[f] / best else NA_real_
    cond_frac <- is.finite(frac) && frac < frac_thresh
    flagged[f] <- cond_spread || cond_frac
  }
  flagged
}
