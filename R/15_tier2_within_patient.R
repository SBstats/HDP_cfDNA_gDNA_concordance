###############################################################################
## R/15_tier2_within_patient.R  —  Tier 2 DEFINITIVE within-patient longitudinal
## sub-clonal concordance (runs on existing fits; no re-fit required).
##
## Upgrades over the Tier-1 provisional (R/tier01_provisional analysis):
##   (1) POSTERIOR-UNCERTAINTY-PROPAGATED subclonal Delta co-movement: the
##       change-in-dominant-weight correlation is computed WITHIN each posterior
##       draw and summarised as a posterior mean + 95% credible interval, instead
##       of a single point estimate from posterior-mean weights.
##   (2) VALID M0 (independent-source) check: within each draw we take gDNA's
##       dominant signature and read cfDNA's loading on the SAME raw component
##       index (omega_g[k] and omega_cf[k] share theta_k in the stored traces),
##       so the comparison is not corrupted by independent per-source relabeling.
##       => co-movement under M0 cannot be a coupling-prior artifact.
##   (3) MODEL-FREE gene-level Delta test reproduced (the primary evidence).
##   (4) Between-chain THETA agreement audit from the per-chain theta_postmean.
##
## Outputs: results/tier2/*.csv
###############################################################################

set.seed(20260724)
suppressPackageStartupMessages({ library(dplyr) })
base_dir <- getwd()  # run from project root
if (interactive() || !nzchar(Sys.getenv("SLURM_JOB_ID"))) setwd(base_dir)
base_dir <- getwd()
res_dir <- file.path(base_dir, "results")
out_dir <- file.path(res_dir, "tier2"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

md   <- readRDS(file.path(res_dir, "model_data.rds"))
Yg   <- md$Y_g; Ycf <- md$Y_cf; N <- md$N_obs; tvec <- md$time; pvec <- md$patient

## consecutive within-patient transitions, strictly increasing time (drop duplicate-timepoint pairs)
build_transitions <- function() {
  do.call(rbind, lapply(unique(pvec), function(pt) {
    idx <- which(pvec == pt); idx <- idx[order(tvec[idx])]
    if (length(idx) < 2) return(NULL)
    keep <- do.call(rbind, lapply(1:(length(idx)-1), function(j) {
      a <- idx[j]; b <- idx[j+1]
      if (tvec[b] == tvec[a]) return(NULL)
      data.frame(patient = pt, a = a, b = b, stratum = paste0(tvec[a], "->", tvec[b]))
    }))
    keep
  }))
}
TR <- build_transitions(); M <- nrow(TR)
cat(sprintf("within-patient transitions: M=%d across %d patients\n", M, length(unique(TR$patient))))

## stratified pairing-shuffle of a length-M vector (permute within transition-type)
strat_shuffle <- function(x) {
  y <- x
  for (s in unique(TR$stratum)) { ii <- which(TR$stratum == s); if (length(ii) > 1) y[ii] <- sample(x[ii]) }
  y
}

###############################################################################
## (1)+(2) posterior-uncertainty-propagated subclonal Delta co-movement
###############################################################################
subclonal_delta <- function(fit_path, label) {
  fit <- readRDS(fit_path)
  G  <- do.call(rbind, lapply(fit$chains, function(c) c$samples$omega_g_trace))   # S x (K*N)
  CF <- do.call(rbind, lapply(fit$chains, function(c) c$samples$omega_cf_trace))
  K  <- fit$chains[[1]]$K_trace; S <- nrow(G)
  cor_obs <- numeric(S); cor_null <- numeric(S)
  for (s in 1:S) {
    Gm  <- matrix(G[s, ],  nrow = N, ncol = K)   # [obs x k], raw order (g & cf share theta_k)
    CFm <- matrix(CF[s, ], nrow = N, ncol = K)
    ref <- which.max(colMeans(Gm))               # gDNA's dominant signature at this draw
    wg  <- Gm[, ref]; wcf <- CFm[, ref]
    dg  <- wg[TR$b]  - wg[TR$a]
    dcf <- wcf[TR$b] - wcf[TR$a]
    cor_obs[s]  <- suppressWarnings(cor(dg, dcf))
    cor_null[s] <- suppressWarnings(cor(dg, strat_shuffle(dcf)))
  }
  obs_mean <- mean(cor_obs, na.rm = TRUE)
  ci <- quantile(cor_obs, c(.025, .975), na.rm = TRUE)
  p_perm <- (1 + sum(cor_null >= obs_mean, na.rm = TRUE)) / (1 + sum(is.finite(cor_null)))
  z <- (obs_mean - mean(cor_null, na.rm = TRUE)) / sd(cor_null, na.rm = TRUE)
  cat(sprintf("\n[%s]  posterior mean r(Δg,Δcf) = %.3f  (95%% CrI %.3f, %.3f)\n", label, obs_mean, ci[1], ci[2]))
  cat(sprintf("   pairing-shuffle null mean = %.3f | combined p = %.4f | z = %.2f\n",
              mean(cor_null, na.rm = TRUE), p_perm, z))
  cat(sprintf("   P(posterior draw co-movement > 0) = %.3f\n", mean(cor_obs > 0, na.rm = TRUE)))
  data.frame(model = label, post_mean_r = round(obs_mean,4),
             cri_lo = round(ci[1],4), cri_hi = round(ci[2],4),
             null_mean = round(mean(cor_null, na.rm=TRUE),4),
             p_combined = round(p_perm,4), z = round(z,2),
             p_draw_gt0 = round(mean(cor_obs > 0, na.rm=TRUE),3))
}
cat("\n===== (1)+(2) subclonal Delta co-movement with posterior uncertainty =====\n")
sub_M1 <- subclonal_delta(file.path(res_dir, "krd_M1_fit.rds"), "M1 tracking (uncertainty-propagated)")
sub_M0 <- subclonal_delta(file.path(res_dir, "krd_M0_fit.rds"), "M0 independent (valid cross-source ref)")
write.csv(rbind(sub_M1, sub_M0), file.path(out_dir, "tier2_subclonal_delta.csv"), row.names = FALSE)

###############################################################################
## (3) model-free gene-level Delta test (primary evidence, reproduced)
###############################################################################
cat("\n===== (3) model-free gene-level Delta co-movement (primary) =====\n")
dG  <- sapply(1:M, function(m) Yg[, TR$b[m]]  - Yg[, TR$a[m]])
dCF <- sapply(1:M, function(m) Ycf[, TR$b[m]] - Ycf[, TR$a[m]])
C   <- cor(dG, dCF); obs <- mean(diag(C))
n_perm <- 5000
strat_null <- replicate(n_perm, { q <- 1:M; for (s in unique(TR$stratum)) { ii <- which(TR$stratum==s); if (length(ii)>1) q[ii] <- sample(ii) }; mean(C[cbind(1:M, q)]) })
p_s <- (1 + sum(strat_null >= obs)) / (1 + n_perm)
## Save the permutation null draws + observed value so the applied paper figS3
## generator (R/26) can draw panel (a)'s histogram from the exact same run.
saveRDS(list(observed = obs, null_draws = strat_null, n_perm = n_perm,
             level_corr = suppressWarnings(mean(diag(cor(Yg[, TR$a], Ycf[, TR$a]))))),
        file.path(base_dir, "results", "tier2", "genelevel_increment_null.rds"))
cat(sprintf("mean Δ-corr (true) = %.3f | stratified null = %.3f | p = %.4f | z = %.2f | frac>0 = %.2f\n",
            obs, mean(strat_null), p_s, (obs-mean(strat_null))/sd(strat_null), mean(diag(C) > 0)))
write.csv(data.frame(metric="gene-level Delta (model-free)", M=M, delta_true=round(obs,4),
                     null=round(mean(strat_null),4), p=round(p_s,4),
                     z=round((obs-mean(strat_null))/sd(strat_null),2), frac_pos=round(mean(diag(C)>0),3)),
          file.path(out_dir, "tier2_genelevel_delta.csv"), row.names = FALSE)

###############################################################################
## (4) between-chain theta agreement audit (from per-chain theta_postmean)
###############################################################################
cat("\n===== (4) between-chain theta agreement (occupied components) =====\n")
## Audited on the POPULATION signature theta_bar [p x K], not the subject-level
## theta [p x K x n].
##
## Two reasons. First, theta_bar is the quantity the manuscript interprets: the
## subject-level signatures are informed by only ~T_i observations each and are
## documented as not consistently estimable (writeup Remark rem:not_claimed), so
## a per-subject between-chain correlation would mostly measure prior noise.
## Second, an earlier version of this function read dim(theta_postmean)[3] as a
## TIMEPOINT index; that dimension is now the SUBJECT, so the old loop both
## mislabelled its output and produced n x Kocc rows instead of T x Kocc.
theta_audit <- function(fit_path, label, Kocc = 2) {
  fit <- readRDS(fit_path)
  Th <- lapply(fit$chains, function(c) c$theta_bar_postmean)   # each [p x K]
  if (any(vapply(Th, is.null, logical(1))))
    stop("theta_bar_postmean missing; refit with the current sampler.")
  nc <- length(Th)
  Kocc <- min(Kocc, ncol(Th[[1]]))
  rows <- list()
  for (k in 1:Kocc) {
    Mk <- sapply(Th, function(a) a[, k])                       # p x nchains
    prs <- combn(nc, 2)
    mincor <- min(sapply(1:ncol(prs), function(j) cor(Mk[, prs[1,j]], Mk[, prs[2,j]])))
    rows[[length(rows)+1]] <- data.frame(model=label, comp=k,
      min_pairwise_cor = round(mincor,4),
      betw_chain_sd = round(mean(apply(Mk,1,sd)),4),
      signal_sd = round(sd(rowMeans(Mk)),4))
  }
  do.call(rbind, rows)
}
th_M1 <- theta_audit(file.path(res_dir, "krd_M1_fit.rds"), "M1")
th_M0 <- theta_audit(file.path(res_dir, "krd_M0_fit.rds"), "M0")
th <- rbind(th_M1, th_M0)
th$ratio_sd <- round(th$betw_chain_sd / th$signal_sd, 3)
print(th, row.names = FALSE)
cat(sprintf("\nθ occupied-component min pairwise between-chain corr: M1 = %.3f, M0 = %.3f\n",
            min(th_M1$min_pairwise_cor), min(th_M0$min_pairwise_cor)))
write.csv(th, file.path(out_dir, "tier2_theta_between_chain.csv"), row.names = FALSE)

cat("\nDONE. Outputs in results/tier2/\n")
