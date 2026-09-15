#!/usr/bin/env Rscript
###############################################################################
## slurm/run_kfold_aggregate.R  --  held-out scoring + M1 vs M0 comparison
##
## Run AFTER all 2*F fold fits from run_kfold_cv.R exist in results/kfold/.
## Short-to-moderate job (held-out lpd is MC over the weight prior).
##
## For each held-out observation (i,t) in fold f and each source, computes the
## cross-validated log predictive density, marginalising the unseen source
## weights over their trained PRIOR:
##    omega^s      ~ Dir(kappa * pi)   [M1]  or Dir(kappa*pi^s)  [M0]
##    theta_{i*jk} ~ N(theta_bar[j,k], varsigma_k^2)
##    theta_0_{i*j}~ N(theta_0_bar[j], varsigma_0^2)
##    omega_0_{i*t}~ Beta(a_0, b_0)
## via Monte Carlo (MC_DRAWS). Every subject-indexed parameter is marginalised
## because a held-out observation's subject may be absent from training.
##
## elpd_kfold(model) = sum of held-out lpd; delta = M1 - M0 with paired SE.
## Output: results/kfold_cv.csv, results/kfold_pointwise.csv
##
## THETA: uses `theta_bar_postmean` (population level) plus `varsigma_k`, in
## canonical component order. Only POPULATION-level quantities are used; the
## subject-level theta is drawn from its prior rather than plugged in. Fits
## lacking `theta_bar_postmean` are rejected with an error rather than silently
## falling back to a stale field.
##
## Scoring machinery is shared with run_lopo_aggregate.R via slurm/lib_heldout.R
## (previously duplicated byte-for-byte in both files).
##
## A PER-FOLD CONVERGENCE GATE excludes folds whose trained log-likelihood is a
## strong low outlier; the within-fit chain-health check cannot see these
## because all chains agree at the same bad mode.
###############################################################################

base_dir <- getwd()
source(file.path(base_dir, "R", "01_lib_core.R"))
source(file.path(base_dir, "slurm", "lib_heldout.R"))
set.seed(20260718)

md <- readRDS(file.path(base_dir, "results", "model_data.rds"))
Y_g <- md$Y_g; Y_cf <- md$Y_cf; N <- md$N_obs; p <- md$p
patient <- md$patient; tvec <- md$time
out_dir <- file.path(base_dir, "results", "kfold")
MC <- as.integer(Sys.getenv("MC_DRAWS", "400"))

fa <- readRDS(file.path(out_dir, "fold_assignment.rds"))
fold_id <- fa$fold_id; F_folds <- fa$F_folds

## ---------------------------------------------------------------------------
## PER-FOLD CONVERGENCE GATE
## The existing chain-health check compares chains WITHIN a fit, so it cannot
## detect a fold where all chains agree at the same bad mode (exactly what
## happened in fold 3 on 2026-07-18: mean loglik ~32k vs 46k-125k elsewhere,
## consistent across all 4 chains). Here we compare each fold's trained mean
## log-likelihood, per gene per training observation, ACROSS folds and flag
## outliers. Flagged folds are excluded from the elpd comparison and reported.
## ---------------------------------------------------------------------------
## PRIMARY criterion: the project's PRE-SPECIFIED within-fit chain-health check
## (identify_trapped_chains, ll_gap = 10,000), applied per fold. A fold is
## pathological when its own chains disagree — i.e. the sampler failed to
## converge on that training subset — which is scale-free and does not penalise
## a fold merely for being intrinsically harder than the others.
##
## Rationale (2026-07-19): a purely cross-fold rule flagged fold 5 even though
## all its chains were tightly converged (spread 1,319 / 644, none trapped); it
## was only "low" because folds 1-2 happened to be high. Excluding it would
## discard valid data. Fold 3 by contrast had M0 chains at
## 32,329 / -27,207 / -6,557 / 32,401 (spread 59,608, two chains trapped).
fold_trapped <- function(f, model) {
  fp <- file.path(out_dir, sprintf("fit_f%d_%s.rds", f, model))
  if (!file.exists(fp)) return(NA)
  tr <- identify_trapped_chains(readRDS(fp), ll_gap = 10000)
  length(tr$drop) > 0
}
fold_ll_norm <- function(f, model) {
  fp <- file.path(out_dir, sprintf("fit_f%d_%s.rds", f, model))
  if (!file.exists(fp)) return(NA_real_)
  fit <- readRDS(fp)
  n_train <- sum(fold_id != f)
  mean(sapply(fit$chains, function(c) mean(c$samples$loglik))) / (p * n_train)
}
gate <- data.frame(fold = 1:F_folds,
                   ll_M1 = sapply(1:F_folds, fold_ll_norm, model = "M1"),
                   ll_M0 = sapply(1:F_folds, fold_ll_norm, model = "M0"))
gate$worst <- pmin(gate$ll_M1, gate$ll_M0)
# Flag folds whose normalised loglik falls far BELOW the best fold, measured as
# a fraction of the spread across folds. Two earlier attempts failed in local
# testing (2026-07-18): a MAD z-score is degenerate at small F (at F=2 every
# fold scores +/-0.6745), and a ratio-to-best inverts when logliks are POSITIVE
# (as they are on the real data, where the pathological fold 3 scored ratio
# 0.25, not >2). This gap-based rule is sign-agnostic.
gate_frac <- as.numeric(Sys.getenv("KFOLD_GATE_FRAC", "0.7"))
best_ll <- max(gate$worst, na.rm = TRUE)          # healthiest fold
gate$gap_vs_best <- best_ll - gate$worst          # >= 0, larger = worse
# Leave-one-out scale: for each fold, how large is its shortfall relative to
# the spread of the OTHER folds? Using the other folds' own range avoids both
# the small-F degeneracy of a MAD z-score and the sign-inversion of a ratio.
gate$excess <- sapply(seq_len(nrow(gate)), function(i) {
  others <- gate$worst[-i]
  others <- others[is.finite(others)]
  # Need >=2 other folds to estimate a spread; with fewer, this signal is
  # undefined (at F=2 it would otherwise return Inf for the lower fold and
  # false-positive on perfectly healthy fits). Fall back to frac_of_best.
  if (length(others) < 2) return(0)
  spread <- diff(range(others))
  # shortfall below the WORST of the other folds, scaled by their spread
  short <- min(others) - gate$worst[i]
  if (!is.finite(short) || short <= 0) return(0)
  if (spread > 0) short / spread else 0
})
# EXCLUSION DECISION: driven by the pre-specified within-fit chain-health check
# (chains disagreeing => sampler failed on that training subset). The cross-fold
# statistics below are retained as REPORTED DIAGNOSTICS only; they are useful
# context but are not used to exclude, because folds legitimately differ in
# training-set difficulty and a cross-fold rule false-positives on healthy folds.
gate$frac_of_best <- gate$gap_vs_best / pmax(abs(best_ll), .Machine$double.eps)
gate$trapped_M1 <- sapply(1:F_folds, fold_trapped, model = "M1")
gate$trapped_M0 <- sapply(1:F_folds, fold_trapped, model = "M0")
gate$flagged <- (!is.na(gate$trapped_M1) & gate$trapped_M1) |
                (!is.na(gate$trapped_M0) & gate$trapped_M0)
# Advisory only: a fold that is a strong cross-fold outlier but whose chains all
# converged is reported for inspection, NOT excluded.
gate$advisory_outlier <- !is.na(gate$excess) &
                         (gate$excess > 1 | gate$frac_of_best > gate_frac) &
                         !gate$flagged
cat("\n=== Per-fold convergence gate (loglik per gene per training obs) ===\n")
print(gate, digits = 4, row.names = FALSE)
if (any(gate$flagged)) {
  cat("  EXCLUDED folds (chains trapped; pre-specified check):",
      paste(gate$fold[gate$flagged], collapse = ", "), "\n")
} else cat("  No folds excluded (all chains converged).\n")
if (any(gate$advisory_outlier)) {
  cat("  ADVISORY (cross-fold outlier but chains converged; RETAINED):",
      paste(gate$fold[gate$advisory_outlier], collapse = ", "), "\n")
}
write.csv(gate, file.path(base_dir, "results", "kfold_convergence_gate.csv"), row.names = FALSE)

pw <- list()
for (f in 1:F_folds) {
  if (isTRUE(gate$flagged[gate$fold == f])) {
    cat("fold", f, "FLAGGED by convergence gate; skipping\n"); next
  }
  te <- which(fold_id == f)
  f1 <- file.path(out_dir, sprintf("fit_f%d_M1.rds", f))
  f0 <- file.path(out_dir, sprintf("fit_f%d_M0.rds", f))
  if (!file.exists(f1) || !file.exists(f0)) { cat("fold", f, "missing fits; skip\n"); next }
  P1 <- trained_params(readRDS(f1), "M1"); P0 <- trained_params(readRDS(f0), "M0")
  for (o in te) {
    s1 <- score_obs(P1, Y_g, Y_cf, o, M = MC, p = p)
    s0 <- score_obs(P0, Y_g, Y_cf, o, M = MC, p = p)
    pw[[length(pw)+1]] <- data.frame(obs = o, fold = f, model = "M1",
                                     lpd_g = s1["lpd_g"], lpd_cf = s1["lpd_cf"], lpd = s1["lpd"],
                                     mcse_g = s1["mcse_g"], mcse_cf = s1["mcse_cf"])
    pw[[length(pw)+1]] <- data.frame(obs = o, fold = f, model = "M0",
                                     lpd_g = s0["lpd_g"], lpd_cf = s0["lpd_cf"], lpd = s0["lpd"],
                                     mcse_g = s0["mcse_g"], mcse_cf = s0["mcse_cf"])
  }
  cat("fold", f, "scored (", length(te), "held-out obs)\n")
}
pw_all <- do.call(rbind, pw)
write.csv(pw_all, file.path(base_dir, "results", "kfold_pointwise.csv"), row.names = FALSE)

m1 <- pw_all$lpd[pw_all$model == "M1"]; o1 <- pw_all$obs[pw_all$model == "M1"]
m0 <- pw_all$lpd[pw_all$model == "M0"]; o0 <- pw_all$obs[pw_all$model == "M0"]
d <- m1 - m0[match(o1, o0)]                        # paired by obs
summary_df <- data.frame(elpd_M1 = sum(m1), elpd_M0 = sum(m0),
                         delta_elpd = sum(d), se_delta = sqrt(length(d)) * sd(d),
                         n_obs = length(d),
                         median_delta_per_obs = median(d),
                         n_obs_favouring_M1 = sum(d > 0),
                         folds_used = sum(!gate$flagged, na.rm = TRUE),
                         folds_flagged = paste(gate$fold[gate$flagged], collapse = ";"),
                         favours = ifelse(sum(d) > 0, "M1", "M0"))
# Guard: a mean/median disagreement signals outlier-driven results (the failure
# mode of the invalid 2026-07-18 run). Surface it rather than reporting a
# headline number that a handful of observations produced.
stopifnot("held-out lpd contains NA - inspect kfold_pointwise.csv" = !anyNA(d))
if (!is.na(sum(d)) && !is.na(median(d)) && sign(sum(d)) != sign(median(d)))
  cat("\n*** WARNING: total and median delta_elpd disagree in sign ***\n",
      "    Result is outlier-driven; inspect kfold_pointwise.csv per fold\n",
      "    BEFORE using this number.\n", sep = "")
write.csv(summary_df, file.path(base_dir, "results", "kfold_cv.csv"), row.names = FALSE)
cat("\n==== K-fold CV ====\n"); print(summary_df)
