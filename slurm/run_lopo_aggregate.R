#!/usr/bin/env Rscript
###############################################################################
## slurm/run_lopo_aggregate.R -- Tier 4: held-out scoring + M1 vs M0 comparison
## for the leave-one-PATIENT-out fits produced by run_lopo_cv.R.
##
## Scoring machinery now lives in slurm/lib_heldout.R, shared with
## run_kfold_aggregate.R (previously duplicated byte-for-byte in both files).
##
## Held-out log predictive density marginalises over the ENTIRE subject-level
## prior -- weights, signatures, background mean and contamination fraction --
## because a held-out patient has none of those parameters. See lib_heldout.R
## and writeup Sec 7 (eq:lopo_omega--eq:lopo_w0). This is noisier than the
## earlier plug-in scoring by design; Monte Carlo SEs are reported.
##
## Tier-4 additions:
##   - patient-grouped folds (fold_assignment.rds written by run_lopo_cv.R);
##   - SOURCE-SPECIFIC delta elpd: delta_g and delta_cf isolate whether SHARING
##     the composition across sources (M1) improves held-out prediction of each
##     source. delta_cf > 0 is the pairing-specific statement that knowing gDNA
##     improves prediction of the paired cfDNA (and vice versa for delta_g);
##   - per-patient delta elpd (results/lopo_per_patient.csv) for a forest plot.
##
## Output: results/lopo_cv.csv, results/lopo_pointwise.csv, results/lopo_per_patient.csv
###############################################################################

base_dir <- getwd()
source(file.path(base_dir, "R", "01_lib_core.R"))
source(file.path(base_dir, "slurm", "lib_heldout.R"))
set.seed(20260724)

res_dir <- Sys.getenv("LOPO_RES_DIR", file.path(base_dir, "results"))  # override for local testing
md <- readRDS(file.path(res_dir, "model_data.rds"))
Y_g <- md$Y_g; Y_cf <- md$Y_cf; N <- md$N_obs; p <- md$p
patient <- md$patient; tvec <- md$time
out_dir <- file.path(res_dir, "lopo")
MC <- as.integer(Sys.getenv("MC_DRAWS", "400"))

fa <- readRDS(file.path(out_dir, "fold_assignment.rds"))
fold_id <- fa$fold_id; F_folds <- fa$F_folds

## Per-fold convergence gate: pre-specified within-fit chain-health check.
fold_trapped <- function(f, model) {
  fp <- file.path(out_dir, sprintf("fit_f%d_%s.rds", f, model))
  if (!file.exists(fp)) return(NA)
  length(identify_trapped_chains(readRDS(fp), ll_gap = 10000)$drop) > 0
}
gate <- data.frame(fold = 1:F_folds,
                   trapped_M1 = sapply(1:F_folds, fold_trapped, model = "M1"),
                   trapped_M0 = sapply(1:F_folds, fold_trapped, model = "M0"))
gate$flagged <- (!is.na(gate$trapped_M1) & gate$trapped_M1) |
                (!is.na(gate$trapped_M0) & gate$trapped_M0)
cat("\n=== Per-fold convergence gate (patient-grouped LOPO) ===\n")
print(gate, row.names = FALSE)
if (any(gate$flagged)) cat("  EXCLUDED folds (chains trapped):",
                           paste(gate$fold[gate$flagged], collapse = ", "), "\n")
write.csv(gate, file.path(res_dir, "lopo_convergence_gate.csv"), row.names = FALSE)

pw <- list()
for (f in 1:F_folds) {
  if (isTRUE(gate$flagged[gate$fold == f])) { cat("fold", f, "flagged; skip\n"); next }
  te <- which(fold_id == f)
  f1 <- file.path(out_dir, sprintf("fit_f%d_M1.rds", f))
  f0 <- file.path(out_dir, sprintf("fit_f%d_M0.rds", f))
  if (!file.exists(f1) || !file.exists(f0)) { cat("fold", f, "missing fits; skip\n"); next }
  P1 <- trained_params(readRDS(f1), "M1"); P0 <- trained_params(readRDS(f0), "M0")
  for (o in te) {
    s1 <- score_obs(P1, Y_g, Y_cf, o, M = MC, p = p)
    s0 <- score_obs(P0, Y_g, Y_cf, o, M = MC, p = p)
    pw[[length(pw)+1]] <- data.frame(obs=o, patient=patient[o], fold=f, model="M1",
                                     lpd_g=s1["lpd_g"], lpd_cf=s1["lpd_cf"], lpd=s1["lpd"],
                                     mcse_g=s1["mcse_g"], mcse_cf=s1["mcse_cf"])
    pw[[length(pw)+1]] <- data.frame(obs=o, patient=patient[o], fold=f, model="M0",
                                     lpd_g=s0["lpd_g"], lpd_cf=s0["lpd_cf"], lpd=s0["lpd"],
                                     mcse_g=s0["mcse_g"], mcse_cf=s0["mcse_cf"])
  }
  cat("fold", f, "scored (", length(te), "held-out obs)\n")
}
pw_all <- do.call(rbind, pw)
write.csv(pw_all, file.path(res_dir, "lopo_pointwise.csv"), row.names = FALSE)

## paired-by-observation deltas (M1 - M0), total and per source
w1 <- pw_all[pw_all$model=="M1", ]; w0 <- pw_all[pw_all$model=="M0", ]
w0 <- w0[match(w1$obs, w0$obs), ]
paired_se <- function(d) sqrt(length(d)) * sd(d)
d   <- w1$lpd    - w0$lpd
dg  <- w1$lpd_g  - w0$lpd_g
dcf <- w1$lpd_cf - w0$lpd_cf

summary_df <- data.frame(
  elpd_M1 = sum(w1$lpd), elpd_M0 = sum(w0$lpd),
  delta_elpd = sum(d), se_delta = paired_se(d),
  delta_elpd_g = sum(dg), se_delta_g = paired_se(dg),        # sharing helps predict gDNA
  delta_elpd_cf = sum(dcf), se_delta_cf = paired_se(dcf),    # ... and paired cfDNA (pairing-specific)
  n_obs = length(d), n_obs_favouring_M1 = sum(d > 0),
  median_delta_per_obs = median(d),
  folds_used = sum(!gate$flagged, na.rm = TRUE),
  folds_flagged = paste(gate$fold[gate$flagged], collapse = ";"),
  favours = ifelse(sum(d) > 0, "M1", "M0"))
stopifnot("held-out lpd contains NA" = !anyNA(d))
if (sign(sum(d)) != sign(median(d)))
  cat("\n*** WARNING: total and median delta_elpd disagree in sign (outlier-driven) ***\n")
write.csv(summary_df, file.path(res_dir, "lopo_cv.csv"), row.names = FALSE)

## per-patient delta (for a forest plot / heterogeneity)
pp <- aggregate(cbind(d = d, dg = dg, dcf = dcf), by = list(patient = w1$patient), sum)
write.csv(pp, file.path(res_dir, "lopo_per_patient.csv"), row.names = FALSE)

cat("\n==== Leave-one-patient-out CV (Tier 4) ====\n"); print(summary_df)
cat(sprintf("\nPairing-specific readout: delta_elpd_cf = %.1f (SE %.1f) => sharing gDNA %s held-out cfDNA prediction.\n",
            summary_df$delta_elpd_cf, summary_df$se_delta_cf,
            ifelse(summary_df$delta_elpd_cf > 0, "IMPROVES", "does NOT improve")))
