#!/usr/bin/env Rscript
###############################################################################
## slurm/run_krd_combine.R
## Combine M1 and M0 fits, compute diagnostics and model comparison.
## Run AFTER both 04a_krd_M1.sbatch and 04b_krd_M0.sbatch have completed.
##
## Reads: results/krd_M1_fit.rds, results/krd_M0_fit.rds
## Saves: all CSV result files
###############################################################################

base_dir <- getwd()
source(file.path(base_dir, "R", "00_bootstrap.R"))
source(file.path(base_dir, "R", "01_lib_core.R"))

out_dir <- file.path(base_dir, "results")

# Observation metadata: needed to map omega_0's per-observation columns back to
# (patient, timepoint). omega_0 is indexed by observation, not by patient.
md <- readRDS(file.path(out_dir, "model_data.rds"))

cat("Loading M1 fit...\n")
multi_M1 <- readRDS(file.path(out_dir, "krd_M1_fit.rds"))
cat("Loading M0 fit...\n")
multi_M0 <- readRDS(file.path(out_dir, "krd_M0_fit.rds"))

# ===========================================================================
# Helper: generate per-chain overlaid traceplots for key parameters
# ===========================================================================
generate_traceplots <- function(chains_M1, out_path, title_prefix = "") {
  n_ch <- length(chains_M1)
  N_obs <- ncol(chains_M1[[1]]$final_omega_g)
  K_trace <- chains_M1[[1]]$K_trace
  chain_cols <- c("steelblue", "firebrick", "forestgreen", "darkorange")

  plot_chain_scalar <- function(param_name, ylab_expr, main_expr) {
    vals <- lapply(1:n_ch, function(i) chains_M1[[i]]$samples[[param_name]])
    y_range <- range(unlist(vals))
    plot(vals[[1]], type = "l", col = chain_cols[1], ylim = y_range,
         xlab = "Iteration (post burn-in, thinned)", ylab = ylab_expr, main = main_expr)
    for (ch in 2:n_ch) lines(vals[[ch]], col = chain_cols[ch])
    legend("topright", paste("Chain", 1:n_ch), col = chain_cols[1:n_ch],
           lty = 1, cex = 0.6, bg = "white")
  }

  plot_chain_omega <- function(col_idx, ylab_expr, main_expr, trace_name = "omega_g_trace") {
    vals <- lapply(1:n_ch, function(i) chains_M1[[i]]$samples[[trace_name]][, col_idx])
    y_range <- range(unlist(vals))
    plot(vals[[1]], type = "l", col = chain_cols[1], ylim = y_range,
         xlab = "Iteration (post burn-in, thinned)", ylab = ylab_expr, main = main_expr)
    for (ch in 2:n_ch) lines(vals[[ch]], col = chain_cols[ch])
    legend("topright", paste("Chain", 1:n_ch), col = chain_cols[1:n_ch],
           lty = 1, cex = 0.6, bg = "white")
  }

  pdf(out_path, width = 12, height = 17)
  par(mfrow = c(5, 2), mar = c(4, 4, 3, 1))

  pfx <- if (nchar(title_prefix) > 0) paste0(title_prefix, ": ") else ""

  # Row 1: K+ and loglik
  plot_chain_scalar("K_plus", expression(K^"+"),
                    bquote(.(pfx) * "Trace: " * K^"+"))
  plot_chain_scalar("loglik", "Log-likelihood",
                    bquote(.(pfx) * "Log-likelihood (" * M[1] * ")"))

  # Row 2: sigma_0 and omega_0 (observation 1)
  plot_chain_scalar("sigma_0", expression(sigma[0]),
                    bquote(.(pfx) * "Trace: " * sigma[0]))
  wN_vals <- lapply(1:n_ch, function(i) chains_M1[[i]]$samples$omega_0[, 1])
  y_range_wN <- range(unlist(wN_vals))
  plot(wN_vals[[1]], type = "l", col = chain_cols[1], ylim = y_range_wN,
       xlab = "Iteration (post burn-in, thinned)", ylab = expression(omega[0]),
       main = bquote(.(pfx) * omega[0] * " (observation 1)"))
  for (ch in 2:n_ch) lines(wN_vals[[ch]], col = chain_cols[ch])
  legend("topright", paste("Chain", 1:n_ch), col = chain_cols[1:n_ch],
         lty = 1, cex = 0.6, bg = "white")

  # Row 3: omega_g and omega_cf component 1, observation 1
  if (!is.null(chains_M1[[1]]$samples$omega_g_trace)) {
    plot_chain_omega(1, expression(omega[1]^g),
                     bquote(.(pfx) * omega[1]^g * " (obs 1, comp 1)"), "omega_g_trace")
    plot_chain_omega(1, expression(omega[1]^cf),
                     bquote(.(pfx) * omega[1]^cf * " (obs 1, comp 1)"), "omega_cf_trace")
  } else {
    plot.new(); text(0.5, 0.5, "omega traces\nnot stored")
    plot.new(); text(0.5, 0.5, "omega traces\nnot stored")
  }

  # Row 4: omega_g and omega_cf component 2, observation 1
  if (!is.null(chains_M1[[1]]$samples$omega_g_trace) && K_trace >= 2) {
    col2 <- N_obs + 1
    plot_chain_omega(col2, expression(omega[2]^g),
                     bquote(.(pfx) * omega[2]^g * " (obs 1, comp 2)"), "omega_g_trace")
    plot_chain_omega(col2, expression(omega[2]^cf),
                     bquote(.(pfx) * omega[2]^cf * " (obs 1, comp 2)"), "omega_cf_trace")
  } else {
    plot.new(); text(0.5, 0.5, "comp 2 traces\nnot available")
    plot.new(); text(0.5, 0.5, "comp 2 traces\nnot available")
  }

  # Row 5: log-likelihood
  plot_chain_scalar("loglik", "Log-likelihood",
                    bquote(.(pfx) * "Log-likelihood (" * M[1] * ")"))
  plot.new()
  text(0.5, 0.5, paste0(pfx, "\n", n_ch, " chains"))

  dev.off()
}

# ===========================================================================
# (i) BEFORE label switching: generate traceplots from raw chains
# ===========================================================================
cat("\n=== Generating traceplots BEFORE label switching ===\n")
fig_dir <- file.path(base_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE)
generate_traceplots(multi_M1$chains,
                    file.path(fig_dir, "krd_traceplots_before_relabeling.pdf"),
                    "Before relabeling")
cat("  Saved: krd_traceplots_before_relabeling.pdf\n")

# --- Label switching resolution ---
# Reorder components by decreasing pi at each saved iteration (per source under M0).
cat("\n=== Label Switching Resolution ===\n")
cat("  M1:\n")
for (i in seq_along(multi_M1$chains)) {
  multi_M1$chains[[i]] <- relabel_by_weight(multi_M1$chains[[i]])
  cat(sprintf("    Chain %d relabeled\n", i))
}
cat("  Aligning M1 component labels across chains...\n")
multi_M1$chains <- align_chains_by_signature(multi_M1$chains)

cat("  M0:\n")
for (i in seq_along(multi_M0$chains)) {
  multi_M0$chains[[i]] <- relabel_by_weight(multi_M0$chains[[i]])
  cat(sprintf("    Chain %d relabeled\n", i))
}
cat("  Aligning M0 component labels across chains...\n")
multi_M0$chains <- align_chains_by_signature(multi_M0$chains)

# ===========================================================================
# (ii) AFTER label switching: generate traceplots from relabeled chains
# ===========================================================================
cat("\n=== Generating traceplots AFTER label switching ===\n")
generate_traceplots(multi_M1$chains,
                    file.path(fig_dir, "krd_traceplots_after_relabeling.pdf"),
                    "After relabeling")
cat("  Saved: krd_traceplots_after_relabeling.pdf\n")

# Re-merge after relabeling + alignment
multi_M1$merged <- merge_chain_samples(multi_M1$chains)
multi_M0$merged <- merge_chain_samples(multi_M0$chains)
cat("  Merged samples updated\n")

# ===========================================================================
# Signature-hierarchy diagnostics
#
# RETIRED: rho^theta_{kt} = cor(theta_g[,k], theta_cf[,k]) previously lived here.
# Under the current model a single theta is SHARED by both sources, so there is
# no cf-vs-g signature pair to correlate and the quantity does not exist
# (see the retirement note in R/01_lib_core.R). `krd_rho_theta.csv` is no
# longer produced.
#
# It is replaced by the two diagnostics the revised specification requires:
#   theta_separation_check()   -- empirical check of Assumption A1, i.e. whether
#                                 the between-subject dispersion varsigma_k^2 is
#                                 smaller than the population component
#                                 separation. If it fails, component labels lose
#                                 their cross-subject meaning.
#   theta_label_switch_rate()  -- fraction of subjects whose fitted component
#                                 ordering disagrees with the cohort ordering.
# ===========================================================================
cat("\n=== Signature-hierarchy diagnostics ===\n")
tryCatch({
  diag_rows <- list()
  for (lab in c("M1", "M0")) {
    mf <- if (lab == "M1") multi_M1 else multi_M0
    ch <- mf$chains[[1]]
    sep <- theta_separation_check(ch)
    ls_ <- theta_label_switch_rate(ch)
    cat(sprintf("  %s: A1 separation satisfied = %s | label-switch rate = %.3f (%d/%d subjects)\n",
                lab, isTRUE(attr(sep, "all_satisfied")), ls_$rate,
                length(ls_$switched), ls_$n_subj))
    if (nrow(sep)) {
      sep$model <- lab
      sep$label_switch_rate <- ls_$rate
      diag_rows[[length(diag_rows) + 1L]] <- sep
    }
  }
  if (length(diag_rows)) {
    dd <- do.call(rbind, diag_rows)
    write.csv(dd, file.path(out_dir, "krd_theta_separation.csv"), row.names = FALSE)
    cat(sprintf("  Saved: krd_theta_separation.csv (%d rows)\n", nrow(dd)))
  }
}, error = function(e) {
  cat("signature-hierarchy diagnostics FAILED:", conditionMessage(e), "\n")
})
# --- Log-likelihood spread: identify trapped chains ---
trap_M1 <- identify_trapped_chains(multi_M1, ll_gap = 10000)
print_ll_spread(trap_M1, label = "M1")
trap_M0 <- identify_trapped_chains(multi_M0, ll_gap = 10000)
print_ll_spread(trap_M0, label = "M0")

# --- Convergence diagnostics ---
cat("\n=== Convergence Diagnostics ===\n")
# alpha_dp and kappa are both fixed — not monitored
# Monitor sigma_0, K_plus, omega_0 (per-OBSERVATION) and the new hierarchical
# dispersions. varsigma_* are a known locus of poor mixing (funnel geometry),
# so they are monitored explicitly rather than folded into a summary.
# and expand to N_obs*K_trace columns; compute_convergence handles matrices.
conv_params <- c("sigma_0", "K_plus", "omega_0", "varsigma_0")
conv_M1 <- compute_convergence(multi_M1, params = conv_params)
print_convergence(conv_M1)

conv_M0 <- compute_convergence(multi_M0, params = conv_params)
print_convergence(conv_M0)

# --- Model comparison ---
cat("\n=== Model Comparison ===\n")
# Full: all chains (original behavior, for comparison)
waic_M1_full <- NULL; waic_M0_full <- NULL
tryCatch({
  waic_M1_full <- loo::waic(multi_M1$merged$pointwise_ll)
  waic_M0_full <- loo::waic(multi_M0$merged$pointwise_ll)
  dw_full <- waic_M0_full$estimates["waic", "Estimate"] -
             waic_M1_full$estimates["waic", "Estimate"]
  cat(sprintf("[ALL chains]  Delta WAIC (M0-M1): %.2f (positive favors M1)\n", dw_full))
}, error = function(e) {
  cat("WAIC (all chains) failed:", conditionMessage(e), "\n")
})

# Clean: exclude chains flagged as trapped (LL spread check)
waic_M1_clean <- NULL; waic_M0_clean <- NULL
tryCatch({
  if (length(trap_M1$keep) >= 1 && length(trap_M0$keep) >= 1) {
    ss_M1 <- waic_loo_subset(multi_M1, trap_M1$keep)
    ss_M0 <- waic_loo_subset(multi_M0, trap_M0$keep)
    waic_M1_clean <- ss_M1$waic
    waic_M0_clean <- ss_M0$waic
    dw_clean <- waic_M0_clean$estimates["waic", "Estimate"] -
                waic_M1_clean$estimates["waic", "Estimate"]
    cat(sprintf("[CLEAN chains] Delta WAIC (M0-M1): %.2f  (M1 kept: %s; M0 kept: %s)\n",
                dw_clean,
                paste(trap_M1$keep, collapse = ","),
                paste(trap_M0$keep, collapse = ",")))
    if (!is.null(waic_M1_full)) {
      pct <- 100 * abs(dw_clean - dw_full) / max(abs(dw_full), 1)
      cat(sprintf("Shift in Delta WAIC after excluding trapped: %.2f (%.1f%%)\n",
                  dw_clean - dw_full, pct))
    }
  } else {
    cat("Not enough non-trapped chains in one or both models for clean WAIC.\n")
  }
}, error = function(e) {
  cat("WAIC (clean) failed:", conditionMessage(e), "\n")
})

# Assign legacy names for downstream code
waic_M1 <- if (!is.null(waic_M1_clean)) waic_M1_clean else waic_M1_full
waic_M0 <- if (!is.null(waic_M0_clean)) waic_M0_clean else waic_M0_full

# --- Save CSVs ---

# Posterior summary
wN_all_iters <- rowMeans(multi_M1$merged$omega_0)
# alpha_dp=1 and kappa=1 are both fixed — report as constants
summary_df <- data.frame(
  parameter = c("K_plus", "sigma_0", "omega_0_mean", "kappa", "alpha_dp"),
  mean = c(mean(multi_M1$merged$K_plus), mean(multi_M1$merged$sigma_0),
           mean(wN_all_iters), 1.0, 1.0),
  median = c(median(multi_M1$merged$K_plus), median(multi_M1$merged$sigma_0),
             median(wN_all_iters), 1.0, 1.0),
  ci_lo = c(quantile(multi_M1$merged$K_plus, 0.025), quantile(multi_M1$merged$sigma_0, 0.025),
            quantile(wN_all_iters, 0.025), 1.0, 1.0),
  ci_hi = c(quantile(multi_M1$merged$K_plus, 0.975), quantile(multi_M1$merged$sigma_0, 0.975),
            quantile(wN_all_iters, 0.975), 1.0, 1.0)
)
write.csv(summary_df, file.path(out_dir, "krd_posterior_summary.csv"), row.names = FALSE)

# Model comparison: both ALL-chains and CLEAN-chains rows
comparison_df <- data.frame(
  metric = c("delta_WAIC_all", "delta_LOO_all",
             "delta_WAIC_clean", "delta_LOO_clean"),
  value = NA_real_, favors = NA_character_, stringsAsFactors = FALSE
)

# Delta WAIC: all
tryCatch({
  dw <- waic_M0_full$estimates["waic", "Estimate"] -
        waic_M1_full$estimates["waic", "Estimate"]
  comparison_df$value[1] <- dw
  comparison_df$favors[1] <- if (dw > 0) "M1" else "M0"
}, error = function(e) NULL)

# Delta LOO: all
tryCatch({
  loo_M1_full <- loo::loo(multi_M1$merged$pointwise_ll)
  loo_M0_full <- loo::loo(multi_M0$merged$pointwise_ll)
  dl <- loo_M0_full$estimates["looic", "Estimate"] -
        loo_M1_full$estimates["looic", "Estimate"]
  comparison_df$value[2] <- dl
  comparison_df$favors[2] <- if (dl > 0) "M1" else "M0"
}, error = function(e) cat("LOO (all) failed:", conditionMessage(e), "\n"))

# Delta WAIC: clean
tryCatch({
  if (!is.null(waic_M1_clean) && !is.null(waic_M0_clean)) {
    dw <- waic_M0_clean$estimates["waic", "Estimate"] -
          waic_M1_clean$estimates["waic", "Estimate"]
    comparison_df$value[3] <- dw
    comparison_df$favors[3] <- if (dw > 0) "M1" else "M0"
  }
}, error = function(e) NULL)

# Delta LOO: clean
tryCatch({
  if (length(trap_M1$keep) >= 1 && length(trap_M0$keep) >= 1) {
    loo_M1_clean <- waic_loo_subset(multi_M1, trap_M1$keep)$loo
    loo_M0_clean <- waic_loo_subset(multi_M0, trap_M0$keep)$loo
    if (!is.null(loo_M1_clean) && !is.null(loo_M0_clean)) {
      dl <- loo_M0_clean$estimates["looic", "Estimate"] -
            loo_M1_clean$estimates["looic", "Estimate"]
      comparison_df$value[4] <- dl
      comparison_df$favors[4] <- if (dl > 0) "M1" else "M0"
    }
  }
}, error = function(e) cat("LOO (clean) failed:", conditionMessage(e), "\n"))

write.csv(comparison_df, file.path(out_dir, "krd_model_comparison.csv"), row.names = FALSE)

# Chain health CSV (so we have a record of what was trapped)
chain_health_df <- data.frame(
  model = c(rep("M1", length(trap_M1$ll_means)), rep("M0", length(trap_M0$ll_means))),
  chain = c(seq_along(trap_M1$ll_means), seq_along(trap_M0$ll_means)),
  mean_loglik = c(trap_M1$ll_means, trap_M0$ll_means),
  gap_from_best = c(trap_M1$best - trap_M1$ll_means,
                    trap_M0$best - trap_M0$ll_means),
  trapped = c(seq_along(trap_M1$ll_means) %in% trap_M1$drop,
              seq_along(trap_M0$ll_means) %in% trap_M0$drop)
)
write.csv(chain_health_df, file.path(out_dir, "krd_chain_health.csv"), row.names = FALSE)

# Contamination fractions.
#
# omega_0 is now indexed by OBSERVATION, not patient, so the primary table has
# one row per (patient, timepoint). A per-patient aggregate is also written
# because several manuscript tables want one row per patient -- but note it is
# an AVERAGE over that patient's timepoints and therefore discards exactly the
# temporal variation the respecification was introduced to capture. Prefer the
# per-observation file for anything substantive.
om0 <- multi_M1$merged$omega_0
stopifnot(ncol(om0) == md$N_obs)
omega0_obs_df <- data.frame(
  obs       = seq_len(ncol(om0)),
  patient   = md$patient,
  timepoint = md$time,
  post_mean = colMeans(om0),
  ci_lo     = apply(om0, 2, quantile, 0.025),
  ci_hi     = apply(om0, 2, quantile, 0.975)
)
write.csv(omega0_obs_df, file.path(out_dir, "krd_omega0_by_observation.csv"),
          row.names = FALSE)

## Per-patient aggregate. NOTE this discards the temporal variation that
## omega_0 was respecified to capture; krd_omega0_by_observation.csv is the
## primary artifact and this one exists for patient-level reporting.
##
## The credible interval is formed by averaging each patient's observations
## WITHIN every posterior draw and then taking quantiles across draws, so it is
## a genuine posterior interval for the patient mean. (Taking quantiles of the
## per-observation posterior means instead would report between-observation
## spread, which is a different quantity.) Downstream consumers -- Table 4 in
## R/08_tables.R and the per-patient table in R/09_refresh_manuscript.R -- read
## ci_lo/ci_hi, and previously found them only in the retired
## krd_wN_by_patient.csv.
pat_ids <- sort(unique(md$patient))
om0_pat <- vapply(pat_ids, function(pp) {
  cols <- which(md$patient == pp)
  rowMeans(om0[, cols, drop = FALSE])      # [n_save] draws of the patient mean
}, numeric(nrow(om0)))                      # [n_save x n_patients]

omega0_pat_df <- data.frame(
  patient   = pat_ids,
  post_mean = colMeans(om0_pat),
  ci_lo     = apply(om0_pat, 2, quantile, 0.025),
  ci_hi     = apply(om0_pat, 2, quantile, 0.975),
  n_obs     = as.integer(table(md$patient)[as.character(pat_ids)])
)
write.csv(omega0_pat_df, file.path(out_dir, "krd_omega0_by_patient.csv"),
          row.names = FALSE)
cat(sprintf("  Saved: krd_omega0_by_observation.csv (%d rows), krd_omega0_by_patient.csv (%d rows)\n",
            nrow(omega0_obs_df), nrow(omega0_pat_df)))

# Convergence
write.csv(conv_M1$summary, file.path(out_dir, "krd_convergence_M1.csv"), row.names = FALSE)
write.csv(conv_M0$summary, file.path(out_dir, "krd_convergence_M0.csv"), row.names = FALSE)

cat("\nAll core results saved.\n")

# ===========================================================================
# Derived results needed for both manuscripts. Each script reads from results/
# and writes additional artifacts there (and a few PDFs under figures/ or
# nature_manuscript/figures/). They are self-contained sources of CSV outputs
# that downstream table/figure scripts consume.
# ===========================================================================

cat("\n=== Derived results: posterior concordance metrics (R/04_concordance_metrics) ===\n")
tryCatch({
  source(file.path(base_dir, "R", "04_concordance_metrics.R"))
}, error = function(e) {
  cat("R/04_concordance_metrics FAILED:", conditionMessage(e), "\n")
})

cat("\n=== Derived results: Pareto-k diagnostics (R/05_pareto_loo_diagnostic) ===\n")
tryCatch({
  source(file.path(base_dir, "R", "05_pareto_loo_diagnostic.R"))
}, error = function(e) {
  cat("R/05_pareto_loo_diagnostic FAILED:", conditionMessage(e), "\n")
})

cat("\n=== Derived results: per-component variance diagnostic (R/06_sigma_k_diagnostic) ===\n")
tryCatch({
  source(file.path(base_dir, "R", "06_sigma_k_diagnostic.R"))
}, error = function(e) {
  cat("R/06_sigma_k_diagnostic FAILED:", conditionMessage(e), "\n")
})

cat("\nAll core and derived results saved.\n")

sink(file.path(out_dir, "sessionInfo_krd_production.txt"))
cat("KRd production analysis:", format(Sys.time()), "\n\n")
sessionInfo()
sink()
