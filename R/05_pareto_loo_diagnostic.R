###############################################################################
## 24_loo_pareto_diagnostic.R
##
## Round 3 diagnostic (Priority 1a): extract Pareto k-hat values from the PSIS-LOO
## computation on the existing M1 and M0 fits, to assess whether the ΔLOO value
## reported in results/krd_model_comparison.csv is actually trustworthy.
##
## Background:
## The model stores pointwise_ll as [n_save x (2*N_obs)], where each column is a
## sum over p = 11,837 gene-level log-likelihoods for a single (observation,
## source) pair. Leaving out one such aggregate unit in PSIS-LOO may produce
## heavy-tailed importance weights because the posterior for θ_jk(t) and the
## horseshoe hyperparameters τ_k² re-shape substantially when ~12k gene
## contributions are dropped. This diagnostic quantifies how pervasive the
## pathology is via Pareto k-hat.
##
## Inputs:
##   results/krd_M1_fit.rds, results/krd_M0_fit.rds
##
## Outputs:
##   results/krd_loo_M1.rds, results/krd_loo_M0.rds (full loo objects)
##   results/krd_pareto_k_M1.csv, results/krd_pareto_k_M0.csv
##   results/krd_loo_reliability_summary.csv
##
## Pass/fail rule: if more than ~10% of the 2*N_obs = 168 units have k̂ > 0.7,
## PSIS-LOO is formally unreliable and ΔLOO should be de-emphasized.
###############################################################################

suppressPackageStartupMessages({
  library(loo)
})

base_dir <- getwd()
results_dir <- file.path(base_dir, "results")

cat("=== Pareto k-hat diagnostic on existing fits ===\n\n")

# ---------------------------------------------------------------------------
# Load fits and run PSIS-LOO (keep full loo objects)
# ---------------------------------------------------------------------------
multi_M1 <- readRDS(file.path(results_dir, "krd_M1_fit.rds"))
multi_M0 <- readRDS(file.path(results_dir, "krd_M0_fit.rds"))

ll_M1 <- multi_M1$merged$pointwise_ll
ll_M0 <- multi_M0$merged$pointwise_ll

n_obs_pairs <- ncol(ll_M1) / 2  # 84 observations; columns 1..84 = gDNA, 85..168 = cfDNA
cat("pointwise_ll dimensions:\n")
cat("  M1:", nrow(ll_M1), "iterations x", ncol(ll_M1), "units\n")
cat("  M0:", nrow(ll_M0), "iterations x", ncol(ll_M0), "units\n")
cat("  N_obs =", n_obs_pairs, "; columns split gDNA (1..", n_obs_pairs,
    ") and cfDNA (", n_obs_pairs+1, "..", 2*n_obs_pairs, ")\n\n")

cat("Running loo::loo() on M1...\n")
loo_M1 <- loo::loo(ll_M1)
cat("Running loo::loo() on M0...\n")
loo_M0 <- loo::loo(ll_M0)

saveRDS(loo_M1, file.path(results_dir, "krd_loo_M1.rds"))
saveRDS(loo_M0, file.path(results_dir, "krd_loo_M0.rds"))

# ---------------------------------------------------------------------------
# Extract per-unit Pareto k-hat and metadata
# ---------------------------------------------------------------------------
make_k_df <- function(loo_obj, model_label) {
  k_vec <- loo_obj$diagnostics$pareto_k
  N_units <- length(k_vec)
  n_pairs <- N_units / 2
  data.frame(
    model     = model_label,
    unit      = seq_len(N_units),
    source    = c(rep("gDNA", n_pairs), rep("cfDNA", n_pairs)),
    obs_idx   = rep(seq_len(n_pairs), 2),
    pareto_k  = k_vec,
    stringsAsFactors = FALSE
  )
}

k_M1 <- make_k_df(loo_M1, "M1")
k_M0 <- make_k_df(loo_M0, "M0")

write.csv(k_M1, file.path(results_dir, "krd_pareto_k_M1.csv"), row.names = FALSE)
write.csv(k_M0, file.path(results_dir, "krd_pareto_k_M0.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# Summary tables and reliability verdict
# ---------------------------------------------------------------------------
bin_k <- function(k) {
  cut(k,
      breaks = c(-Inf, 0.5, 0.7, 1.0, Inf),
      labels = c("good (<=0.5)", "ok (0.5-0.7)", "bad (0.7-1.0)", "very bad (>1.0)"),
      right  = TRUE)
}

summarize_k <- function(k_df, label) {
  tbl_all <- table(bin_k(k_df$pareto_k))
  tbl_g   <- table(bin_k(k_df$pareto_k[k_df$source == "gDNA"]))
  tbl_cf  <- table(bin_k(k_df$pareto_k[k_df$source == "cfDNA"]))

  n_total      <- nrow(k_df)
  n_bad        <- sum(k_df$pareto_k > 0.7)
  n_very_bad   <- sum(k_df$pareto_k > 1.0)
  pct_bad      <- 100 * n_bad / n_total
  pct_very_bad <- 100 * n_very_bad / n_total

  cat(sprintf("\n--- %s Pareto k-hat summary ---\n", label))
  cat("All units (n =", n_total, "):\n"); print(tbl_all)
  cat("gDNA only (n =", n_total/2, "):\n"); print(tbl_g)
  cat("cfDNA only (n =", n_total/2, "):\n"); print(tbl_cf)
  cat(sprintf("Max k-hat:            %.3f\n", max(k_df$pareto_k)))
  cat(sprintf("Units with k > 0.7:   %d / %d  (%.1f%%)\n", n_bad, n_total, pct_bad))
  cat(sprintf("Units with k > 1.0:   %d / %d  (%.1f%%)\n",
              n_very_bad, n_total, pct_very_bad))

  list(
    model          = label,
    n_total        = n_total,
    n_good         = as.integer(tbl_all["good (<=0.5)"]),
    n_ok           = as.integer(tbl_all["ok (0.5-0.7)"]),
    n_bad          = as.integer(tbl_all["bad (0.7-1.0)"]),
    n_very_bad     = as.integer(tbl_all["very bad (>1.0)"]),
    pct_above_0_7  = pct_bad,
    pct_above_1_0  = pct_very_bad,
    max_k          = max(k_df$pareto_k),
    reliable       = pct_bad < 10
  )
}

s_M1 <- summarize_k(k_M1, "M1")
s_M0 <- summarize_k(k_M0, "M0")

reliability_df <- do.call(rbind, lapply(list(s_M1, s_M0), as.data.frame))
write.csv(reliability_df, file.path(results_dir, "krd_loo_reliability_summary.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# Patient / timepoint cross-tab of very-bad k-hat (if model_data available)
# ---------------------------------------------------------------------------
md_file <- file.path(results_dir, "model_data.rds")
if (file.exists(md_file)) {
  md <- readRDS(md_file)
  k_M1$patient <- md$patient[k_M1$obs_idx]
  k_M1$time    <- md$time[k_M1$obs_idx]
  k_M0$patient <- md$patient[k_M0$obs_idx]
  k_M0$time    <- md$time[k_M0$obs_idx]

  cat("\n--- M1 k > 0.7 by (patient, time, source) ---\n")
  bad_M1 <- subset(k_M1, pareto_k > 0.7)
  if (nrow(bad_M1) > 0) {
    print(bad_M1[order(-bad_M1$pareto_k), c("patient", "time", "source", "pareto_k")],
          row.names = FALSE, digits = 3)
  } else {
    cat("  (none)\n")
  }
  cat("\n--- M0 k > 0.7 by (patient, time, source) ---\n")
  bad_M0 <- subset(k_M0, pareto_k > 0.7)
  if (nrow(bad_M0) > 0) {
    print(bad_M0[order(-bad_M0$pareto_k), c("patient", "time", "source", "pareto_k")],
          row.names = FALSE, digits = 3)
  } else {
    cat("  (none)\n")
  }
}

# ---------------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------------
cat("\n=== Reliability verdict ===\n")
verdict <- function(s) {
  if (s$pct_above_0_7 < 10) {
    sprintf("%s: PSIS-LOO RELIABLE (%.1f%% k > 0.7)", s$model, s$pct_above_0_7)
  } else if (s$pct_above_0_7 < 50) {
    sprintf("%s: PSIS-LOO QUESTIONABLE (%.1f%% k > 0.7)",
            s$model, s$pct_above_0_7)
  } else {
    sprintf("%s: PSIS-LOO UNRELIABLE (%.1f%% k > 0.7) -- de-emphasize LOO, use WAIC",
            s$model, s$pct_above_0_7)
  }
}
cat(verdict(s_M1), "\n")
cat(verdict(s_M0), "\n")

cat("\nFiles written:\n")
cat("  results/krd_loo_M1.rds, krd_loo_M0.rds\n")
cat("  results/krd_pareto_k_M1.csv, krd_pareto_k_M0.csv\n")
cat("  results/krd_loo_reliability_summary.csv\n")
