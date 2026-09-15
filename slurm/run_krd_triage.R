#!/usr/bin/env Rscript
###############################################################################
## slurm/run_krd_triage.R
## Triage: recompute WAIC/LOO on existing M1/M0 fits, including and excluding
## chains flagged by the log-likelihood spread check. Does NOT rerun MCMC.
##
## Reads:  results/krd_M1_fit.rds, results/krd_M0_fit.rds
## Writes: results/krd_triage_waic.csv, results/krd_chain_health.csv
##
## Use this to assess whether the published Delta-WAIC is being driven by
## chains trapped in low-density modes, before investing in a full rerun.
###############################################################################

base_dir <- getwd()
source(file.path(base_dir, "R", "00_bootstrap.R"))
source(file.path(base_dir, "R", "01_lib_core.R"))
suppressPackageStartupMessages(library(loo))

out_dir <- file.path(base_dir, "results")

cat("Loading fits...\n")
multi_M1 <- readRDS(file.path(out_dir, "krd_M1_fit.rds"))
multi_M0 <- readRDS(file.path(out_dir, "krd_M0_fit.rds"))

# --- LL spread ---
trap_M1 <- identify_trapped_chains(multi_M1, ll_gap = 10000)
print_ll_spread(trap_M1, label = "M1")
trap_M0 <- identify_trapped_chains(multi_M0, ll_gap = 10000)
print_ll_spread(trap_M0, label = "M0")

# --- WAIC/LOO: all chains (baseline) ---
cat("\n=== WAIC/LOO: ALL chains ===\n")
waic_M1_all <- loo::waic(multi_M1$merged$pointwise_ll)
waic_M0_all <- loo::waic(multi_M0$merged$pointwise_ll)
loo_M1_all  <- loo::loo(multi_M1$merged$pointwise_ll)
loo_M0_all  <- loo::loo(multi_M0$merged$pointwise_ll)
dw_all <- waic_M0_all$estimates["waic", "Estimate"]  - waic_M1_all$estimates["waic", "Estimate"]
dl_all <- loo_M0_all$estimates["looic", "Estimate"] - loo_M1_all$estimates["looic", "Estimate"]
cat(sprintf("Delta WAIC (all): %.2f | Delta LOO (all): %.2f\n", dw_all, dl_all))

# --- WAIC/LOO: clean chains only ---
cat("\n=== WAIC/LOO: CLEAN chains only ===\n")
if (length(trap_M1$keep) >= 1 && length(trap_M0$keep) >= 1) {
  ss_M1 <- waic_loo_subset(multi_M1, trap_M1$keep)
  ss_M0 <- waic_loo_subset(multi_M0, trap_M0$keep)
  if (!is.null(ss_M1$waic) && !is.null(ss_M0$waic)) {
    dw_clean <- ss_M0$waic$estimates["waic", "Estimate"] -
                ss_M1$waic$estimates["waic", "Estimate"]
    cat(sprintf("Delta WAIC (clean): %.2f\n", dw_clean))
    cat(sprintf("  Shift vs all: %.2f (%.1f%% of all)\n",
                dw_clean - dw_all, 100 * (dw_clean - dw_all) / max(abs(dw_all), 1)))
  } else {
    dw_clean <- NA_real_
  }
  if (!is.null(ss_M1$loo) && !is.null(ss_M0$loo)) {
    dl_clean <- ss_M0$loo$estimates["looic", "Estimate"] -
                ss_M1$loo$estimates["looic", "Estimate"]
    cat(sprintf("Delta LOO (clean):  %.2f\n", dl_clean))
    cat(sprintf("  Shift vs all: %.2f (%.1f%% of all)\n",
                dl_clean - dl_all, 100 * (dl_clean - dl_all) / max(abs(dl_all), 1)))
  } else {
    dl_clean <- NA_real_
  }
} else {
  cat("Not enough non-trapped chains; skipping clean WAIC/LOO.\n")
  dw_clean <- NA_real_
  dl_clean <- NA_real_
}

# --- Save ---
out_df <- data.frame(
  metric = c("delta_WAIC_all", "delta_LOO_all",
             "delta_WAIC_clean", "delta_LOO_clean"),
  value = c(dw_all, dl_all, dw_clean, dl_clean)
)
write.csv(out_df, file.path(out_dir, "krd_triage_waic.csv"), row.names = FALSE)

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

cat("\nTriage outputs written:\n")
cat("  ", file.path(out_dir, "krd_triage_waic.csv"), "\n")
cat("  ", file.path(out_dir, "krd_chain_health.csv"), "\n")
cat("\nDone at:", format(Sys.time()), "\n")
