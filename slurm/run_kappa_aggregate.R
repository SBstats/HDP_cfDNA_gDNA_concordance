#!/usr/bin/env Rscript
###############################################################################
## slurm/run_kappa_aggregate.R  --  build Supplementary Table SX from fits
##
## Run AFTER all 8 array tasks of run_kappa_sensitivity.R have produced
## results/kappa/fit_<model>_kappa<kappa>.rds. Short job (minutes).
##
## For each kappa: Delta WAIC (M0 - M1; +ve favours M1) on clean chains,
## K+ median under M1, mean L1 subclonal concordance under M1, rho exact.
## Output: results/kappa_sensitivity.csv
###############################################################################

base_dir <- getwd()
source(file.path(base_dir, "R", "01_lib_core.R"))
suppressPackageStartupMessages(library(loo))

md <- readRDS(file.path(base_dir, "results", "model_data.rds")); N <- md$N_obs
out_dir <- file.path(base_dir, "results", "kappa")
alpha <- 1.0; kappas <- c(0.5, 1.0, 2.0, 5.0)

waic_clean <- function(fit) {
  trap <- identify_trapped_chains(fit, ll_gap = 10000)
  keep <- if (length(trap$keep)) trap$keep else seq_along(fit$chains)
  ll <- do.call(rbind, lapply(fit$chains[keep], function(c) c$samples$pointwise_ll))
  suppressWarnings(loo::waic(ll)$estimates["waic", "Estimate"])
}
mean_L1 <- function(fit) {
  for (i in seq_along(fit$chains)) fit$chains[[i]] <- relabel_by_weight(fit$chains[[i]])
  fit$chains <- align_chains_by_signature(fit$chains)
  m <- merge_chain_samples(fit$chains); K <- ncol(m$omega_g_trace) / N
  Wg  <- sapply(1:K, function(k) sapply(1:N, function(o) mean(m$omega_g_trace[, (k-1)*N+o])))
  Wcf <- sapply(1:K, function(k) sapply(1:N, function(o) mean(m$omega_cf_trace[, (k-1)*N+o])))
  mean(sapply(1:N, function(o) 1 - 0.5 * sum(abs(Wg[o, ] - Wcf[o, ]))))
}
med_Kplus <- function(fit) median(unlist(lapply(fit$chains, function(c) c$samples$K_plus)))

rows <- list()
for (kap in kappas) {
  f1p <- file.path(out_dir, sprintf("fit_M1_kappa%.2f.rds", kap))
  f0p <- file.path(out_dir, sprintf("fit_M0_kappa%.2f.rds", kap))
  if (!file.exists(f1p) || !file.exists(f0p)) {
    cat(sprintf("  kappa=%.2f: missing fit(s), skipping\n", kap)); next
  }
  f1 <- readRDS(f1p); f0 <- readRDS(f0p)
  rows[[length(rows)+1]] <- data.frame(
    kappa = kap,
    delta_WAIC = waic_clean(f0) - waic_clean(f1),
    Kplus_median = med_Kplus(f1),
    mean_L1_concordance = mean_L1(f1),
    rho_exact = (kap + 1) / (alpha + kap + 1))
  cat(sprintf("  kappa=%.2f done\n", kap))
}
res <- do.call(rbind, rows)
write.csv(res, file.path(base_dir, "results", "kappa_sensitivity.csv"), row.names = FALSE)
cat("\nSaved results/kappa_sensitivity.csv\n"); print(res)
