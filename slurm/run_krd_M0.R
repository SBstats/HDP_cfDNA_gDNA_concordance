#!/usr/bin/env Rscript
###############################################################################
## slurm/run_krd_M0.R
## Fit M0 (non-tracking model) on KRd trial data — 4 parallel chains
## Saves: results/krd_M0_fit.rds
###############################################################################

base_dir <- getwd()
cat("base_dir:", base_dir, "\n")
cat("Loading sampler...\n")
source(file.path(base_dir, "R", "00_bootstrap.R"))
source(file.path(base_dir, "R", "01_lib_core.R"))
if (!exists(".rcpp_available") || !.rcpp_available) {
  stop("FATAL: Rcpp compilation failed. Cannot run production analysis without C++ acceleration.\n",
       "Try: module load gcc; module load R/4.3.1; then resubmit.")
}
cat("Rcpp OK\n")

model_data <- readRDS(file.path(base_dir, "results", "model_data.rds"))
Y_g  <- model_data$Y_g
Y_cf <- model_data$Y_cf

cat(sprintf("KRd M0: p=%d, n=%d, N_obs=%d, T=%d\n",
            nrow(Y_g), model_data$n, model_data$N_obs, max(model_data$time)))

cat("\n=== Fitting M0 (Non-tracking, 4 chains) ===\n")
multi_M0 <- fit_multi_chain(
  Y_g, Y_cf, model_data$patient, model_data$time,
  n_chains = 4, seed = 100,
  K = 10, n_iter = 22500, n_burn = 20000, thin = 10,
  model = "M0",
  anneal = TRUE, T_anneal = 20.0, n_cool_buffer = 10000,
  store_pointwise_ll = TRUE
)

out_dir <- file.path(base_dir, "results")
dir.create(out_dir, showWarnings = FALSE)
saveRDS(multi_M0, file.path(out_dir, "krd_M0_fit.rds"))
cat("\nM0 fit saved to results/krd_M0_fit.rds\n")
cat("Done at:", format(Sys.time()), "\n")
