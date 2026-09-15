#!/usr/bin/env Rscript
###############################################################################
## slurm/run_krd_theta_trace.R
## Tier 2 confirmatory re-fit that STORES A THINNED THETA TRACE so the shared
## subclonal signatures theta can be formally audited with between-chain
## split-Rhat (R/16_theta_convergence_audit.R).
##
## Identical production configuration to run_krd_M1.R / run_krd_M0.R, plus
## store_theta_trace = TRUE. Model chosen by the first command-line arg.
##   Rscript slurm/run_krd_theta_trace.R M1   # or M0  (default M1)
## Saves: results/krd_<model>_fit_thetatrace.rds
##
## NOTE: the local Tier 2 between-chain theta agreement audit already shows the
## occupied-component signatures agree across chains at corr >= 0.994, so this
## re-fit is CONFIRMATORY (produces the formal per-gene Rhat).
###############################################################################

args  <- commandArgs(trailingOnly = TRUE)
model <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "M1"
stopifnot(model %in% c("M1", "M0"))

base_dir <- getwd()
cat("base_dir:", base_dir, "| model:", model, "\n")
source(file.path(base_dir, "R", "00_bootstrap.R"))
source(file.path(base_dir, "R", "01_lib_core.R"))
if (!exists(".rcpp_available") || !.rcpp_available) {
  stop("FATAL: Rcpp compilation failed. Try: module load gcc; module load R; then resubmit.")
}
cat("Rcpp OK\n")

model_data <- readRDS(file.path(base_dir, "results", "model_data.rds"))
cat(sprintf("KRd %s theta-trace: p=%d, n=%d, N_obs=%d, T=%d\n",
            model, nrow(model_data$Y_g), model_data$n, model_data$N_obs, max(model_data$time)))

fit <- fit_multi_chain(
  model_data$Y_g, model_data$Y_cf, model_data$patient, model_data$time,
  n_chains = 4, seed = 42,
  K = 10, n_iter = 22500, n_burn = 20000, thin = 10,
  model = model, kappa_fixed = 1.0,
  anneal = TRUE, T_anneal = 20.0, n_cool_buffer = 10000,
  store_pointwise_ll = TRUE,
  store_theta_trace = TRUE, theta_trace_thin = 5     # 250 saved / 5 = 50 theta draws/chain
)

out_dir <- file.path(base_dir, "results")
dir.create(out_dir, showWarnings = FALSE)
saveRDS(fit, file.path(out_dir, sprintf("krd_%s_fit_thetatrace.rds", model)))
cat(sprintf("\nSaved results/krd_%s_fit_thetatrace.rds\n", model))
cat("Done at:", format(Sys.time()), "\n")
