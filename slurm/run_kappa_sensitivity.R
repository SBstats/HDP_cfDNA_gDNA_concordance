#!/usr/bin/env Rscript
###############################################################################
## slurm/run_kappa_sensitivity.R  --  FIT ONE (model, kappa) PER ARRAY TASK
##
## Supplementary Table SX: robustness of the concordance verdict to the fixed
## tracking precision kappa. caslake caps wall clock at 36h and each production
## fit is ~18-26h, so we fit exactly ONE model at ONE kappa per array task.
##
## Task grid (SLURM_ARRAY_TASK_ID = 0..7):
##   idx  kappa  model
##    0   0.5    M1        4   2.0   M1
##    1   0.5    M0        5   2.0   M0
##    2   1.0    M1        6   5.0   M1
##    3   1.0    M0        7   5.0   M0
##
## Each task writes results/kappa/fit_<model>_kappa<kappa>.rds.
## After all 8 tasks finish, run slurm/run_kappa_aggregate.R (06b_*.sbatch).
###############################################################################

base_dir <- getwd()
source(file.path(base_dir, "R", "00_bootstrap.R"))
source(file.path(base_dir, "R", "01_lib_core.R"))
if (!exists(".rcpp_available") || !.rcpp_available)
  stop("FATAL: Rcpp required for production fits.")

md <- readRDS(file.path(base_dir, "results", "model_data.rds"))
out_dir <- file.path(base_dir, "results", "kappa")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## build the 8-row task grid
grid <- expand.grid(model = c("M1", "M0"), kappa = c(0.5, 1.0, 2.0, 5.0),
                    stringsAsFactors = FALSE)
grid <- grid[order(grid$kappa, grid$model == "M0"), ]   # M1 before M0 within kappa
rownames(grid) <- NULL

tid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "-1"))
tasks <- if (tid >= 0) (tid + 1) else seq_len(nrow(grid))   # no array => run all (local test)

for (r in tasks) {
  kap <- grid$kappa[r]; model <- grid$model[r]
  cat(sprintf("\n=== Task %d: fit %s at kappa=%.2f ===\n", r - 1, model, kap))
  fit <- fit_multi_chain(md$Y_g, md$Y_cf, md$patient, md$time,
                         n_chains = 4, seed = 42, K = 10,
                         n_iter = 22500, n_burn = 20000, thin = 10,
                         model = model, kappa_fixed = kap,
                         anneal = TRUE, T_anneal = 20.0, n_cool_buffer = 10000,
                         store_pointwise_ll = TRUE)
  f <- file.path(out_dir, sprintf("fit_%s_kappa%.2f.rds", model, kap))
  saveRDS(fit, f); cat("  saved", f, "\n")
}
cat("\nTask(s) complete. After all 8 finish, run slurm/run_kappa_aggregate.R\n")
