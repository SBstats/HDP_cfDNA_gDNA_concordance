#!/usr/bin/env Rscript
###############################################################################
## slurm/run_kfold_cv.R  --  FIT ONE (fold, model) PER ARRAY TASK
##
## PSIS-free K-fold CV for M1 vs M0 (replaces the non-interpretable
## Delta WAIC = 3.28e6 / Pareto k-hat > 1 comparison). caslake caps wall clock
## at 36h and each fit is ~15-20h, so we fit exactly ONE model on ONE fold's
## TRAINING observations per array task.
##
## Task grid (SLURM_ARRAY_TASK_ID = 0 .. 2*F-1), F folds:
##   idx = 2*(fold-1) + {0:M1, 1:M0}
##   F=5 -> tasks 0..9.
##
## Fold assignment is DETERMINISTIC (fixed seed, stratified by timepoint) and
## identical here and in run_kfold_aggregate.R, so held-out scoring uses the
## same partition. Each task writes results/kfold/fit_f<fold>_<model>.rds.
## After all tasks finish, run slurm/run_kfold_aggregate.R (07b_*.sbatch).
##
## Env: N_FOLDS (default 5). For a fast validation run set N_FOLDS=2 and
## optionally KFOLD_NITER/KFOLD_NBURN (defaults 22500/20000).
###############################################################################

base_dir <- getwd()
source(file.path(base_dir, "R", "00_bootstrap.R"))
source(file.path(base_dir, "R", "01_lib_core.R"))
if (!exists(".rcpp_available") || !.rcpp_available)
  stop("FATAL: Rcpp required for production fits.")

md <- readRDS(file.path(base_dir, "results", "model_data.rds"))
N <- md$N_obs; tvec <- md$time
F_folds <- as.integer(Sys.getenv("N_FOLDS", "5"))
n_iter  <- as.integer(Sys.getenv("KFOLD_NITER", "22500"))
n_burn  <- as.integer(Sys.getenv("KFOLD_NBURN", "20000"))
out_dir <- file.path(base_dir, "results", "kfold")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## DETERMINISTIC stratified fold assignment (MUST match aggregate script)
make_folds <- function(tvec, F_folds) {
  set.seed(20260718)
  fid <- integer(length(tvec))
  for (tp in sort(unique(tvec))) {
    idx <- which(tvec == tp)
    fid[idx] <- sample(rep_len(1:F_folds, length(idx)))
  }
  fid
}
fold_id <- make_folds(tvec, F_folds)
saveRDS(list(fold_id = fold_id, F_folds = F_folds),
        file.path(out_dir, "fold_assignment.rds"))   # for provenance

tid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "-1"))
tasks <- if (tid >= 0) (tid + 1) else seq_len(2 * F_folds)   # no array => all (local)

for (r in tasks) {
  fold  <- ((r - 1) %/% 2) + 1
  model <- c("M1", "M0")[((r - 1) %% 2) + 1]
  tr <- which(fold_id != fold)
  cat(sprintf("\n=== Task %d: fold %d/%d, %s, %d train obs ===\n",
              r - 1, fold, F_folds, model, length(tr)))
  fit <- fit_multi_chain(md$Y_g[, tr, drop = FALSE], md$Y_cf[, tr, drop = FALSE],
                         md$patient[tr], tvec[tr],
                         n_chains = 4, seed = 42, K = 10,
                         n_iter = n_iter, n_burn = n_burn, thin = 10,
                         model = model, kappa_fixed = 1.0,
                         anneal = TRUE, T_anneal = 20.0, n_cool_buffer = 10000,
                         store_pointwise_ll = FALSE)
  f <- file.path(out_dir, sprintf("fit_f%d_%s.rds", fold, model))
  saveRDS(fit, f); cat("  saved", f, "\n")
}
cat("\nTask(s) complete. After all", 2 * F_folds,
    "finish, run slurm/run_kfold_aggregate.R\n")
