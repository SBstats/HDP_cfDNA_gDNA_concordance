#!/usr/bin/env Rscript
###############################################################################
## slurm/run_lopo_cv.R  --  Tier 4: leave-one-PATIENT-out predictive comparison
## FIT ONE (fold, model) PER ARRAY TASK.
##
## Improves on the timepoint-stratified k-fold (run_kfold_cv.R) by grouping folds
## by PATIENT: a held-out patient's cfDNA AND gDNA observations are entirely
## absent from training, so held-out prediction is a genuine "new patient" test
## with no within-patient leakage. With N_FOLDS = (#patients) this is true
## leave-one-patient-out; with fewer folds it is patient-grouped K-fold.
##
## Task grid (SLURM_ARRAY_TASK_ID = 0 .. 2*F-1):
##   fold  = (task %/% 2) + 1 ;  model = c("M1","M0")[task %% 2 + 1]
## Writes results/lopo/fit_f<fold>_<model>.rds. Then run run_lopo_aggregate.R.
##
## Env: N_FOLDS (default = number of patients => true LOPO).
##      LOPO_NITER / LOPO_NBURN (defaults 22500 / 20000).
## Fast local validation: N_FOLDS=2 LOPO_NITER=80 LOPO_NBURN=50 Rscript ...
###############################################################################

base_dir <- getwd()
source(file.path(base_dir, "R", "00_bootstrap.R"))
source(file.path(base_dir, "R", "01_lib_core.R"))
if (!exists(".rcpp_available") || !.rcpp_available)
  stop("FATAL: Rcpp required for production fits.")

res_dir <- Sys.getenv("LOPO_RES_DIR", file.path(base_dir, "results"))  # override for local testing
md <- readRDS(file.path(res_dir, "model_data.rds"))
patient <- md$patient; tvec <- md$time
n_pat   <- length(unique(patient))
F_folds <- as.integer(Sys.getenv("N_FOLDS", as.character(n_pat)))
n_iter  <- as.integer(Sys.getenv("LOPO_NITER", "22500"))
n_burn  <- as.integer(Sys.getenv("LOPO_NBURN", "20000"))
out_dir <- file.path(res_dir, "lopo")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## DETERMINISTIC patient-grouped fold assignment (each patient wholly in one
## fold). MUST match run_lopo_aggregate.R -> persisted to fold_assignment.rds.
make_folds_patient <- function(patient, F_folds) {
  set.seed(20260724)
  ups <- sort(unique(patient))
  stopifnot(F_folds <= length(ups))
  pf <- sample(rep_len(1:F_folds, length(ups)))       # patient -> fold (balanced)
  names(pf) <- as.character(ups)
  as.integer(pf[as.character(patient)])               # fold per observation
}
fold_id <- make_folds_patient(patient, F_folds)
saveRDS(list(fold_id = fold_id, F_folds = F_folds, by = "patient"),
        file.path(out_dir, "fold_assignment.rds"))

tid   <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "-1"))
tasks <- if (tid >= 0) (tid + 1) else seq_len(2 * F_folds)   # no array => all (local)

for (r in tasks) {
  fold  <- ((r - 1) %/% 2) + 1
  model <- c("M1", "M0")[((r - 1) %% 2) + 1]
  tr <- which(fold_id != fold)
  cat(sprintf("\n=== Task %d: fold %d/%d, %s, %d train obs (%d held-out patient(s)) ===\n",
              r - 1, fold, F_folds, model, length(tr),
              length(unique(patient[fold_id == fold]))))
  fit <- fit_multi_chain(md$Y_g[, tr, drop = FALSE], md$Y_cf[, tr, drop = FALSE],
                         patient[tr], tvec[tr],
                         n_chains = 4, seed = 42, K = 10,
                         n_iter = n_iter, n_burn = n_burn, thin = 10,
                         model = model, kappa_fixed = 1.0,
                         anneal = TRUE, T_anneal = 20.0, n_cool_buffer = 10000,
                         store_pointwise_ll = FALSE)
  f <- file.path(out_dir, sprintf("fit_f%d_%s.rds", fold, model))
  saveRDS(fit, f); cat("  saved", f, "\n")
}
cat("\nTask(s) complete. After all", 2 * F_folds,
    "finish, run slurm/run_lopo_aggregate.R\n")
