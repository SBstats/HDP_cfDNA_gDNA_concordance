#!/bin/bash
#SBATCH --job-name=krd_lopo
#SBATCH --output=logs/krd/lopo_%A_%a.out
#SBATCH --error=logs/krd/lopo_%A_%a.err
#SBATCH --time=36:00:00                  # caslake max is 36h; if SLURM rejects this as
                                         # over-limit, drop to 35:59:00 or split
                                         # the run (see README, HPC notes).
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=pi-jiyuan
#SBATCH --partition=caslake
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sbhandari52@uchicago.edu
#SBATCH --cpus-per-task=4
#SBATCH --array=0-25

################################################################################
# 09_lopo_cv.sbatch — Tier 4: leave-one-PATIENT-out predictive comparison, FIT phase.
#
# ONE (fold, model) fit per array task. Patient-grouped folds (a held-out
# patient's cfDNA + gDNA are entirely out of training). idx = 2*(fold-1)+{0:M1,1:M0}.
#
# DEFAULT: N_FOLDS=13 (2 patients/fold) -> 2*13 = 26 tasks (array 0-25). This is
# a compute-balanced patient-grouped K-fold: each model trains on 24 patients,
# and all 84 observations are scored once. Precision of delta-elpd is set by the
# 84 held-out obs, so 13-fold is nearly as informative as true LOPO at half cost.
#
# TRUE leave-one-patient-out (26 folds, trains on 25 patients each):
#   N_FOLDS=26 sbatch --array=0-51 slurm/09_lopo_cv.sbatch     # 52 fits (~2x)
#
# VALIDATE FIRST (cheap) — held-out scoring in the aggregate is the new code:
#   N_FOLDS=2 LOPO_NITER=3000 LOPO_NBURN=2000 \
#     sbatch --array=0-3 --time=04:00:00 slurm/09_lopo_cv.sbatch
#   sbatch slurm/09b_lopo_aggregate.sbatch     # -> results/lopo_cv.csv
#
# FULL RUN (default 13-fold):
#   mkdir -p logs/krd
#   AID=$(sbatch --parsable slurm/09_lopo_cv.sbatch)
#   sbatch --dependency=afterok:$AID slurm/09b_lopo_aggregate.sbatch
#
# Keep N_FOLDS and --array in sync: N_FOLDS=F -> --array=0-(2F-1).
#
################################################################################

echo "lopo task $SLURM_ARRAY_TASK_ID on $(hostname) @ $(date)"
# UChicago RCC (Midway3) toolchain. Rcpp requires the matching gcc.
module load gcc/13.2.0 R/4.4.2+gcc-13.2.0
cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs/krd

export N_FOLDS=${N_FOLDS:-13}

Rscript slurm/run_lopo_cv.R

echo "End @ $(date), exit $?"
