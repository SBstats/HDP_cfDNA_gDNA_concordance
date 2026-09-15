#!/bin/bash
#SBATCH --job-name=krd_kfold
#SBATCH --output=logs/krd/kfold_%A_%a.out
#SBATCH --error=logs/krd/kfold_%A_%a.err
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
#SBATCH --array=0-9

################################################################################
# 07_kfold_cv.sbatch — PSIS-free K-fold CV (M1 vs M0), FIT phase
#
# ONE (fold, model) fit per array task so each stays under the caslake 36h cap.
# Default N_FOLDS=5 -> 2*5 = 10 tasks (array 0-9). idx = 2*(fold-1)+{0:M1,1:M0}.
#
# VALIDATE FIRST (cheap) before the full run — the held-out-lpd marginalisation
# in the aggregate step is new code:
#   N_FOLDS=2 KFOLD_NITER=3000 KFOLD_NBURN=2000 \
#     sbatch --array=0-3 --time=04:00:00 slurm/07_kfold_cv.sbatch
#   sbatch slurm/07b_kfold_aggregate.sbatch      # check results/kfold_cv.csv
#
# FULL RUN:
#   mkdir -p logs/krd
#   AID=$(sbatch --parsable slurm/07_kfold_cv.sbatch)
#   sbatch --dependency=afterok:$AID slurm/07b_kfold_aggregate.sbatch
#
# To change fold count, set N_FOLDS and the array range together:
#   N_FOLDS=5 -> --array=0-9 ; N_FOLDS=10 -> --array=0-19
################################################################################

echo "kfold task $SLURM_ARRAY_TASK_ID on $(hostname) @ $(date)"
# UChicago RCC (Midway3) toolchain. Rcpp requires the matching gcc.
module load gcc/13.2.0 R/4.4.2+gcc-13.2.0
cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs/krd

export N_FOLDS=${N_FOLDS:-5}

Rscript slurm/run_kfold_cv.R

echo "End @ $(date), exit $?"
