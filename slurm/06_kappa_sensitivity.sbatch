#!/bin/bash
#SBATCH --job-name=krd_kappa
#SBATCH --output=logs/krd/kappa_%A_%a.out
#SBATCH --error=logs/krd/kappa_%A_%a.err
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
#SBATCH --array=0-7

################################################################################
# 06_kappa_sensitivity.sbatch — Supplementary Table SX (kappa robustness)
#
# ONE (model, kappa) fit per array task (8 tasks) so each stays under the
# caslake 36h cap (each fit ~18-26h).
#   idx: 0=(k0.5,M1) 1=(k0.5,M0) 2=(k1,M1) 3=(k1,M0)
#        4=(k2,M1)   5=(k2,M0)   6=(k5,M1) 7=(k5,M0)
#
# USAGE:
#   mkdir -p logs/krd
#   AID=$(sbatch --parsable slurm/06_kappa_sensitivity.sbatch)
#   sbatch --dependency=afterok:$AID slurm/06b_kappa_aggregate.sbatch
################################################################################

echo "kappa task $SLURM_ARRAY_TASK_ID on $(hostname) @ $(date)"
# UChicago RCC (Midway3) toolchain. Rcpp requires the matching gcc.
module load gcc/13.2.0 R/4.4.2+gcc-13.2.0
cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs/krd

Rscript slurm/run_kappa_sensitivity.R

echo "End @ $(date), exit $?"
