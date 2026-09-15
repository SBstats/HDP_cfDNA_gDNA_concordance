#!/bin/bash
#SBATCH --job-name=krd_M0
#SBATCH --output=logs/krd/krd_M0_%j.out
#SBATCH --error=logs/krd/krd_M0_%j.err
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

################################################################################
# 04b_krd_M0.sbatch — Fit M0 (non-tracking model), 4 parallel chains
#
# USAGE:
#   cd /path/to/project
#   sbatch slurm/04b_krd_M0.sbatch
#
# Can be submitted in parallel with 04a_krd_M1.sbatch
#
# ESTIMATED TIME: ~2-4 hours (p=310 DE genes, K=10, T=5, 22500 iters, thin=10, Rcpp, annealing)
# OUTPUT: results/krd_M0_fit.rds
################################################################################

echo "=============================================="
echo "KRd M0 Analysis"
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "Running on host: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Start time: $(date)"
echo "=============================================="

# UChicago RCC (Midway3) toolchain. Rcpp requires the matching gcc.
module load gcc/13.2.0 R/4.4.2+gcc-13.2.0

cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs/krd

Rscript slurm/run_krd_M0.R

echo "=============================================="
echo "End time: $(date)"
echo "Exit status: $?"
echo "=============================================="
