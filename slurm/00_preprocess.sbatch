#!/bin/bash
#SBATCH --job-name=krd_preprocess
#SBATCH --output=logs/krd/krd_preprocess_%j.out
#SBATCH --error=logs/krd/krd_preprocess_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=pi-jiyuan
#SBATCH --partition=caslake
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sbhandari52@uchicago.edu
#SBATCH --cpus-per-task=1

################################################################################
# 00_preprocess.sbatch — Preprocess KRd 5-hMC data
#
# Run ONCE before any KRd analysis jobs (04a, 04b, 04c).
# Produces results/model_data.rds (p=310 DE genes, n=26, T=5).
#
# USAGE:
#   cd /path/to/project
#   bash slurm/setup_logs.sh
#   sbatch slurm/00_preprocess.sbatch
#
# OUTPUT: results/model_data.rds
################################################################################

echo "=============================================="
echo "KRd Preprocessing"
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "=============================================="

# UChicago RCC (Midway3) toolchain. Rcpp requires the matching gcc.
module load gcc/13.2.0 R/4.4.2+gcc-13.2.0

cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs/krd results

Rscript R/02_preprocess.R

echo "=============================================="
echo "End time: $(date)"
echo "Exit status: $?"
echo "=============================================="
