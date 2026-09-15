#!/bin/bash
#SBATCH --job-name=krd_kappa_agg
#SBATCH --output=logs/krd/kappa_agg_%j.out
#SBATCH --error=logs/krd/kappa_agg_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=pi-jiyuan
#SBATCH --partition=caslake
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sbhandari52@uchicago.edu
#SBATCH --cpus-per-task=2

################################################################################
# 06b_kappa_aggregate.sbatch — build results/kappa_sensitivity.csv from the 8
# fits produced by 06_kappa_sensitivity.sbatch. Submit with an afterok
# dependency on the array job (see 06 header).
################################################################################

echo "kappa aggregate on $(hostname) @ $(date)"
# UChicago RCC (Midway3) toolchain. Rcpp requires the matching gcc.
module load gcc/13.2.0 R/4.4.2+gcc-13.2.0
cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs/krd

Rscript slurm/run_kappa_aggregate.R

echo "End @ $(date), exit $?"
