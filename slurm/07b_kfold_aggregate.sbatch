#!/bin/bash
#SBATCH --job-name=krd_kfold_agg
#SBATCH --output=logs/krd/kfold_agg_%j.out
#SBATCH --error=logs/krd/kfold_agg_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=pi-jiyuan
#SBATCH --partition=caslake
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sbhandari52@uchicago.edu
#SBATCH --cpus-per-task=4

################################################################################
# 07b_kfold_aggregate.sbatch — held-out scoring + M1-vs-M0 elpd comparison.
# Run after all fold fits from 07_kfold_cv.sbatch exist (afterok dependency).
# Set N_FOLDS to match the fit phase; MC_DRAWS controls held-out MC precision.
################################################################################

echo "kfold aggregate on $(hostname) @ $(date)"
# UChicago RCC (Midway3) toolchain. Rcpp requires the matching gcc.
module load gcc/13.2.0 R/4.4.2+gcc-13.2.0
cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs/krd

export N_FOLDS=${N_FOLDS:-5}
export MC_DRAWS=${MC_DRAWS:-400}

Rscript slurm/run_kfold_aggregate.R

echo "End @ $(date), exit $?"
