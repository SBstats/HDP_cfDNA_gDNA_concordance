#!/bin/bash
#SBATCH --job-name=krd_lopo_agg
#SBATCH --output=logs/krd/lopo_agg_%j.out
#SBATCH --error=logs/krd/lopo_agg_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=pi-jiyuan
#SBATCH --partition=caslake
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sbhandari52@uchicago.edu
#SBATCH --cpus-per-task=1

################################################################################
# 09b_lopo_aggregate.sbatch — Tier 4 scoring/aggregation phase.
#
# Run AFTER all fits from 09_lopo_cv.sbatch exist in results/lopo/. Computes
# held-out log predictive density (MC over the trained weight prior, using
# theta_postmean), the total and SOURCE-SPECIFIC delta-elpd (M1 vs M0), a
# per-fold convergence gate, and a per-patient breakdown.
#   Outputs: results/lopo_cv.csv, lopo_pointwise.csv, lopo_per_patient.csv,
#            lopo_convergence_gate.csv
# Env: MC_DRAWS (default 400).
################################################################################

echo "lopo aggregate on $(hostname) @ $(date)"
# UChicago RCC (Midway3) toolchain. Rcpp requires the matching gcc.
module load gcc/13.2.0 R/4.4.2+gcc-13.2.0
cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs/krd

Rscript slurm/run_lopo_aggregate.R

echo "End @ $(date), exit $?"
