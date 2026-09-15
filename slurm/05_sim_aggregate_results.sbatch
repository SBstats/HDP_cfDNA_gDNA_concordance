#!/bin/bash
#SBATCH --job-name=aggregate
#SBATCH --output=logs/sim/aggregate_%j.out
#SBATCH --error=logs/sim/aggregate_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=pi-jiyuan
#SBATCH --partition=caslake
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sbhandari52@uchicago.edu
#SBATCH --cpus-per-task=1


################################################################################
# 05_aggregate_results.sbatch
#
# SLURM batch script for aggregating simulation results.
# Run AFTER ALL simulation array jobs (01-05) have completed.
#
# USAGE:
#   cd /path/to/project
#   sbatch slurm/05_sim_aggregate_results.sbatch
#
# OUTPUTS (in results/sim_aggregated/):
#   - scenario_1_summary.csv  (signature estimation, 15 configs)
#   - scenario_2_summary.csv  (signature recovery vs n, 6 configs)
#   - scenario_3_summary.csv  (sample-size power curve, 7 configs)
#   - scenario_4_summary.csv  (misspecification robustness, 4 configs)
#   - scenario_5_summary.csv  (compartment-specific signatures, 4 configs)
#
################################################################################

echo "=============================================="
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "Running on host: $(hostname)"
echo "Start time: $(date)"
echo "=============================================="

# UChicago RCC (Midway3) toolchain. Rcpp requires the matching gcc.
module load gcc/13.2.0 R/4.4.2+gcc-13.2.0

cd "${SLURM_SUBMIT_DIR}"

mkdir -p logs/sim

Rscript slurm/aggregate_results.R

echo "=============================================="
echo "End time: $(date)"
echo "Exit status: $?"
echo "=============================================="
