#!/bin/bash
#SBATCH --job-name=krd_thetatrace
#SBATCH --output=logs/krd/krd_thetatrace_%j.out
#SBATCH --error=logs/krd/krd_thetatrace_%j.err
#SBATCH --time=36:00:00                  # caslake max is 36h; if SLURM rejects this as
                                         # over-limit, drop to 35:59:00 or split
                                         # the run (see README, HPC notes).
#SBATCH --mem=64G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=pi-jiyuan
#SBATCH --partition=caslake
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sbhandari52@uchicago.edu
#SBATCH --cpus-per-task=4

################################################################################
# 08_theta_trace_refit.sbatch — Tier 2 confirmatory re-fit storing a thinned
# theta trace, so the shared subclonal signatures can be audited with formal
# between-chain split-Rhat (R/16_theta_convergence_audit.R).
#
# USAGE (model defaults to M1):
#   cd /path/to/project
#   mkdir -p logs/krd
#   sbatch slurm/08_theta_trace_refit.sbatch M1
#   sbatch slurm/08_theta_trace_refit.sbatch M0
#
# Mirrors 04a_krd_M1.sbatch. Mem bumped 32G -> 64G for the theta-trace arrays
# (50 draws x 11,837 x 10 x 5 doubles ~ 236 MB/chain).
# OUTPUT: results/krd_<model>_fit_thetatrace.rds  (~1 GB)
#
# Everything below the directives is cluster-agnostic.
################################################################################

echo "=============================================="
echo "KRd theta-trace re-fit | model=${1:-M1}"
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "Running on host: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK | Start: $(date)"
echo "=============================================="

# UChicago RCC (Midway3) toolchain. Rcpp requires the matching gcc.
module load gcc/13.2.0 R/4.4.2+gcc-13.2.0

cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs/krd

MODEL=${1:-M1}
Rscript slurm/run_krd_theta_trace.R "${MODEL}"

echo "=============================================="
echo "End time: $(date) | Exit status: $?"
echo "=============================================="
