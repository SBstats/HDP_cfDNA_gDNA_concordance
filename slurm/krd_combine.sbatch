#!/bin/bash
#SBATCH --job-name=krd_combine
#SBATCH --output=logs/krd/krd_combine_%j.out
#SBATCH --error=logs/krd/krd_combine_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=pi-jiyuan
#SBATCH --partition=caslake
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sbhandari52@uchicago.edu
#SBATCH --cpus-per-task=1

################################################################################
# 04c_krd_combine.sbatch — Combine M1+M0 fits, compute diagnostics & comparison
#
# Run AFTER both 04a_krd_M1.sbatch and 04b_krd_M0.sbatch have completed.
#
# USAGE:
#   cd /path/to/project
#   sbatch slurm/04c_krd_combine.sbatch
#
# OUTPUTS (in results/):
#   Core (from run_krd_combine.R itself):
#     - krd_posterior_summary.csv
#     - krd_model_comparison.csv
#     - krd_chain_health.csv
#     - krd_wN_by_patient.csv
#     - krd_convergence_M1.csv, krd_convergence_M0.csv
#     - krd_rho_theta.csv          (rho^theta_{kt}: signature concordance per component/timepoint)
#     - sessionInfo_krd_production.txt
#   Derived (via sourcing R/22, R/24, R/25 at the end):
#     - krd_L1_concordance_posterior.csv  (R/22; feeds Figure 3b)
#     - krd_prediction_correlation.csv    (R/22; feeds Figure 3c)
#     - krd_loo_M1.rds, krd_loo_M0.rds    (R/24)
#     - krd_pareto_k_M1.csv, krd_pareto_k_M0.csv (R/24)
#     - krd_loo_reliability_summary.csv   (R/24)
#     - krd_sigma_k_by_chain_M1.csv, krd_sigma_k_by_chain_M0.csv (R/25)
#     - krd_sigma_k_locking_verdict.csv   (R/25)
#   PDF figures (also produced by the derived scripts):
#     - nature_manuscript/figures/fig_L1_concordance_forest.pdf, fig3_model_concordance_heatmap.pdf, fig_prediction_correlation_heatmap.pdf  (R/22)
#     - figures/krd_sigma_k_boxplot_by_chain.pdf  (R/25)
#     - figures/krd_traceplots_*.pdf  (this script: relabeling traceplots)
################################################################################

echo "=============================================="
echo "KRd Combine Results"
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "=============================================="

# UChicago RCC (Midway3) toolchain. Rcpp requires the matching gcc.
module load gcc/13.2.0 R/4.4.2+gcc-13.2.0

cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs/krd

Rscript slurm/run_krd_combine.R

echo "=============================================="
echo "End time: $(date)"
echo "Exit status: $?"
echo "=============================================="
