#!/bin/bash
#SBATCH --job-name=sim_figs
#SBATCH --output=logs/sim/sim_figs_%j.out
#SBATCH --error=logs/sim/sim_figs_%j.err
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
# 05b_sim_figures.sbatch
# Build simulation figures/tables from results/sim_aggregated/*.csv.
# Run AFTER 05_aggregate_results.sbatch.
#   figures/sim_scenario1_contraction.pdf, sim_scenario2a_robustness.pdf,
#   sim_scenario2b_pi_consistency.pdf, sim_scenario3_roc.pdf
#   tables/sim_scenario{1,2,3,4}_*.tex
################################################################################

echo "Start: $(date)"
# UChicago RCC (Midway3) toolchain. Rcpp requires the matching gcc.
module load gcc/13.2.0 R/4.4.2+gcc-13.2.0
cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs/sim
Rscript R/25_sim_figures_tables.R
echo "End: $(date)  Exit: $?"
