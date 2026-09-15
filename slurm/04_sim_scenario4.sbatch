#!/bin/bash
#SBATCH --job-name=sim_S4
#SBATCH --output=logs/scenario_4/sim_S4_%a_%A.out
#SBATCH --error=logs/scenario_4/sim_S4_%a_%A.err
#SBATCH --array=0-399
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=pi-jiyuan
#SBATCH --partition=caslake
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sbhandari52@uchicago.edu
#SBATCH --cpus-per-task=8


################################################################################
# 04_sim_scenario4.sbatch  --- Simulation Scenario 4 (misspecification robustness)
#
# ONE array task == ONE (config, replicate) pair (fits BOTH M1 and M0 once).
# Flattened from the old "loop all configs inside one task" layout because the
# hardened Scenario-4 sampler (8 chains x 30000 iters) makes a single fit-pair
# ~8 h; four serial configs would sit dangerously close to Midway's 36 h wall
# clock. One task per (config, rep) keeps every task to a single fit-pair (~8 h)
# and is trivially parallel.
#
#   array index a in [0, N_CONFIGS*R - 1]
#     CONFIG = a / R + 1        (1-based; R = reps per config = 100)
#     REP    = a % R + 1
#
# If N_CONFIGS or R changes, update --array=0-(N_CONFIGS*R - 1) accordingly.
# Current grid: {t3, corr, hk, all} misspecification departures = 4 configs;
# R = 100 => 400 tasks (0-399). Uses the hardened sampler (8 chains) =>
# cpus-per-task=8.
#
# USAGE:   cd /path/to/project ; sbatch slurm/04_sim_scenario4.sbatch
# AFTER:   sbatch slurm/05_sim_aggregate_results.sbatch ; sbatch slurm/05b_sim_figures.sbatch
################################################################################

echo "=============================================="
echo "SLURM Job ID: $SLURM_JOB_ID  Array Task: $SLURM_ARRAY_TASK_ID"
echo "Host: $(hostname)  Start: $(date)"
echo "=============================================="

# UChicago RCC (Midway3) toolchain. Rcpp requires the matching gcc.
module load gcc/13.2.0 R/4.4.2+gcc-13.2.0

cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs/scenario_4

SCEN=4
R_PER_CONFIG=100
N_CONFIGS=$(Rscript -e 'source("R/simulation/lib_simulation.R"); cat(nrow(get_scenario_grid(4,"production")))')

A=${SLURM_ARRAY_TASK_ID}
CONFIG=$(( A / R_PER_CONFIG + 1 ))
REP=$(( A % R_PER_CONFIG + 1 ))

if [ "${CONFIG}" -gt "${N_CONFIGS}" ]; then
  echo "Array index ${A} -> config ${CONFIG} exceeds N_CONFIGS=${N_CONFIGS}; nothing to do."
  exit 0
fi

echo "Scenario ${SCEN}: config ${CONFIG}/${N_CONFIGS}, replicate ${REP}"
Rscript slurm/run_one_sim.R ${SCEN} ${CONFIG} ${REP} production

echo "=============================================="
echo "End: $(date)  Exit: $?"
echo "=============================================="
