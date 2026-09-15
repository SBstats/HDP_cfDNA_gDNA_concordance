#!/bin/bash
#SBATCH --job-name=sim_S5
#SBATCH --output=logs/scenario_5/sim_S5_%a_%A.out
#SBATCH --error=logs/scenario_5/sim_S5_%a_%A.err
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
# 05_sim_scenario5.sbatch  --- Simulation Scenario 5 (compartment-specific signatures)
#
# REDESIGNED: tests whether the model correctly recovers rho^theta (signature
# concordance) under graded inter-source divergence. Extended delta_theta grid
# {0, 1, 2, 4} spans rho_true from 1.00 down to ~0.30 -- well above the bias
# floor -- replacing the old {0, 0.1, 0.25, 0.5} grid which had rho_true
# variation of only 0.036, smaller than the estimator's ~0.04 bias.
#
# ONE array task == ONE (config, replicate) pair.
#
#   array index a in [0, N_CONFIGS*R - 1]
#     CONFIG = a / R + 1        (1-based; R = reps per config = 100)
#     REP    = a % R + 1
#
# Current grid: 4 delta_theta values = 4 configs; R = 100 => 400 tasks (0-399).
# Uses the hardened sampler (8 chains) => cpus-per-task=8.
#
# USAGE:   cd /path/to/project ; sbatch slurm/05_sim_scenario5.sbatch
# AFTER:   sbatch slurm/05_sim_aggregate_results.sbatch ; sbatch slurm/05b_sim_figures.sbatch
################################################################################

echo "=============================================="
echo "SLURM Job ID: $SLURM_JOB_ID  Array Task: $SLURM_ARRAY_TASK_ID"
echo "Host: $(hostname)  Start: $(date)"
echo "=============================================="

# UChicago RCC (Midway3) toolchain. Rcpp requires the matching gcc.
module load gcc/13.2.0 R/4.4.2+gcc-13.2.0

cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs/scenario_5

SCEN=5
R_PER_CONFIG=100
N_CONFIGS=$(Rscript -e 'source("R/simulation/lib_simulation.R"); cat(nrow(get_scenario_grid(5,"production")))')

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
