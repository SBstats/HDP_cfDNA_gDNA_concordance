#!/bin/bash
#SBATCH --job-name=sim_S2
#SBATCH --output=logs/scenario_2/sim_S2_%a_%A.out
#SBATCH --error=logs/scenario_2/sim_S2_%a_%A.err
#SBATCH --array=0-599
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=pi-jiyuan
#SBATCH --partition=caslake
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sbhandari52@uchicago.edu
#SBATCH --cpus-per-task=4


################################################################################
# 02_sim_scenario2.sbatch  --- Simulation Scenario 2 (signature recovery vs n)
#
# REDESIGNED: tests whether signature recovery rate and rho^theta improve
# monotonically with cohort size n. patient_pi=FALSE so a single shared
# cohort composition pi* exists as the estimand.
#
# ONE array task == ONE (config, replicate) pair.
#
#   array index a in [0, N_CONFIGS*R - 1]
#     CONFIG = a / R + 1        (1-based; R = reps per config = 100)
#     REP    = a % R + 1
#
# Grid: n in {5,10,20,50,100,200} x p=300 = 6 configs; R=100 => 600 tasks (0-599).
# Uses PRODUCTION sampler (4 chains) => cpus-per-task=4.
#
# USAGE:   cd /path/to/project ; sbatch slurm/02_sim_scenario2.sbatch
# AFTER:   sbatch slurm/05_sim_aggregate_results.sbatch ; sbatch slurm/05b_sim_figures.sbatch
################################################################################

echo "=============================================="
echo "SLURM Job ID: $SLURM_JOB_ID  Array Task: $SLURM_ARRAY_TASK_ID"
echo "Host: $(hostname)  Start: $(date)"
echo "=============================================="

# UChicago RCC (Midway3) toolchain. Rcpp requires the matching gcc.
module load gcc/13.2.0 R/4.4.2+gcc-13.2.0

cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs/scenario_2

SCEN=2
R_PER_CONFIG=100
N_CONFIGS=$(Rscript -e 'source("R/simulation/lib_simulation.R"); cat(nrow(get_scenario_grid(2,"production")))')

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
