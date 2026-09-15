#!/bin/bash
#SBATCH --job-name=sim_S3
#SBATCH --output=logs/scenario_3/sim_S3_%a_%A.out
#SBATCH --error=logs/scenario_3/sim_S3_%a_%A.err
#SBATCH --array=0-699
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=pi-jiyuan
#SBATCH --partition=caslake
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sbhandari52@uchicago.edu
#SBATCH --cpus-per-task=8


################################################################################
# 03_sim_scenario3.sbatch  --- Simulation Scenario 3 (sample-size power curve)
#
# REDESIGNED: tests whether WAIC detection rate and increment detection rate
# increase monotonically with cohort size n under M1, and are controlled at
# nominal 5% under M0. kappa*=1 (diffuse compositions) throughout.
#
# ONE array task == ONE (config, replicate) pair (fits BOTH M1 and M0 once).
#
#   array index a in [0, N_CONFIGS*R - 1]
#     CONFIG = a / R + 1        (1-based; R = reps per config = 100)
#     REP    = a % R + 1
#
# If N_CONFIGS or R changes, update --array=0-(N_CONFIGS*R - 1) accordingly.
# Current grid: M0-null (n=26) + M1 x {n=5,10,20,26,50,100} = 7 configs; R=100
# => 700 tasks (0-699). p=300 (avoids WAIC DOF imbalance at p>>N_obs).
# Uses the hardened sampler (8 chains) => cpus-per-task=8.
#
# USAGE:   cd /path/to/project ; sbatch slurm/03_sim_scenario3.sbatch
# AFTER:   sbatch slurm/05_sim_aggregate_results.sbatch ; sbatch slurm/05b_sim_figures.sbatch
################################################################################

echo "=============================================="
echo "SLURM Job ID: $SLURM_JOB_ID  Array Task: $SLURM_ARRAY_TASK_ID"
echo "Host: $(hostname)  Start: $(date)"
echo "=============================================="

# UChicago RCC (Midway3) toolchain. Rcpp requires the matching gcc.
module load gcc/13.2.0 R/4.4.2+gcc-13.2.0

cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs/scenario_3

SCEN=3
R_PER_CONFIG=100
N_CONFIGS=$(Rscript -e 'source("R/simulation/lib_simulation.R"); cat(nrow(get_scenario_grid(3,"production")))')

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
