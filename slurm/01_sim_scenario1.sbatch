#!/bin/bash
#SBATCH --job-name=sim_S1
#SBATCH --output=logs/scenario_1/sim_S1_%a_%A.out
#SBATCH --error=logs/scenario_1/sim_S1_%a_%A.err
#SBATCH --array=0-1499
#SBATCH --time=36:00:00                  # caslake max is 36h; if SLURM rejects this as
                                         # over-limit, drop to 35:59:00 or split
                                         # the run (see README, HPC notes).
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=pi-jiyuan
#SBATCH --partition=caslake
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sbhandari52@uchicago.edu
#SBATCH --cpus-per-task=4


################################################################################
# 01_sim_scenario1.sbatch  --- Simulation Scenario 1 (signature estimation)
#
# Tests near-minimax contraction rate: normalized risk = MSE / (s_k * log(p/s_k))
# should be approximately flat in p (within each sparsity stratum), validating
# the horseshoe-regularized estimator. Also benchmarks against unregularized,
# lasso, and spike-and-slab competitors, and computes 95% CI coverage.
#
# ONE array task == ONE (config, replicate) pair (fits M1 + M0 once).
#
#   array index a in [0, N_CONFIGS*R - 1]
#     CONFIG = a / R + 1        (1-based; R = reps per config = 100)
#     REP    = a % R + 1
#
# If N_CONFIGS or R changes, update --array=0-(N_CONFIGS*R - 1) accordingly.
# Current grid: p in {150,300,500,1000,2000} x s_k/p in {0.01,0.05,0.10} = 15 configs;
# R = 100 => 1500 tasks (0-1499). Uses the PRODUCTION sampler (4 chains) => cpus-per-task=4.
#
# USAGE:   cd /path/to/project ; sbatch slurm/01_sim_scenario1.sbatch
# AFTER:   sbatch slurm/05_sim_aggregate_results.sbatch ; sbatch slurm/05b_sim_figures.sbatch
################################################################################

echo "=============================================="
echo "SLURM Job ID: $SLURM_JOB_ID  Array Task: $SLURM_ARRAY_TASK_ID"
echo "Host: $(hostname)  Start: $(date)"
echo "=============================================="

# UChicago RCC (Midway3) toolchain. Rcpp requires the matching gcc.
module load gcc/13.2.0 R/4.4.2+gcc-13.2.0

cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs/scenario_1

SCEN=1
R_PER_CONFIG=100
N_CONFIGS=$(Rscript -e 'source("R/simulation/lib_simulation.R"); cat(nrow(get_scenario_grid(1,"production")))')

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
