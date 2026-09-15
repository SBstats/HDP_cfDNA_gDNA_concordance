#!/bin/bash
###############################################################################
## slurm/submit_all.sh
## Submit all simulation scenarios to SLURM using per-scenario sbatch files.
##
## Usage: bash submit_all.sh [scale]
##        scale = "production" (default) or "local"
##
## Each scenario has its own sbatch file with the correct account, QOS,
## memory, cpus-per-task, and array size for that scenario.
###############################################################################

set -euo pipefail

SCALE=${1:-production}
PROJDIR="$(cd "$(dirname "$0")/.." && pwd)"

echo "======================================"
echo "Submitting all simulation scenarios"
echo "Scale: ${SCALE}"
echo "Project: ${PROJDIR}"
echo "======================================"

cd "${PROJDIR}"

echo ""
echo "Scenario 1 (signature estimation): 15 configs x 100 reps = 1500 tasks"
sbatch slurm/01_sim_scenario1.sbatch

echo ""
echo "Scenario 2 (signature recovery vs n): 6 configs x 100 reps = 600 tasks"
sbatch slurm/02_sim_scenario2.sbatch

echo ""
echo "Scenario 3 (sample-size power curve): 7 configs x 100 reps = 700 tasks"
sbatch slurm/03_sim_scenario3.sbatch

echo ""
echo "Scenario 4 (misspecification robustness): 4 configs x 100 reps = 400 tasks"
sbatch slurm/04_sim_scenario4.sbatch

echo ""
echo "Scenario 5 (compartment-specific signatures): 4 configs x 100 reps = 400 tasks"
sbatch slurm/05_sim_scenario5.sbatch

echo ""
echo "All simulation jobs submitted."
