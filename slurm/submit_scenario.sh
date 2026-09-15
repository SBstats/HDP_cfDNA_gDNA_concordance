#!/bin/bash
###############################################################################
## slurm/submit_scenario.sh
## Submit a SLURM array job for one scenario + config combination
##
## Usage: bash submit_scenario.sh <scenario> <config_id> <n_reps> [scale]
##
## Example: bash submit_scenario.sh 1 3 100 production
###############################################################################

set -euo pipefail

SCENARIO=${1:?Usage: submit_scenario.sh <scenario> <config_id> <n_reps> [scale]}
CONFIG=${2:?Missing config_id}
NREPS=${3:-100}
SCALE=${4:-production}

# Project directory (parent of slurm/)
PROJDIR="$(cd "$(dirname "$0")/.." && pwd)"
LOGDIR="${PROJDIR}/logs/scenario_${SCENARIO}/config_${CONFIG}"
mkdir -p "${LOGDIR}"

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=sim_S${SCENARIO}_C${CONFIG}
#SBATCH --array=1-${NREPS}
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=${LOGDIR}/rep_%a_%A.out
#SBATCH --error=${LOGDIR}/rep_%a_%A.err

module load R/4.3.1 2>/dev/null || true

cd "${PROJDIR}"
Rscript slurm/run_one_sim.R ${SCENARIO} ${CONFIG} \${SLURM_ARRAY_TASK_ID} ${SCALE}
EOF

echo "Submitted: Scenario ${SCENARIO}, Config ${CONFIG}, ${NREPS} reps"
