#!/bin/bash
# Create all log directories the sbatch scripts write to.
# SLURM opens --output/--error paths at job-START (before the job body runs),
# so these dirs MUST exist BEFORE `sbatch`, or the job fails instantly with
# NonZeroExitCode and produces NO .out/.err file. .gitkeep placeholders exist
# locally but are often dropped by file-transfer tools, so run this once on the
# cluster after uploading, from the project root:
#
#   bash slurm/setup_logs.sh
#
mkdir -p logs/krd \
         logs/scenario_1 logs/scenario_2 logs/scenario_3 logs/scenario_4 logs/scenario_5 \
         logs/sim
echo "Log directories ready:"
ls -ld logs/krd logs/scenario_1 logs/scenario_2 logs/scenario_3 logs/scenario_4 logs/scenario_5 logs/sim
