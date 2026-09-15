#!/usr/bin/env Rscript
###############################################################################
## slurm/run_one_sim.R
## CLI wrapper for simulation runner. Called by SLURM array jobs.
##
## Usage: Rscript slurm/run_one_sim.R <scenario> <config_id> <rep_id> [scale]
##
## NOTE: Must be called from the project root directory (the SLURM scripts
## do `cd "${PROJDIR}"` before invoking this).
###############################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript run_one_sim.R <scenario> <config_id> <rep_id> [scale]")
}

scenario  <- as.integer(args[1])
config_id <- as.integer(args[2])
rep_id    <- as.integer(args[3])
scale     <- if (length(args) >= 4) args[4] else "production"

base_dir <- getwd()

source(file.path(base_dir, "R", "simulation", "lib_simulation.R"))

run_one_sim(scenario = scenario,
            config_id = config_id,
            rep_id = rep_id,
            scale = scale,
            base_dir = base_dir)
