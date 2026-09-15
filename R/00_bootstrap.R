#!/usr/bin/env Rscript
###############################################################################
## R/00_bootstrap.R
##
## Tiny bootstrap script that compiles the Rcpp inner loops at the TOP LEVEL
## of a fresh R session. Call this once at the start of any analysis script
## (interactive or batch) BEFORE sourcing R/01_lib_core.R.
##
## Background: on some platforms (notably macOS) `sourceCpp()` accumulates
## enough C-stack depth that, when invoked from inside `source()`, it
## triggers "C stack usage too close to the limit" errors. Running
## `sourceCpp()` here at the top level avoids that.
##
## Usage:
##   source("R/00_bootstrap.R")
##   source("R/01_lib_core.R")
##
## After bootstrap, .rcpp_available will be TRUE if the C++ code compiled,
## FALSE otherwise (and 01_lib_core.R will then use the pure-R fallback).
###############################################################################

suppressPackageStartupMessages({
  library(Rcpp)
})

.cpp_file <- file.path(getwd(), "R", "cpp", "sampler_core.cpp")
if (!file.exists(.cpp_file)) {
  message("R/cpp/sampler_core.cpp not found; pure-R fallback will be used.")
  .rcpp_available <- FALSE
} else {
  ok <- tryCatch({
    sourceCpp(.cpp_file)
    TRUE
  }, error = function(e) {
    message("Rcpp compilation failed: ", conditionMessage(e))
    FALSE
  })
  .rcpp_available <- ok
  if (ok) {
    message("[bootstrap] Rcpp compiled OK; .rcpp_available = TRUE")
  } else {
    message("[bootstrap] Rcpp failed; will use pure-R fallback. ",
            "Production fits should be run on a machine where Rcpp compiles.")
  }
}
rm(.cpp_file)
