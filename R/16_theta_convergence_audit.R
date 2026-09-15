#!/usr/bin/env Rscript
###############################################################################
## R/16_theta_convergence_audit.R
## Formal between-chain convergence audit of the shared subclonal signatures
## theta, using the thinned theta trace produced by run_krd_theta_trace.R.
## Complements the Tier 2 local between-chain agreement check (R/15) with a
## per-gene split-Rhat (vectorised, rank-free basic Gelman-Rubin).
##
##   Rscript R/16_theta_convergence_audit.R M1   # or M0  (default M1)
## Requires: results/krd_<model>_fit_thetatrace.rds
## Outputs:  results/tier2/theta_rhat_<model>.csv
###############################################################################

args  <- commandArgs(trailingOnly = TRUE)
model <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "M1"
base_dir <- getwd()
## The label-switch and Assumption-A1 diagnostics reported at the foot of this
## script live in R/01_lib_core.R. Without this source() they are undefined,
## and because they were called inside tryCatch(error = NULL) the failure was
## silent: the section header printed with nothing under it and the output CSV
## was never written. Load explicitly.
suppressWarnings(suppressMessages(source(file.path(base_dir, "R", "01_lib_core.R"))))
stopifnot(exists("theta_label_switch_rate"), exists("theta_separation_check"))
res_dir  <- file.path(base_dir, "results")
out_dir  <- file.path(res_dir, "tier2"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fp <- file.path(res_dir, sprintf("krd_%s_fit_thetatrace.rds", model))
if (!file.exists(fp)) stop("Missing ", fp, " — run: sbatch slurm/08_theta_trace_refit.sbatch ", model)
fit <- readRDS(fp)
if (is.null(fit$chains[[1]]$theta_bar_trace))
  stop("Fit has no theta_bar_trace; re-fit with store_theta_trace=TRUE (run_krd_theta_trace.R).")

## Audited on the POPULATION signature trace theta_bar_trace [n_tr x p x K].
##
## An earlier version used theta_trace and read its 4th dimension as a
## TIMEPOINT; that dimension is now the SUBJECT. More importantly, the
## population signature is the quantity whose convergence matters: the
## subject-level signatures are weakly informed by design and not consistently
## estimable (writeup Remark rem:not_claimed), so their split-Rhat is dominated
## by prior noise rather than being a convergence signal.
##
## theta_label_switch_rate() is the companion diagnostic and is reported below:
## with subject-indexed signatures a chain can settle on different component
## orderings for different subjects, which no per-parameter Rhat would reveal.
nc  <- length(fit$chains)
dd  <- dim(fit$chains[[1]]$theta_bar_trace)   # [n_tr x p x K]
ntr <- dd[1]; p <- dd[2]; K <- dd[3]
cat(sprintf("model=%s  chains=%d  theta_bar draws/chain=%d  p=%d  K=%d\n", model, nc, ntr, p, K))

## vectorised basic split-Rhat across genes: A is [n_draws x n_chains x p]
split_rhat_vec <- function(A) {
  n <- dim(A)[1]; m <- dim(A)[2]; h <- floor(n / 2)
  halves <- array(0, dim = c(h, 2 * m, dim(A)[3]))
  halves[, 1:m, ]           <- A[1:h, , , drop = FALSE]
  halves[, (m + 1):(2 * m), ] <- A[(n - h + 1):n, , , drop = FALSE]
  cmean <- apply(halves, c(2, 3), mean)
  cvar  <- apply(halves, c(2, 3), var)
  B <- h * apply(cmean, 2, var)
  W <- apply(cvar, 2, mean)
  ok <- W > 1e-12
  rh <- rep(NA_real_, length(W))
  rh[ok] <- sqrt(((h - 1) / h * W[ok] + B[ok] / h) / W[ok])
  rh
}

Kocc <- min(3, K)
rows <- list()
for (k in 1:Kocc) {
  A <- array(0, dim = c(ntr, nc, p))
  for (ci in 1:nc) A[, ci, ] <- fit$chains[[ci]]$theta_bar_trace[, , k]
  rh <- split_rhat_vec(A)
  rows[[length(rows) + 1]] <- data.frame(
    model = model, comp = k, n_genes = sum(is.finite(rh)),
    rhat_median = round(median(rh, na.rm = TRUE), 4),
    rhat_p95    = round(quantile(rh, 0.95, na.rm = TRUE), 4),
    rhat_max    = round(max(rh, na.rm = TRUE), 4),
    pct_gt_1.01 = round(100 * mean(rh > 1.01, na.rm = TRUE), 2),
    pct_gt_1.1  = round(100 * mean(rh > 1.10, na.rm = TRUE), 2))
}
tab <- do.call(rbind, rows)
print(tab, row.names = FALSE)
write.csv(tab, file.path(out_dir, sprintf("theta_rhat_%s.csv", model)), row.names = FALSE)
cat(sprintf("\nSaved results/tier2/theta_rhat_%s.csv\n", model))
cat(sprintf("VERDICT: occupied-component theta_bar %s converged (max Rhat = %.3f over %d genes x comps).\n",
            if (max(tab$rhat_max[tab$comp <= 2], na.rm = TRUE) < 1.1) "APPEARS" else "MAY NOT HAVE",
            max(tab$rhat_max[tab$comp <= 2], na.rm = TRUE), p))

## ---------------------------------------------------------------------------
## Companion diagnostics required by the subject-indexed specification.
##
## A per-parameter Rhat cannot detect a chain that has settled on DIFFERENT
## component orderings for different subjects: every marginal looks healthy
## while the component label has lost its cross-subject meaning, invalidating
## cohort-level statements about pi or theta_bar. We therefore also report the
## within-subject label-switch rate and the empirical Assumption A1 check.
## ---------------------------------------------------------------------------
cat("\n===== within-subject label switching and Assumption A1 =====\n")
ls_rows <- list()
for (ci in seq_len(nc)) {
  ch <- fit$chains[[ci]]
  ## Capture the error MESSAGE rather than discarding it. theta_separation_check
  ## deliberately errors when the fit carries no varsigma_k (a fit with no
  ## dispersion estimate cannot be checked against A1); swallowing that into NA
  ## turns a loud, actionable failure into a blank cell that reads as "nothing
  ## to report". A1_status distinguishes the three outcomes explicitly.
  ls_err <- NULL; sep_err <- NULL
  ls_ <- tryCatch(theta_label_switch_rate(ch),
                  error = function(e) { ls_err <<- conditionMessage(e); NULL })
  sep <- tryCatch(theta_separation_check(ch),
                  error = function(e) { sep_err <<- conditionMessage(e); NULL })
  if (!is.null(ls_err))
    cat(sprintf("  chain %d: label-switch diagnostic ERRORED: %s\n", ci, ls_err))
  if (!is.null(sep_err))
    cat(sprintf("  chain %d: Assumption A1 check ERRORED: %s\n", ci, sep_err))

  a1_status <- if (!is.null(sep_err)) "ERROR"
               else if (is.null(sep)) "UNAVAILABLE"
               else if (isTRUE(attr(sep, "all_satisfied"))) "OK"
               else "FAIL"

  if (!is.null(ls_)) {
    cat(sprintf("  chain %d: label-switch rate = %.3f (%d/%d subjects) | A1 = %s\n",
                ci, ls_$rate, length(ls_$switched), ls_$n_subj, a1_status))
    ls_rows[[length(ls_rows) + 1L]] <- data.frame(
      model = model, chain = ci,
      label_switch_rate = ls_$rate,
      n_switched = length(ls_$switched), n_subj = ls_$n_subj,
      A1_status = a1_status,
      A1_m_req = if (!is.null(sep) && !is.null(sep$m_req)) sep$m_req[1] else NA_integer_,
      A1_satisfied = if (is.null(sep)) NA else isTRUE(attr(sep, "all_satisfied")))
  } else {
    ## Even with no label-switch result, record that the A1 check was attempted,
    ## so an errored diagnostic is visible in the output rather than absent.
    ls_rows[[length(ls_rows) + 1L]] <- data.frame(
      model = model, chain = ci,
      label_switch_rate = NA_real_, n_switched = NA_integer_, n_subj = NA_integer_,
      A1_status = a1_status, A1_m_req = NA_integer_, A1_satisfied = NA)
  }
}
if (length(ls_rows)) {
  ls_tab <- do.call(rbind, ls_rows)
  write.csv(ls_tab, file.path(out_dir, sprintf("theta_label_switch_%s.csv", model)),
            row.names = FALSE)
  cat(sprintf("Saved results/tier2/theta_label_switch_%s.csv\n", model))
  if (any(ls_tab$label_switch_rate > 0, na.rm = TRUE))
    cat("WARNING: nonzero label-switch rate -- component labels are NOT\n",
        "         cohort-consistent; interpret the decomposition per subject.\n")
  if (any(ls_tab$A1_status == "FAIL", na.rm = TRUE))
    cat("WARNING: Assumption A1 separation condition FAILS in at least one chain\n",
        "         (between-subject dispersion exceeds component separation).\n")
  ## A diagnostic that could not run is NOT a passing diagnostic. Warn
  ## separately so an unavailable check is never mistaken for a clean one.
  if (any(ls_tab$A1_status %in% c("ERROR", "UNAVAILABLE"), na.rm = TRUE))
    cat("WARNING: Assumption A1 check did not run for at least one chain\n",
        "         (status ERROR/UNAVAILABLE). This is NOT a pass -- the fit may\n",
        "         carry no varsigma_k. Refit with the current sampler.\n")
  if (all(is.na(ls_tab$label_switch_rate)))
    cat("WARNING: label-switch rate unavailable for every chain.\n")
}
