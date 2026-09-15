###############################################################################
## 23_refresh_manuscript_outputs.R
##
## Regenerate manuscript-ready figures and tables from the post-remediation
## .rds outputs in results/. Produces both LaTeX tables (tables/) and PDF
## figures (figures/) for both the methods paper and the applied paper.
##
## Inputs (expected in results/):
##   krd_M1_fit.rds, krd_M0_fit.rds
##   krd_chain_health.csv, krd_model_comparison.csv,
##   krd_posterior_summary.csv, krd_convergence_M1.csv, krd_convergence_M0.csv,
##   krd_omega0_by_patient.csv
##
## Outputs:
##   figures/krd_chain_health.pdf           — LL per chain with 10k-gap band
##   figures/krd_loglik_traces.pdf          — M1 and M0 LL traces, 4 chains
##   figures/krd_traceplots_remediated.pdf  — K+, sigma_0, omega_0, omega traces
##   tables/krd_chain_health.tex            — chain-level LL gap table
##   tables/krd_model_comparison.tex        — dWAIC/dLOO with clean/all split
##   tables/krd_posterior_summary.tex       — K+, sigma_0, omega_0_mean
##   tables/krd_convergence_summary.tex     — per-parameter Rhat + ESS summary
##   tables/krd_omega0_by_patient.tex       — per-patient contamination table
###############################################################################

base_dir <- getwd()
results_dir <- file.path(base_dir, "results")
figures_dir <- file.path(base_dir, "figures")
tables_dir  <- file.path(base_dir, "tables")
dir.create(figures_dir, showWarnings = FALSE)
dir.create(tables_dir,  showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
})

# ---------------------------------------------------------------------------
# Load fresh results
# ---------------------------------------------------------------------------
multi_M1 <- readRDS(file.path(results_dir, "krd_M1_fit.rds"))
multi_M0 <- readRDS(file.path(results_dir, "krd_M0_fit.rds"))

ch_health   <- read.csv(file.path(results_dir, "krd_chain_health.csv"))
cmp         <- read.csv(file.path(results_dir, "krd_model_comparison.csv"))
post_summ   <- read.csv(file.path(results_dir, "krd_posterior_summary.csv"))
conv_M1     <- read.csv(file.path(results_dir, "krd_convergence_M1.csv"))
conv_M0     <- read.csv(file.path(results_dir, "krd_convergence_M0.csv"))
# Per-patient contamination AGGREGATE. Note omega_0 is estimated per
# OBSERVATION; this file averages within patient, which discards the temporal
# variation the respecification introduced. Use
# krd_omega0_by_observation.csv for anything substantive.
wN_patient  <- read.csv(file.path(results_dir, "krd_omega0_by_patient.csv"))
## Fail loudly on a schema drift rather than emitting character(0) into the
## table (which surfaces as an opaque "differing number of rows" error).
stopifnot(all(c("patient", "post_mean", "ci_lo", "ci_hi") %in% names(wN_patient)))

# ---------------------------------------------------------------------------
# Figure: chain health (mean LL per chain, 10k-gap threshold annotated)
# ---------------------------------------------------------------------------
ch_health$model <- factor(ch_health$model, levels = c("M1", "M0"))
ch_health$label <- ifelse(ch_health$trapped, "trapped", "retained")

p_health <- ggplot(ch_health,
                   aes(x = factor(chain), y = mean_loglik,
                       fill = label)) +
  geom_col(width = 0.7) +
  geom_hline(data = do.call(rbind, lapply(split(ch_health, ch_health$model),
                                          function(d) data.frame(
                                            model = d$model[1],
                                            y = max(d$mean_loglik) - 10000))),
             aes(yintercept = y), linetype = "dashed", colour = "grey40") +
  facet_wrap(~ model, scales = "free_y") +
  scale_fill_manual(values = c(retained = "#4C9F70", trapped = "#D1495B")) +
  labs(x = "Chain", y = "Mean post-burn-in log-likelihood",
       fill = NULL,
       title = "Chain health check (dashed line = best chain − 10,000)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(file.path(figures_dir, "krd_chain_health.pdf"),
       p_health, width = 7, height = 4)

# ---------------------------------------------------------------------------
# Figure: LL traces (M1 and M0), 4 chains overlaid
# ---------------------------------------------------------------------------
extract_ll_traces <- function(multi_fit, model_label) {
  do.call(rbind, lapply(seq_along(multi_fit$chains), function(i) {
    data.frame(
      iter   = seq_along(multi_fit$chains[[i]]$samples$loglik),
      loglik = multi_fit$chains[[i]]$samples$loglik,
      chain  = factor(i),
      model  = model_label
    )
  }))
}
ll_df <- rbind(extract_ll_traces(multi_M1, "M1"),
               extract_ll_traces(multi_M0, "M0"))
ll_df$model <- factor(ll_df$model, levels = c("M1", "M0"))

p_ll <- ggplot(ll_df, aes(x = iter, y = loglik, colour = chain)) +
  geom_line(linewidth = 0.3) +
  facet_wrap(~ model, scales = "free_y", ncol = 1) +
  labs(x = "Iteration (post burn-in, thinned)",
       y = "Log-likelihood",
       colour = "Chain",
       title = "Post-remediation: log-likelihood traces across 4 chains") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(file.path(figures_dir, "krd_loglik_traces.pdf"),
       p_ll, width = 8, height = 6)

# ---------------------------------------------------------------------------
# Figure: remediated traceplots (K+, sigma_0, omega_0 obs 1, log-lik)
# ---------------------------------------------------------------------------
extract_scalar <- function(multi_fit, param, model_label) {
  do.call(rbind, lapply(seq_along(multi_fit$chains), function(i) {
    val <- multi_fit$chains[[i]]$samples[[param]]
    if (is.matrix(val)) val <- val[, 1]
    data.frame(iter = seq_along(val), value = val,
               chain = factor(i), model = model_label, param = param)
  }))
}

trace_panels <- rbind(
  extract_scalar(multi_M1, "K_plus",  "M1"),
  extract_scalar(multi_M1, "sigma_0", "M1"),
  extract_scalar(multi_M1, "omega_0",      "M1"),
  extract_scalar(multi_M1, "loglik",  "M1")
)
trace_panels$param <- factor(trace_panels$param,
                             levels = c("K_plus","loglik","sigma_0","omega_0"),
                             labels = c("K+", "log-likelihood",
                                        "sigma_0", "omega_0 (obs 1)"))

p_trace <- ggplot(trace_panels, aes(x = iter, y = value, colour = chain)) +
  geom_line(linewidth = 0.3) +
  facet_wrap(~ param, scales = "free_y", ncol = 2) +
  labs(x = "Iteration (post burn-in, thinned)", y = NULL,
       colour = "Chain",
       title = "M1 trace plots after remediation (4 chains, 250 samples each)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(file.path(figures_dir, "krd_traceplots_remediated.pdf"),
       p_trace, width = 8, height = 6)

# ---------------------------------------------------------------------------
# LaTeX table: chain health
# ---------------------------------------------------------------------------
esc <- function(x) x  # placeholder for LaTeX escaping if needed

chain_health_tex <- c(
  "\\begin{tabular}{llrrl}",
  "\\toprule",
  "Model & Chain & Mean log-lik & Gap from best & Status \\\\",
  "\\midrule",
  paste0(
    ch_health$model, " & ", ch_health$chain, " & ",
    formatC(ch_health$mean_loglik, format = "f", big.mark = ",", digits = 1),
    " & ",
    formatC(ch_health$gap_from_best, format = "f", big.mark = ",", digits = 1),
    " & ", ifelse(ch_health$trapped, "trapped", "retained"),
    " \\\\"),
  "\\bottomrule",
  "\\end{tabular}"
)
writeLines(chain_health_tex, file.path(tables_dir, "krd_chain_health.tex"))

# ---------------------------------------------------------------------------
# LaTeX table: model comparison (all vs clean)
# ---------------------------------------------------------------------------
dw_all   <- cmp$value[cmp$metric == "delta_WAIC_all"]
dl_all   <- cmp$value[cmp$metric == "delta_LOO_all"]
dw_clean <- cmp$value[cmp$metric == "delta_WAIC_clean"]
dl_clean <- cmp$value[cmp$metric == "delta_LOO_clean"]

model_cmp_tex <- c(
  "\\begin{tabular}{lrrl}",
  "\\toprule",
  "Criterion & All chains & Clean chains & Favours \\\\",
  "\\midrule",
  sprintf("$\\Delta$WAIC ($\\text{WAIC}_{\\mathcal{M}_0} - \\text{WAIC}_{\\mathcal{M}_1}$) & %s & %s & %s \\\\",
          formatC(dw_all,   format = "f", big.mark = ",", digits = 0),
          formatC(dw_clean, format = "f", big.mark = ",", digits = 0),
          if (dw_clean > 0) "$\\mathcal{M}_1$" else "$\\mathcal{M}_0$"),
  sprintf("$\\Delta$LOO  ($\\text{LOO}_{\\mathcal{M}_0} - \\text{LOO}_{\\mathcal{M}_1}$) & %s & %s & %s \\\\",
          formatC(dl_all,   format = "f", big.mark = ",", digits = 0),
          formatC(dl_clean, format = "f", big.mark = ",", digits = 0),
          if (dl_clean > 0) "$\\mathcal{M}_1$" else "$\\mathcal{M}_0$"),
  "\\bottomrule",
  "\\end{tabular}"
)
writeLines(model_cmp_tex, file.path(tables_dir, "krd_model_comparison.tex"))

# ---------------------------------------------------------------------------
# LaTeX table: posterior summary (K+, sigma_N, wN_mean, kappa, alpha)
# ---------------------------------------------------------------------------
fmt <- function(x, digits = 3) formatC(x, format = "g", digits = digits)
fmt_int <- function(x) formatC(x, format = "d")

ps <- post_summ
post_tex <- c(
  "\\begin{tabular}{lcccl}",
  "\\toprule",
  "Parameter & Mean & Median & 95\\% CI & Note \\\\",
  "\\midrule",
  sprintf("$K^+$ & %s & %s & (%s, %s) & Occupied components \\\\",
          fmt_int(ps$mean[ps$parameter == "K_plus"]),
          fmt_int(ps$median[ps$parameter == "K_plus"]),
          fmt_int(ps$ci_lo[ps$parameter == "K_plus"]),
          fmt_int(ps$ci_hi[ps$parameter == "K_plus"])),
  sprintf("$\\sigma_N$ & %.2f & %.2f & (%.2f, %.2f) & Background SD \\\\",
          ps$mean[ps$parameter == "sigma_0"],
          ps$median[ps$parameter == "sigma_0"],
          ps$ci_lo[ps$parameter == "sigma_0"],
          ps$ci_hi[ps$parameter == "sigma_0"]),
  sprintf("$\\bar{\\omega}^N$ & %.2e & %.2e & (%.2e, %.2e) & Mean contamination \\\\",
          ps$mean[ps$parameter == "wN_mean"],
          ps$median[ps$parameter == "wN_mean"],
          ps$ci_lo[ps$parameter == "wN_mean"],
          ps$ci_hi[ps$parameter == "wN_mean"]),
  "$\\kappa$ & \\multicolumn{3}{c}{Fixed at 1} & Tracking precision \\\\",
  "$\\alpha$ & \\multicolumn{3}{c}{Fixed at 1} & DP concentration \\\\",
  "Tracking corr.\\ & \\multicolumn{3}{c}{$(\\kappa+1)/(\\alpha+\\kappa+1) = 2/3$} & Structural \\\\",
  "\\bottomrule",
  "\\end{tabular}"
)
writeLines(post_tex, file.path(tables_dir, "krd_posterior_summary.tex"))

# ---------------------------------------------------------------------------
# LaTeX table: convergence summary (sigma_N + wN summary)
# ---------------------------------------------------------------------------
wN_rows_M1 <- conv_M1[grepl("^omega_0", conv_M1$parameter), ]
wN_rows_M0 <- conv_M0[grepl("^omega_0", conv_M0$parameter), ]

conv_tex <- c(
  "\\begin{tabular}{llrrr}",
  "\\toprule",
  "Model & Parameter set & Median $\\hat{R}$ (max) & Median bulk ESS (min) & N params \\\\",
  "\\midrule",
  sprintf("$\\mathcal{M}_1$ & $\\sigma_0$ & %.3f & %.0f & 1 \\\\",
          conv_M1$rhat[conv_M1$parameter == "sigma_0"],
          conv_M1$bulk_ess[conv_M1$parameter == "sigma_0"]),
  sprintf("$\\mathcal{M}_1$ & $\\omega_{0it}$ (per-observation) & %.3f (%.3f) & %.0f (%.0f) & %d \\\\",
          median(wN_rows_M1$rhat, na.rm = TRUE),
          max(wN_rows_M1$rhat, na.rm = TRUE),
          median(wN_rows_M1$bulk_ess, na.rm = TRUE),
          min(wN_rows_M1$bulk_ess, na.rm = TRUE),
          nrow(wN_rows_M1)),
  sprintf("$\\mathcal{M}_0$ & $\\sigma_0$ & %.3f & %.0f & 1 \\\\",
          conv_M0$rhat[conv_M0$parameter == "sigma_0"],
          conv_M0$bulk_ess[conv_M0$parameter == "sigma_0"]),
  sprintf("$\\mathcal{M}_0$ & $\\omega_{0it}$ (per-observation) & %.3f (%.3f) & %.0f (%.0f) & %d \\\\",
          median(wN_rows_M0$rhat, na.rm = TRUE),
          max(wN_rows_M0$rhat, na.rm = TRUE),
          median(wN_rows_M0$bulk_ess, na.rm = TRUE),
          min(wN_rows_M0$bulk_ess, na.rm = TRUE),
          nrow(wN_rows_M0)),
  "\\bottomrule",
  "\\end{tabular}"
)
writeLines(conv_tex, file.path(tables_dir, "krd_convergence_summary.tex"))

# Also write the outliers (patients with R-hat > 1.05) to a separate table
outlier_M1 <- wN_rows_M1[wN_rows_M1$rhat > 1.05 & !is.na(wN_rows_M1$rhat), ]
outlier_M0 <- wN_rows_M0[wN_rows_M0$rhat > 1.05 & !is.na(wN_rows_M0$rhat), ]

if (nrow(outlier_M1) > 0 || nrow(outlier_M0) > 0) {
  outlier_tex <- c(
    "\\begin{tabular}{llrr}",
    "\\toprule",
    "Model & Observation & $\\hat{R}$ & Bulk ESS \\\\",
    "\\midrule"
  )
  if (nrow(outlier_M1) > 0) {
    outlier_tex <- c(outlier_tex,
      sprintf("$\\mathcal{M}_1$ & %s & %.3f & %.0f \\\\",
              gsub("omega_0\\[|\\]", "", outlier_M1$parameter),
              outlier_M1$rhat, outlier_M1$bulk_ess))
  }
  if (nrow(outlier_M0) > 0) {
    outlier_tex <- c(outlier_tex,
      sprintf("$\\mathcal{M}_0$ & %s & %.3f & %.0f \\\\",
              gsub("omega_0\\[|\\]", "", outlier_M0$parameter),
              outlier_M0$rhat, outlier_M0$bulk_ess))
  }
  outlier_tex <- c(outlier_tex, "\\bottomrule", "\\end{tabular}")
  writeLines(outlier_tex, file.path(tables_dir, "krd_convergence_outliers.tex"))
}

# ---------------------------------------------------------------------------
# LaTeX table: per-patient wN
# ---------------------------------------------------------------------------
wN_tex <- c(
  "\\begin{tabular}{rrrr}",
  "\\toprule",
  "Patient & Post.\\ mean & 2.5\\% & 97.5\\% \\\\",
  "\\midrule",
  sprintf("%d & %.2e & %.2e & %.2e \\\\",
          wN_patient$patient,
          wN_patient$post_mean,
          wN_patient$ci_lo,
          wN_patient$ci_hi),
  "\\bottomrule",
  "\\end{tabular}"
)
writeLines(wN_tex, file.path(tables_dir, "krd_omega0_by_patient.tex"))

# ---------------------------------------------------------------------------
# Summary report to console
# ---------------------------------------------------------------------------
cat("\n=== Manuscript outputs regenerated from post-remediation results ===\n")
cat("\nModel comparison:\n")
cat(sprintf("  Delta WAIC (clean):  %+10.1f  [favours %s]\n",
            dw_clean, if (dw_clean > 0) "M1" else "M0"))
cat(sprintf("  Delta  LOO (clean):  %+10.1f  [favours %s]\n",
            dl_clean, if (dl_clean > 0) "M1" else "M0"))

cat("\nChain health (M1):\n")
m1_gap <- diff(range(ch_health$mean_loglik[ch_health$model == "M1"]))
cat(sprintf("  Spread max-min:  %.1f log-units  (chains trapped: %d)\n",
            m1_gap, sum(ch_health$trapped[ch_health$model == "M1"])))

cat("\nChain health (M0):\n")
m0_gap <- diff(range(ch_health$mean_loglik[ch_health$model == "M0"]))
cat(sprintf("  Spread max-min:  %.1f log-units  (chains trapped: %d)\n",
            m0_gap, sum(ch_health$trapped[ch_health$model == "M0"])))

cat("\nPosterior summary:\n")
cat(sprintf("  K+ (median)      = %d\n", post_summ$median[post_summ$parameter == "K_plus"]))
cat(sprintf("  sigma_N  (mean)  = %.2f\n", post_summ$mean[post_summ$parameter == "sigma_0"]))
cat(sprintf("  bar(omega_N)     = %.2e\n", post_summ$mean[post_summ$parameter == "wN_mean"]))

cat("\nFiles written:\n")
cat("  figures/krd_chain_health.pdf\n")
cat("  figures/krd_loglik_traces.pdf\n")
cat("  figures/krd_traceplots_remediated.pdf\n")
cat("  tables/krd_chain_health.tex\n")
cat("  tables/krd_model_comparison.tex\n")
cat("  tables/krd_posterior_summary.tex\n")
cat("  tables/krd_convergence_summary.tex\n")
if (nrow(outlier_M1) > 0 || nrow(outlier_M0) > 0) {
  cat("  tables/krd_convergence_outliers.tex\n")
}
cat("  tables/krd_omega0_by_patient.tex\n\n")
