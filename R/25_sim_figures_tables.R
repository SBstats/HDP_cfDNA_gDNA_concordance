###############################################################################
## R/25_sim_figures_tables.R
##
## Build the simulation figures and LaTeX tables from the aggregated summaries
## (results/sim_aggregated/scenario_{1,2,3,4}_summary.csv), matching the exact
## file names the methods paper (Writeup Draft2.tex) references as
## commented \input / \includegraphics placeholders:
##
##   figures/sim_scenario1_contraction.pdf      (S1 log-log normalized-risk rate)
##   figures/sim_scenario2b_pi_consistency.pdf  (S2 pi L1 error vs p; consistency)
##   figures/sim_scenario3_roc.pdf              (S3 detection vs kappa* + FPR)
##   tables/sim_scenario1_theta.tex
##   tables/sim_scenario2_consistency.tex
##   tables/sim_scenario3_operating.tex
##   tables/sim_scenario4_misspecification.tex
##
## Usage: Rscript R/25_sim_figures_tables.R [base_dir]
###############################################################################

suppressPackageStartupMessages({ library(ggplot2) })
args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[1] else getwd()
agg_dir <- file.path(base_dir, "results", "sim_aggregated")
fig_dir <- file.path(base_dir, "figures");  dir.create(fig_dir, showWarnings = FALSE)
tab_dir <- file.path(base_dir, "tables");   dir.create(tab_dir, showWarnings = FALSE)

rd <- function(s) {
  f <- file.path(agg_dir, sprintf("scenario_%d_summary.csv", s))
  if (!file.exists(f)) { message("missing: ", f); return(NULL) }
  read.csv(f, stringsAsFactors = FALSE)
}
fmt <- function(x, d = 3) ifelse(is.na(x), "--", formatC(x, format = "f", digits = d))

## ---------------------------------------------------------------------------
## Scenario 1 --- contraction-rate figure (log-log) + signature/competitor table
## ---------------------------------------------------------------------------
s1 <- rd(1)
if (!is.null(s1)) {
  s1$sk <- factor(s1$s_k_frac)
  p1 <- ggplot(s1, aes(p, norm_mse_mean, color = sk, group = sk)) +
    geom_line() + geom_point() +
    scale_x_log10() + scale_y_log10() +
    labs(x = "p (log scale)", y = "Normalized risk MSE / [s_k log(p/s_k)] (log scale)",
         color = expression(s[k]/p),
         title = "Scenario 1: posterior contraction rate") +
    theme_bw(base_size = 10)
  ggsave(file.path(fig_dir, "sim_scenario1_contraction.pdf"), p1, width = 6.5, height = 4)

  ## table: horseshoe vs competitors (per-gene MSE, FPR) + recovery + coverage.
  ## MSE columns are PER-GENE means (÷p), comparable across p; "Rec." is the
  ## fraction of the K_true true signatures recovered (junk-robust matcher).
  ## Competitor MSE is rescaled to per-gene as well when a total was stored.
  ord <- order(s1$s_k_frac, s1$p)
  s1o <- s1[ord, ]
  has_rec <- "recovery_rate_mean" %in% names(s1o)
  lines <- c(
    "\\begin{tabular}{rrrrrrrr}", "\\toprule",
    "$p$ & $s_k/p$ & MSE (HS) & Norm. & Rec. & TPR & FPR & MSE (lasso) \\\\",
    "\\midrule")
  for (i in seq_len(nrow(s1o))) {
    r <- s1o[i, ]
    rec <- if (has_rec) fmt(r$recovery_rate_mean,2) else "--"
    lines <- c(lines, sprintf("%d & %.2f & %s & %s & %s & %s & %s & %s \\\\",
      r$p, r$s_k_frac, fmt(r$mse_mean,3), fmt(r$norm_mse_mean,3), rec, fmt(r$tpr_mean,2),
      fmt(r$fpr_mean,3), fmt(r$lasso_mse_mean,2)))
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}")
  writeLines(lines, file.path(tab_dir, "sim_scenario1_theta.tex"))
}

## ---------------------------------------------------------------------------
## Scenario 2 --- posterior consistency of pi (L1 error vs p) + table
## ---------------------------------------------------------------------------
s2 <- rd(2)
if (!is.null(s2)) {
  s2b <- s2[order(s2$p), ]
  p2b <- ggplot(s2b, aes(p, pi_l1_mean)) +
    geom_line() + geom_point() +
    geom_errorbar(aes(ymin = pi_l1_mean - pi_l1_mcse, ymax = pi_l1_mean + pi_l1_mcse), width = 0.03) +
    scale_x_log10() +
    labs(x = "p (log scale)", y = expression("Posterior "*L[1]*" error "*"||"*hat(pi)-pi^"*"*"||"[1]),
         title = "Scenario 2: posterior consistency of pi") +
    theme_bw(base_size = 10)
  ggsave(file.path(fig_dir, "sim_scenario2b_pi_consistency.pdf"), p2b, width = 6, height = 4)

  ## table: pi L1 (and L2) error vs p, should decrease with p
  lines <- c("\\begin{tabular}{rrrr}", "\\toprule",
    "$p$ & $L_1$ error & (MCSE) & $L_2$ error \\\\", "\\midrule")
  for (i in seq_len(nrow(s2b))) { r <- s2b[i, ]
    lines <- c(lines, sprintf("%g & %s & %s & %s \\\\", r$p,
      fmt(r$pi_l1_mean,3), fmt(r$pi_l1_mcse,4), fmt(r$pi_l2_mean,3))) }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}")
  writeLines(lines, file.path(tab_dir, "sim_scenario2_consistency.tex"))
}

## ---------------------------------------------------------------------------
## Scenario 3 --- detection-vs-kappa figure (+FPR marker) + operating table
## ---------------------------------------------------------------------------
s3 <- rd(3)
if (!is.null(s3)) {
  alt  <- s3[s3$sub == "3_alt", ]
  null <- s3[s3$sub == "3_null", ]
  ## PRIMARY detection = the change-based increment test (R5); the increment
  ## EXCESS is the continuous tracking score plotted against kappa*.
  fpr0 <- if (nrow(null)) null$incr_detection_rate[1] else NA
  if (nrow(alt)) {
    p3 <- ggplot(alt, aes(kappa_true, incr_detection_rate)) +
      geom_line() + geom_point() +
      { if (!is.na(fpr0)) geom_point(data = data.frame(kappa_true = 0, incr_detection_rate = fpr0),
          aes(kappa_true, incr_detection_rate), color = "red", size = 3) } +
      annotate("text", x = 0, y = ifelse(is.na(fpr0), 0, fpr0), label = "FPR (M0 true)",
               hjust = -0.1, size = 3, color = "red") +
      ylim(0, 1) +
      labs(x = expression(kappa^"*"~"(0 = M0 truth)"),
           y = "Increment-test detection rate",
           title = "Scenario 3: tracking-detection operating characteristics (increment estimand)") +
      theme_bw(base_size = 10)
    ggsave(file.path(fig_dir, "sim_scenario3_roc.pdf"), p3, width = 6, height = 4)
  }
  ## table: increment (primary) + correlation (secondary) detection/excess by truth
  lines <- c("\\begin{tabular}{lrrr}", "\\toprule",
    "Truth & incr.\\ detect/FPR & incr.\\ excess & occ.\\ corr.\\ excess \\\\", "\\midrule")
  if (nrow(null)) lines <- c(lines, sprintf("$M_0$ (null) & %s & %s & %s \\\\",
    fmt(null$incr_detection_rate[1],3), fmt(null$incr_excess_mean[1],3), fmt(null$cor_excess_mean[1],3)))
  for (i in seq_len(nrow(alt))) { r <- alt[order(alt$kappa_true), ][i, ]
    lines <- c(lines, sprintf("$M_1,\\ \\kappa^*=%g$ & %s & %s & %s \\\\", r$kappa_true,
      fmt(r$incr_detection_rate,3), fmt(r$incr_excess_mean,3), fmt(r$cor_excess_mean,3))) }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}")
  writeLines(lines, file.path(tab_dir, "sim_scenario3_operating.tex"))
}

## ---------------------------------------------------------------------------
## Scenario 4 --- misspecification robustness table (no figure in the writeup)
## ---------------------------------------------------------------------------
s4 <- rd(4)
if (!is.null(s4)) {
  lab <- c("4_t3" = "$t_3$ noise", "4_corr" = "correlated genes",
           "4_hk" = "housekeeping bg", "4_all" = "all combined")
  lines <- c("\\begin{tabular}{lrrrr}", "\\toprule",
    "Departure & occ.\\ corr.\\ excess & incr.\\ excess & median $\\hat K^+$ & favors $M_1$ \\\\", "\\midrule")
  for (i in seq_len(nrow(s4))) { r <- s4[i, ]
    lines <- c(lines, sprintf("%s & %s & %s & %s & %s \\\\",
      ifelse(is.na(lab[r$sub]), r$sub, lab[r$sub]),
      fmt(r$cor_excess_mean,3), fmt(r$incr_excess_mean,3),
      fmt(r$Kplus_median_mean,2), fmt(r$favors_M1_rate,2))) }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}")
  writeLines(lines, file.path(tab_dir, "sim_scenario4_misspecification.tex"))
}

cat("Simulation figures and tables written to figures/ and tables/.\n")
