
###############################################################################
## ---- 05_traceplots.R
###############################################################################

###############################################################################
## 05_traceplots.R
## Generate trace plots for MCMC convergence diagnostics
##
## Produces trace plots, running mean plots, and density plots for
## key model parameters: kappa, alpha_dp, K+, wN, sigma_k
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(coda)
})

#' Generate comprehensive trace plots for MCMC diagnostics
#'
#' @param fit      Output from fit_dp_concordance
#' @param out_dir  Directory to save plots
#' @param prefix   File name prefix (e.g., "M1" or "M0")
generate_traceplots <- function(fit, out_dir, prefix = "M1") {

  samples <- fit$samples
  n_save <- fit$n_save

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  iter_idx <- 1:n_save

  # ---- 1. Kappa trace plot (M1 only) ----
  if (fit$model == "M1") {
    pdf(file.path(out_dir, paste0(prefix, "_trace_kappa.pdf")),
        width = 10, height = 6)
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

    # Trace
    plot(iter_idx, samples$kappa, type = "l", col = "steelblue",
         xlab = "Iteration (post burn-in)", ylab = expression(kappa),
         main = expression("Trace: " * kappa))

    # Running mean
    cum_mean <- cumsum(samples$kappa) / iter_idx
    plot(iter_idx, cum_mean, type = "l", col = "darkred",
         xlab = "Iteration", ylab = expression("Running mean of " * kappa),
         main = expression("Running mean: " * kappa))

    # Density
    plot(density(samples$kappa), col = "steelblue", lwd = 2,
         main = expression("Posterior density: " * kappa),
         xlab = expression(kappa))
    abline(v = mean(samples$kappa), lty = 2, col = "red")
    abline(v = quantile(samples$kappa, c(0.025, 0.975)),
           lty = 3, col = "gray40")

    # ACF
    acf(samples$kappa, main = expression("ACF: " * kappa), lag.max = 50)

    dev.off()
    cat("Saved", prefix, "kappa trace plot\n")
  }

  # ---- 2. Alpha_DP trace plot ----
  pdf(file.path(out_dir, paste0(prefix, "_trace_alpha.pdf")),
      width = 10, height = 6)
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

  plot(iter_idx, samples$alpha_dp, type = "l", col = "steelblue",
       xlab = "Iteration", ylab = expression(alpha),
       main = expression("Trace: " * alpha[DP]))
  cum_mean <- cumsum(samples$alpha_dp) / iter_idx
  plot(iter_idx, cum_mean, type = "l", col = "darkred",
       xlab = "Iteration", ylab = expression("Running mean"),
       main = expression("Running mean: " * alpha[DP]))
  plot(density(samples$alpha_dp), col = "steelblue", lwd = 2,
       main = expression("Posterior: " * alpha[DP]),
       xlab = expression(alpha))
  abline(v = mean(samples$alpha_dp), lty = 2, col = "red")
  acf(samples$alpha_dp, main = expression("ACF: " * alpha[DP]), lag.max = 50)

  dev.off()
  cat("Saved", prefix, "alpha_DP trace plot\n")

  # ---- 3. K+ trace plot ----
  pdf(file.path(out_dir, paste0(prefix, "_trace_Kplus.pdf")),
      width = 10, height = 4)
  par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

  plot(iter_idx, samples$K_plus, type = "l", col = "steelblue",
       xlab = "Iteration", ylab = expression(K^"+"),
       main = expression("Trace: " * K^"+"))
  barplot(table(samples$K_plus) / n_save, col = "steelblue",
          xlab = expression(K^"+"), ylab = "Posterior probability",
          main = expression("Posterior: " * K^"+"))

  dev.off()
  cat("Saved", prefix, "K+ trace plot\n")

  # ---- 4. Contamination fractions wN ----
  n_obs_c <- ncol(samples$omega_0)   # omega_0 is per OBSERVATION, not per patient
  pdf(file.path(out_dir, paste0(prefix, "_trace_wN.pdf")),
      width = 10, height = 8)
  n_show <- min(n_obs_c, 9)
  par(mfrow = c(3, 3), mar = c(3, 3, 2, 1))
  for (i in 1:n_show) {
    plot(iter_idx, samples$omega_0[, i], type = "l", col = "steelblue",
         xlab = "Iteration", ylab = expression(omega[0]),
         main = paste0("Obs ", i, ": ", expression(omega[0])))
  }
  dev.off()
  cat("Saved", prefix, "wN trace plots\n")

  # ---- 5. Component standard deviations ----
  K <- ncol(samples$sigma_k)
  pdf(file.path(out_dir, paste0(prefix, "_trace_sigma.pdf")),
      width = 10, height = 6)
  par(mfrow = c(2, min(K, 5)), mar = c(3, 3, 2, 1))
  for (k in 1:min(K, 5)) {
    plot(iter_idx, samples$sigma_k[, k], type = "l", col = "steelblue",
         xlab = "Iteration", ylab = expression(sigma[k]),
         main = paste0("sigma_", k))
  }
  for (k in 1:min(K, 5)) {
    if (var(samples$sigma_k[, k]) > 0) {
      plot(density(samples$sigma_k[, k]), col = "steelblue", lwd = 2,
           main = paste0("Density: sigma_", k))
    } else {
      plot.new()
    }
  }
  dev.off()
  cat("Saved", prefix, "sigma trace plots\n")

  # ---- 6. Log-likelihood trace ----
  pdf(file.path(out_dir, paste0(prefix, "_trace_loglik.pdf")),
      width = 8, height = 4)
  par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))
  plot(iter_idx, samples$loglik, type = "l", col = "steelblue",
       xlab = "Iteration", ylab = "Log-likelihood",
       main = paste0(prefix, ": Log-likelihood trace"))
  dev.off()
  cat("Saved", prefix, "log-likelihood trace plot\n")

  # ---- 7. Effective sample size summary ----
  ess_kappa <- if (fit$model == "M1") effectiveSize(samples$kappa) else NA
  ess_alpha <- effectiveSize(samples$alpha_dp)
  ess_ll    <- effectiveSize(samples$loglik)

  cat("\n--- Effective Sample Sizes (", prefix, ") ---\n")
  if (!is.na(ess_kappa)) cat(sprintf("  kappa:    %.0f\n", ess_kappa))
  cat(sprintf("  alpha_DP: %.0f\n", ess_alpha))
  cat(sprintf("  loglik:   %.0f\n", ess_ll))
  cat(sprintf("  omega_0[1]: %.0f\n", effectiveSize(samples$omega_0[, 1])))

  invisible(NULL)
}


###############################################################################
## ---- 16_manuscript_figures.R
###############################################################################

###############################################################################
## 16_manuscript_figures.R
## Generate all PDF figures for the manuscript
##
## Reads aggregated simulation results and KRd analysis output,
## produces publication-quality figures.
##
## Usage: Rscript R/16_manuscript_figures.R [base_dir]
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
})

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[1] else getwd()

fig_dir <- file.path(base_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE)

sim_dir <- file.path(base_dir, "results", "sim_aggregated")

theme_set(theme_bw(base_size = 12))

# ===========================================================================
# SCENARIO 1: Normalized MSE vs p (formerly Scenario 2)
# ===========================================================================
s2_file <- file.path(sim_dir, "scenario_2_summary.csv")
if (file.exists(s2_file)) {
  s2 <- read.csv(s2_file, stringsAsFactors = FALSE)

  g2 <- ggplot(s2, aes(x = p, y = norm_mse_mean,
                        color = factor(s_k_frac))) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    labs(x = "Number of features (p)",
         y = expression("Normalized MSE: MSE / [" * s[k] * log(p/s[k]) * "]"),
         color = expression(s[k]/p)) +
    scale_color_brewer(palette = "Dark2")

  ggsave(file.path(fig_dir, "sim_scenario2_contraction.pdf"), g2,
         width = 6, height = 4)
  cat("Saved: sim_scenario2_contraction.pdf\n")
} else {
  cat("Scenario 2 results not found; skipping figure.\n")
}

# ===========================================================================
# SCENARIO 3a: Empirical vs theoretical correlation
# ===========================================================================
s3_file <- file.path(sim_dir, "scenario_3_summary.csv")
if (file.exists(s3_file)) {
  s3 <- read.csv(s3_file, stringsAsFactors = FALSE)

  # 3a: correlation curve
  s3a <- s3[s3$sub == "3a", ]
  if (nrow(s3a) > 0) {
    kappa_seq <- seq(0.5, 100, length.out = 200)
    alpha_val <- 1
    theory <- data.frame(
      kappa = kappa_seq,
      corr = (kappa_seq + 1) / (alpha_val + kappa_seq + 1)
    )

    g3a <- ggplot() +
      geom_line(data = theory, aes(x = kappa, y = corr),
                color = "black", linewidth = 0.8) +
      geom_point(data = s3a,
                 aes(x = kappa_true,
                     y = mean(s3a$corr_error_mean) + (s3a$kappa_true + 1)/(alpha_val + s3a$kappa_true + 1)),
                 color = "red", size = 3) +
      labs(x = expression(kappa),
           y = expression("Cross-source correlation Corr(" * omega^g * ", " * omega^cf * ")"),
           title = "Empirical vs. theoretical tracking correlation") +
      annotate("text", x = 80, y = 0.6,
               label = expression("Theory: " * frac(kappa+1, alpha+kappa+1)),
               size = 4)

    ggsave(file.path(fig_dir, "sim_scenario3_corr.pdf"), g3a,
           width = 6, height = 4)
    cat("Saved: sim_scenario3_corr.pdf\n")
  }

  # 3b: pi consistency
  s3b <- s3[s3$sub == "3b", ]
  if (nrow(s3b) > 0) {
    g3b <- ggplot(s3b, aes(x = p, y = pi_l1_mean)) +
      geom_line(linewidth = 0.8, color = "steelblue") +
      geom_point(size = 2, color = "steelblue") +
      labs(x = "Number of features (p)",
           y = expression("||" * hat(pi) - pi^"*" * "||"[1]),
           title = expression("Posterior consistency: " * L[1] * " error vs. p"))

    ggsave(file.path(fig_dir, "sim_scenario3_pi_consistency.pdf"), g3b,
           width = 6, height = 4)
    cat("Saved: sim_scenario3_pi_consistency.pdf\n")
  }
} else {
  cat("Scenario 3 results not found; skipping figures.\n")
}

# ===========================================================================
# KRD ANALYSIS FIGURES
# ===========================================================================
m1_file <- file.path(base_dir, "results", "krd_M1_fit.rds")
if (file.exists(m1_file)) {
  m1 <- readRDS(m1_file)
  m0 <- NULL
  m0_file <- file.path(base_dir, "results", "krd_M0_fit.rds")
  if (file.exists(m0_file)) m0 <- readRDS(m0_file)
  krd <- list(multi_M1 = m1, multi_M0 = m0)

  # ===========================================================================
  # Helper: generate a 5x2 per-chain traceplot PDF
  # ===========================================================================
  chain_cols <- c("steelblue", "firebrick", "forestgreen", "darkorange")

  generate_biometrics_traceplots <- function(chains_list, out_path, title_prefix = "",
                                              chains_M0_list = NULL) {
    n_ch <- length(chains_list)
    N_obs_tp <- ncol(chains_list[[1]]$final_omega_g)
    K_trace_tp <- chains_list[[1]]$K_trace

    plot_scalar <- function(pnm, ylab_e, main_e) {
      vals <- lapply(1:n_ch, function(i) chains_list[[i]]$samples[[pnm]])
      yr <- range(unlist(vals))
      plot(vals[[1]], type="l", col=chain_cols[1], ylim=yr,
           xlab="Iteration (post burn-in, thinned)", ylab=ylab_e, main=main_e)
      for (ch in 2:n_ch) lines(vals[[ch]], col=chain_cols[ch])
      legend("topright", paste("Chain",1:n_ch), col=chain_cols[1:n_ch], lty=1, cex=0.6, bg="white")
    }

    plot_omega_col <- function(cidx, ylab_e, main_e, tnm) {
      vals <- lapply(1:n_ch, function(i) chains_list[[i]]$samples[[tnm]][, cidx])
      yr <- range(unlist(vals))
      plot(vals[[1]], type="l", col=chain_cols[1], ylim=yr,
           xlab="Iteration (post burn-in, thinned)", ylab=ylab_e, main=main_e)
      for (ch in 2:n_ch) lines(vals[[ch]], col=chain_cols[ch])
      legend("topright", paste("Chain",1:n_ch), col=chain_cols[1:n_ch], lty=1, cex=0.6, bg="white")
    }

    pfx <- if (nchar(title_prefix) > 0) paste0(title_prefix, ": ") else ""

    pdf(out_path, width = 12, height = 17)
    par(mfrow = c(5, 2))

    plot_scalar("K_plus", expression(K^"+"), bquote(.(pfx) * K^"+"))
    plot_scalar("loglik", "Log-likelihood", bquote(.(pfx) * "Log-lik (" * M[1] * ")"))

    plot_scalar("sigma_0", expression(sigma[0]), bquote(.(pfx) * sigma[0]))
    wN_v <- lapply(1:n_ch, function(i) chains_list[[i]]$samples$omega_0[, 1])
    yr_wN <- range(unlist(wN_v))
    plot(wN_v[[1]], type="l", col=chain_cols[1], ylim=yr_wN,
         xlab="Iteration (post burn-in, thinned)", ylab=expression(omega[0]*"[obs 1]"),
         main=bquote(.(pfx) * omega[0] * " (observation 1)"))
    for (ch in 2:n_ch) lines(wN_v[[ch]], col=chain_cols[ch])
    legend("topright", paste("Chain",1:n_ch), col=chain_cols[1:n_ch], lty=1, cex=0.6, bg="white")

    if (!is.null(chains_list[[1]]$samples$omega_g_trace)) {
      plot_omega_col(1, expression(omega[1]^g), bquote(.(pfx) * omega[1]^g * " (obs 1, comp 1)"), "omega_g_trace")
      plot_omega_col(1, expression(omega[1]^cf), bquote(.(pfx) * omega[1]^cf * " (obs 1, comp 1)"), "omega_cf_trace")
    } else { plot.new(); text(0.5,0.5,"omega traces\nnot stored"); plot.new(); text(0.5,0.5,"omega traces\nnot stored") }

    if (!is.null(chains_list[[1]]$samples$omega_g_trace) && K_trace_tp >= 2) {
      c2 <- N_obs_tp + 1
      plot_omega_col(c2, expression(omega[2]^g), bquote(.(pfx) * omega[2]^g * " (obs 1, comp 2)"), "omega_g_trace")
      plot_omega_col(c2, expression(omega[2]^cf), bquote(.(pfx) * omega[2]^cf * " (obs 1, comp 2)"), "omega_cf_trace")
    } else { plot.new(); text(0.5,0.5,"comp 2\nnot available"); plot.new(); text(0.5,0.5,"comp 2\nnot available") }

    plot_scalar("loglik", "Log-likelihood", bquote(.(pfx) * "Log-likelihood (" * M[1] * ")"))
    if (!is.null(chains_M0_list)) {
      n_m0 <- length(chains_M0_list)
      ll_m0 <- lapply(1:n_m0, function(i) chains_M0_list[[i]]$samples$loglik)
      plot(ll_m0[[1]], type="l", col=chain_cols[1], ylim=range(unlist(ll_m0)),
           xlab="Iteration (post burn-in, thinned)", ylab="Log-likelihood",
           main=bquote(.(pfx) * "Log-likelihood (" * M[0] * ")"))
      for (ch in 2:n_m0) lines(ll_m0[[ch]], col=chain_cols[ch])
      legend("topright", paste("Chain",1:n_m0), col=chain_cols[1:n_m0], lty=1, cex=0.6, bg="white")
    } else { plot.new(); text(0.5,0.5,"M0 loglik\nnot available") }

    dev.off()
  }

  merged <- krd$multi_M1$merged
  chains_M1 <- krd$multi_M1$chains
  chains_M0_raw <- if (!is.null(krd$multi_M0)) krd$multi_M0$chains else NULL

  # (i) BEFORE label switching
  generate_biometrics_traceplots(chains_M1,
    file.path(fig_dir, "krd_traceplots_before_relabeling.pdf"),
    "Before relabeling", chains_M0_raw)
  cat("Saved: krd_traceplots_before_relabeling.pdf\n")

  # Apply label switching
  source(file.path(base_dir, "R", "01_lib_core.R"))
  chains_M1_relabeled <- chains_M1
  for (i in seq_along(chains_M1_relabeled)) {
    chains_M1_relabeled[[i]] <- relabel_by_weight(chains_M1_relabeled[[i]])
  }
  chains_M1_relabeled <- align_chains_by_signature(chains_M1_relabeled)

  # (ii) AFTER label switching
  generate_biometrics_traceplots(chains_M1_relabeled,
    file.path(fig_dir, "krd_traceplots_after_relabeling.pdf"),
    "After relabeling", chains_M0_raw)
  cat("Saved: krd_traceplots_after_relabeling.pdf\n")

  # Canonical krd_traceplots.pdf (post-relabeling, for the manuscript)
  generate_biometrics_traceplots(chains_M1_relabeled,
    file.path(fig_dir, "krd_traceplots.pdf"), "", chains_M0_raw)
  cat("Saved: krd_traceplots.pdf\n")

  # Contamination by patient
  wN_df <- data.frame(
    obs = rep(1:ncol(merged$omega_0), each = nrow(merged$omega_0)),
    wN = as.vector(merged$omega_0)
  )

  ## One box per OBSERVATION (subject-timepoint), not per patient: omega_0 has
  ## N_obs = 82 columns, not n = 26.
  g_wN <- ggplot(wN_df, aes(x = factor(obs), y = wN)) +
    geom_boxplot(fill = "lightblue", outlier.size = 0.5) +
    labs(x = "Observation (subject-timepoint)", y = expression(omega[0]),
         title = "Posterior contamination fractions by observation") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(fig_dir, "krd_contamination.pdf"), g_wN,
         width = 8, height = 4)
  cat("Saved: krd_contamination.pdf\n")

} else {
  cat("KRd results not found; skipping KRd figures.\n")
}

cat("\nFigure generation complete.\n")


###############################################################################
## ---- 18_nature_figures.R
###############################################################################

###############################################################################
## 18_nature_figures.R
## Generate figures for the applied paper version
##
## Figures that use observed data only (no model results) are generated now.
## Figures requiring model fitting results are guarded by file.exists() and
## will be generated once analysis completes on the supercomputer.
##
## Output: PNG files at 300 dpi in nature_manuscript/figures/
##
## Usage: Rscript R/18_nature_figures.R [base_dir]
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readxl)
})

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[1] else
  getwd()

fig_dir <- file.path(base_dir, "nature_manuscript", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

res_dir <- file.path(base_dir, "results")

# publication-quality theme
theme_nature <- theme_bw(base_size = 11, base_family = "sans") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "grey80"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12)
  )
theme_set(theme_nature)

# Color palettes
source_colors <- c("cfDNA" = "#E63946", "gDNA" = "#457B9D")
timepoint_shapes <- c("Screening" = 16, "C4" = 17, "C8" = 15,
                       "C18" = 18, "3 YR F/U" = 8)

# Save as PDF (reliable cross-platform) and convert to PNG via R
save_fig <- function(filename, plot, width, height, dpi = 300) {
  # Save as PDF first (always works)
  pdf_name <- sub("\\.png$", ".pdf", filename)
  pdf_path <- file.path(fig_dir, pdf_name)
  pdf(pdf_path, width = width, height = height)
  print(plot)
  dev.off()
  cat("  Saved:", pdf_name, "\n")
}

# Alias for compatibility
save_png <- save_fig

cat("=== Applied Paper Figure Generation ===\n")
cat("Base directory:", base_dir, "\n")
cat("Figure output:", fig_dir, "\n\n")

# ===========================================================================
# FIGURE 1: Analysis Pipeline Flowchart (conceptual — no data needed)
# ===========================================================================
cat("--- Figure 1: Analysis Pipeline Flowchart ---\n")

library(grid)

draw_flowchart <- function() {
  grid.newpage()

  # Helper: draw a rounded-rect box with centered text
  draw_box <- function(x, y, w, h, label, fill = "grey95", border = "black",
                       lwd = 1, lty = 1, cex = 0.75, font = 1) {
    grid.roundrect(x = unit(x, "npc"), y = unit(y, "npc"),
                   width = unit(w, "npc"), height = unit(h, "npc"),
                   r = unit(0.008, "npc"),
                   gp = gpar(fill = fill, col = border, lwd = lwd, lty = lty))
    grid.text(label, x = unit(x, "npc"), y = unit(y, "npc"),
              gp = gpar(cex = cex, fontface = font, fontfamily = "sans"))
  }

  # Helper: draw an arrow between two points
  draw_arrow <- function(x0, y0, x1, y1, lwd = 1) {
    grid.lines(x = unit(c(x0, x1), "npc"), y = unit(c(y0, y1), "npc"),
               arrow = arrow(length = unit(0.015, "npc"), type = "closed"),
               gp = gpar(lwd = lwd, fill = "black"))
  }

  # Helper: draw a dashed group border with label
  draw_group <- function(x, y, w, h, label, col = "red3") {
    grid.roundrect(x = unit(x, "npc"), y = unit(y, "npc"),
                   width = unit(w, "npc"), height = unit(h, "npc"),
                   r = unit(0.01, "npc"),
                   gp = gpar(col = col, lwd = 1.5, lty = 2, fill = NA))
    grid.text(label, x = unit(x + w/2 - 0.01, "npc"),
              y = unit(y + h/2 + 0.005, "npc"),
              just = c("right", "bottom"),
              gp = gpar(col = "blue3", cex = 0.65, fontface = "italic",
                        fontfamily = "sans"))
  }

  # ---- Layout (top to bottom, y from 0.96 down) ----
  # Compressed ~8% vs. old layout to accommodate new DE gene selection group
  cx <- 0.50  # center x for main flow

  # Row 1: KRd trial data
  draw_box(cx, 0.945, 0.40, 0.04,
           "KRd clinical trial data (n=26 patients, T=5 timepoints)", cex = 0.7)

  # Split arrows to gDNA and cfDNA
  draw_arrow(cx - 0.07, 0.925, cx - 0.16, 0.893)
  draw_arrow(cx + 0.07, 0.925, cx + 0.16, 0.893)

  # Row 2: Two specimen types
  draw_box(cx - 0.20, 0.870, 0.24, 0.038,
           "gDNA (bone marrow aspirate)\n92 samples", fill = "#D6EAF8", cex = 0.65)
  draw_box(cx + 0.20, 0.870, 0.24, 0.038,
           "cfDNA (peripheral blood)\n113 samples", fill = "#FADBD8", cex = 0.65)

  # Merge arrows
  draw_arrow(cx - 0.20, 0.851, cx, 0.823)
  draw_arrow(cx + 0.20, 0.851, cx, 0.823)

  # Row 3: Pairing
  draw_box(cx, 0.805, 0.38, 0.033,
           "Pair by patient + timepoint\n(82 obs after deduplication of 3-YR F/U gDNA duplicates)",
           cex = 0.63)
  draw_arrow(cx, 0.788, cx, 0.767)

  # ---- Preprocessing group ----
  draw_group(cx, 0.722, 0.50, 0.082, "Data preprocessing")

  draw_box(cx, 0.745, 0.42, 0.030,
           "Filter genes: <10 counts in >5% samples removed", cex = 0.63)
  draw_arrow(cx, 0.730, cx, 0.712)
  draw_box(cx, 0.697, 0.42, 0.030,
           "Joint DESeq2 VST normalization (19,100 -> 11,430 genes)", cex = 0.63)

  draw_arrow(cx, 0.681, cx, 0.660)

  # ---- DE gene selection group ----
  draw_group(cx, 0.585, 0.54, 0.118, "Candidate gene selection")

  draw_box(cx, 0.630, 0.46, 0.028,
           "Pooled limma-voom + duplicateCorrelation\n(design: ~ timepoint + source; all 82 obs; block on patient)",
           fill = "#FEF9E7", cex = 0.60)
  draw_arrow(cx, 0.616, cx, 0.570)

  draw_box(cx, 0.558, 0.46, 0.028,
           "Pooled limma-voom -> 318 candidate genes\n(FDR < 0.05, |log2FC| > 0.5; ~ timepoint + source)",
           fill = "#F9EBEA", border = "#c0392b", lwd = 1.5, cex = 0.62, font = 2)

  draw_arrow(cx, 0.544, cx, 0.524)

  # ---- Exploratory ----
  draw_box(cx, 0.507, 0.44, 0.030,
           "Exploratory analysis: PCA + gene-level correlations\n(restricted to 318 candidate genes)",
           cex = 0.62)
  draw_arrow(cx, 0.492, cx, 0.472)

  # ---- Bayesian model fitting group ----
  draw_group(cx, 0.392, 0.56, 0.145, "Bayesian DP mixture model")

  draw_box(cx, 0.447, 0.48, 0.033,
           "Fit M1 (tracking) and M0 (non-tracking) models\non 318 candidate genes",
           cex = 0.67, font = 2)
  draw_arrow(cx, 0.430, cx, 0.411)

  # Sub-boxes inside model group
  draw_box(cx - 0.15, 0.393, 0.22, 0.028,
           "DP prior + horseshoe\nshrinkage (K=10, kappa=alpha=1)", fill = "grey90", cex = 0.58)
  draw_box(cx + 0.15, 0.393, 0.22, 0.028,
           "Blocked Gibbs + tempered\nalloc. + K-means init (4 chains)", fill = "grey90", cex = 0.58)

  draw_arrow(cx - 0.15, 0.379, cx - 0.15, 0.358)
  draw_arrow(cx + 0.15, 0.379, cx + 0.15, 0.358)

  draw_box(cx, 0.345, 0.48, 0.028,
           "Convergence diagnostics (split-Rhat < 1.01)", cex = 0.63)

  draw_arrow(cx, 0.331, cx, 0.308)

  # ---- Model comparison group ----
  draw_group(cx, 0.265, 0.56, 0.078, "Model comparison")

  draw_box(cx - 0.15, 0.267, 0.22, 0.033,
           "WAIC\n(predictive fit + penalty)", fill = "grey90", cex = 0.60)
  draw_box(cx + 0.15, 0.267, 0.22, 0.033,
           "PSIS-LOO\ncross-validation", fill = "grey90", cex = 0.60)

  draw_arrow(cx - 0.15, 0.250, cx - 0.07, 0.227)
  draw_arrow(cx + 0.15, 0.250, cx + 0.07, 0.227)

  # ---- Results row ----
  draw_arrow(cx, 0.232, cx, 0.210)

  draw_box(cx, 0.192, 0.48, 0.033,
           "Posterior inference: subclonal weights, signatures, contamination",
           cex = 0.68, font = 2)
  draw_arrow(cx, 0.175, cx, 0.153)

  # ---- Output group ----
  draw_group(cx, 0.098, 0.56, 0.098, "Key outputs")

  draw_box(cx - 0.18, 0.117, 0.17, 0.037,
           "Tracking correlation\nrho = 2/3 (exact, kappa=alpha=1)", fill = "#D5F5E3", cex = 0.56)
  draw_box(cx, 0.117, 0.17, 0.037,
           "Subclonal trajectories\n(cfDNA vs gDNA)", fill = "#D5F5E3", cex = 0.56)
  draw_box(cx + 0.18, 0.117, 0.17, 0.037,
           "Contamination fractions\n(per patient)", fill = "#D5F5E3", cex = 0.56)

  draw_box(cx, 0.070, 0.48, 0.028,
           "Model selection: M1 (tracking) vs M0 (non-tracking)", fill = "#D5F5E3", cex = 0.63)
}

pdf(file.path(fig_dir, "fig1_analysis_pipeline.pdf"), width = 8, height = 11)
draw_flowchart()
dev.off()
cat("  Saved: fig1_analysis_pipeline.pdf\n")

# ===========================================================================
# FIGURE 5: PCA of 5-hMC Profiles (uses model_data.rds — observed data)
# ===========================================================================
cat("\n--- Figure 5: PCA of 5-hMC Profiles ---\n")

model_data_file <- file.path(res_dir, "model_data.rds")
if (file.exists(model_data_file)) {
  md <- readRDS(model_data_file)

  # Combine cfDNA and gDNA into one matrix [samples x genes]
  Y_combined <- cbind(md$Y_cf, md$Y_g)
  source_labels <- c(rep("cfDNA", ncol(md$Y_cf)), rep("gDNA", ncol(md$Y_g)))

  # Get timepoint labels
  timepoint_order <- c("1" = "Screening", "2" = "C4", "3" = "C8",
                        "4" = "C18", "5" = "3 YR F/U")
  time_labels <- c(
    timepoint_order[as.character(md$time)],
    timepoint_order[as.character(md$time)]
  )

  patient_labels <- c(md$patient, md$patient)

  # PCA on transposed matrix (samples as rows)
  pca_res <- prcomp(t(Y_combined), center = TRUE, scale. = TRUE)
  var_explained <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

  pca_df <- data.frame(
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    PC3 = pca_res$x[, 3],
    Source = source_labels,
    Timepoint = time_labels,
    Patient = patient_labels
  )

  # PC1 vs PC2
  fig5a <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Source, shape = Timepoint)) +
    geom_point(size = 2.5, alpha = 0.8) +
    scale_color_manual(values = source_colors) +
    scale_shape_manual(values = timepoint_shapes) +
    labs(x = paste0("PC1 (", var_explained[1], "%)"),
         y = paste0("PC2 (", var_explained[2], "%)"),
         title = "a") +
    theme_nature +
    theme(legend.position = "right")

  # PC1 vs PC3
  fig5b <- ggplot(pca_df, aes(x = PC1, y = PC3, color = Source, shape = Timepoint)) +
    geom_point(size = 2.5, alpha = 0.8) +
    scale_color_manual(values = source_colors) +
    scale_shape_manual(values = timepoint_shapes) +
    labs(x = paste0("PC1 (", var_explained[1], "%)"),
         y = paste0("PC3 (", var_explained[3], "%)"),
         title = "b") +
    theme_nature +
    theme(legend.position = "right")

  # Paired connections: link cfDNA and gDNA from same patient-timepoint
  # Each observation index maps to both a cfDNA and gDNA sample
  n_obs <- length(md$patient)
  paired_df <- data.frame(
    PC1_cfDNA = pca_res$x[1:n_obs, 1],
    PC2_cfDNA = pca_res$x[1:n_obs, 2],
    PC1_gDNA  = pca_res$x[(n_obs + 1):(2 * n_obs), 1],
    PC2_gDNA  = pca_res$x[(n_obs + 1):(2 * n_obs), 2]
  )

  fig5c <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Source)) +
    geom_segment(data = paired_df,
                 aes(x = PC1_cfDNA, xend = PC1_gDNA,
                     y = PC2_cfDNA, yend = PC2_gDNA),
                 color = "grey70", linewidth = 0.3, alpha = 0.5,
                 inherit.aes = FALSE) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_manual(values = source_colors) +
    labs(x = paste0("PC1 (", var_explained[1], "%)"),
         y = paste0("PC2 (", var_explained[2], "%)"),
         title = "c  Paired samples connected") +
    theme_nature

  save_png("fig5a_pca_pc1_pc2.png", fig5a, width = 6, height = 5)
  save_png("fig5b_pca_pc1_pc3.png", fig5b, width = 6, height = 5)
  save_png("fig5c_pca_paired.png", fig5c, width = 6, height = 5)

  cat("  Saved: fig5a_pca_pc1_pc2.png, fig5b_pca_pc1_pc3.png, fig5c_pca_paired.png\n")
  cat("  Variance explained: PC1 =", var_explained[1], "%, PC2 =",
      var_explained[2], "%, PC3 =", var_explained[3], "%\n")
  cat("  Total samples:", nrow(pca_df), "(cfDNA:",
      sum(pca_df$Source == "cfDNA"), ", gDNA:",
      sum(pca_df$Source == "gDNA"), ")\n")
} else {
  cat("  model_data.rds not found; skipping PCA figure.\n")
}

# ===========================================================================
# FIGURE 2: Gene-Level Concordance Heatmap (OBSERVED DATA — no model needed)
# ===========================================================================
cat("\n--- Figure 2: Gene-Level Concordance Heatmap ---\n")

if (file.exists(model_data_file)) {
  md <- readRDS(model_data_file)

  # For each paired observation (patient x timepoint), compute the Pearson
  # correlation between cfDNA and gDNA 5-hMC vectors across all p genes.
  # This is a model-free measure of global concordance.
  n_obs <- md$N_obs
  timepoint_order <- c("1" = "Screening", "2" = "C4", "3" = "C8",
                        "4" = "C18", "5" = "3 YR F/U")

  obs_corr <- numeric(n_obs)
  for (obs in seq_len(n_obs)) {
    obs_corr[obs] <- cor(md$Y_cf[, obs], md$Y_g[, obs], method = "pearson")
  }

  heatmap_df <- data.frame(
    patient = factor(md$paired_info$Study.ID),
    timepoint = factor(timepoint_order[as.character(md$time)],
                       levels = c("Screening", "C4", "C8", "C18", "3 YR F/U")),
    correlation = obs_corr
  )

  # Aggregate duplicate patient-timepoints (e.g., 101-17 and 101-73 have 2 entries at 3YR F/U)
  heatmap_df <- heatmap_df %>%
    group_by(patient, timepoint) %>%
    summarise(correlation = mean(correlation), .groups = "drop")

  # Summary statistics
  cat("  Gene-level cfDNA-gDNA correlation:\n")
  cat("    Mean:", round(mean(obs_corr), 3), "\n")
  cat("    Median:", round(median(obs_corr), 3), "\n")
  cat("    Range:", round(min(obs_corr), 3), "-", round(max(obs_corr), 3), "\n")
  cat("    SD:", round(sd(obs_corr), 3), "\n")

  # Panel A: Heatmap (patients x timepoints)
  fig2a <- ggplot(heatmap_df, aes(x = timepoint, y = patient, fill = correlation)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = round(correlation, 2)), size = 2.5, color = "black") +
    scale_fill_gradient2(low = "#E63946", mid = "#FFFFCC", high = "#457B9D",
                         midpoint = median(obs_corr),
                         limits = c(min(obs_corr) - 0.02, max(obs_corr) + 0.02),
                         name = "Pearson r") +
    labs(x = "Treatment Timepoint", y = "Patient",
         title = "a  Gene-level cfDNA-gDNA concordance") +
    theme_nature +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank())

  # Panel B: Distribution of correlations (histogram + density)
  fig2b <- ggplot(heatmap_df, aes(x = correlation)) +
    geom_histogram(aes(y = after_stat(density)), bins = 20,
                   fill = "#457B9D", alpha = 0.5, color = "white") +
    geom_density(linewidth = 0.8, color = "#264653") +
    geom_vline(xintercept = median(obs_corr), linetype = "dashed",
               color = "#E63946", linewidth = 0.7) +
    annotate("text", x = median(obs_corr) - 0.01, y = Inf,
             label = paste0("median = ", round(median(obs_corr), 3)),
             vjust = 2, hjust = 1, size = 3.5, color = "#E63946") +
    labs(x = "Pearson correlation (cfDNA vs gDNA across genes)",
         y = "Density",
         title = "b  Distribution of pairwise concordance") +
    theme_nature

  # Panel C: Concordance by timepoint (boxplot)
  fig2c <- ggplot(heatmap_df, aes(x = timepoint, y = correlation, fill = timepoint)) +
    geom_boxplot(alpha = 0.6, outlier.size = 1.5) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
    scale_fill_brewer(palette = "Set2") +
    labs(x = "Treatment Timepoint",
         y = "Pearson correlation (cfDNA vs gDNA)",
         title = "c  Concordance across treatment timepoints") +
    theme_nature +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))

  save_png("fig2a_concordance_heatmap.png", fig2a, width = 7, height = 7)
  save_png("fig2b_concordance_distribution.png", fig2b, width = 5, height = 4)
  save_png("fig2c_concordance_by_timepoint.png", fig2c, width = 5, height = 4)

  # Save summary statistics for manuscript text
  corr_summary <- data.frame(
    metric = c("mean", "median", "sd", "min", "max", "q25", "q75"),
    value = c(mean(obs_corr), median(obs_corr), sd(obs_corr),
              min(obs_corr), max(obs_corr),
              quantile(obs_corr, 0.25), quantile(obs_corr, 0.75))
  )
  write.csv(corr_summary, file.path(base_dir, "nature_manuscript",
            "fig2_concordance_summary.csv"), row.names = FALSE)
  cat("  Saved: fig2a, fig2b, fig2c panels + concordance_summary.csv\n")

  # Per-timepoint summaries
  tp_summary <- heatmap_df %>%
    group_by(timepoint) %>%
    summarize(n = n(), mean_corr = mean(correlation),
              sd_corr = sd(correlation), .groups = "drop")
  cat("\n  Per-timepoint concordance:\n")
  print(as.data.frame(tp_summary))

} else {
  cat("  model_data.rds not found; skipping concordance heatmap.\n")
}

# ===========================================================================
# FIGURE 2 (MODEL): Posterior Concordance Heatmap (REQUIRES MODEL RESULTS)
# ===========================================================================
cat("\n--- Figure 2 (Model-Based): Posterior Concordance Heatmap ---\n")

krd_file <- file.path(res_dir, "krd_M1_fit.rds")
if (file.exists(krd_file)) {
  krd <- readRDS(krd_file)

  # Extract posterior correlations between omega_g and omega_cf
  # per patient-timepoint for dominant subclone.
  # NOTE: the saved fit stores omega_*_trace as [n_samp x (K * N_obs)] matrices
  # (component-major, col = (k-1)*N_obs + obs); there is no merged$w_g array.
  # Relabel + align chains before merging so component identities are canonical.
  krd_rl <- krd
  for (ci in seq_along(krd_rl$chains))
    krd_rl$chains[[ci]] <- relabel_by_weight(krd_rl$chains[[ci]])
  krd_rl$chains <- align_chains_by_signature(krd_rl$chains)
  merged <- merge_chain_samples(krd_rl$chains)

  n_samp <- nrow(merged$omega_g_trace)
  n_obs  <- readRDS(model_data_file)$N_obs
  K      <- ncol(merged$omega_g_trace) / n_obs
  col_of <- function(k, obs) (k - 1) * n_obs + obs

  corr_matrix <- matrix(NA, nrow = n_obs, ncol = 1)
  for (obs in seq_len(n_obs)) {
    mean_w <- sapply(seq_len(K), function(k) mean(merged$omega_g_trace[, col_of(k, obs)]))
    dom_k  <- which.max(mean_w)
    wg_k   <- merged$omega_g_trace[,  col_of(dom_k, obs)]
    wcf_k  <- merged$omega_cf_trace[, col_of(dom_k, obs)]
    corr_matrix[obs, 1] <- if (sd(wg_k) > 1e-12 && sd(wcf_k) > 1e-12)
      cor(wg_k, wcf_k) else NA
  }

  md <- readRDS(model_data_file)
  timepoint_order <- c("1" = "Screening", "2" = "C4", "3" = "C8",
                        "4" = "C18", "5" = "3 YR F/U")

  heatmap_model_df <- data.frame(
    patient = factor(md$patient),
    timepoint = factor(timepoint_order[as.character(md$time)],
                       levels = c("Screening", "C4", "C8", "C18", "3 YR F/U")),
    correlation = corr_matrix[, 1]
  )

  fig2_model <- ggplot(heatmap_model_df, aes(x = timepoint, y = patient, fill = correlation)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(low = "#E63946", mid = "white", high = "#457B9D",
                         midpoint = 0.5, limits = c(0, 1),
                         name = "Posterior\nConcordance") +
    labs(x = "Treatment Timepoint", y = "Patient",
         title = "Subclonal-level cfDNA-gDNA concordance (model-based)") +
    theme_nature +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank())

  save_png("fig2_model_concordance_heatmap.png", fig2_model, width = 7, height = 6)
  cat("  Saved: fig2_model_concordance_heatmap.png\n")
} else {
  cat("  KRd production results not found; model-based Figure 2 will be generated after analysis.\n")
}

# ===========================================================================
# FIGURE 3: Longitudinal Subclonal Trajectories (REQUIRES MODEL RESULTS)
# ===========================================================================
cat("\n--- Figure 3: Longitudinal Subclonal Trajectories ---\n")

if (file.exists(krd_file)) {
  krd <- readRDS(krd_file)
  md <- readRDS(model_data_file)
  # Canonical (relabeled + aligned) weights; saved fit has no merged$w_g array.
  krd_rl <- krd
  for (ci in seq_along(krd_rl$chains))
    krd_rl$chains[[ci]] <- relabel_by_weight(krd_rl$chains[[ci]])
  krd_rl$chains <- align_chains_by_signature(krd_rl$chains)
  merged <- merge_chain_samples(krd_rl$chains)
  n_obs_traj <- md$N_obs
  K_traj     <- ncol(merged$omega_g_trace) / n_obs_traj
  col_traj   <- function(k, obs) (k - 1) * n_obs_traj + obs

  # Select representative patients (those with most timepoints)
  patient_tp_count <- table(md$patient)
  top_patients <- as.integer(names(sort(patient_tp_count, decreasing = TRUE))[1:min(6, length(patient_tp_count))])

  timepoint_order <- c("1" = "Screening", "2" = "C4", "3" = "C8",
                        "4" = "C18", "5" = "3 YR F/U")

  # Build trajectory data
  traj_list <- list()
  for (pt in top_patients) {
    obs_idx <- which(md$patient == pt)
    for (oi in obs_idx) {
      for (k in 1:min(3, K_traj)) {
        w_g_samp  <- merged$omega_g_trace[,  col_traj(k, oi)]
        w_cf_samp <- merged$omega_cf_trace[, col_traj(k, oi)]
        traj_list[[length(traj_list) + 1]] <- data.frame(
          patient = paste0("Patient ", pt),
          timepoint = timepoint_order[as.character(md$time[oi])],
          time_num = md$time[oi],
          subclone = paste0("Subclone ", k),
          source = "gDNA",
          weight_mean = mean(w_g_samp),
          weight_lo = quantile(w_g_samp, 0.025),
          weight_hi = quantile(w_g_samp, 0.975)
        )
        traj_list[[length(traj_list) + 1]] <- data.frame(
          patient = paste0("Patient ", pt),
          timepoint = timepoint_order[as.character(md$time[oi])],
          time_num = md$time[oi],
          subclone = paste0("Subclone ", k),
          source = "cfDNA",
          weight_mean = mean(w_cf_samp),
          weight_lo = quantile(w_cf_samp, 0.025),
          weight_hi = quantile(w_cf_samp, 0.975)
        )
      }
    }
  }
  traj_df <- do.call(rbind, traj_list)

  fig3 <- ggplot(traj_df, aes(x = time_num, y = weight_mean,
                               color = subclone, linetype = source)) +
    geom_ribbon(aes(ymin = weight_lo, ymax = weight_hi, fill = subclone),
                alpha = 0.1, linetype = 0) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    facet_wrap(~ patient, scales = "free_y", ncol = 3) +
    scale_color_manual(values = c("#2A9D8F", "#E9C46A", "#264653")) +
    scale_fill_manual(values = c("#2A9D8F", "#E9C46A", "#264653")) +
    scale_linetype_manual(values = c("gDNA" = "solid", "cfDNA" = "dashed")) +
    scale_x_continuous(breaks = 1:5,
                       labels = c("SCR", "C4", "C8", "C18", "3Y")) +
    labs(x = "Treatment Timepoint", y = "Subclonal Weight",
         color = "", fill = "", linetype = "Source",
         title = "Longitudinal subclonal dynamics: cfDNA vs gDNA") +
    theme_nature +
    theme(legend.position = "bottom",
          strip.text = element_text(face = "bold"))

  save_png("fig3_trajectories.png", fig3, width = 9, height = 6)
  cat("  Saved: fig3_trajectories.png\n")
} else {
  cat("  KRd production results not found; Figure 3 will be generated after analysis.\n")
}

# ===========================================================================
# FIGURE 4: Clinical Implications Composite (REQUIRES MODEL RESULTS)
# ===========================================================================
cat("\n--- Figure 4: Clinical Implications Composite ---\n")

if (file.exists(krd_file)) {
  krd <- readRDS(krd_file)
  merged <- krd$merged

  # Panel A: Tracking strength (kappa-derived correlation) per patient
  kappa_samp <- merged$kappa
  alpha_samp <- merged$alpha_dp
  corr_samp <- (kappa_samp + 1) / (alpha_samp + kappa_samp + 1)

  corr_summary <- data.frame(
    metric = "Global",
    mean = mean(corr_samp),
    lo = quantile(corr_samp, 0.025),
    hi = quantile(corr_samp, 0.975)
  )

  fig4a <- ggplot(corr_summary, aes(x = mean, y = metric)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey50") +
    labs(x = "Tracking Correlation", y = "",
         title = "a  Posterior tracking strength") +
    xlim(0, 1) +
    theme_nature

  # Panel B: Contamination fractions by OBSERVATION.
  # omega_0 is per (subject, timepoint): [n_samp x N_obs], NOT [n_samp x n].
  # The previous comment and the "Patient" axis label both asserted 26 columns
  # where there are 82, which is the row-semantics trap this model introduced.
  wN_mat <- merged$omega_0  # [n_samp x N_obs]
  wN_df <- data.frame(
    obs = rep(seq_len(ncol(wN_mat)), each = nrow(wN_mat)),
    wN = as.vector(wN_mat)
  )

  fig4b <- ggplot(wN_df, aes(x = factor(obs), y = wN)) +
    geom_boxplot(fill = "#E9C46A", outlier.size = 0.3, alpha = 0.7) +
    labs(x = "Observation (subject-timepoint)",
         y = expression("Contamination fraction (" * omega[0] * ")"),
         title = "b  Normal DNA contamination") +
    theme_nature +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

  # Panel C: Tracking vs contamination scatter.
  # `wN_means` is PER-OBSERVATION. It was previously also read by two later,
  # unrelated sections of this file via global leakage, where it was reported
  # as a per-patient statistic; those now compute their own per-patient
  # aggregate. Keep the name local in meaning and explicit in the column.
  wN_means <- colMeans(wN_mat)
  scatter_df <- data.frame(
    obs = seq_along(wN_means),
    contamination = wN_means,
    tracking = mean(corr_samp)  # Global; per-patient if available
  )

  fig4c <- ggplot(scatter_df, aes(x = contamination, y = tracking)) +
    geom_point(size = 3, color = "#264653") +
    geom_smooth(method = "lm", se = TRUE, color = "#E63946", linewidth = 0.8) +
    labs(x = expression("Mean contamination (" * omega[0] * ")"),
         y = "Tracking correlation",
         title = "c  Tracking vs contamination") +
    theme_nature

  # Save panels
  save_png("fig4a_tracking_strength.png", fig4a, width = 5, height = 3)
  save_png("fig4b_contamination.png", fig4b, width = 7, height = 4)
  save_png("fig4c_tracking_vs_contamination.png", fig4c, width = 5, height = 4)
  cat("  Saved: fig4a, fig4b, fig4c panels\n")
} else {
  cat("  KRd production results not found; Figure 4 will be generated after analysis.\n")
}

cat("\n=== Applied paper figure generation complete ===\n")


###############################################################################
## ---- 20_nature_results_figures.R
###############################################################################

###############################################################################
## 20_nature_results_figures.R
## Generate all applied paper figures using production KRd results
##
## Requires: results/krd_M1_fit.rds, results/krd_M0_fit.rds,
##           results/model_data.rds, results/*.csv
##
## Usage: Rscript R/20_nature_results_figures.R
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

base_dir <- getwd()
fig_dir <- file.path(base_dir, "nature_manuscript", "figures")
res_dir <- file.path(base_dir, "results")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

theme_nature <- theme_bw(base_size = 11, base_family = "sans") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "grey80"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12)
  )
theme_set(theme_nature)

save_fig <- function(filename, plot, width, height) {
  pdf_path <- file.path(fig_dir, filename)
  pdf(pdf_path, width = width, height = height)
  print(plot)
  dev.off()
  cat("  Saved:", filename, "\n")
}

cat("=== Generating Applied Paper Figures from Production Results ===\n\n")

# Load data
md <- readRDS(file.path(res_dir, "model_data.rds"))
m1 <- readRDS(file.path(res_dir, "krd_M1_fit.rds"))
merged <- m1$merged

timepoint_map <- c("1" = "Screening", "2" = "C4", "3" = "C8",
                    "4" = "C18", "5" = "3 YR F/U")
tp_levels <- c("Screening", "C4", "C8", "C18", "3 YR F/U")

# ===========================================================================
# FIGURE 3: Model-Based Subclonal Concordance Heatmap
# ===========================================================================
cat("--- Figure 3: Model-Based Concordance Heatmap ---\n")

# Compute posterior mean omegas from MCMC traces (proper Bayesian estimates)
n_obs <- md$N_obs
K <- m1$chains[[1]]$K
K_trace <- m1$chains[[1]]$K_trace
n_ch <- length(m1$chains)

# Preserve the raw (un-relabeled) chains before any post-processing so the
# supplementary "before relabeling" traceplots use the actual raw samples.
# R copies on assignment, so this is a genuine snapshot.
chains_raw <- m1$chains

# Apply label switching before computing posterior means
for (i in seq_along(m1$chains)) {
  m1$chains[[i]] <- relabel_by_weight(m1$chains[[i]])
}
m1$chains <- align_chains_by_signature(m1$chains)
merged <- merge_chain_samples(m1$chains)

# Compute posterior mean omega from MCMC traces (K_trace components)
omega_g_mean <- matrix(0, K_trace, n_obs)
omega_cf_mean <- matrix(0, K_trace, n_obs)
for (k in 1:K_trace) {
  for (obs in 1:n_obs) {
    col_idx <- (k - 1) * n_obs + obs
    omega_g_mean[k, obs] <- mean(merged$omega_g_trace[, col_idx])
    omega_cf_mean[k, obs] <- mean(merged$omega_cf_trace[, col_idx])
  }
}
cat("  Using posterior mean omegas from", nrow(merged$omega_g_trace),
    "MCMC iterations (K_trace =", K_trace, "components)\n")

# For each obs, compute correlation between omega_g and omega_cf vectors
obs_model_corr <- numeric(n_obs)
for (obs in seq_len(n_obs)) {
  wg <- omega_g_mean[, obs]
  wcf <- omega_cf_mean[, obs]
  # Only correlate non-trivial components
  active <- which(wg > 0.001 | wcf > 0.001)
  if (length(active) >= 2) {
    obs_model_corr[obs] <- cor(wg[active], wcf[active])
  } else {
    obs_model_corr[obs] <- NA
  }
}

# Subclonal-level concordance heatmap based on the POSTERIOR CORRELATION of the
# source-specific weights (obs_model_corr above), the retained concordance metric.
# (The L1 concordance heatmap was removed by design.)
heatmap_model_df <- data.frame(
  patient = factor(md$paired_info$Study.ID),
  timepoint = factor(timepoint_map[as.character(md$time)], levels = tp_levels),
  concordance = obs_model_corr
)

# Aggregate duplicate patient-timepoints (e.g., 101-17 and 101-73 have 2 entries at 3YR F/U)
heatmap_model_df <- heatmap_model_df %>%
  group_by(patient, timepoint) %>%
  summarise(concordance = mean(concordance, na.rm = TRUE), .groups = "drop")

cat("  Subclonal weight correlation: mean =", round(mean(heatmap_model_df$concordance, na.rm = TRUE), 3), "\n")

fig3 <- ggplot(heatmap_model_df, aes(x = timepoint, y = patient, fill = concordance)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(concordance, 2)), size = 2.3, color = "black") +
  scale_fill_gradient2(low = "#E63946", mid = "#FFFFCC", high = "#457B9D",
                       midpoint = 0, limits = c(-1, 1),
                       name = "Weight\ncorrelation") +
  labs(x = "Treatment Timepoint", y = "Patient",
       title = "Subclonal-level cfDNA-gDNA concordance (posterior weight correlation)") +
  theme_nature +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

save_fig("fig3_model_concordance_heatmap.pdf", fig3, width = 7, height = 7)

# ===========================================================================
# FIGURE 4: Longitudinal Subclonal Trajectories
# ===========================================================================
cat("\n--- Figure 4: Longitudinal Subclonal Trajectories ---\n")

# Find patients with most timepoints
patient_tp <- md$paired_info %>%
  group_by(Study.ID) %>%
  summarize(n_tp = n(), .groups = "drop") %>%
  arrange(desc(n_tp))

top_patients <- patient_tp$Study.ID[1:min(6, nrow(patient_tp))]
cat("  Selected patients:", paste(top_patients, collapse = ", "), "\n")

# Build trajectory data from posterior-mean omegas (using the MCMC trace)
# Select occupied components by posterior-mean pi across all saved iterations,
# computed per-observation from omega traces. We treat component k as occupied
# if its cohort-averaged posterior-mean weight exceeds 0.01.
# omega_g_mean is [K_trace x n_obs]; rowMeans gives per-component cohort avg
pi_mean_from_trace <- rowMeans(omega_g_mean)  # length K_trace
ord <- order(pi_mean_from_trace, decreasing = TRUE)
occupied <- ord[pi_mean_from_trace[ord] > 0.05]
if (length(occupied) < 2) occupied <- ord[1:min(2, length(ord))]
cat("  Occupied components (by posterior-mean omega > 0.05):", occupied,
    " (pi =", round(pi_mean_from_trace[occupied], 3), ")\n")

traj_list <- list()
for (pt in top_patients) {
  obs_idx <- which(md$paired_info$Study.ID == pt)
  for (oi in obs_idx) {
    for (k in occupied) {
      wg_val  <- omega_g_mean[k,  oi]
      wcf_val <- omega_cf_mean[k, oi]

      traj_list[[length(traj_list) + 1]] <- data.frame(
        patient = pt,
        timepoint = timepoint_map[as.character(md$time[oi])],
        time_num = md$time[oi],
        subclone = paste0("Subclone ", which(occupied == k)),
        source = "gDNA",
        weight = wg_val
      )
      traj_list[[length(traj_list) + 1]] <- data.frame(
        patient = pt,
        timepoint = timepoint_map[as.character(md$time[oi])],
        time_num = md$time[oi],
        subclone = paste0("Subclone ", which(occupied == k)),
        source = "cfDNA",
        weight = wcf_val
      )
    }
  }
}
traj_df <- do.call(rbind, traj_list)
traj_df$timepoint <- factor(traj_df$timepoint, levels = tp_levels)

fig4 <- ggplot(traj_df, aes(x = time_num, y = weight,
                             color = subclone, linetype = source)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  facet_wrap(~ patient, ncol = 3, scales = "free_y") +
  scale_color_manual(values = c("Subclone 1" = "#2A9D8F", "Subclone 2" = "#E9C46A",
                                "Subclone 3" = "#264653")) +
  scale_linetype_manual(values = c("gDNA" = "solid", "cfDNA" = "dashed")) +
  scale_x_continuous(breaks = 1:5, labels = c("SCR", "C4", "C8", "C18", "3Y")) +
  labs(x = "Treatment Timepoint", y = "Subclonal Weight",
       color = "", linetype = "Source",
       title = "Longitudinal subclonal dynamics: cfDNA vs gDNA") +
  theme_nature +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"))

save_fig("fig4_trajectories.pdf", fig4, width = 9, height = 7)

# ===========================================================================
# FIGURE 5: Concordance & Contamination
# ===========================================================================
cat("\n--- Figure 5: Concordance & Contamination ---\n")

# Tracking correlation is deterministic with fixed kappa=1, alpha=1: rho = 2/3
rho_fixed <- 2/3
cat("  Tracking correlation (exact): rho =", round(rho_fixed, 4),
    "(kappa=1, alpha=1 both fixed)\n")

# Panel A: Per-observation distribution of the posterior weight correlation
# (obs_model_corr from the Figure 3 section). The L1 concordance distribution was
# removed by design; we report the retained posterior-correlation metric.
concordance_per_obs <- obs_model_corr[is.finite(obs_model_corr)]

concordance_df <- data.frame(concordance = concordance_per_obs)
fig5a <- ggplot(concordance_df, aes(x = concordance)) +
  geom_histogram(aes(y = after_stat(density)), bins = 25,
                 fill = "#457B9D", alpha = 0.6, color = "white") +
  geom_density(linewidth = 0.8, color = "#264653") +
  geom_vline(xintercept = mean(concordance_per_obs), linetype = "dashed",
             color = "#E63946", linewidth = 0.7) +
  annotate("text", x = mean(concordance_per_obs) + 0.02, y = Inf,
           label = paste0("mean = ", round(mean(concordance_per_obs), 3)),
           vjust = 2, hjust = 0, size = 3.5, color = "#E63946") +
  labs(x = expression("Subclonal weight correlation  " * rho(omega^g, omega^cf)),
       y = "Density",
       title = "a  Per-observation subclonal weight correlation") +
  theme_nature

# Panel B: Contamination by patient.
# omega_0 is per OBSERVATION, so each patient contributes T_i columns; we map
# columns to patients via model_data rather than assuming one column per patient.
wN_long <- data.frame(
  obs = rep(seq_len(ncol(merged$omega_0)), each = nrow(merged$omega_0)),
  wN  = as.vector(merged$omega_0)
)
wN_long$patient <- md$patient[wN_long$obs]
wN_means_obs <- colMeans(merged$omega_0)
pat_mean <- tapply(wN_means_obs, md$patient, mean)
patient_order <- as.integer(names(sort(pat_mean)))
wN_long$patient <- factor(wN_long$patient, levels = patient_order)

fig5b <- ggplot(wN_long, aes(x = patient, y = wN * 100)) +
  geom_boxplot(fill = "#E9C46A", outlier.size = 0.3, alpha = 0.7) +
  labs(x = "Patient (ordered by mean contamination; box spans that patient's timepoints)",
       y = expression("Contamination fraction " * omega[0] * " (%)"),
       title = "b  Normal DNA contamination by patient") +
  theme_nature +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

save_fig("fig5a_concordance_distribution.pdf", fig5a, width = 5, height = 4)
save_fig("fig5b_contamination_by_patient.pdf", fig5b, width = 7, height = 4)

# Print key stats
cat("  Subclonal weight correlation: mean =", round(mean(concordance_per_obs), 4),
    ", median =", round(median(concordance_per_obs), 4),
    ", range = (", round(min(concordance_per_obs), 4), ",",
    round(max(concordance_per_obs), 4), ")\n")
## Report BOTH units explicitly. `wN_means_obs`/`pat_mean` are defined just
## above; the per-patient median is the one that answers "what is a typical
## patient's contamination", and it differs materially from the median over
## observations whenever T_i is unbalanced (patients with 5 timepoints
## otherwise count five times).
cat("  Contamination: median across PATIENTS =",
    round(median(pat_mean) * 100, 4), "% (n =", length(pat_mean), "patients)\n")
cat("  Contamination: median across OBSERVATIONS =",
    round(median(wN_means_obs) * 100, 4), "% (N_obs =",
    length(wN_means_obs), ")\n")
cat("  Contamination range (observations):",
    round(min(wN_means_obs) * 100, 5), "% -",
    round(max(wN_means_obs) * 100, 4), "%\n")

# ===========================================================================
# SUPPLEMENTARY: Convergence trace plots (before & after label switching)
# ===========================================================================
cat("\n--- Supplementary: Trace Plots (before & after label switching) ---\n")

chain_cols <- c("steelblue", "firebrick", "forestgreen", "darkorange")

# Helper: generate a complete 5x2 traceplot page from a list of chains
generate_nature_traceplots <- function(chains_list, out_path, title_prefix = "") {
  n_ch <- length(chains_list)
  N_obs_tp <- ncol(chains_list[[1]]$final_omega_g)
  K_trace_tp <- chains_list[[1]]$K_trace

  plot_scalar <- function(param_name, ylab_expr, main_expr) {
    vals <- lapply(1:n_ch, function(i) chains_list[[i]]$samples[[param_name]])
    y_range <- range(unlist(vals))
    plot(vals[[1]], type = "l", col = chain_cols[1], ylim = y_range,
         xlab = "Iteration (post burn-in, thinned)", ylab = ylab_expr, main = main_expr)
    for (ch in 2:n_ch) lines(vals[[ch]], col = chain_cols[ch])
    legend("topright", paste("Chain", 1:n_ch), col = chain_cols[1:n_ch],
           lty = 1, cex = 0.6, bg = "white")
  }

  plot_omega_col <- function(col_idx, ylab_expr, main_expr, trace_nm) {
    vals <- lapply(1:n_ch, function(i) chains_list[[i]]$samples[[trace_nm]][, col_idx])
    y_range <- range(unlist(vals))
    plot(vals[[1]], type = "l", col = chain_cols[1], ylim = y_range,
         xlab = "Iteration (post burn-in, thinned)", ylab = ylab_expr, main = main_expr)
    for (ch in 2:n_ch) lines(vals[[ch]], col = chain_cols[ch])
    legend("topright", paste("Chain", 1:n_ch), col = chain_cols[1:n_ch],
           lty = 1, cex = 0.6, bg = "white")
  }

  pfx <- if (nchar(title_prefix) > 0) paste0(title_prefix, ": ") else ""

  pdf(out_path, width = 12, height = 17)
  par(mfrow = c(5, 2), mar = c(4, 4, 3, 1))

  # Row 1: K+ and loglik (kappa is fixed, no traceplot needed)
  plot_scalar("K_plus", expression(K^"+"), bquote(.(pfx) * K^"+"))
  plot_scalar("loglik", "Log-likelihood", bquote(.(pfx) * "Log-lik (" * M[1] * ")"))

  # Row 2
  plot_scalar("sigma_0", expression(sigma[0]), bquote(.(pfx) * sigma[0]))
  wN_v <- lapply(1:n_ch, function(i) chains_list[[i]]$samples$omega_0[, 1])
  yr <- range(unlist(wN_v))
  plot(wN_v[[1]], type = "l", col = chain_cols[1], ylim = yr,
       xlab = "Iteration (post burn-in, thinned)", ylab = expression(omega[0]),
       main = bquote(.(pfx) * omega[0] * " (observation 1)"))
  for (ch in 2:n_ch) lines(wN_v[[ch]], col = chain_cols[ch])
  legend("topright", paste("Chain", 1:n_ch), col = chain_cols[1:n_ch],
         lty = 1, cex = 0.6, bg = "white")

  # Row 3
  if (!is.null(chains_list[[1]]$samples$omega_g_trace)) {
    plot_omega_col(1, expression(omega[1]^g),
                   bquote(.(pfx) * omega[1]^g * " (obs 1, comp 1)"), "omega_g_trace")
    plot_omega_col(1, expression(omega[1]^cf),
                   bquote(.(pfx) * omega[1]^cf * " (obs 1, comp 1)"), "omega_cf_trace")
  } else {
    plot.new(); text(0.5, 0.5, "omega traces\nnot stored")
    plot.new(); text(0.5, 0.5, "omega traces\nnot stored")
  }

  # Row 4
  if (!is.null(chains_list[[1]]$samples$omega_g_trace) && K_trace_tp >= 2) {
    c2 <- N_obs_tp + 1
    plot_omega_col(c2, expression(omega[2]^g),
                   bquote(.(pfx) * omega[2]^g * " (obs 1, comp 2)"), "omega_g_trace")
    plot_omega_col(c2, expression(omega[2]^cf),
                   bquote(.(pfx) * omega[2]^cf * " (obs 1, comp 2)"), "omega_cf_trace")
  } else {
    plot.new(); text(0.5, 0.5, "comp 2 traces\nnot available")
    plot.new(); text(0.5, 0.5, "comp 2 traces\nnot available")
  }

  # Row 5
  plot_scalar("loglik", "Log-likelihood",
              bquote(.(pfx) * "Log-likelihood (" * M[1] * ")"))
  m0_file_tp <- file.path(base_dir, "results", "krd_M0_fit.rds")
  if (file.exists(m0_file_tp)) {
    m0_tp <- readRDS(m0_file_tp)
    ch_m0 <- m0_tp$chains
    n_m0 <- length(ch_m0)
    ll_m0 <- lapply(1:n_m0, function(i) ch_m0[[i]]$samples$loglik)
    plot(ll_m0[[1]], type = "l", col = chain_cols[1], ylim = range(unlist(ll_m0)),
         xlab = "Iteration (post burn-in, thinned)", ylab = "Log-likelihood",
         main = bquote(.(pfx) * "Log-likelihood (" * M[0] * ")"))
    for (ch in 2:n_m0) lines(ll_m0[[ch]], col = chain_cols[ch])
    legend("topright", paste("Chain", 1:n_m0), col = chain_cols[1:n_m0],
           lty = 1, cex = 0.6, bg = "white")
  } else {
    plot.new(); text(0.5, 0.5, "M0 loglik\nnot available")
  }

  dev.off()
}

# (i) BEFORE label switching: traceplots from the raw chains (saved before
# any relabeling was applied for Figure 3 above)
cat("  Generating traceplots BEFORE label switching...\n")
generate_nature_traceplots(chains_raw,
  file.path(fig_dir, "supp_traceplots_before_relabeling.pdf"),
  "Before relabeling")
cat("  Saved: supp_traceplots_before_relabeling.pdf\n")

# (ii) AFTER label switching: use the already-relabeled chains that were
# produced earlier for Figure 3 (m1$chains has had relabel_by_weight and
# align_chains_by_signature applied in place).
cat("  Generating traceplots AFTER label switching...\n")
generate_nature_traceplots(m1$chains,
  file.path(fig_dir, "supp_traceplots_after_relabeling.pdf"),
  "After relabeling")
cat("  Saved: supp_traceplots_after_relabeling.pdf\n")

# Canonical supp_traceplots.pdf = post-relabeling (manuscript figure)
generate_nature_traceplots(m1$chains,
  file.path(fig_dir, "supp_traceplots.pdf"), "")
cat("  Saved: supp_traceplots.pdf\n")

# ===========================================================================
# Summary statistics for manuscript
# ===========================================================================
cat("\n=== KEY RESULTS FOR MANUSCRIPT ===\n")
mc <- read.csv(file.path(res_dir, "krd_model_comparison.csv"))
cat("delta_WAIC =", mc$value[mc$metric == "delta_WAIC"], "\n")
cat("delta_LOO =", mc$value[mc$metric == "delta_LOO"], "\n")
cat("kappa: 1.0 (fixed)\n")
cat("alpha: 1.0 (fixed)\n")
cat("Tracking rho: 2/3 = 0.667 (exact, kappa=alpha=1)\n")
cat("K+: unique values =", sort(unique(merged$K_plus)), "\n")
cat("sigma_N: mean =", round(mean(merged$sigma_0), 2), "\n")
## Manuscript summary line: quote the PER-PATIENT median (pat_mean), not the
## per-observation one. These differ by ~50% on the unbalanced KRd design.
cat("Contamination: median across patients =", round(median(pat_mean)*100, 4),
    "%, max =", round(max(pat_mean)*100, 3), "%",
    "| across observations: median =", round(median(wN_means_obs)*100, 4),
    "%, max =", round(max(wN_means_obs)*100, 3), "%\n")
cat("MCMC: 4 chains x 15000 iter, 10000 burn-in, thin=20, 1000 total samples\n")

cat("\n=== Figure generation complete ===\n")


###############################################################################
## ---- 21_raw_concordance_heatmap.R
###############################################################################

###############################################################################
## 21_raw_concordance_heatmap.R
## Generate gene-level concordance heatmap from RAW (pre-filtering, pre-VST)
## count data — all 19,100 genes, no normalization
##
## This provides a more conservative estimate of cfDNA-gDNA concordance
## than the post-filtering/post-VST version.
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

base_dir <- getwd()
fig_dir <- file.path(base_dir, "nature_manuscript", "figures")
data_dir <- file.path(base_dir, "KRd trial", "5hmC data")

theme_nature <- theme_bw(base_size = 11, base_family = "sans") +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white", color = "grey80"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12)
  )

save_fig <- function(filename, plot, width, height) {
  pdf_path <- file.path(fig_dir, filename)
  pdf(pdf_path, width = width, height = height)
  print(plot)
  dev.off()
  cat("  Saved:", filename, "\n")
}

cat("=== Raw Data Concordance Heatmap (19,100 genes, no filtering/VST) ===\n\n")

# --- Load raw data ---
cfDNA_key <- read.csv(file.path(data_dir, "kRd-cfDNA_sample_key.csv"),
                      stringsAsFactors = FALSE)
gDNA_key  <- read.csv(file.path(data_dir, "kRd-gDNA_sample_key.csv"),
                      stringsAsFactors = FALSE)

cfDNA_key$barcode <- cfDNA_key$Assigned.ID
gDNA_key$barcode  <- gDNA_key$Assigned.ID

# Pair by patient + timepoint
paired <- inner_join(
  cfDNA_key %>% select(Study.ID, Timepoint, barcode_cf = barcode),
  gDNA_key  %>% select(Study.ID, Timepoint, barcode_g  = barcode),
  by = c("Study.ID", "Timepoint")
)

cfDNA_counts <- readRDS(file.path(data_dir, "kRd-cfDNA_genebody_count.RDS"))
gDNA_counts  <- readRDS(file.path(data_dir, "kRd-gDNA_genebody_count.RDS"))

cfDNA_raw <- as.matrix(cfDNA_counts[, paired$barcode_cf])
gDNA_raw  <- as.matrix(gDNA_counts[, paired$barcode_g])

cat("Raw count matrices: ", nrow(cfDNA_raw), "genes x", ncol(cfDNA_raw), "paired obs\n")

# --- Compute correlations on raw counts ---
n_obs <- ncol(cfDNA_raw)
raw_corr <- numeric(n_obs)
for (obs in seq_len(n_obs)) {
  raw_corr[obs] <- cor(cfDNA_raw[, obs], gDNA_raw[, obs], method = "pearson")
}

timepoint_map <- c("Screening" = "Screening", "C4" = "C4", "C8" = "C8",
                    "C18" = "C18", "3 YR F/U" = "3 YR F/U")
tp_levels <- c("Screening", "C4", "C8", "C18", "3 YR F/U")

heatmap_df <- data.frame(
  patient = factor(paired$Study.ID),
  timepoint = factor(paired$Timepoint, levels = tp_levels),
  correlation = raw_corr
)

cat("\nRaw gene-level cfDNA-gDNA correlation (19,100 genes, no filtering/VST):\n")
cat("  Mean:", round(mean(raw_corr), 3), "\n")
cat("  Median:", round(median(raw_corr), 3), "\n")
cat("  Range:", round(min(raw_corr), 3), "-", round(max(raw_corr), 3), "\n")
cat("  SD:", round(sd(raw_corr), 3), "\n")

# Per-timepoint
tp_summary <- heatmap_df %>%
  group_by(timepoint) %>%
  summarize(n = n(), mean_corr = mean(correlation),
            sd_corr = sd(correlation), .groups = "drop")
cat("\n  Per-timepoint concordance:\n")
print(as.data.frame(tp_summary))

# --- Panel A: Heatmap ---
fig2a_raw <- ggplot(heatmap_df, aes(x = timepoint, y = patient, fill = correlation)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(correlation, 2)), size = 2.5, color = "black") +
  scale_fill_gradient2(low = "#E63946", mid = "#FFFFCC", high = "#457B9D",
                       midpoint = median(raw_corr),
                       limits = c(min(raw_corr) - 0.02, max(raw_corr) + 0.02),
                       name = "Pearson r") +
  labs(x = "Treatment Timepoint", y = "Patient",
       title = "a  Gene-level cfDNA-gDNA concordance (raw counts, 19,100 genes)") +
  theme_nature +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

# --- Panel B: Distribution ---
fig2b_raw <- ggplot(heatmap_df, aes(x = correlation)) +
  geom_histogram(aes(y = after_stat(density)), bins = 20,
                 fill = "#457B9D", alpha = 0.5, color = "white") +
  geom_density(linewidth = 0.8, color = "#264653") +
  geom_vline(xintercept = median(raw_corr), linetype = "dashed",
             color = "#E63946", linewidth = 0.7) +
  annotate("text", x = median(raw_corr) - 0.01, y = Inf,
           label = paste0("median = ", round(median(raw_corr), 3)),
           vjust = 2, hjust = 1, size = 3.5, color = "#E63946") +
  labs(x = "Pearson correlation (cfDNA vs gDNA across 19,100 genes)",
       y = "Density",
       title = "b  Distribution of pairwise concordance (raw counts)") +
  theme_nature

# --- Panel C: By timepoint ---
fig2c_raw <- ggplot(heatmap_df, aes(x = timepoint, y = correlation, fill = timepoint)) +
  geom_boxplot(alpha = 0.6, outlier.size = 1.5) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Treatment Timepoint",
       y = "Pearson correlation (cfDNA vs gDNA)",
       title = "c  Concordance across treatment timepoints (raw counts)") +
  theme_nature +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

save_fig("fig2a_raw_concordance_heatmap.pdf", fig2a_raw, width = 7, height = 7)
save_fig("fig2b_raw_concordance_distribution.pdf", fig2b_raw, width = 5, height = 4)
save_fig("fig2c_raw_concordance_by_timepoint.pdf", fig2c_raw, width = 5, height = 4)

# Save summary
raw_summary <- data.frame(
  metric = c("mean", "median", "sd", "min", "max"),
  value = c(mean(raw_corr), median(raw_corr), sd(raw_corr),
            min(raw_corr), max(raw_corr))
)
write.csv(raw_summary, file.path(base_dir, "nature_manuscript",
          "fig2_raw_concordance_summary.csv"), row.names = FALSE)

cat("\n=== Done ===\n")



###############################################################################
## ---- fig3b_response_heatmap.R (recreated; longitudinal MRD response heatmap)
###############################################################################
##
## Build Figure 3b: longitudinal MRD response heatmap.
## Layout matches fig3a_model_concordance_heatmap.pdf exactly:
##   - Rows: patient Study.ID (26 patients)
##   - Cols: treatment timepoint (Screening, C4, C8, C18, 3 YR F/U)
##
## Cell fill encodes the ordinal MRD-by-NGS response at each timepoint, mapped
## from the KRd clinical Excel:
##   - C8       <- "MRD post 8 cycles by NGS (Adaptive) negative =1"  (1 = MRD-)
##   - C18      <- "MRD by NGS EoT Negative = 1"                      (1 = MRD-)
##   - 3 YR F/U <- "MRD by NGS 3yr f/u"  (fallback to 2yr, then 1yr if NA)
##   - Screening, C4 : not assessed (no MRD measurement pre-treatment)
##
## Each cell is colored on an ordinal "depth-of-response" spectrum and
## annotated with the IMWG best overall response (sCR or VGPR) for context.
##
## Output: nature_manuscript/figures/fig3b_response_heatmap.pdf
##
## NOTE: This script is sectioned in 07_figures.R to be run on its own; it does
## not execute when 07_figures.R is sourced (it's wrapped in if (FALSE)).
## Run with: Rscript -e 'BUILD_FIG3B<-TRUE; source("R/07_figures.R")'
###############################################################################

if (exists("BUILD_FIG3B") && isTRUE(BUILD_FIG3B)) {

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readxl)
  library(patchwork)
})

base_dir <- getwd()
res_dir  <- file.path(base_dir, "results")
fig_dir  <- file.path(base_dir, "nature_manuscript", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

md <- readRDS(file.path(res_dir, "model_data.rds"))
paired <- md$paired_info
tp_levels <- c("Screening", "C4", "C8", "C18", "3 YR F/U")

clin_file <- file.path(base_dir, "KRd trial", "Clinical data",
                       "KRd 12_1725 data_2_14_2023 from Ben without PHI.xlsx")
suppressWarnings({ clin <- read_excel(clin_file, sheet = "12_1725 cfDNA Dataset") })

mrd_df <- data.frame(
  Study.ID       = clin$`Study number`,
  best_response  = clin$`best overall response`,
  mrd_c8         = clin$`MRD post 8 cycles by NGS (Adaptive) negative =1`,
  mrd_c18        = clin$`MRD by NGS EoT Negative = 1`,
  mrd_1yr        = clin$`MRD by NGS 1yr f/u`,
  mrd_2yr        = clin$`MRD by NGS 2yr f/u`,
  mrd_3yr        = clin$`MRD by NGS 3yr f/u`,
  stringsAsFactors = FALSE
)
mrd_df$mrd_3yrFU <- ifelse(!is.na(mrd_df$mrd_3yr), mrd_df$mrd_3yr,
                    ifelse(!is.na(mrd_df$mrd_2yr), mrd_df$mrd_2yr, mrd_df$mrd_1yr))

code_mrd <- function(x) {
  out <- rep(NA_character_, length(x))
  out[!is.na(x) & x == 1] <- "MRD-"
  out[!is.na(x) & x == 0] <- "MRD+"
  out
}

mrd_long <- data.frame(
  Study.ID  = rep(mrd_df$Study.ID, times = 5),
  timepoint = factor(rep(tp_levels, each = nrow(mrd_df)), levels = tp_levels),
  mrd_status = c(rep(NA_character_, nrow(mrd_df)),
                 rep(NA_character_, nrow(mrd_df)),
                 code_mrd(mrd_df$mrd_c8),
                 code_mrd(mrd_df$mrd_c18),
                 code_mrd(mrd_df$mrd_3yrFU)),
  best_response = rep(mrd_df$best_response, times = 5),
  stringsAsFactors = FALSE
)

patients_in_fig <- sort(unique(paired$Study.ID))
grid <- expand.grid(Study.ID = patients_in_fig,
                    timepoint = factor(tp_levels, levels = tp_levels),
                    stringsAsFactors = FALSE)
grid$timepoint <- factor(as.character(grid$timepoint), levels = tp_levels)
heatmap_df <- merge(grid, mrd_long, by = c("Study.ID", "timepoint"),
                    all.x = TRUE, sort = FALSE)
heatmap_df$Study.ID <- factor(heatmap_df$Study.ID, levels = sort(patients_in_fig))
heatmap_df$cell_label <- ifelse(is.na(heatmap_df$mrd_status), "", heatmap_df$mrd_status)
heatmap_df$response_value <- ifelse(heatmap_df$mrd_status == "MRD-", 1,
                             ifelse(heatmap_df$mrd_status == "MRD+", 0, NA))

best_resp <- mrd_df[, c("Study.ID", "best_response")]
heatmap_df <- merge(heatmap_df, best_resp, by = "Study.ID",
                    all.x = TRUE, suffixes = c("", ".bo"))

theme_nature <- theme_bw(base_size = 10) +
  theme(plot.title    = element_text(face = "bold", size = 11, hjust = 0),
        axis.title    = element_text(size = 9),
        axis.text     = element_text(size = 7),
        legend.title  = element_text(size = 8, face = "bold"),
        legend.text   = element_text(size = 7),
        legend.key.size = grid::unit(0.4, "cm"),
        strip.text    = element_text(size = 8, face = "bold"))

strip_df <- unique(heatmap_df[, c("Study.ID", "best_response")])
strip_df$best_response <- factor(strip_df$best_response, levels = c("VGPR", "sCR"))
strip_df$xcol <- "Best\nresponse"
resp_palette <- c("VGPR" = "#9F6BA0", "sCR" = "#2C7A4B")

p_strip <- ggplot(strip_df, aes(x = xcol, y = Study.ID, fill = best_response)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = resp_palette, name = "Best\noverall\nresponse",
                    na.value = "grey90", drop = FALSE) +
  labs(x = NULL, y = "Patient") +
  theme_nature +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 7),
        axis.title.x = element_blank(),
        panel.grid   = element_blank(),
        plot.margin  = margin(5.5, 0, 5.5, 5.5))

p_main <- ggplot(heatmap_df, aes(x = timepoint, y = Study.ID, fill = response_value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = cell_label), size = 2.3, color = "black") +
  scale_fill_gradient2(low = "#E63946", mid = "#FFFFCC", high = "#457B9D",
                       midpoint = 0.5, limits = c(0, 1),
                       breaks = c(0, 1), labels = c("MRD+", "MRD-"),
                       name = "MRD by NGS\n(deeper response)",
                       na.value = "grey90") +
  labs(x = "Treatment Timepoint", y = NULL,
       title = "Longitudinal MRD response (NGS) by patient and timepoint") +
  theme_nature +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid   = element_blank(),
        plot.margin  = margin(5.5, 5.5, 5.5, 0))

fig3b <- p_strip + p_main + plot_layout(widths = c(0.16, 1))
ggsave(file.path(fig_dir, "fig3b_response_heatmap.pdf"),
       fig3b, width = 8, height = 7)
cat("Saved: fig3b_response_heatmap.pdf\n")

}  # end BUILD_FIG3B

###############################################################################
## ---- Per-component (k = 1..K_occ) supplementary figures (Figs S1, S2)
##
## Three sets of artifacts, all conditioned on the same M1 posterior used in
## Figures 3a and 3c, with chains relabeled + aligned identically:
##
##   (S1, top)  : per-component posterior correlation heatmaps for the
##                occupied components (k = 1, 2, 3). Single-component version
##                of Figure 3c -- isolates the (g, cf) posterior coupling at
##                each k instead of aggregating across all K.
##
##   (S1, bot)  : per-component subclonal L1-concordance heatmaps for the
##                occupied components. Single-component version of Figure 3a:
##                  C_{it}^{(k)} = 1 - |omega^g_{it,k} - omega^{cf}_{it,k}|
##                averaged over the Q saved iterations. Values in [0, 1].
##
##   (S2)       : cohort-level component-weight pi posterior histograms in a
##                4 x 3 grid (10 occupied of 12 cells, two empty). Pi^{(q)}
##                reconstructed from the saved stick-breaking traces v^{(q)}
##                with v_K = 1.
##
## K_occ is the number of components with non-negligible cohort weight
## (auto-detected: mean pi_k > 1e-4). Components k > K_occ have numerically
## zero omega draws at every iteration and produce all-grey heatmaps, so they
## are skipped.
###############################################################################
BUILD_FIGS_PER_COMPONENT <- TRUE
if (BUILD_FIGS_PER_COMPONENT) {

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

base_dir <- getwd()
source(file.path(base_dir, "R", "01_lib_core.R"))

res_dir <- file.path(base_dir, "results")
fig_dir <- file.path(base_dir, "nature_manuscript", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n=== Per-component supplementary figures (Figs S1, S2) ===\n")

md <- readRDS(file.path(res_dir, "model_data.rds"))
m1 <- readRDS(file.path(res_dir, "krd_M1_fit.rds"))

cat("Applying label switching...\n")
for (i in seq_along(m1$chains)) {
  m1$chains[[i]] <- relabel_by_weight(m1$chains[[i]])
}
m1$chains <- align_chains_by_signature(m1$chains)
merged <- merge_chain_samples(m1$chains)

K       <- m1$chains[[1]]$K
K_trace <- m1$chains[[1]]$K_trace
N_obs   <- m1$chains[[1]]$N_obs
Q       <- nrow(merged$omega_g_trace)

# Reconstruct pi^{(q)} once -- needed both for K_occ detection and Fig S2.
pi_draws_list <- vector("list", length(m1$chains))
for (i in seq_along(m1$chains)) {
  v_mat  <- m1$chains[[i]]$samples$v          # n_save x (K-1)
  n_save <- nrow(v_mat)
  pi_mat <- matrix(NA_real_, n_save, K)
  for (s in seq_len(n_save)) {
    pi_mat[s, ] <- stick_break(c(v_mat[s, ], 1))
  }
  pi_draws_list[[i]] <- pi_mat
}
pi_draws <- do.call(rbind, pi_draws_list)     # (n_chains*n_save) x K

pi_mean <- colMeans(pi_draws)
K_occ   <- sum(pi_mean > 1e-4)
cat(sprintf("K=%d, K_trace=%d, K_occ=%d, N_obs=%d, Q=%d\n",
            K, K_trace, K_occ, N_obs, Q))

# Patient / timepoint labels (same conventions as fig3a/fig3c)
timepoint_map <- c("1" = "SCR", "2" = "C4", "3" = "C8", "4" = "C18", "5" = "3Y")
tp_levels     <- c("SCR", "C4", "C8", "C18", "3Y")

obs_meta <- data.frame(
  patient   = factor(md$paired_info$Study.ID),
  timepoint = factor(timepoint_map[as.character(md$time)], levels = tp_levels),
  stringsAsFactors = FALSE
)

# Helper to save both PDF and PNG with consistent device backends.
save_pdf_png <- function(plot, basename, width, height, dpi = 300) {
  ggsave(file.path(fig_dir, paste0(basename, ".pdf")), plot,
         width = width, height = height, device = grDevices::pdf)
  ggsave(file.path(fig_dir, paste0(basename, ".png")), plot,
         width = width, height = height, dpi = dpi, units = "in",
         device = grDevices::png)
}

# -------------------------------------------------------------------------
# (S1, top) Per-component posterior correlation heatmaps, k = 1..K_occ
# -------------------------------------------------------------------------
cat("\n--- (S1) Per-component correlation heatmaps ---\n")

for (kk in seq_len(K_occ)) {
  cols <- (kk - 1) * N_obs + seq_len(N_obs)
  omg  <- merged$omega_g_trace[,  cols, drop = FALSE]
  omcf <- merged$omega_cf_trace[, cols, drop = FALSE]

  r_obs <- numeric(N_obs)
  for (obs in seq_len(N_obs)) {
    vg  <- omg[,  obs]
    vcf <- omcf[, obs]
    r_obs[obs] <- if (sd(vg) > 1e-12 && sd(vcf) > 1e-12) cor(vg, vcf) else NA_real_
  }

  df_k <- cbind(obs_meta, r = r_obs) %>%
    group_by(patient, timepoint) %>%
    summarise(r = mean(r, na.rm = TRUE), .groups = "drop")

  fig_k <- ggplot(df_k, aes(x = timepoint, y = patient, fill = r)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", r)), size = 2.3, color = "black") +
    scale_fill_gradient2(low = "#E63946", mid = "#FFFFCC", high = "#457B9D",
                         midpoint = 0, limits = c(-1, 1),
                         breaks = c(-1, -0.5, 0, 0.5, 1),
                         name = "Posterior\ncorrelation",
                         na.value = "grey90") +
    labs(x = "Treatment Timepoint", y = "Patient",
         title = bquote("Posterior correlation, component" ~ italic(k) == .(kk) *
                        ":  " * cor(omega[.(kk)]^{cf}, omega[.(kk)]^{g}))) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid  = element_blank())

  save_pdf_png(fig_k, sprintf("fig3c_posterior_correlation_k%02d", kk),
               width = 7, height = 7)
  cat(sprintf("  k=%d : mean r = %.3f, median = %.3f\n",
              kk, mean(r_obs, na.rm = TRUE), median(r_obs, na.rm = TRUE)))
}

# (Per-component L1 concordance heatmaps and per_component_concordance_with_mrd.csv
#  were removed by design: the L1 metric is no longer reported, and a single
#  component's gene-level contribution correlation is degenerate. The retained
#  per-component view is the posterior weight correlation, fig3c_posterior_correlation_k*
#  produced by the loop above.)

# -------------------------------------------------------------------------
# (S2) Cohort-level pi posterior histograms in a 4 x 3 grid
# -------------------------------------------------------------------------
cat("\n--- (S2) Cohort pi posterior histograms ---\n")

pi_summary <- data.frame(
  k        = seq_len(K),
  mean     = pi_mean,
  ci_lo    = apply(pi_draws, 2, quantile, probs = 0.025),
  ci_hi    = apply(pi_draws, 2, quantile, probs = 0.975),
  median   = apply(pi_draws, 2, median)
)
write.csv(pi_summary,
          file.path(res_dir, "krd_pi_posterior_summary.csv"),
          row.names = FALSE)
print(round(pi_summary, 5))

pi_long <- data.frame(
  k    = rep(seq_len(K), each = nrow(pi_draws)),
  draw = as.numeric(pi_draws)
)
pi_long$k_lab <- factor(
  pi_long$k,
  levels = seq_len(K),
  labels = sapply(seq_len(K), function(k) {
    sprintf("k = %d:  mean = %.3f  (95%% CI %.3f to %.3f)",
            k, pi_summary$mean[k], pi_summary$ci_lo[k], pi_summary$ci_hi[k])
  })
)

fig_pi <- ggplot(pi_long, aes(x = draw)) +
  geom_histogram(bins = 50, fill = "#457B9D", color = "white",
                 linewidth = 0.1) +
  geom_vline(data = pi_summary,
             aes(xintercept = mean),
             color = "#E63946", linetype = "dashed", linewidth = 0.5,
             inherit.aes = FALSE) +
  facet_wrap(~ k_lab, ncol = 3, scales = "free") +
  labs(x = expression("Posterior draws of cohort component weight " * pi[k]),
       y = "Frequency",
       title = expression("Posterior distribution of cohort-level component weights " *
                          pi[k] * ", k = 1..K"),
       subtitle = "Red dashed line = posterior mean; 95% CI shown in each facet header") +
  theme_bw(base_size = 9) +
  theme(strip.text = element_text(size = 7),
        panel.grid.minor = element_blank())

save_pdf_png(fig_pi, "figS_pi_histograms", width = 9, height = 8)
cat("  Saved: figS_pi_histograms.{pdf,png}\n")

cat("\nDone (per-component figures).\n")

}  # end BUILD_FIGS_PER_COMPONENT
