###############################################################################
## 22_posterior_concordance.R
## Compute posterior-based concordance metrics from MCMC output
##
## Metrics (the L1 subclonal concordance metric was REMOVED by design; the paper
## reports only the posterior correlation of omega_cf and omega_g):
##   (c) Component-level correlation for occupied components
##   (d) Model-based gene-level prediction correlation (Corr of Yhat_g, Yhat_cf)
##   (e) Posterior iteration-level correlation of omega_g, omega_cf  [PRIMARY]
##
## Requires: results/krd_M1_fit.rds, results/model_data.rds
## Usage: source("R/22_posterior_concordance.R")
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

base_dir <- getwd()
source(file.path(base_dir, "R", "01_lib_core.R"))

res_dir <- file.path(base_dir, "results")
fig_dir <- file.path(base_dir, "nature_manuscript", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Posterior Concordance Metrics ===\n\n")

# Load data
md <- readRDS(file.path(res_dir, "model_data.rds"))
m1 <- readRDS(file.path(res_dir, "krd_M1_fit.rds"))

# Apply label switching
cat("Applying label switching...\n")
for (i in seq_along(m1$chains)) {
  m1$chains[[i]] <- relabel_by_weight(m1$chains[[i]])
}
m1$chains <- align_chains_by_signature(m1$chains)
merged <- merge_chain_samples(m1$chains)

K <- m1$chains[[1]]$K
K_trace <- m1$chains[[1]]$K_trace
N_obs <- m1$chains[[1]]$N_obs
n_ch <- length(m1$chains)
n_total <- nrow(merged$omega_g_trace)  # total saved iterations across chains

cat(sprintf("K=%d, K_trace=%d, N_obs=%d, n_total=%d\n", K, K_trace, N_obs, n_total))

# Helper: extract omega vector for iteration s, observation obs from merged traces
extract_omega <- function(trace_matrix, s, obs, K_trace, N_obs) {
  sapply(1:K_trace, function(k) trace_matrix[s, (k - 1) * N_obs + obs])
}

# ===========================================================================
# Posterior-mean source-specific weights (shared by metrics (d) and (e)).
# NOTE: the L1 subclonal concordance metric was removed by design --- the paper
# reports ONLY the posterior correlation of omega_cf and omega_g (metric (e))
# plus the model-based prediction correlation (metric (d)). We retain the
# posterior-mean omega matrices below because metric (d) needs them.
# ===========================================================================
cat("\n--- Posterior-mean source-specific weights ---\n")

omega_g_post <- matrix(NA, K_trace, N_obs)
omega_cf_post <- matrix(NA, K_trace, N_obs)
for (k in 1:K_trace) {
  for (obs in 1:N_obs) {
    col_idx <- (k - 1) * N_obs + obs
    omega_g_post[k, obs] <- mean(merged$omega_g_trace[, col_idx])
    omega_cf_post[k, obs] <- mean(merged$omega_cf_trace[, col_idx])
  }
}

# Timepoint label maps (relocated from the removed L1 block; used by metric (e)).
# Adapt to the analysis coding: 5-timepoint clinical schedule or 2-phase collapse.
if (max(md$time, na.rm = TRUE) <= 2) {
  timepoint_map <- c("1" = "early", "2" = "late")
  tp_levels     <- c("early", "late")
} else {
  timepoint_map <- c("1" = "SCR", "2" = "C4", "3" = "C8", "4" = "C18", "5" = "3Y")
  tp_levels     <- c("SCR", "C4", "C8", "C18", "3Y")
}

# ===========================================================================
# METRIC (c): Component-level correlation (occupied components only)
# ===========================================================================
cat("\n--- Metric (c): Component-level correlation ---\n")

comp_corr_matrix <- matrix(NA, n_total, N_obs)
for (s in 1:n_total) {
  # Reconstruct pi from v at this iteration
  v_s <- c(merged$v[s, ], 1)
  pi_s <- stick_break(v_s)
  occupied <- which(pi_s > 0.01)

  for (obs in 1:N_obs) {
    og <- extract_omega(merged$omega_g_trace, s, obs, K_trace, N_obs)
    oc <- extract_omega(merged$omega_cf_trace, s, obs, K_trace, N_obs)
    if (length(occupied) >= 2) {
      og_occ <- og[occupied]
      oc_occ <- oc[occupied]
      sd_g <- sd(og_occ, na.rm = TRUE)
      sd_c <- sd(oc_occ, na.rm = TRUE)
      if (!is.na(sd_g) && !is.na(sd_c) && sd_g > 1e-10 && sd_c > 1e-10) {
        comp_corr_matrix[s, obs] <- cor(og_occ, oc_occ, use = "complete.obs")
      }
    }
  }
  if (s %% 200 == 0) cat("    iteration", s, "/", n_total, "\n")
}

comp_corr_mean <- colMeans(comp_corr_matrix, na.rm = TRUE)
comp_corr_lo <- apply(comp_corr_matrix, 2, quantile, 0.025, na.rm = TRUE)
comp_corr_hi <- apply(comp_corr_matrix, 2, quantile, 0.975, na.rm = TRUE)

cat("  Component-level correlation (occupied only):\n")
cat("    overall mean =", round(mean(comp_corr_mean, na.rm = TRUE), 4), "\n")
cat("    overall median =", round(median(comp_corr_mean, na.rm = TRUE), 4), "\n")

# ===========================================================================
# METRIC (d): Model-based gene-level prediction correlation
# ===========================================================================
cat("\n--- Metric (d): Model-based gene-level prediction correlation ---\n")

# Average the posterior-mean signatures across chains. Use theta_postmean, NOT
# final_theta: final_theta is a single last iteration and is invalid under the
# posterior multimodality of this sampler (see R/01_lib_core.R). theta_postmean
# is the online posterior mean in canonical component order.
#
# NOTE: theta_postmean is [p x K x n] -- the third dimension indexes the
# SUBJECT, not the timepoint. (Signatures are time-invariant but
# subject-specific; see writeup Sec 4.1.) An earlier version of this script
# allocated [p x K x T_max] and sliced by md$time[obs], which is wrong on both
# counts and errors outright whenever n != T_max.
p_genes <- m1$chains[[1]]$p
n_subj  <- m1$chains[[1]]$n
theta_avg <- array(0, dim = c(p_genes, K, n_subj))
for (ch in 1:n_ch) {
  theta_avg <- theta_avg + m1$chains[[ch]]$theta_postmean
}
theta_avg <- theta_avg / n_ch

# Compute model-based prediction correlation per observation
pred_corr <- numeric(N_obs)
for (obs in 1:N_obs) {
  i <- md$patient[obs]
  # yhat_s[j] = sum_k omega_s[k,obs] * theta[j,k,i]  -- subject i's signatures
  yhat_g  <- as.vector(theta_avg[, , i] %*% omega_g_post[, obs])
  yhat_cf <- as.vector(theta_avg[, , i] %*% omega_cf_post[, obs])
  pred_corr[obs] <- cor(yhat_g, yhat_cf)
}

cat("  Model-based gene-level prediction correlation:\n")
cat("    mean =", round(mean(pred_corr), 4),
    ", median =", round(median(pred_corr), 4),
    ", range = (", round(min(pred_corr), 4), ",",
    round(max(pred_corr), 4), ")\n")

# Save prediction correlation
pred_corr_df <- data.frame(
  obs = 1:N_obs, patient = md$patient, time = md$time,
  pred_correlation = round(pred_corr, 4)
)
write.csv(pred_corr_df, file.path(res_dir, "krd_prediction_correlation.csv"),
          row.names = FALSE)
cat("  Saved: krd_prediction_correlation.csv\n")

# Prediction correlation heatmap (same format as naive Pearson r)
pred_df <- data.frame(
  patient = factor(md$paired_info$Study.ID),
  timepoint = factor(md$paired_info$Timepoint,
                     levels = c("Screening", "C4", "C8", "C18", "3 YR F/U")),
  correlation = pred_corr
)

# Aggregate duplicate patient-timepoints (e.g., 101-17 and 101-73 have 2 entries at 3YR F/U)
pred_df <- pred_df %>%
  group_by(patient, timepoint) %>%
  summarise(correlation = mean(correlation), .groups = "drop")

fig_pred_heatmap <- ggplot(pred_df, aes(x = timepoint, y = patient, fill = correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(correlation, 2)), size = 2.5) +
  scale_fill_gradient2(low = "#FEE090", mid = "#4575B4", high = "#313695",
                       midpoint = median(pred_corr), limits = c(min(pred_corr), 1),
                       name = "Model r") +
  labs(x = "Treatment Timepoint", y = "Patient",
       title = "Model-based gene-level prediction correlation") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file.path(fig_dir, "fig_prediction_correlation_heatmap.pdf"), width = 7, height = 8)
print(fig_pred_heatmap)
dev.off()
cat("  Saved: fig_prediction_correlation_heatmap.pdf\n")

# ===========================================================================
# SUMMARY
# ===========================================================================
cat("\n=== SUMMARY OF POSTERIOR CONCORDANCE METRICS ===\n")
cat("(c) Component-level correlation (occupied): mean =", round(mean(comp_corr_mean, na.rm=TRUE), 4), "\n")
cat("(d) Model-based gene-level prediction corr: mean =", round(mean(pred_corr), 4), "\n")
cat("(e) Posterior correlation of omega_g, omega_cf: reported below.\n")

# ===========================================================================
# METRIC (e): Posterior iteration-level correlation between source-specific
#             weight draws (Fig 3c).
#
# For each paired observation (i, t), compute the Pearson correlation across
# posterior iterations and components q,k between the source-specific weight
# vectors:
#   r_{it} = cor( vec(omega^{g,q}_{it,k}) , vec(omega^{cf,q}_{it,k}) )
# where the vectors have length Q*K (Q = saved iterations, K = K_trace).
# Captures whether gDNA and cfDNA mixture weight posteriors move together.
# ===========================================================================
cat("\n--- Metric (e): Posterior iteration-level correlation (Fig 3c) ---\n")

post_corr_obs <- numeric(N_obs)
for (obs in 1:N_obs) {
  cols <- ((seq_len(K_trace) - 1) * N_obs) + obs
  omega_g_qk  <- merged$omega_g_trace[,  cols, drop = FALSE]   # Q x K
  omega_cf_qk <- merged$omega_cf_trace[, cols, drop = FALSE]   # Q x K
  vg  <- as.numeric(omega_g_qk)
  vcf <- as.numeric(omega_cf_qk)
  if (sd(vg) > 1e-12 && sd(vcf) > 1e-12) {
    post_corr_obs[obs] <- cor(vg, vcf)
  } else {
    post_corr_obs[obs] <- NA_real_
  }
}

cat(sprintf("  Posterior correlation across iterations x components:\n"))
cat(sprintf("    mean   = %.4f\n",   mean(post_corr_obs,  na.rm = TRUE)))
cat(sprintf("    median = %.4f\n", median(post_corr_obs,  na.rm = TRUE)))
cat(sprintf("    range  = (%.4f, %.4f)\n",
            min(post_corr_obs, na.rm = TRUE), max(post_corr_obs, na.rm = TRUE)))

# Heatmap data frame on the same patient x timepoint grid.
# timepoint_map / tp_levels are defined once near the top (adaptive to 5-cycle vs collapse).
post_corr_df <- data.frame(
  patient    = factor(md$paired_info$Study.ID),
  timepoint  = factor(timepoint_map[as.character(md$time)], levels = tp_levels),
  post_corr  = post_corr_obs,
  stringsAsFactors = FALSE
) %>%
  group_by(patient, timepoint) %>%
  summarise(post_corr = mean(post_corr, na.rm = TRUE), .groups = "drop")

write.csv(post_corr_df, file.path(res_dir, "krd_posterior_correlation.csv"), row.names = FALSE)
cat("  Saved: krd_posterior_correlation.csv\n")

fig_post_corr <- ggplot(post_corr_df, aes(x = timepoint, y = patient, fill = post_corr)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", post_corr)), size = 2.3, color = "black") +
  scale_fill_gradient2(low = "#E63946", mid = "#FFFFCC", high = "#457B9D",
                       midpoint = 0, limits = c(-1, 1),
                       breaks = c(-1, -0.5, 0, 0.5, 1),
                       name = "Posterior\ncorrelation",
                       na.value = "grey90") +
  labs(x = "Treatment Timepoint", y = "Patient",
       title = expression("Posterior correlation across iterations: " *
                          cor(omega^{cf}, omega^{g}))) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

ggsave(file.path(fig_dir, "fig3c_posterior_correlation_heatmap.pdf"),
       fig_post_corr, width = 7, height = 7)
cat("  Saved: fig3c_posterior_correlation_heatmap.pdf\n")
