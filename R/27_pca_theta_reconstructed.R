###############################################################################
## R/27_pca_theta_reconstructed.R
##
## Model-based counterpart of the observed per-cycle paired PCA (R/19). The
## latent signatures theta are SHARED across sources by construction (a single
## theta_{ijk} serves both, writeup Sec 4.1); the only source-specific model
## quantity is the reconstructed profile
##       Yhat^s_it = sum_k omega^s_it,k * theta_{i,.,k} = theta(.,.,i) %*% omega^s_it,
## where the third index of theta is the SUBJECT i (signatures are
## time-invariant but subject-specific). This uses the posterior-mean
## signatures (theta_postmean) and posterior-mean source-specific weights
## under the tracking model M1. We run a SINGLE PCA on
## the combined reconstructed matrix [Yhat_cf | Yhat_g] and facet the paired
## points/segments by treatment timepoint, exactly mirroring R/19 so the two
## PCAs (observed = fig5c_pca_paired_bycycle.pdf; model = this script) form the
## two panels of the lead figure requested by the PI.
##
## Output: nature_manuscript/figures/fig_pca_model_bycycle.pdf
###############################################################################

suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

base_dir <- getwd()  # run from project root
setwd(base_dir)
res_dir <- file.path(base_dir, "results")
fig_dir <- file.path(base_dir, "nature_manuscript", "figures")
source(file.path(base_dir, "R", "01_lib_core.R"))

## --- Load data + fit, canonicalize (same pipeline as R/04) --------------------
md <- readRDS(file.path(res_dir, "model_data.rds"))
m1 <- readRDS(file.path(res_dir, "krd_M1_fit.rds"))
for (i in seq_along(m1$chains)) m1$chains[[i]] <- relabel_by_weight(m1$chains[[i]])
m1$chains <- align_chains_by_signature(m1$chains)
merged <- merge_chain_samples(m1$chains)

K       <- m1$chains[[1]]$K
K_trace <- m1$chains[[1]]$K_trace
N_obs   <- m1$chains[[1]]$N_obs
n_ch    <- length(m1$chains)
p_genes <- nrow(md$Y_g)
n_subj  <- md$n          # third dim of theta_postmean is the SUBJECT
T_max   <- max(md$time)

## --- Posterior-mean source-specific weights [K_trace x N_obs] -----------------
omega_g_post  <- matrix(NA_real_, K_trace, N_obs)
omega_cf_post <- matrix(NA_real_, K_trace, N_obs)
for (k in 1:K_trace) {
  for (obs in 1:N_obs) {
    col_idx <- (k - 1) * N_obs + obs
    omega_g_post[k, obs]  <- mean(merged$omega_g_trace[, col_idx])
    omega_cf_post[k, obs] <- mean(merged$omega_cf_trace[, col_idx])
  }
}

## --- Posterior-mean signatures theta_avg [p x K x n] (canonical, NOT final) ---
## The third dimension indexes the SUBJECT: signatures are time-invariant but
## subject-specific (writeup Sec 4.1). An earlier version allocated
## [p x K x T_max] and sliced by timepoint, which errors whenever n != T_max.
theta_avg <- array(0, dim = c(p_genes, K, n_subj))
for (ch in seq_len(n_ch)) theta_avg <- theta_avg + m1$chains[[ch]]$theta_postmean
theta_avg <- theta_avg / n_ch

## --- Reconstruct per-observation profiles Yhat^s [p x N_obs] ------------------
Yhat_cf <- matrix(NA_real_, p_genes, N_obs)
Yhat_g  <- matrix(NA_real_, p_genes, N_obs)
for (obs in 1:N_obs) {
  i <- md$patient[obs]
  Yhat_cf[, obs] <- as.vector(theta_avg[, , i] %*% omega_cf_post[, obs])
  Yhat_g[,  obs] <- as.vector(theta_avg[, , i] %*% omega_g_post[,  obs])
}

## --- SINGLE PCA on combined reconstructed matrix (mirror R/19) ----------------
tp_levels <- c("Screening", "C4", "C8", "C18", "3 YR F/U")
tp <- factor(tp_levels[md$time], levels = tp_levels)

Y_combined <- cbind(Yhat_cf, Yhat_g)          # cols 1..n = cfDNA, n+1..2n = gDNA
pca <- prcomp(t(Y_combined), center = TRUE, scale. = TRUE)
ve  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

pts <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                  Source = c(rep("cfDNA", N_obs), rep("gDNA", N_obs)),
                  Timepoint = rep(tp, 2))
seg <- data.frame(PC1_cf = pca$x[1:N_obs, 1],              PC2_cf = pca$x[1:N_obs, 2],
                  PC1_g  = pca$x[(N_obs+1):(2*N_obs), 1],  PC2_g  = pca$x[(N_obs+1):(2*N_obs), 2],
                  Timepoint = tp)

p <- ggplot() +
  geom_segment(data = seg, aes(x = PC1_cf, y = PC2_cf, xend = PC1_g, yend = PC2_g),
               color = "grey70", linewidth = 0.3) +
  geom_point(data = pts, aes(PC1, PC2, color = Source), size = 1) +
  facet_grid(. ~ Timepoint) +
  scale_color_manual(values = c(cfDNA = "#E63946", gDNA = "#457B9D")) +
  labs(x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)"),
       title = "Model-reconstructed paired cfDNA-gDNA profiles in PC space, by treatment timepoint") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "fig_pca_model_bycycle.pdf"), p, width = 11, height = 2.9)
cat("Saved: nature_manuscript/figures/fig_pca_model_bycycle.pdf\n")
cat(sprintf("Reconstructed PCA: PC1 %.1f%%, PC2 %.1f%%; per-cycle n: %s\n", ve[1], ve[2],
            paste(sprintf("%s=%d", tp_levels, as.integer(table(tp)[tp_levels])), collapse=", ")))
