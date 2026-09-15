###############################################################################
## R/14_permutation_negative_control.R
##
## GOAL: prove the observed cfDNA-gDNA concordance is REAL (pairing-specific),
## not an artifact of 5hmC profiles being broadly similar across all samples.
## Strategy: break the true patient-timepoint pairing by permuting which cfDNA
## sample is matched to which gDNA sample, recompute every concordance metric,
## and show the observed values exceed the permutation null.
##
## Two levels, both from EXISTING fits (no re-fit needed):
##   (1) Model-free gene-level Pearson r  (the r = 0.956 headline)
##   (2) Model-based subclonal metric: posterior correlation of the source-specific
##       omega weights (Fig 3c), from canonical (relabeled/aligned) omega weights.
##       (The L1 concordance metric was removed by design.)
##
## Null = random pairing within the SAME timepoint stratum (so the null is not
## trivially inflated by timepoint structure) AND, as a stricter alternative,
## fully random pairing across all observations. Derangements only (no self-pair).
##
## Outputs (results/): permutation_negative_control.csv
## Outputs (figures/): fig_permutation_null.pdf
###############################################################################

set.seed(20260718)
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

base_dir <- getwd(); res_dir <- file.path(base_dir, "results")
fig_dir  <- file.path(base_dir, "nature_manuscript", "figures")
source(file.path(base_dir, "R", "01_lib_core.R"))

md  <- readRDS(file.path(res_dir, "model_data.rds"))
Yg  <- md$Y_g; Ycf <- md$Y_cf                 # [p x N_obs] VST signals
N   <- md$N_obs; tvec <- md$time

## canonical posterior-mean weights per obs (relabel/align/merge)
krd <- readRDS(file.path(res_dir, "krd_M1_fit.rds")); m1 <- krd
for (i in seq_along(m1$chains)) m1$chains[[i]] <- relabel_by_weight(m1$chains[[i]])
m1$chains <- align_chains_by_signature(m1$chains)
merged <- merge_chain_samples(m1$chains)
K <- ncol(merged$omega_g_trace) / N
Wg  <- sapply(1:K, function(k) sapply(1:N, function(o) mean(merged$omega_g_trace[, (k-1)*N+o])))
Wcf <- sapply(1:K, function(k) sapply(1:N, function(o) mean(merged$omega_cf_trace[, (k-1)*N+o])))
# rows = obs, cols = component

## ---- metric functions (per matched index vector `match_cf`) --------------
gene_r <- function(match_cf)
  mean(sapply(1:N, function(i) cor(Yg[, i], Ycf[, match_cf[i]])))
post_r <- function(match_cf)
  mean(sapply(1:N, function(i) {
    a <- Wg[i, ]; b <- Wcf[match_cf[i], ]
    if (sd(a) < 1e-9 || sd(b) < 1e-9) NA else cor(a, b) }), na.rm = TRUE)

identity_match <- 1:N
obs <- c(gene_r = gene_r(identity_match),
         post_r  = post_r(identity_match))

## ---- permutation nulls ---------------------------------------------------
n_perm <- 2000
## stratified derangement: permute cf indices WITHIN timepoint stratum
strat_derange <- function() {
  m <- 1:N
  for (tp in unique(tvec)) {
    idx <- which(tvec == tp)
    if (length(idx) > 1) {
      repeat { p <- sample(idx); if (all(p != idx)) break }
      m[idx] <- p
    }
  }
  m
}
## global derangement across all observations
glob_derange <- function() { repeat { p <- sample(1:N); if (all(p != 1:N)) break }; p }

run_null <- function(derange_fn) {
  out <- matrix(NA, n_perm, 2, dimnames = list(NULL, c("gene_r","post_r")))
  for (b in 1:n_perm) {
    m <- derange_fn()
    out[b, ] <- c(gene_r(m), post_r(m))
  }
  out
}
cat("Running stratified null (", n_perm, "perms)...\n"); null_s <- run_null(strat_derange)
cat("Running global null...\n");                          null_g <- run_null(glob_derange)

## ---- summarize -----------------------------------------------------------
summ <- function(null_mat, label) {
  do.call(rbind, lapply(colnames(null_mat), function(mn) {
    nv <- null_mat[, mn]; ov <- obs[[mn]]
    data.frame(null = label, metric = mn, observed = ov,
               null_mean = mean(nv, na.rm = TRUE), null_sd = sd(nv, na.rm = TRUE),
               null_q975 = quantile(nv, 0.975, na.rm = TRUE),
               p_perm = (1 + sum(nv >= ov, na.rm = TRUE)) / (1 + sum(!is.na(nv))),
               z = (ov - mean(nv, na.rm = TRUE)) / sd(nv, na.rm = TRUE))
  }))
}
res <- rbind(summ(null_s, "stratified_within_timepoint"), summ(null_g, "global"))
res$metric <- recode(res$metric, gene_r = "Gene-level Pearson r",
                     post_r = "Posterior correlation")
cat("\n==== Permutation negative control ====\n")
print(res, row.names = FALSE, digits = 4)
write.csv(res, file.path(res_dir, "permutation_negative_control.csv"), row.names = FALSE)
cat("\nSaved: results/permutation_negative_control.csv\n")

## ---- figure: pairing-specific gain on a COMMON (null-centred) x-axis ------
## The earlier version used per-panel free x-axes (scales="free"), each zoomed
## to its own null width; this made the tiny gene-level gain (~0.009) look as
## separated from its null as the genuine subclonal gains. Here we centre every
## metric's null at its own null mean and plot the deviation from that null
## (the "pairing-specific gain" = value - null mean) on a SHARED x-axis. The
## dashed line marks the null (zero gain) and the red line the observed gain, so
## the size of the gain is directly comparable: the gene-level red sits almost
## on top of its null, whereas the subclonal red lines are clearly to the right.
nmean <- c(gene_r  = mean(null_s[, "gene_r"]),
           post_r  = mean(null_s[, "post_r"]))
mlab  <- c(gene_r  = "Gene-level Pearson r",
           post_r  = "Posterior correlation")
plot_df <- bind_rows(lapply(names(nmean), function(m)
  data.frame(metric = factor(mlab[m], levels = mlab), value = null_s[, m] - nmean[m])))
gains  <- c(obs[["gene_r"]], obs[["post_r"]]) - nmean
obs_df <- data.frame(metric = factor(mlab, levels = mlab), gain = as.numeric(gains),
                     lab = sprintf("observed\ngain +%.3f", as.numeric(gains)))
xr  <- range(c(plot_df$value, obs_df$gain)); xpad <- diff(xr) * 0.10
p <- ggplot(plot_df, aes(value)) +
  geom_histogram(bins = 40, fill = "grey75", color = "white") +
  geom_vline(xintercept = 0, color = "grey35", linetype = "dashed") +
  geom_vline(data = obs_df, aes(xintercept = gain), color = "#E63946", linewidth = 1) +
  geom_text(data = obs_df, aes(x = gain, y = Inf, label = lab),
            color = "#E63946", size = 2.6, hjust = 1.1, vjust = 1.4) +
  facet_wrap(~ metric, ncol = 2, scales = "free_y") +
  coord_cartesian(xlim = c(xr[1] - xpad, xr[2] + xpad)) +
  labs(x = "Concordance minus its own random-pairing null (pairing-specific gain)",
       y = "Permutations (stratified null)",
       title = "Pairing-specific gain on a common scale: observed (red) vs random-pairing null (grey; dashed = null)") +
  theme_bw(base_size = 10)
ggsave(file.path(fig_dir, "fig_permutation_null.pdf"), p, width = 6.5, height = 3.4)
cat("Saved: nature_manuscript/figures/fig_permutation_null.pdf\n\nDone.\n")
