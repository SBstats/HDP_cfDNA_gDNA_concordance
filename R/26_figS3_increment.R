###############################################################################
## R/26_figS3_increment.R
## Applied paper Supplementary Figure S3 — within-patient longitudinal increment
## concordance between cfDNA and gDNA.
##
## (a) permutation-null histogram of the model-free gene-level increment
##     correlation with the observed value marked;
## (b) three-estimator comparison (model-free gene-level; model-based M1 and M0
##     dominant-subclone weight) with 95% CrIs and pairing-shuffle nulls.
##
## Inputs (all produced by R/15_tier2_within_patient.R):
##   results/tier2/genelevel_increment_null.rds   (observed + 5000 null draws)
##   results/tier2/tier2_genelevel_delta.csv       (gene-level summary)
##   results/tier2/tier2_subclonal_delta.csv       (M1/M0 summary)
## Output: nature_manuscript/figures/figS3_longitudinal_increment.pdf
## Usage: Rscript R/26_figS3_increment.R
###############################################################################

suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
base_dir <- getwd()
res <- file.path(base_dir, "results")
fig <- file.path(base_dir, "nature_manuscript", "figures"); dir.create(fig, showWarnings = FALSE, recursive = TRUE)

nullobj <- readRDS(file.path(res, "tier2", "genelevel_increment_null.rds"))
gl <- read.csv(file.path(res, "tier2", "tier2_genelevel_delta.csv"), stringsAsFactors = FALSE)
sc <- read.csv(file.path(res, "tier2", "tier2_subclonal_delta.csv"), stringsAsFactors = FALSE)

obs      <- nullobj$observed
null_dr  <- nullobj$null_draws
p_gl     <- gl$p[1]
lvl_corr <- if (!is.null(nullobj$level_corr) && is.finite(nullobj$level_corr)) nullobj$level_corr else NA_real_

## ---- Panel (a): permutation null histogram ----
pa <- ggplot(data.frame(x = null_dr), aes(x)) +
  geom_histogram(bins = 40, fill = "grey75", colour = "grey55", linewidth = 0.2) +
  geom_vline(xintercept = obs, colour = "firebrick", linewidth = 1) +
  annotate("text", x = obs, y = Inf, vjust = 1.5, hjust = -0.05,
           label = sprintf("observed = %.3f\n(P = %s)", obs, format(p_gl, scientific = TRUE)),
           colour = "firebrick", size = 3) +
  annotate("text", x = mean(null_dr), y = Inf, vjust = 1.5, hjust = 1.05,
           label = sprintf("null mean = %.3f", mean(null_dr)), colour = "grey30", size = 3) +
  labs(x = "Mean within-patient increment correlation",
       y = "Permutations",
       subtitle = if (is.finite(lvl_corr))
         sprintf("(a) First-differencing collapses the level correlation (~%.2f) to the null (~%.2f)", lvl_corr, mean(null_dr))
       else "(a) Model-free gene-level increment vs pairing-shuffle null") +
  theme_bw(base_size = 10)

## ---- Panel (b): three-estimator comparison ----
est <- data.frame(
  estimator = factor(c("Model-free\n(gene-level)", "Model-based\nM1 (tracking)", "Model-based\nM0 (independent)"),
                     levels = c("Model-free\n(gene-level)", "Model-based\nM1 (tracking)", "Model-based\nM0 (independent)")),
  value = c(gl$delta_true[1], sc$post_mean_r[1], sc$post_mean_r[2]),
  lo    = c(NA, sc$cri_lo[1], sc$cri_lo[2]),
  hi    = c(NA, sc$cri_hi[1], sc$cri_hi[2]),
  null  = c(gl$null[1], sc$null_mean[1], sc$null_mean[2])
)
pb <- ggplot(est, aes(value, estimator)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15, na.rm = TRUE, colour = "steelblue4") +
  geom_point(size = 3, colour = "steelblue4") +
  geom_point(aes(x = null), shape = 4, size = 2.5, colour = "grey40") +
  labs(x = "Increment co-movement correlation (cfDNA vs gDNA)", y = NULL,
       subtitle = "(b) Estimates (points), 95% CrIs (bars), pairing-shuffle nulls (crosses)") +
  theme_bw(base_size = 10)

g <- pa / pb + plot_layout(heights = c(1, 0.9))
ggsave(file.path(fig, "figS3_longitudinal_increment.pdf"), g, width = 7, height = 7)
cat("Wrote nature_manuscript/figures/figS3_longitudinal_increment.pdf\n")
cat(sprintf("  panel a: observed %.3f vs null %.3f (P=%s); panel b: MF %.3f / M1 %.3f / M0 %.3f\n",
            obs, mean(null_dr), format(p_gl), gl$delta_true[1], sc$post_mean_r[1], sc$post_mean_r[2]))
