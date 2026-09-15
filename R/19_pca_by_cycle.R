###############################################################################
## R/19_pca_by_cycle.R
## Per-cycle version of the paired-sample PCA (previously fig5c, all cycles in
## one panel). A SINGLE PCA is computed on the combined cfDNA+gDNA VST matrix
## (same as R/07_figures.R), then the paired points/lines are faceted into one
## row with a separate column per treatment timepoint, so within-pair proximity
## can be read cycle by cycle in a common coordinate system.
##
## Output: nature_manuscript/figures/fig5c_pca_paired_bycycle.pdf
###############################################################################

suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
base_dir <- getwd()  # run from project root
setwd(base_dir)
res_dir <- file.path(base_dir, "results")
fig_dir <- file.path(base_dir, "nature_manuscript", "figures")

md <- readRDS(file.path(res_dir, "model_data.rds"))
Ycf <- md$Y_cf; Yg <- md$Y_g; n <- md$N_obs
tp_levels <- c("Screening", "C4", "C8", "C18", "3 YR F/U")
tp <- factor(tp_levels[md$time], levels = tp_levels)

## ONE PCA on the combined matrix (cfDNA columns 1..n, gDNA columns n+1..2n)
Y_combined <- cbind(Ycf, Yg)
pca <- prcomp(t(Y_combined), center = TRUE, scale. = TRUE)
ve  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

## points (both sources) and paired segments, each tagged with its timepoint
pts <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                  Source = c(rep("cfDNA", n), rep("gDNA", n)),
                  Timepoint = rep(tp, 2))
seg <- data.frame(PC1_cf = pca$x[1:n, 1],           PC2_cf = pca$x[1:n, 2],
                  PC1_g  = pca$x[(n+1):(2*n), 1],   PC2_g  = pca$x[(n+1):(2*n), 2],
                  Timepoint = tp)

p <- ggplot() +
  geom_segment(data = seg, aes(x = PC1_cf, y = PC2_cf, xend = PC1_g, yend = PC2_g),
               color = "grey70", linewidth = 0.3) +
  geom_point(data = pts, aes(PC1, PC2, color = Source), size = 1) +
  facet_grid(. ~ Timepoint) +
  scale_color_manual(values = c(cfDNA = "#E63946", gDNA = "#457B9D")) +
  labs(x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)"),
       title = "Paired cfDNA-gDNA samples in PC space, by treatment timepoint") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "fig5c_pca_paired_bycycle.pdf"), p, width = 11, height = 2.9)
cat("Saved: nature_manuscript/figures/fig5c_pca_paired_bycycle.pdf\n")
cat(sprintf("PC1 %.1f%%, PC2 %.1f%%; per-cycle n: %s\n", ve[1], ve[2],
            paste(sprintf("%s=%d", tp_levels, as.integer(table(tp)[tp_levels])), collapse=", ")))
