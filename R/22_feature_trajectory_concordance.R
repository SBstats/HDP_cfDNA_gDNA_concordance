###############################################################################
## R/22_feature_trajectory_concordance.R
##
## FEATURE-LEVEL trajectory concordance (storyline item 4), DESCRIPTIVE.
##
## For each patient and each of the 11,837 genes, correlate the gene's cfDNA
## time-trajectory with its gDNA time-trajectory across that patient's shared
## timepoints:  r_ij = cor( Ycf[j, patient i's timepoints], Yg[j, same] ).
## We then summarise the DISTRIBUTION of r_ij across genes (pooled over patients)
## and compare it to a pairing-shuffle null (cf trajectory of a gene matched to a
## DIFFERENT patient's g trajectory of the same gene).
##
## IMPORTANT CAVEAT (stated in the figure and the manuscript): with only T<=5
## timepoints per patient, each individual per-gene correlation is very noisy
## (3 d.f. at T=5). This analysis is therefore DESCRIPTIVE / hypothesis-
## generating about the *shape of the distribution*, NOT a precise per-gene
## estimate. Restricted to patients with >=3 shared timepoints.
##
## Outputs: results/feature_trajectory_concordance.csv          (summary stats)
##          results/feature_trajectory_per_patient.csv          (per-patient means)
##          figures/fig_feature_trajectory_distribution.pdf
##          nature_manuscript/figures/fig_feature_trajectory_distribution.pdf
###############################################################################

set.seed(20260808)
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
base_dir <- getwd()  # run from project root
setwd(base_dir)
res_dir <- file.path(base_dir, "results")
fig_dirs <- c(file.path(base_dir, "figures"), file.path(base_dir, "nature_manuscript", "figures"))
for (d in fig_dirs) dir.create(d, showWarnings = FALSE, recursive = TRUE)

md  <- readRDS(file.path(res_dir, "model_data.rds"))
Yg  <- md$Y_g; Ycf <- md$Y_cf; tvec <- md$time; pvec <- md$patient
MIN_TP <- 3

## vectorised per-gene correlation of two [p x T] matrices (over columns = time)
rowwise_cor <- function(A, B) {
  Ac <- A - rowMeans(A); Bc <- B - rowMeans(B)
  num <- rowSums(Ac * Bc)
  den <- sqrt(rowSums(Ac^2) * rowSums(Bc^2))
  r <- num / den
  r[!is.finite(r)] <- NA           # genes with zero variance in a source
  r
}

pts <- sort(unique(pvec))
usable <- pts[sapply(pts, function(pt) length(unique(tvec[pvec == pt])) >= MIN_TP)]
cat(sprintf("patients with >=%d shared timepoints: %d of %d\n", MIN_TP, length(usable), length(pts)))

## observed: per-patient [p x T] slices (order columns by time)
slices_cf <- lapply(usable, function(pt) { ii <- which(pvec == pt); ii <- ii[order(tvec[ii])]; Ycf[, ii, drop = FALSE] })
slices_g  <- lapply(usable, function(pt) { ii <- which(pvec == pt); ii <- ii[order(tvec[ii])]; Yg[,  ii, drop = FALSE] })
names(slices_cf) <- names(slices_g) <- usable

r_obs <- do.call(cbind, lapply(seq_along(usable), function(i) rowwise_cor(slices_cf[[i]], slices_g[[i]])))  # [p x P]
colnames(r_obs) <- usable
r_obs_vec <- as.vector(r_obs)

## null: match each patient's cf trajectory to a DIFFERENT patient's g trajectory
## (only when the two have the same number of timepoints, so the T dimension lines up)
Pn <- length(usable); ncols <- sapply(slices_cf, ncol)
null_means <- replicate(200, {
  perm <- sample(Pn); ok <- perm != (1:Pn) & ncols[perm] == ncols
  if (!any(ok)) return(NA)
  rr <- unlist(lapply(which(ok), function(i) rowwise_cor(slices_cf[[i]], slices_g[[perm[i]]])))
  mean(rr, na.rm = TRUE)
})

## per-patient observed mean r
per_pt <- data.frame(patient = usable,
                     n_timepoints = ncols,
                     mean_r = round(colMeans(r_obs, na.rm = TRUE), 4),
                     frac_gt0 = round(colMeans(r_obs > 0, na.rm = TRUE), 4))
write.csv(per_pt, file.path(res_dir, "feature_trajectory_per_patient.csv"), row.names = FALSE)

## summary
summ <- data.frame(
  n_patients = Pn, min_timepoints = MIN_TP, n_genes = nrow(Yg),
  mean_r        = round(mean(r_obs_vec, na.rm = TRUE), 4),
  median_r      = round(median(r_obs_vec, na.rm = TRUE), 4),
  frac_gt_0     = round(mean(r_obs_vec > 0, na.rm = TRUE), 4),
  frac_gt_0.5   = round(mean(r_obs_vec > 0.5, na.rm = TRUE), 4),
  frac_gt_0.8   = round(mean(r_obs_vec > 0.8, na.rm = TRUE), 4),
  null_mean_r   = round(mean(null_means, na.rm = TRUE), 4),
  obs_minus_null = round(mean(r_obs_vec, na.rm = TRUE) - mean(null_means, na.rm = TRUE), 4)
)
write.csv(summ, file.path(res_dir, "feature_trajectory_concordance.csv"), row.names = FALSE)
cat("\n=== Feature-level trajectory concordance (descriptive) ===\n")
print(summ, row.names = FALSE)

## figure: observed distribution of per-gene r vs the null distribution
df <- rbind(
  data.frame(r = r_obs_vec,                       what = "Observed (true pairing)"),
  data.frame(r = as.vector(sapply(1:Pn, function(i) {  # a matched-null sample for the histogram
      j <- ((i %% Pn) + 1); if (ncols[j] != ncols[i]) return(rep(NA, nrow(Yg))); rowwise_cor(slices_cf[[i]], slices_g[[j]]) })),
    what = "Random-pairing null")
)
df <- df[is.finite(df$r), ]
p <- ggplot(df, aes(r, fill = what, color = what)) +
  geom_density(alpha = 0.35, linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "grey40", linetype = "dashed") +
  scale_fill_manual(values = c("Observed (true pairing)" = "#E63946", "Random-pairing null" = "grey55")) +
  scale_color_manual(values = c("Observed (true pairing)" = "#E63946", "Random-pairing null" = "grey40")) +
  labs(x = "Per-gene cfDNA-gDNA trajectory correlation (across a patient's timepoints)",
       y = "Density", fill = NULL, color = NULL,
       title = "Distribution of gene-level trajectory concordance across ~11,837 genes",
       subtitle = sprintf("Descriptive only: T<=5 timepoints => each per-gene r is noisy (%d patients with >=%d timepoints)",
                          Pn, MIN_TP)) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 8))
for (d in fig_dirs) ggsave(file.path(d, "fig_feature_trajectory_distribution.pdf"), p, width = 8, height = 4.2)

cat("\nSaved: results/feature_trajectory_concordance.csv, feature_trajectory_per_patient.csv,\n")
cat("       fig_feature_trajectory_distribution.pdf (figures/ + nature_manuscript/figures/)\n")
