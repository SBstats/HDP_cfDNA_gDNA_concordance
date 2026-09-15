###############################################################################
## R/21_concordance_hierarchy.R
##
## Assembles the cfDNA-gDNA concordance HIERARCHY into a single ordered table and
## a single "ladder" figure. The storyline (from simplest/least-specific to
## strongest/most pairing-specific) is:
##
##   (1) LEVEL          cor of per-patient MEAN profiles           [looks alike]
##   (2) SNAPSHOT       gene-level r at each obs (cross-sectional) [looks alike]
##   (3) TRAJECTORY     cor of the mean-removed DEVIATIONS         [tracks in time]
##   (4) CHANGE (gene)  cor of first-differences DeltaY            [tracks in time]
##   (5) LATENT LEVEL   subclonal weight cross-source correlation   [shared biology]
##   (6) LATENT CHANGE  cor of subclonal weight increments Dw      [tracks biology]
##
## Every rung is reported as OBSERVED vs its own PAIRING-SHUFFLE / random-pairing
## NULL (permutation negative control), so the reader can see directly that
## "level agreement" is largely non-specific (observed ~= null) while the
## change/tracking rungs are pairing-specific (observed clears its null).
##
## Runs entirely on EXISTING fits (no re-fit). Reuses helpers from R/01 and the
## within-patient transition / stratified-shuffle logic from R/15.
##
## Outputs: results/concordance_hierarchy.csv
##          figures/fig_concordance_hierarchy.pdf
##          nature_manuscript/figures/fig_concordance_hierarchy.pdf
###############################################################################

set.seed(20260808)
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
base_dir <- getwd()  # run from project root
setwd(base_dir)
suppressWarnings(suppressMessages(source(file.path(base_dir, "R", "01_lib_core.R"))))
res_dir <- file.path(base_dir, "results")
fig_dirs <- c(file.path(base_dir, "figures"), file.path(base_dir, "nature_manuscript", "figures"))
for (d in fig_dirs) dir.create(d, showWarnings = FALSE, recursive = TRUE)

md  <- readRDS(file.path(res_dir, "model_data.rds"))
Yg  <- md$Y_g; Ycf <- md$Y_cf; N <- md$N_obs; tvec <- md$time; pvec <- md$patient
n_perm <- 2000

## ---------------------------------------------------------------------------
## within-patient transitions + stratified pairing-shuffle  (as in R/15:33-53)
## ---------------------------------------------------------------------------
build_transitions <- function() {
  do.call(rbind, lapply(unique(pvec), function(pt) {
    idx <- which(pvec == pt); idx <- idx[order(tvec[idx])]
    if (length(idx) < 2) return(NULL)
    do.call(rbind, lapply(1:(length(idx) - 1), function(j) {
      a <- idx[j]; b <- idx[j + 1]
      if (tvec[b] == tvec[a]) return(NULL)
      data.frame(patient = pt, a = a, b = b, stratum = paste0(tvec[a], "->", tvec[b]))
    }))
  }))
}
TR <- build_transitions(); Mtr <- nrow(TR)
cat(sprintf("within-patient transitions: M=%d across %d patients\n",
            Mtr, length(unique(TR$patient))))

## ---------------------------------------------------------------------------
## RUNGS 1 & 3 -- level vs trajectory (deviation) decomposition
## For each patient with >=2 timepoints:
##   level_i = cor( mean_t Ycf_it , mean_t Yg_it )          (patient-mean profiles)
##   traj_i  = cor( vec(Ycf - patient-mean), vec(Yg - patient-mean) )  (deviations)
## Null = shuffle which patient's cf mean/deviation is matched to which g.
## ---------------------------------------------------------------------------
pts_long <- unique(pvec[ave(pvec, pvec, FUN = length) >= 2])  # patients with >=2 obs
pts_long <- sort(unique(TR$patient))                          # (same set: >=2 timepoints)
cf_mean <- sapply(pts_long, function(pt) rowMeans(Ycf[, which(pvec == pt), drop = FALSE]))
g_mean  <- sapply(pts_long, function(pt) rowMeans(Yg[,  which(pvec == pt), drop = FALSE]))  # [p x P]

level_i <- sapply(seq_along(pts_long), function(i) cor(cf_mean[, i], g_mean[, i]))

# deviations: within each patient subtract that patient's mean, then vectorise per patient
dev_cor <- sapply(pts_long, function(pt) {
  ii <- which(pvec == pt)
  Dcf <- Ycf[, ii, drop = FALSE] - rowMeans(Ycf[, ii, drop = FALSE])
  Dg  <- Yg[,  ii, drop = FALSE] - rowMeans(Yg[,  ii, drop = FALSE])
  cor(as.vector(Dcf), as.vector(Dg))
})

# per-patient deviation matrices kept for the shuffle null
dev_list_cf <- lapply(pts_long, function(pt) { ii <- which(pvec == pt); Ycf[, ii, drop = FALSE] - rowMeans(Ycf[, ii, drop = FALSE]) })
dev_list_g  <- lapply(pts_long, function(pt) { ii <- which(pvec == pt); Yg[,  ii, drop = FALSE] - rowMeans(Yg[,  ii, drop = FALSE]) })

P <- length(pts_long)
level_null <- replicate(n_perm, {
  perm <- sample(P); if (all(perm == 1:P)) perm <- rev(perm)
  mean(sapply(1:P, function(i) cor(cf_mean[, i], g_mean[, perm[i]])))
})
traj_null <- replicate(n_perm, {
  perm <- sample(P); if (all(perm == 1:P)) perm <- rev(perm)
  # only patients whose shuffled partner has the SAME number of timepoints can be vectorised
  vals <- sapply(1:P, function(i) {
    if (ncol(dev_list_cf[[i]]) != ncol(dev_list_g[[perm[i]]])) return(NA)
    cor(as.vector(dev_list_cf[[i]]), as.vector(dev_list_g[[perm[i]]]))
  })
  mean(vals, na.rm = TRUE)
})

## ---------------------------------------------------------------------------
## RUNG 2 -- snapshot gene-level r (cross-sectional, all 84 obs)  [item 1]
## Null = random within-timepoint derangement of the cf<->g pairing.
## ---------------------------------------------------------------------------
snap_obs <- mean(sapply(1:N, function(i) cor(Yg[, i], Ycf[, i])))
strat_derange_idx <- function() {
  m <- 1:N
  for (s in unique(tvec)) {
    ii <- which(tvec == s)
    if (length(ii) > 1) { repeat { p <- sample(ii); if (all(p != ii)) break }; m[ii] <- p }
  }
  m
}
snap_null <- replicate(n_perm, { mm <- strat_derange_idx(); mean(sapply(1:N, function(i) cor(Yg[, i], Ycf[, mm[i]]))) })

## ---------------------------------------------------------------------------
## RUNG 4 -- gene-level CHANGE (first-difference) concordance  [item 3, as R/15]
## ---------------------------------------------------------------------------
dG  <- sapply(1:Mtr, function(m) Yg[, TR$b[m]]  - Yg[, TR$a[m]])
dCF <- sapply(1:Mtr, function(m) Ycf[, TR$b[m]] - Ycf[, TR$a[m]])
Cd  <- cor(dG, dCF)
change_obs <- mean(diag(Cd))
change_null <- replicate(n_perm, {
  q <- 1:Mtr
  for (s in unique(TR$stratum)) { ii <- which(TR$stratum == s); if (length(ii) > 1) q[ii] <- sample(ii) }
  mean(Cd[cbind(1:Mtr, q)])
})

## ---------------------------------------------------------------------------
## RUNGS 5 & 6 -- latent (subclonal weight) level & change  [items 5,6]
## canonical posterior-mean weights (relabel/align/merge) as in R/14:34-41
## ---------------------------------------------------------------------------
m1 <- readRDS(file.path(res_dir, "krd_M1_fit.rds"))
for (i in seq_along(m1$chains)) m1$chains[[i]] <- relabel_by_weight(m1$chains[[i]])
m1$chains <- align_chains_by_signature(m1$chains)
merged <- merge_chain_samples(m1$chains); K <- ncol(merged$omega_g_trace) / N
Wg  <- sapply(1:K, function(k) colMeans(merged$omega_g_trace[,  ((k-1)*N + 1):(k*N), drop = FALSE]))  # [N x K]
Wcf <- sapply(1:K, function(k) colMeans(merged$omega_cf_trace[, ((k-1)*N + 1):(k*N), drop = FALSE]))

# latent LEVEL = cross-source correlation of the posterior-mean subclonal weights
# across observations, averaged over occupied components. (The L1 concordance
# metric was removed by design; the latent-level rung now uses the retained
# posterior-correlation measure of omega.)
occ <- which((colMeans(Wg) + colMeans(Wcf)) / 2 > 0.01)         # occupied components
if (length(occ) < 2) occ <- order(colMeans(Wg), decreasing = TRUE)[1:2]
comp_corr <- function(g, cf) mean(sapply(occ, function(k)
  if (sd(g[, k]) > 1e-9 && sd(cf[, k]) > 1e-9) cor(g[, k], cf[, k]) else NA_real_), na.rm = TRUE)
latlevel_obs <- comp_corr(Wg, Wcf)
# null: random within-timepoint re-pairing of the cfDNA weight rows
latlevel_null <- replicate(n_perm, {
  mm <- strat_derange_idx()
  comp_corr(Wg, Wcf[mm, , drop = FALSE])
})

# latent CHANGE = cor of dominant-subclone weight increments (posterior-mean weights)
ref <- which.max(colMeans(Wg))                       # gDNA's dominant signature
dwg  <- Wg[TR$b, ref]  - Wg[TR$a, ref]
dwcf <- Wcf[TR$b, ref] - Wcf[TR$a, ref]
latchange_obs <- cor(dwg, dwcf)
latchange_null <- replicate(n_perm, {
  y <- dwcf
  for (s in unique(TR$stratum)) { ii <- which(TR$stratum == s); if (length(ii) > 1) y[ii] <- sample(dwcf[ii]) }
  cor(dwg, y)
})

## ---------------------------------------------------------------------------
## assemble the ladder
## ---------------------------------------------------------------------------
rung <- function(name, kind, obs, nullvec) {
  nm <- mean(nullvec, na.rm = TRUE); sdn <- sd(nullvec, na.rm = TRUE)
  data.frame(metric = name, kind = kind, observed = obs, null_mean = nm,
             gain = obs - nm, null_sd = sdn,
             z = (obs - nm) / sdn,
             p_perm = (1 + sum(nullvec >= obs, na.rm = TRUE)) / (1 + sum(is.finite(nullvec))))
}
H <- rbind(
  rung("1. Level (patient-mean profiles)",        "looks alike", mean(level_i),  level_null),
  rung("2. Snapshot (gene-level, per obs)",        "looks alike", snap_obs,      snap_null),
  rung("3. Trajectory (mean-removed deviations)",  "tracks",      mean(dev_cor), traj_null),
  rung("4. Change (gene-level first differences)", "tracks",      change_obs,    change_null),
  rung("5. Latent level (subclonal weight corr.)", "shared biology", latlevel_obs, latlevel_null),
  rung("6. Latent change (dominant-weight incr.)", "tracks biology", latchange_obs, latchange_null)
)
H[ , c("observed","null_mean","gain","null_sd","z","p_perm")] <-
  round(H[ , c("observed","null_mean","gain","null_sd","z","p_perm")], 4)
H$metric <- factor(H$metric, levels = rev(H$metric))   # top rung at top of plot

write.csv(H, file.path(res_dir, "concordance_hierarchy.csv"), row.names = FALSE)
cat("\n=== Concordance hierarchy (observed vs pairing-shuffle null) ===\n")
print(H[order(H$metric, decreasing = TRUE), c("metric","observed","null_mean","gain","p_perm","z")],
      row.names = FALSE)

## ---------------------------------------------------------------------------
## ladder figure: observed (filled) vs null (open) per rung, with gain arrow
## ---------------------------------------------------------------------------
plotd <- rbind(
  data.frame(metric = H$metric, value = H$observed, what = "Observed (true pairing)"),
  data.frame(metric = H$metric, value = H$null_mean, what = "Random-pairing null")
)
p <- ggplot() +
  geom_segment(data = H, aes(x = null_mean, xend = observed, y = metric, yend = metric),
               color = "grey60", linewidth = 0.5) +
  geom_point(data = plotd, aes(value, metric, color = what, shape = what), size = 3) +
  scale_color_manual(values = c("Observed (true pairing)" = "#E63946",
                                "Random-pairing null" = "grey45")) +
  scale_shape_manual(values = c("Observed (true pairing)" = 16, "Random-pairing null" = 1)) +
  labs(x = "Concordance (observed vs its own random-pairing null)", y = NULL,
       color = NULL, shape = NULL,
       title = "A hierarchy of cfDNA-gDNA concordance",
       subtitle = "Level agreement is largely non-specific; tracking (change) rungs clear their null") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 8.5))
for (d in fig_dirs) ggsave(file.path(d, "fig_concordance_hierarchy.pdf"), p, width = 8.5, height = 4)

cat("\nSaved: results/concordance_hierarchy.csv and fig_concordance_hierarchy.pdf (figures/ + nature_manuscript/figures/)\n")
cat(sprintf("\nReproduction check -- gene-level change rung: obs=%.3f null=%.3f (expect ~0.138 / ~0.028)\n",
            change_obs, mean(change_null)))
cat(sprintf("Latent level (subclonal weight correlation): %.3f\n", latlevel_obs))
