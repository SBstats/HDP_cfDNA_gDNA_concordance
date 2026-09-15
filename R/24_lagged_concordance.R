###############################################################################
## R/24_lagged_concordance.R
##
## LAGGED concordance (storyline item 8): does cfDNA at cycle t agree best with
## gDNA at the SAME cycle (lag 0), or does it ANTICIPATE (lag +1: cfDNA_t vs
## gDNA_{t+1}) or LAG BEHIND (lag -1: cfDNA_t vs gDNA_{t-1}) the marrow state?
##
## For each within-patient ordered pair of timepoints (t, t') with t' the next
## sampled cycle, we can compare:
##   lag  0 : cor(Ycf_it , Yg_it)                 [contemporaneous, per obs]
##   lag +1 : cor(Ycf_it , Yg_it')   (cf leads)   [cfDNA anticipates next marrow]
##   lag -1 : cor(Ycf_it', Yg_it)    (cf trails)  [cfDNA reflects previous marrow]
## Reported at the GENE level (model-free). A pairing-shuffle null (break the
## patient matching) is built per lag so each lag is judged against its own null.
##
## EXPLORATORY: with only 5 cycles the number of usable lagged pairs is small,
## so this is a directional check, not a powered test.
##
## Outputs: results/lagged_concordance.csv
##          figures/fig_lagged_concordance.pdf (+ nature_manuscript/figures/)
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
n_perm <- 2000

## ordered consecutive within-patient timepoint pairs (obs index a earlier, b later)
pairs <- do.call(rbind, lapply(sort(unique(pvec)), function(pt) {
  idx <- which(pvec == pt); idx <- idx[order(tvec[idx])]
  if (length(idx) < 2) return(NULL)
  do.call(rbind, lapply(1:(length(idx) - 1), function(j) {
    a <- idx[j]; b <- idx[j + 1]
    if (tvec[b] == tvec[a]) return(NULL)
    data.frame(patient = pt, a = a, b = b)
  }))
}))
Mpair <- nrow(pairs)
cat(sprintf("consecutive within-patient timepoint pairs: %d\n", Mpair))

## lag definitions as (cf index, g index) per pair
lag_index <- list(
  `-1 (cfDNA trails marrow)` = data.frame(cf = pairs$b, g = pairs$a),  # cf_t' vs g_t
  ` 0 (contemporaneous)`     = data.frame(cf = pairs$a, g = pairs$a),  # cf_t  vs g_t
  `+1 (cfDNA anticipates)`   = data.frame(cf = pairs$a, g = pairs$b)   # cf_t  vs g_t'
)

## observed mean gene-level r per lag + pairing-shuffle null (break patient match)
lag_result <- lapply(names(lag_index), function(lag) {
  ix <- lag_index[[lag]]
  obs <- mean(sapply(1:Mpair, function(m) cor(Ycf[, ix$cf[m]], Yg[, ix$g[m]])))
  # precompute the Mpair x Mpair cross matrix once: C[i,j] = cor(cf of pair i, g of pair j)
  C <- cor(Ycf[, ix$cf, drop = FALSE], Yg[, ix$g, drop = FALSE])
  nullv <- replicate(n_perm, { repeat { q <- sample(Mpair); if (all(q != 1:Mpair)) break }; mean(C[cbind(1:Mpair, q)]) })
  nm <- mean(nullv); sdn <- sd(nullv)
  data.frame(lag = lag, n_pairs = Mpair, observed = obs, null_mean = nm,
             gain = obs - nm, z = (obs - nm) / sdn,
             p_perm = (1 + sum(nullv >= obs)) / (1 + n_perm))
})
R <- do.call(rbind, lag_result)
R[ , c("observed","null_mean","gain","z","p_perm")] <- round(R[ , c("observed","null_mean","gain","z","p_perm")], 4)
R$lag <- factor(R$lag, levels = c("-1 (cfDNA trails marrow)", " 0 (contemporaneous)", "+1 (cfDNA anticipates)"))
write.csv(R, file.path(res_dir, "lagged_concordance.csv"), row.names = FALSE)
cat("\n=== Lagged gene-level concordance (exploratory) ===\n")
print(R, row.names = FALSE)

best <- R$lag[which.max(R$gain)]
cat(sprintf("\nLargest pairing-specific gain at lag: %s\n", as.character(best)))

## figure: pairing-specific gain (observed - null) per lag
p <- ggplot(R, aes(lag, gain, fill = lag)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_hline(yintercept = 0, color = "grey50") +
  geom_text(aes(label = sprintf("obs %.3f\nnull %.3f\np=%.3f", observed, null_mean, p_perm)),
            vjust = -0.2, size = 2.6) +
  scale_fill_manual(values = c("-1 (cfDNA trails marrow)" = "#A8DADC",
                               " 0 (contemporaneous)" = "#457B9D",
                               "+1 (cfDNA anticipates)" = "#E63946")) +
  labs(x = "cfDNA-to-gDNA temporal lag", y = "Pairing-specific gain (observed - random-pairing null)",
       title = "Lagged cfDNA-gDNA gene-level concordance",
       subtitle = "Does cfDNA agree best with same-cycle, previous-cycle, or next-cycle marrow? (exploratory, T<=5)") +
  expand_limits(y = max(R$gain) * 1.25) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 8))
for (d in fig_dirs) ggsave(file.path(d, "fig_lagged_concordance.pdf"), p, width = 7.5, height = 4.2)

cat("Saved: results/lagged_concordance.csv, fig_lagged_concordance.pdf (figures/ + nature_manuscript/figures/)\n")
