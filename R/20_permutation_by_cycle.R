###############################################################################
## R/20_permutation_by_cycle.R
## Per-cycle version of the permutation negative control (R/14 pools all
## timepoints; here the pairing is broken and the null built SEPARATELY within
## each treatment timepoint). For every (metric, cycle) we report the observed
## pairing-specific gain (observed value minus that cycle's own random-pairing
## null mean) against the within-cycle null distribution.
##
## NOTE: each cycle has only 13-20 paired observations, so per-cycle nulls are
## noisier (wider) and less powerful than the pooled 84-observation null.
##
## Outputs: results/permutation_by_cycle.csv
##          nature_manuscript/figures/fig_permutation_null_bycycle.pdf
###############################################################################

set.seed(20260718)
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
base_dir <- getwd()  # run from project root
setwd(base_dir)
source(file.path(base_dir, "R", "01_lib_core.R"))
res_dir <- file.path(base_dir, "results"); fig_dir <- file.path(base_dir, "nature_manuscript", "figures")

md <- readRDS(file.path(res_dir, "model_data.rds"))
Yg <- md$Y_g; Ycf <- md$Y_cf; N <- md$N_obs; tvec <- md$time
tp_levels <- c("Screening", "C4", "C8", "C18", "3 YR F/U")

## canonical posterior-mean subclonal weights per obs (relabel/align/merge) -- as in R/14
m1 <- readRDS(file.path(res_dir, "krd_M1_fit.rds"))
for (i in seq_along(m1$chains)) m1$chains[[i]] <- relabel_by_weight(m1$chains[[i]])
m1$chains <- align_chains_by_signature(m1$chains)
merged <- merge_chain_samples(m1$chains); K <- ncol(merged$omega_g_trace) / N
Wg  <- sapply(1:K, function(k) colMeans(merged$omega_g_trace[,  ((k-1)*N + 1):(k*N), drop=FALSE]))
Wcf <- sapply(1:K, function(k) colMeans(merged$omega_cf_trace[, ((k-1)*N + 1):(k*N), drop=FALSE]))  # [N x K]

metrics <- c(gene_r = "Gene-level Pearson r",
             post_r = "Posterior correlation")   # L1 removed by design
n_perm <- 2000
derange <- function(m) { if (m < 2) return(1L); repeat { p <- sample(m); if (all(p != 1:m)) break }; p }

## For a cycle's observation set S, precompute the |S| x |S| cross-source
## matrices so permutations are fast index lookups; C[a,b] pairs gDNA obs S[a]
## with cfDNA obs S[b]. Observed = mean(diag); null perm = mean(C[cbind(1:s, q)]).
per_cycle <- function(S) {
  s <- length(S)
  Cg <- suppressWarnings(cor(Yg[, S, drop=FALSE], Ycf[, S, drop=FALSE]))                 # gene-level r
  Cp <- suppressWarnings(cor(t(Wg[S, , drop=FALSE]), t(Wcf[S, , drop=FALSE])))            # posterior corr
  mats <- list(gene_r = Cg, post_r = Cp)
  out <- lapply(names(mats), function(mn) {
    Cm <- mats[[mn]]; obs <- mean(diag(Cm), na.rm=TRUE)
    nulls <- replicate(n_perm, { q <- derange(s); mean(Cm[cbind(1:s, q)], na.rm=TRUE) })
    list(obs = obs, nulls = nulls, mn = mn)
  })
  names(out) <- names(mats); out
}

summ <- list(); plot_df <- list(); obs_df <- list()
for (ti in seq_along(tp_levels)) {
  S <- which(tvec == ti); if (length(S) < 2) next
  pc <- per_cycle(S)
  for (mn in names(metrics)) {
    o <- pc[[mn]]; nm <- mean(o$nulls, na.rm=TRUE); sdn <- sd(o$nulls, na.rm=TRUE)
    p_perm <- (1 + sum(o$nulls >= o$obs, na.rm=TRUE)) / (1 + sum(is.finite(o$nulls)))
    summ[[length(summ)+1]] <- data.frame(cycle = tp_levels[ti], n_pairs = length(S),
      metric = metrics[mn], observed = o$obs, null_mean = nm, null_sd = sdn,
      gain = o$obs - nm, p_perm = p_perm, z = (o$obs - nm)/sdn)
    plot_df[[length(plot_df)+1]] <- data.frame(cycle = factor(tp_levels[ti], levels=tp_levels),
      metric = factor(metrics[mn], levels=metrics), value = o$nulls - nm)
    obs_df[[length(obs_df)+1]] <- data.frame(cycle = factor(tp_levels[ti], levels=tp_levels),
      metric = factor(metrics[mn], levels=metrics), gain = o$obs - nm)
  }
}
summ <- do.call(rbind, summ); plot_df <- do.call(rbind, plot_df); obs_df <- do.call(rbind, obs_df)
write.csv(summ, file.path(res_dir, "permutation_by_cycle.csv"), row.names = FALSE)
cat("\n=== Per-cycle permutation negative control ===\n"); print(summ, row.names=FALSE, digits=3)

## grid figure: rows = metric, columns = cycle; shared x (pairing-specific gain), free y
p <- ggplot(plot_df, aes(value)) +
  geom_histogram(bins = 30, fill = "grey75", color = "white") +
  geom_vline(xintercept = 0, color = "grey35", linetype = "dashed") +
  geom_vline(data = obs_df, aes(xintercept = gain), color = "#E63946", linewidth = 0.8) +
  facet_grid(metric ~ cycle, scales = "free_y") +
  labs(x = "Concordance minus its own within-cycle random-pairing null (pairing-specific gain)",
       y = "Permutations",
       title = "Per-cycle permutation negative control: observed gain (red) vs within-cycle null (grey; dashed = null)") +
  theme_bw(base_size = 8) + theme(strip.text.y = element_text(angle = 0))
ggsave(file.path(fig_dir, "fig_permutation_null_bycycle.pdf"), p, width = 11, height = 5)
cat("\nSaved: results/permutation_by_cycle.csv and nature_manuscript/figures/fig_permutation_null_bycycle.pdf\n")
