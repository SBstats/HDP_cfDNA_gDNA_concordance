###############################################################################
## R/30_generate_all_manuscript_outputs.R
##
## Single consolidated figure/table generation script for BOTH manuscript
## tracks (methods paper and applied paper).
##
## Figure numbering (applied paper):
##   Fig 1  — Analysis pipeline (manual/Illustrator; not generated here)
##   Fig 2  — Raw gene-level concordance heatmap (2a only; 2b/2c removed)
##   Fig 3  — Per-patient x per-cycle theta concordance heatmap
##             corr_j(theta_k^{cf}(t), theta_k^g(t)), one panel per occ. subclone
##   Fig 4  — Per-cycle permutation negative control
##   Fig 5  — Omega posterior concordance heatmap (secondary estimand)
##   Fig 6  — Clinical outcomes heatmap (best response + longitudinal MRD)
##   Fig 7  — Longitudinal subclonal weight trajectories, all occupied components
##   Sup    — Chain health, traceplots, kappa sensitivity, LOPO forest,
##             theta heatmaps (k_dom, k_2nd), theta scatter grid, theta/timepoint
##
## Usage: Rscript R/30_generate_all_manuscript_outputs.R
##        (run from project root)
###############################################################################

set.seed(20260819)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(patchwork)
  library(ggrepel)
})

base_dir <- getwd()
res_dir  <- file.path(base_dir, "results")
fig_b    <- file.path(base_dir, "figures");                  dir.create(fig_b, showWarnings=FALSE, recursive=TRUE)
tab_b    <- file.path(base_dir, "tables");                   dir.create(tab_b, showWarnings=FALSE, recursive=TRUE)
nm_dir   <- file.path(base_dir, "nature_manuscript");        dir.create(nm_dir, showWarnings=FALSE, recursive=TRUE)
fig_n    <- file.path(nm_dir, "figures");                    dir.create(fig_n, showWarnings=FALSE, recursive=TRUE)
tab_n    <- file.path(nm_dir, "tables_csv");                 dir.create(tab_n, showWarnings=FALSE, recursive=TRUE)

source(file.path(base_dir, "R", "01_lib_core.R"))

clin_file <- file.path(base_dir, "KRd trial", "Clinical data",
                       "KRd 12_1725 data_2_14_2023 from Ben without PHI.xlsx")

cat("=======================================================\n")
cat("Loading results...\n")
md   <- readRDS(file.path(res_dir, "model_data.rds"))
fit1 <- readRDS(file.path(res_dir, "krd_M1_fit.rds"))
fit0 <- readRDS(file.path(res_dir, "krd_M0_fit.rds"))
tt   <- readRDS(file.path(res_dir, "krd_M1_fit_thetatrace.rds"))
cat("Loaded: model_data, krd_M1_fit, krd_M0_fit, krd_M1_fit_thetatrace\n")

## ---- shared helpers ----------------------------------------------------------
pi_info  <- md$paired_info
tvec     <- md$time
pvec     <- md$patient
Yg       <- md$Y_g;  Ycf <- md$Y_cf
n_obs    <- md$N_obs; p <- md$p; n_pat <- md$n; T_max <- md$T_max
K_trace  <- fit1$chains[[1]]$K
n_ch_M1  <- length(fit1$chains)
n_ch_M0  <- length(fit0$chains)

tp_map    <- c("1"="SCR","2"="C4","3"="C8","4"="C18","5"="3Y")
tp_levels <- c("SCR","C4","C8","C18","3Y")

tp_map_long    <- c("1"="Screening","2"="C4","3"="C8","4"="C18","5"="3 YR F/U")
tp_levels_long <- c("Screening","C4","C8","C18","3 YR F/U")

## Clean-chain identification.
## M1: trapped chains have LOWER loglik → keep near MAX.
## M0: trapped chains have HIGHER loglik → keep near MIN.
chain_mean_ll <- function(fit) {
  sapply(seq_along(fit$chains), function(i) mean(fit$chains[[i]]$samples$loglik))
}
m1_ll  <- chain_mean_ll(fit1)
m0_ll  <- chain_mean_ll(fit0)
clean1 <- which(m1_ll >= max(m1_ll) - 10000)
clean0 <- which(m0_ll <= min(m0_ll) + 10000)
cat(sprintf("M1 clean chains: %s  (mean loglik: %s)\n",
            paste(clean1,collapse=","), paste(round(m1_ll[clean1]),collapse=",")))
cat(sprintf("M0 clean chains: %s  (mean loglik: %s)\n",
            paste(clean0,collapse=","), paste(round(m0_ll[clean0]),collapse=",")))

pool_clean <- function(fit, clean_idx) {
  do.call(Map, c(f=rbind, lapply(clean_idx, function(i) fit$chains[[i]]$samples)))
}
m1c <- pool_clean(fit1, clean1)
m0c <- pool_clean(fit0, clean0)

## Pool theta postmeans across all thetatrace chains.
## A single POPULATION signature theta_bar [p x K] now serves BOTH sources
## (writeup Sec 4.1), so the previous theta_g / theta_cf pair collapses to one
## array. theta_bar is the quantity to interpret; the subject-level
## theta [p x K x n] is weakly informed by design (Remark rem:not_claimed).
tt_thbar <- array(0, dim=c(p, K_trace))
tt_th_i  <- array(0, dim=c(p, K_trace, n_pat))
for(i in seq_along(tt$chains)) {
  tt_thbar <- tt_thbar + tt$chains[[i]]$theta_bar_postmean
  tt_th_i  <- tt_th_i  + tt$chains[[i]]$theta_postmean
}
tt_thbar <- tt_thbar / length(tt$chains)
tt_th_i  <- tt_th_i  / length(tt$chains)

## Significant components via posterior-mean pi (averaged over all saved iterations, all chains)
pi_draws_all <- do.call(rbind, lapply(tt$chains, function(ch) {
  v <- ch$samples$v   # [n_save x (K-1)]
  t(apply(v, 1, function(vi) stick_break(c(vi, 1))))
}))
final_pi <- colMeans(pi_draws_all)   # posterior mean pi_k
k_rank   <- order(final_pi, decreasing=TRUE)
k_dom    <- k_rank[1]
k_2nd    <- k_rank[2]
k_occ    <- which(final_pi >= 0.01)   # significant subclones: posterior mean pi >= 1%
K_plus_mode <- as.integer(names(which.max(table(m1c$K_plus))))
cat(sprintf("K+ mode: %d  Dominant k=%d (pi=%.4f)  2nd k=%d (pi=%.4f)\n",
            K_plus_mode, k_dom, final_pi[k_dom], k_2nd, final_pi[k_2nd]))

smry <- function(x) {
  q <- quantile(x, c(0.025, 0.975))
  c(mean=mean(x), median=median(x), lo=unname(q[1]), hi=unname(q[2]))
}

## Posterior-mean omega matrices [K_trace x n_obs]
omega_g_pm  <- matrix(NA, K_trace, n_obs)
omega_cf_pm <- matrix(NA, K_trace, n_obs)
for(k in 1:K_trace) {
  for(obs in 1:n_obs) {
    col_idx <- (k-1)*n_obs + obs
    omega_g_pm[k, obs]  <- mean(m1c$omega_g_trace[, col_idx])
    omega_cf_pm[k, obs] <- mean(m1c$omega_cf_trace[, col_idx])
  }
}

## Posterior-mean non-tumour contamination fraction omega_0.
## NOTE: omega_0 is indexed by OBSERVATION (length N_obs), not by patient --
## the fraction is time-resolved under the present model. The per-patient
## vectors below are AVERAGES over each patient's timepoints, retained for
## back-compatibility with manuscript tables that want one row per patient;
## they discard exactly the temporal variation the respecification captures.
omega0_pm_obs <- colMeans(m1c$omega_0)
omega0_lo_obs <- apply(m1c$omega_0, 2, quantile, 0.025)
omega0_hi_obs <- apply(m1c$omega_0, 2, quantile, 0.975)
wN_pm  <- as.numeric(tapply(omega0_pm_obs, pi_info$patient_id, mean))
wN_lo  <- as.numeric(tapply(omega0_lo_obs, pi_info$patient_id, mean))
wN_hi  <- as.numeric(tapply(omega0_hi_obs, pi_info$patient_id, mean))

## Save posterior-mean omega for downstream figure generation (new Fig 4 trajectory grid).
## obs_meta maps each of the N_obs observations to (Study.ID, Timepoint string, t_num).
## omega_g_pm[k, obs] / omega_cf_pm[k, obs] are the posterior mean weights for
## component k at the obs-th paired observation (rows = components 1..K_trace,
## cols = observations 1..N_obs in the same order as obs_meta).
omega_pm_list <- list(
  omega_g_pm  = omega_g_pm,    # [K_trace x N_obs] gDNA posterior mean subclonal weights
  omega_cf_pm = omega_cf_pm,   # [K_trace x N_obs] cfDNA posterior mean subclonal weights
  ## Per-OBSERVATION contamination (primary; aligned with omega_*_pm columns)
  omega0_pm   = omega0_pm_obs, # [N_obs] posterior mean contamination fraction
  omega0_lo   = omega0_lo_obs, # [N_obs] 2.5th percentile
  omega0_hi   = omega0_hi_obs, # [N_obs] 97.5th percentile
  ## Per-patient AVERAGES (back-compat; discards temporal variation)
  wN_pm       = wN_pm,         # [n_pat] mean over that patient's timepoints
  wN_lo       = wN_lo,         # [n_pat]
  wN_hi       = wN_hi,         # [n_pat]
  K_trace     = K_trace,       # number of components (truncation level)
  N_obs       = n_obs,         # number of paired observations (82)
  k_occ       = k_occ,         # integer vector: significant component indices (pi >= 0.01)
  final_pi    = final_pi,      # length-K_trace vector: posterior mean pi_k for all k
  # Observation metadata — one row per observation, same ordering as omega columns:
  Study.ID    = pi_info$Study.ID,     # patient ID for each observation
  Timepoint   = pi_info$Timepoint,    # full timepoint string: "Screening","C4","C8","C18","3 YR F/U"
  t_num       = pi_info$t_num,        # integer timepoint index 1..5
  patient_id  = pi_info$patient_id    # integer patient index 1..26
)
saveRDS(omega_pm_list, file.path(res_dir, "omega_pm.rds"))
cat(sprintf("  Saved: results/omega_pm.rds  [omega_g_pm: %dx%d, %d significant subclones, omega0_pm: %d obs / %d patients]\n",
            K_trace, n_obs, length(k_occ), length(omega0_pm_obs), length(wN_pm)))

genes <- if(!is.null(md$genes)) md$genes else paste0("Gene_", seq_len(p))


###############################################################################
## SUPPLEMENTARY: CHAIN HEALTH
###############################################################################
cat("\n[S1] Chain health figure...\n")
ll_df <- rbind(
  data.frame(chain=seq_along(fit1$chains), mean_ll=m1_ll, model="M1",
             clean=seq_along(fit1$chains) %in% clean1),
  data.frame(chain=seq_along(fit0$chains), mean_ll=m0_ll, model="M0",
             clean=seq_along(fit0$chains) %in% clean0)
)
ll_df$status <- ifelse(ll_df$clean, "retained","trapped")
ll_df$model  <- factor(ll_df$model, levels=c("M1","M0"))

thresh_df <- ll_df %>% group_by(model) %>%
  summarise(y=ifelse(model[1]=="M1", max(mean_ll)-10000, min(mean_ll)+10000), .groups="drop")

p_chain <- ggplot(ll_df, aes(x=factor(chain), y=mean_ll, fill=status)) +
  geom_col(width=0.6) +
  geom_hline(data=thresh_df, aes(yintercept=y), linetype="dashed", colour="grey40") +
  facet_wrap(~model, scales="free_y") +
  scale_fill_manual(values=c(retained="#4C9F70", trapped="#D1495B")) +
  labs(x="Chain", y="Mean post-burn-in log-likelihood",
       fill=NULL, title="Chain health (dashed = reference +/- 10,000)") +
  theme_bw(base_size=11) + theme(legend.position="bottom")
ggsave(file.path(fig_b, "krd_chain_health.pdf"), p_chain, width=7, height=4)
ggsave(file.path(fig_n, "fig_chain_health.pdf"), p_chain, width=7, height=4)
cat("  Saved krd_chain_health.pdf\n")


###############################################################################
## SUPPLEMENTARY: LOG-LIKELIHOOD TRACES
###############################################################################
cat("\n[S2] Log-likelihood trace figure...\n")
ll_trace_df <- do.call(rbind, lapply(c("M1","M0"), function(m) {
  fit <- if(m=="M1") fit1 else fit0
  cl  <- if(m=="M1") clean1 else clean0
  do.call(rbind, lapply(cl, function(i) {
    s <- fit$chains[[i]]$samples$loglik
    data.frame(iter=seq_along(s), loglik=s, chain=factor(i), model=m)
  }))
}))
ll_trace_df$model <- factor(ll_trace_df$model, levels=c("M1","M0"))

p_trace <- ggplot(ll_trace_df, aes(x=iter, y=loglik, colour=chain)) +
  geom_line(linewidth=0.3) +
  facet_wrap(~model, scales="free_y", ncol=1) +
  labs(x="Iteration (post burn-in, thinned)", y="Log-likelihood", colour="Chain",
       title="Log-likelihood traces (clean chains only)") +
  theme_bw(base_size=11) + theme(legend.position="bottom")
ggsave(file.path(fig_b, "krd_loglik_traces.pdf"), p_trace, width=8, height=6)
ggsave(file.path(fig_n, "fig_loglik_traces.pdf"), p_trace, width=8, height=6)
cat("  Saved krd_loglik_traces.pdf\n")


###############################################################################
## MODEL COMPARISON (WAIC)
###############################################################################
cat("\n[MC] Model comparison WAIC...\n")
lse <- function(x) { m <- max(x); m + log(sum(exp(x-m))) }
waic_clean <- function(ll_mat) {
  lppd  <- sum(apply(ll_mat, 2, function(x) lse(x) - log(nrow(ll_mat))))
  p_w   <- sum(apply(ll_mat, 2, var))
  c(lppd=lppd, p_waic=p_w, waic=-2*(lppd-p_w))
}
w1 <- waic_clean(m1c$pointwise_ll)
w0 <- waic_clean(m0c$pointwise_ll)
delta_waic <- w0["waic"] - w1["waic"]
favours    <- if(delta_waic > 0) "$\\mathcal{M}_1$" else "$\\mathcal{M}_0$"
cat(sprintf("  WAIC M1=%.1f  M0=%.1f  delta=%.1f  Favours: %s\n",
            w1["waic"], w0["waic"], delta_waic, favours))

cmp_tex <- c("\\begin{tabular}{lrrl}","\\toprule",
  "Criterion & $\\mathcal{M}_1$ & $\\mathcal{M}_0$ & Favours \\\\","\\midrule",
  sprintf("WAIC & %s & %s & %s \\\\",
          formatC(w1["waic"],format="f",big.mark=",",digits=0),
          formatC(w0["waic"],format="f",big.mark=",",digits=0), favours),
  sprintf("lppd & %s & %s & -- \\\\",
          formatC(w1["lppd"],format="f",big.mark=",",digits=0),
          formatC(w0["lppd"],format="f",big.mark=",",digits=0)),
  sprintf("$p_{\\mathrm{WAIC}}$ & %s & %s & -- \\\\",
          formatC(w1["p_waic"],format="f",big.mark=",",digits=0),
          formatC(w0["p_waic"],format="f",big.mark=",",digits=0)),
  "\\bottomrule","\\end{tabular}")
writeLines(cmp_tex, file.path(tab_b, "krd_model_comparison.tex"))
write.csv(data.frame(metric=c("WAIC_M1","WAIC_M0","delta_WAIC","lppd_M1","lppd_M0"),
                     value=c(w1["waic"],w0["waic"],delta_waic,w1["lppd"],w0["lppd"])),
          file.path(tab_n, "model_comparison.csv"), row.names=FALSE)
cat("  Saved krd_model_comparison.tex\n")


###############################################################################
## POSTERIOR SUMMARY TABLE
###############################################################################
cat("\n[PT] Posterior summary table...\n")
kp_s  <- smry(m1c$K_plus)
kap_s <- smry(m1c$kappa)
alp_s <- smry(m1c$alpha_dp)
sn_s  <- smry(m1c$sigma_0)
wn_s  <- smry(rowMeans(m1c$omega_0))

post_tex <- c("\\begin{tabular}{lcccc}","\\toprule",
  "Parameter & Mean & Median & 2.5\\% & 97.5\\% \\\\","\\midrule",
  sprintf("$K^+$ (occupied components) & %.1f & %.0f & %.0f & %.0f \\\\",
          kp_s["mean"],kp_s["median"],kp_s["lo"],kp_s["hi"]),
  sprintf("$\\kappa$ (tracking precision) & %.2f & %.2f & %.2f & %.2f \\\\",
          kap_s["mean"],kap_s["median"],kap_s["lo"],kap_s["hi"]),
  sprintf("$\\alpha$ (DP concentration) & %.2f & %.2f & %.2f & %.2f \\\\",
          alp_s["mean"],alp_s["median"],alp_s["lo"],alp_s["hi"]),
  sprintf("$\\sigma_0$ (background SD) & %.3f & %.3f & %.3f & %.3f \\\\",
          sn_s["mean"],sn_s["median"],sn_s["lo"],sn_s["hi"]),
  sprintf("$\\bar{\\omega}^N$ (mean contam.) & %.2e & %.2e & %.2e & %.2e \\\\",
          wn_s["mean"],wn_s["median"],wn_s["lo"],wn_s["hi"]),
  "\\bottomrule","\\end{tabular}")
writeLines(post_tex, file.path(tab_b, "krd_posterior_summary.tex"))
write.csv(data.frame(
  parameter=c("K_plus","kappa","alpha_dp","sigma_0","omega_0_mean"),
  mean=c(kp_s["mean"],kap_s["mean"],alp_s["mean"],sn_s["mean"],wn_s["mean"]),
  median=c(kp_s["median"],kap_s["median"],alp_s["median"],sn_s["median"],wn_s["median"]),
  ci_lo=c(kp_s["lo"],kap_s["lo"],alp_s["lo"],sn_s["lo"],wn_s["lo"]),
  ci_hi=c(kp_s["hi"],kap_s["hi"],alp_s["hi"],sn_s["hi"],wn_s["hi"])),
  file.path(tab_n, "posterior_summary.csv"), row.names=FALSE)
cat("  Saved krd_posterior_summary.tex\n")


###############################################################################
## SUPPLEMENTARY: TRACEPLOTS (clean M1 chains)
###############################################################################
cat("\n[S3] Trace plots...\n")
scalar_trace <- function(fit, param, chains_idx, model_label) {
  do.call(rbind, lapply(chains_idx, function(i) {
    val <- fit$chains[[i]]$samples[[param]]
    if(is.matrix(val)) val <- val[,1]
    data.frame(iter=seq_along(val), value=val, chain=factor(i), param=param, model=model_label)
  }))
}
tp_df <- rbind(
  scalar_trace(fit1, "K_plus",  clean1, "M1"),
  scalar_trace(fit1, "sigma_0", clean1, "M1"),
  scalar_trace(fit1, "loglik",  clean1, "M1")
)
wn1_df <- do.call(rbind, lapply(clean1, function(i) {
  val <- fit1$chains[[i]]$samples$omega_0[,1]
  data.frame(iter=seq_along(val), value=val, chain=factor(i), param="omega_0[obs 1]", model="M1")
}))
tp_df <- rbind(tp_df, wn1_df)
tp_df$param <- factor(tp_df$param,
  levels=c("K_plus","loglik","sigma_0","omega_0[obs 1]"),
  labels=c(expression(K^"+"), "log-lik", expression(sigma[0]), expression(omega[0]["[obs 1]"])))

p_tp <- ggplot(tp_df, aes(x=iter, y=value, colour=chain)) +
  geom_line(linewidth=0.3) +
  facet_wrap(~param, scales="free_y", ncol=2, labeller=label_parsed) +
  labs(x="Iteration (post burn-in)", y=NULL, colour="Chain",
       title=expression("M"[1]*" trace plots (clean chains, separate-"*theta*" model)")) +
  theme_bw(base_size=11) + theme(legend.position="bottom")
ggsave(file.path(fig_b, "krd_traceplots.pdf"), p_tp, width=8, height=6)
ggsave(file.path(fig_n, "fig_traceplots.pdf"),  p_tp, width=8, height=6)
cat("  Saved krd_traceplots.pdf\n")


###############################################################################
## FIGURE 2 — Raw gene-level concordance heatmap (2a only)
## Overlap fix: ggrepel for the 3Y column labels
###############################################################################
cat("\n[Fig 2] Raw gene-level concordance heatmap...\n")

data_dir <- file.path(base_dir, "KRd trial", "5hmC data")
cfDNA_key <- read.csv(file.path(data_dir, "kRd-cfDNA_sample_key.csv"), stringsAsFactors=FALSE)
gDNA_key  <- read.csv(file.path(data_dir, "kRd-gDNA_sample_key.csv"),  stringsAsFactors=FALSE)
cfDNA_key$barcode <- cfDNA_key$Assigned.ID
gDNA_key$barcode  <- gDNA_key$Assigned.ID

paired_raw <- inner_join(
  cfDNA_key %>% select(Study.ID, Timepoint, barcode_cf=barcode),
  gDNA_key  %>% select(Study.ID, Timepoint, barcode_g=barcode),
  by=c("Study.ID","Timepoint")
)
cfDNA_counts <- readRDS(file.path(data_dir, "kRd-cfDNA_genebody_count.RDS"))
gDNA_counts  <- readRDS(file.path(data_dir, "kRd-gDNA_genebody_count.RDS"))
cfDNA_raw    <- as.matrix(cfDNA_counts[, paired_raw$barcode_cf])
gDNA_raw     <- as.matrix(gDNA_counts[,  paired_raw$barcode_g])

raw_corr <- sapply(seq_len(ncol(cfDNA_raw)), function(obs)
  cor(cfDNA_raw[,obs], gDNA_raw[,obs], method="pearson"))

hm_df <- data.frame(
  patient   = factor(paired_raw$Study.ID),
  timepoint = factor(paired_raw$Timepoint, levels=tp_levels_long),
  r         = raw_corr
)

## Identify cells in 3Y column where labels would overlap
## Strategy: use geom_text for all non-3Y; ggrepel only for 3Y
hm_3y    <- hm_df[hm_df$timepoint == "3 YR F/U", ]
hm_other <- hm_df[hm_df$timepoint != "3 YR F/U", ]

p_fig2a <- ggplot(hm_df, aes(x=timepoint, y=patient, fill=r)) +
  geom_tile(colour="white", linewidth=0.5) +
  ## standard labels for all non-3Y cells
  geom_text(data=hm_other,
            aes(x=timepoint, y=patient, label=round(r,2)),
            size=2.4, colour="black") +
  ## ggrepel for 3Y column to resolve overlaps
  geom_text_repel(data=hm_3y,
                  aes(x=timepoint, y=patient, label=round(r,2)),
                  size=2.4, colour="black",
                  box.padding=0.08, point.padding=0,
                  segment.size=0.2, segment.colour="grey50",
                  min.segment.length=0,
                  max.overlaps=Inf,
                  direction="y") +
  scale_fill_gradient2(low="#E63946", mid="#FFFFCC", high="#457B9D",
                       midpoint=median(raw_corr, na.rm=TRUE),
                       name="Pearson r", na.value="grey90") +
  labs(x="Treatment Timepoint", y="Patient",
       title="a  Gene-level cfDNA-gDNA concordance (raw counts)") +
  theme_bw(base_size=10) +
  theme(axis.text.x=element_text(angle=45,hjust=1), panel.grid=element_blank())

ggsave(file.path(fig_b, "krd_raw_concordance_heatmap.pdf"), p_fig2a, width=7, height=8)
ggsave(file.path(fig_n, "fig2a_raw_concordance_heatmap.pdf"), p_fig2a, width=7, height=8)
cat("  Saved fig2a_raw_concordance_heatmap.pdf (2b/2c removed)\n")

write.csv(data.frame(
  metric=c("mean","median","sd","min","max"),
  value=c(mean(raw_corr,na.rm=T), median(raw_corr,na.rm=T),
          sd(raw_corr,na.rm=T), min(raw_corr,na.rm=T), max(raw_corr,na.rm=T))),
  file.path(tab_n, "fig2_raw_concordance_summary.csv"), row.names=FALSE)


###############################################################################
## FIGURE 3 — Combined theta-concordance summary bar + omega posterior heatmap
##
## Because r_k(t) = cor_j(theta_cf_k(t), theta_g_k(t)) is a model-level scalar
## (same value for every patient observed at timepoint t), showing it as a
## patient x timepoint heatmap duplicates information: every row is identical
## within a timepoint.  Instead we show it as a compact summary strip (one row
## per subclone, one column per timepoint), then immediately below it show the
## per-observation omega posterior concordance heatmap (which IS patient-specific).
## The two panels are assembled with patchwork into a single combined figure.
##
## Fig 3  top:  theta signature concordance bar — one value per (k, t)
## Fig 3  bot:  omega posterior concordance heatmap — one value per (patient, t)
##             (separately for each occupied subclone k — see Fig 5 below)
###############################################################################
cat("\n[Fig 3] Theta concordance summary bar...\n")

## Build obs_meta (used by later sections too)
obs_meta <- data.frame(
  obs       = seq_len(n_obs),
  Study.ID  = pi_info$Study.ID,
  t_num     = tvec,
  Timepoint = factor(tp_map[as.character(tvec)], levels=tp_levels)
)

## REPLACED ESTIMAND.
## This figure previously showed r_k = cor(theta_g[,k], theta_cf[,k]), the
## cross-SOURCE signature correlation. A single theta now serves both sources,
## so that quantity is identically 1 by construction and carries no
## information. The meaningful analogue under the present model is the
## cross-SUBJECT agreement: how closely each subject's signature tracks the
## population signature for component k. Low values indicate large
## between-subject heterogeneity (large varsigma_k) and hence a component whose
## label is only weakly shared across the cohort.
r_kt_df <- do.call(rbind, lapply(k_occ, function(k) {
  tb <- tt_thbar[, k]
  r_k <- if (sd(tb, na.rm=TRUE) > 1e-12) {
    mean(vapply(seq_len(n_pat), function(i) {
      ti <- tt_th_i[, k, i]
      if (sd(ti, na.rm=TRUE) > 1e-12) cor(ti, tb) else NA_real_
    }, numeric(1)), na.rm = TRUE)
  } else NA_real_
  data.frame(
    k        = k,
    pi_k     = final_pi[k],
    t_num    = NA_integer_,
    Timepoint = NA_character_,
    r_theta  = r_k,
    k_label  = sprintf("k=%d (pi=%.3f)", k, final_pi[k])
  )
}))
r_kt_df$k_label <- factor(r_kt_df$k_label,
  levels = sprintf("k=%d (pi=%.3f)", k_occ, final_pi[k_occ]))

## Summary bar: one bar per subclone (static theta — no timepoint dimension)
p_fig3_top <- ggplot(r_kt_df, aes(x=r_theta, y=k_label, fill=r_theta)) +
  geom_col(width=0.7) +
  geom_text(aes(label=ifelse(!is.na(r_theta), sprintf("r=%.3f", r_theta), "NA")),
            hjust=-0.1, size=3.5, colour="black") +
  scale_fill_gradient2(low="#E63946", mid="#FFFFCC", high="#457B9D",
                       midpoint=0, limits=c(-1, 1),
                       name=expression(italic(r)*"("*theta[i]*","*bar(theta)*")"),
                       na.value="grey90") +
  scale_x_continuous(limits=c(NA, 1.15), expand=expansion(mult=c(0.05, 0.05))) +
  labs(x=expression(italic(r)(theta[i], bar(theta))), y="Subclone",
       title=expression("Subject-to-population signature agreement  "*italic(r)[k])) +
  theme_bw(base_size=10) +
  theme(panel.grid.minor=element_blank(),
        plot.title=element_text(size=10))

## Broadcast static r_k to patient level for CSV (same value across all timepoints)
theta_r_long <- do.call(rbind, lapply(k_occ, function(k) {
  r_k <- r_kt_df$r_theta[r_kt_df$k == k]  # single scalar for static theta
  data.frame(
    obs       = seq_len(n_obs),
    Study.ID  = obs_meta$Study.ID,
    Timepoint = obs_meta$Timepoint,
    k         = k,
    pi_k      = final_pi[k],
    r_theta   = r_k  # same value for all observations (time-invariant theta)
  )
}))
write.csv(theta_r_long, file.path(tab_n, "theta_concordance_patient_cycle.csv"), row.names=FALSE)
write.csv(r_kt_df,      file.path(tab_n, "theta_concordance_by_component_timepoint.csv"), row.names=FALSE)
cat(sprintf("  Built subject-to-population signature agreement bar (%d subclones)\n",
            length(k_occ)))


###############################################################################
## FIGURE 4 — Omega trajectory grid
## For each patient (row) and each significant subclone (column), plot the
## posterior mean omega_{i,t,k} across treatment cycles for cfDNA (blue) and
## gDNA (red). Grid: 26 rows x K_occ columns, each cell a two-line trajectory.
###############################################################################
cat("\n[Fig 4] Omega trajectory grid...\n")

## Build long-format data frame: (Study.ID, Timepoint, k_label, source, omega_pm)
## Timepoint uses full strings from pi_info$Timepoint — already "Screening","C4",...,"3 YR F/U"
tp_levels_traj <- c("Screening","C4","C8","C18","3 YR F/U")
tp_short       <- c("Screening"="Scr","C4"="C4","C8"="C8","C18"="C18","3 YR F/U"="3Y")

## k_occ sorted by descending final_pi for left-to-right column ordering
k_occ_sorted <- k_occ[order(final_pi[k_occ], decreasing=TRUE)]

traj_rows <- vector("list", length(k_occ_sorted) * n_obs * 2)
row_idx <- 0L
for(ki in seq_along(k_occ_sorted)) {
  k     <- k_occ_sorted[ki]
  k_lbl <- sprintf("Subclone %d (pi=%.3f)", ki, final_pi[k])
  for(obs in seq_len(n_obs)) {
    # Use pi_info$Timepoint directly — full strings, correct ordering
    tp_str <- pi_info$Timepoint[obs]
    for(src in c("gDNA","cfDNA")) {
      row_idx <- row_idx + 1L
      # omega_g_pm[k, obs]: component k is the raw index into the K_trace-row matrix
      omega_val <- if(src == "gDNA") omega_g_pm[k, obs] else omega_cf_pm[k, obs]
      traj_rows[[row_idx]] <- data.frame(
        Study.ID  = pi_info$Study.ID[obs],
        Timepoint = tp_str,
        k_label   = k_lbl,
        ki        = ki,
        source    = src,
        omega_pm  = omega_val,
        stringsAsFactors = FALSE
      )
    }
  }
}
traj_df <- do.call(rbind, traj_rows[seq_len(row_idx)])

## Aggregate duplicates (patients 101-17, 101-73 have 2 entries at 3YR F/U)
traj_df <- traj_df %>%
  group_by(Study.ID, Timepoint, k_label, ki, source) %>%
  summarise(omega_pm = mean(omega_pm, na.rm=TRUE), .groups="drop")

## Ordered factors for clean axis/facet display
traj_df$Timepoint <- factor(traj_df$Timepoint, levels=tp_levels_traj)
traj_df$k_label   <- factor(traj_df$k_label,
                             levels=unique(traj_df$k_label[order(traj_df$ki)]))
traj_df$source    <- factor(traj_df$source, levels=c("cfDNA","gDNA"))
patient_order     <- sort(unique(traj_df$Study.ID))
traj_df$Study.ID  <- factor(traj_df$Study.ID, levels=patient_order)

cat(sprintf("  Trajectory df: %d rows, %d patients, %d subclones, %d sources\n",
            nrow(traj_df), n_distinct(traj_df$Study.ID),
            length(k_occ_sorted), n_distinct(traj_df$source)))

p_fig4_traj <- ggplot(traj_df,
    aes(x=Timepoint, y=omega_pm, color=source, group=source)) +
  geom_line(linewidth=0.7, na.rm=TRUE) +
  geom_point(size=1.2, na.rm=TRUE) +
  scale_color_manual(values=c("cfDNA"="#2196F3","gDNA"="#F44336"),
                     labels=c("cfDNA (peripheral blood)","gDNA (bone marrow)")) +
  scale_x_discrete(labels=tp_short, drop=FALSE) +
  scale_y_continuous(labels=scales::number_format(accuracy=0.01)) +
  facet_grid(Study.ID ~ k_label, scales="free_y", switch="y") +
  labs(x=NULL,
       y=expression(paste("Posterior mean ", omega[italic(itk)])),
       color=NULL,
       title="Fig 4: Longitudinal subclonal weight trajectories") +
  theme_bw(base_size=6.5) +
  theme(strip.text.y.left = element_text(size=5, angle=0, hjust=1),
        strip.text.x      = element_text(size=6.5, face="bold"),
        strip.placement   = "outside",
        axis.text.x       = element_text(angle=45, hjust=1, size=5),
        axis.text.y       = element_text(size=4.5),
        axis.title.y      = element_text(size=6),
        panel.spacing     = unit(0.12, "lines"),
        legend.position   = "bottom",
        legend.text       = element_text(size=6),
        plot.title        = element_text(size=8, face="bold"))

h_traj <- 0.34 * n_distinct(traj_df$Study.ID) + 1.5
w_traj <- 2.8  * length(k_occ_sorted) + 1.2
ggsave(file.path(fig_b, "krd_omega_trajectory_grid.pdf"),  p_fig4_traj, width=w_traj, height=h_traj, limitsize=FALSE)
ggsave(file.path(fig_n, "fig4_omega_trajectory_grid.pdf"), p_fig4_traj, width=w_traj, height=h_traj, limitsize=FALSE)
cat(sprintf("  Saved fig4_omega_trajectory_grid.pdf (%d patients x %d subclones, %.1fx%.1f in)\n",
            n_distinct(traj_df$Study.ID), length(k_occ_sorted), w_traj, h_traj))


###############################################################################
## FIGURE 5 — Per-cycle permutation negative control
## (replaces old global permutation test; was Figure 4 before omega grid added)
###############################################################################
cat("\n[Fig 5] Per-cycle permutation negative control...\n")

## Canonical relabeled/aligned weights (used for omega row)
for (i in seq_along(fit1$chains)) fit1$chains[[i]] <- relabel_by_weight(fit1$chains[[i]])
fit1$chains <- align_chains_by_signature(fit1$chains)
merged <- merge_chain_samples(fit1$chains)
K_merged <- ncol(merged$omega_g_trace) / n_obs
Wg_pm  <- sapply(1:K_merged, function(k)
             colMeans(merged$omega_g_trace[, ((k-1)*n_obs+1):(k*n_obs), drop=FALSE]))
Wcf_pm <- sapply(1:K_merged, function(k)
             colMeans(merged$omega_cf_trace[, ((k-1)*n_obs+1):(k*n_obs), drop=FALSE]))

derange <- function(m) {
  if(m < 2) return(1L)
  repeat { pp <- sample(m); if(all(pp != seq_len(m))) break }
  pp
}

n_perm <- 2000
tp_long_levels <- tp_levels_long

## ---- ROW 1: Gene-level Pearson r (patient-level derangement within cycle) ----
summ_list <- list(); plot_list <- list(); obs_list <- list()

for(ti in seq_along(tp_long_levels)) {
  S <- which(tvec == ti)
  if(length(S) < 2) next
  s <- length(S)
  Cg <- suppressWarnings(cor(Yg[, S, drop=FALSE], Ycf[, S, drop=FALSE]))
  obs_val <- mean(diag(Cg), na.rm=TRUE)
  nulls <- replicate(n_perm, {
    q <- derange(s)
    mean(Cg[cbind(seq_len(s), q)], na.rm=TRUE)
  })
  nm <- mean(nulls, na.rm=TRUE); sdn <- sd(nulls, na.rm=TRUE)
  p_perm <- (1 + sum(nulls >= obs_val, na.rm=TRUE)) / (1 + sum(is.finite(nulls)))
  summ_list[[length(summ_list)+1]] <- data.frame(
    cycle=tp_long_levels[ti], n_pairs=s, metric="Gene-level Pearson r",
    observed=obs_val, null_mean=nm, null_sd=sdn,
    gain=obs_val-nm, p_perm=p_perm, z=(obs_val-nm)/max(sdn,1e-12))
  plot_list[[length(plot_list)+1]] <- data.frame(
    cycle=factor(tp_long_levels[ti], levels=tp_long_levels),
    row_label=factor("Gene-level Pearson r"),
    value=nulls - nm)
  obs_list[[length(obs_list)+1]] <- data.frame(
    cycle=factor(tp_long_levels[ti], levels=tp_long_levels),
    row_label=factor("Gene-level Pearson r"),
    gain=obs_val-nm)
}

## ---- ROWS 2+: Theta posterior-mean concordance — one row per occupied subclone ----
## For each k in k_occ, the permutation shuffles the T-cycle cfDNA–gDNA pairing.
## With T_max=5 timepoints, a derangement of 1:T_max swaps which cycle's cfDNA
## theta is matched to which cycle's gDNA theta. The metric per permutation is the
## mean cor_j(theta_g[j,k,t], theta_cf[j,k,sigma(t)]) averaged over t.
## The observed value at each t is r_k(t); the null is that timepoint matching
## is arbitrary (cfDNA at cycle t could just as well match gDNA at any other cycle).

theta_perm_summ <- list()

## Static theta: r_k is a single value per component (no time variation).
## Permutation null: shuffle gene labels to break the gDNA-cfDNA gene pairing.
for(k in k_occ) {
  tg <- tt_thbar[, k]; tcf <- rowMeans(tt_th_i[, k, , drop=FALSE][, 1, , drop=TRUE])
  r_obs_k <- if(sd(tg,na.rm=TRUE) > 1e-12 && sd(tcf,na.rm=TRUE) > 1e-12) cor(tg, tcf) else NA_real_
  k_label_perm <- sprintf("Theta r — k=%d (pi=%.3f)", k, final_pi[k])
  if(is.na(r_obs_k)) next

  nulls_k <- replicate(n_perm, {
    tcf_null <- sample(tcf)
    if(sd(tg,na.rm=TRUE) > 1e-12 && sd(tcf_null,na.rm=TRUE) > 1e-12) cor(tg, tcf_null) else NA_real_
  })
  nulls_k <- nulls_k[is.finite(nulls_k)]
  nm  <- mean(nulls_k); sdn <- sd(nulls_k)
  p_p <- (1 + sum(nulls_k >= r_obs_k)) / (1 + length(nulls_k))

  ## Add one row per cycle (same observed value — static theta) so the
  ## per-cycle permutation figure still renders one panel per timepoint.
  for(t in 1:T_max) {
    summ_list[[length(summ_list)+1]] <- data.frame(
      cycle=tp_long_levels[t], n_pairs=T_max, metric=k_label_perm,
      observed=r_obs_k, null_mean=nm, null_sd=sdn,
      gain=r_obs_k-nm, p_perm=p_p, z=(r_obs_k-nm)/max(sdn,1e-12))
    plot_list[[length(plot_list)+1]] <- data.frame(
      cycle=factor(tp_long_levels[t], levels=tp_long_levels),
      row_label=factor(k_label_perm),
      value=nulls_k - nm)
    obs_list[[length(obs_list)+1]] <- data.frame(
      cycle=factor(tp_long_levels[t], levels=tp_long_levels),
      row_label=factor(k_label_perm),
      gain=r_obs_k - nm)
  }
  theta_perm_summ[[length(theta_perm_summ)+1]] <- data.frame(
    k=k, pi_k=final_pi[k], t=NA, cycle="pooled",
    r_obs=r_obs_k, null_mean=nm, null_sd=sdn,
    gain=r_obs_k-nm, p_perm=p_p, z=(r_obs_k-nm)/max(sdn,1e-12))
}

perm_summ <- do.call(rbind, summ_list)
perm_plot <- do.call(rbind, plot_list)
perm_obs  <- do.call(rbind, obs_list)

## Row order: gene-level first, then one row per k in k_occ
k_row_labels <- sprintf("Theta r — k=%d (pi=%.3f)", k_occ, final_pi[k_occ])
row_levels <- c("Gene-level Pearson r", k_row_labels)
perm_plot$row_label <- factor(perm_plot$row_label, levels=row_levels)
perm_obs$row_label  <- factor(perm_obs$row_label,  levels=row_levels)

write.csv(perm_summ, file.path(res_dir, "permutation_by_cycle.csv"), row.names=FALSE)
if(length(theta_perm_summ) > 0)
  write.csv(do.call(rbind, theta_perm_summ),
            file.path(res_dir, "theta_permutation_by_cycle.csv"), row.names=FALSE)
cat("\n  Per-cycle permutation summary:\n"); print(perm_summ, row.names=FALSE, digits=3)

n_rows4 <- length(row_levels)   # 1 gene row + n_k_occ theta rows
p_fig4 <- ggplot(perm_plot, aes(value)) +
  geom_histogram(bins=30, fill="grey75", colour="white") +
  geom_vline(xintercept=0, colour="grey35", linetype="dashed") +
  geom_vline(data=perm_obs, aes(xintercept=gain), colour="#E63946", linewidth=0.8) +
  facet_grid(row_label ~ cycle, scales="free_y") +
  labs(x="Concordance gain over null (observed minus null mean)",
       y="Permutations",
       title="Per-cycle permutation negative control: observed gain (red) vs null") +
  theme_bw(base_size=8) +
  theme(strip.text.y=element_text(angle=0, size=7),
        strip.text.x=element_text(size=7))

h4 <- 2.2 * n_rows4 + 0.8
ggsave(file.path(fig_b, "krd_permutation_by_cycle.pdf"),     p_fig4, width=12, height=h4, limitsize=FALSE)
ggsave(file.path(fig_n, "fig5_permutation_null_bycycle.pdf"), p_fig4, width=12, height=h4, limitsize=FALSE)
cat(sprintf("  Saved fig5_permutation_null_bycycle.pdf (%d rows x %d cycles)\n",
            n_rows4, length(tp_long_levels)))


## Save Figure 3 = theta signature concordance bar only (omega discordance removed)
ggsave(file.path(fig_b, "krd_theta_concordance_heatmap.pdf"),  p_fig3_top, width=7, height=3)
ggsave(file.path(fig_n, "fig3_theta_concordance_heatmap.pdf"), p_fig3_top, width=7, height=3)
cat("  Saved fig3_theta_concordance_heatmap.pdf (theta bar only)\n")


###############################################################################
## FIGURE 6 — Clinical outcomes heatmap
## Best overall response (strip) + longitudinal MRD by NGS (main panel)
###############################################################################
cat("\n[Fig 6] Clinical outcomes heatmap...\n")

suppressWarnings({ clin <- read_excel(clin_file, sheet="12_1725 cfDNA Dataset") })

mrd_df_clin <- data.frame(
  Study.ID      = clin$`Study number`,
  best_response = clin$`best overall response`,
  mrd_c8        = clin$`MRD post 8 cycles by NGS (Adaptive) negative =1`,
  mrd_c18       = clin$`MRD by NGS EoT Negative = 1`,
  mrd_1yr       = clin$`MRD by NGS 1yr f/u`,
  mrd_2yr       = clin$`MRD by NGS 2yr f/u`,
  mrd_3yr       = clin$`MRD by NGS 3yr f/u`,
  stringsAsFactors=FALSE
)
mrd_df_clin$mrd_3yrFU <- ifelse(!is.na(mrd_df_clin$mrd_3yr), mrd_df_clin$mrd_3yr,
                         ifelse(!is.na(mrd_df_clin$mrd_2yr), mrd_df_clin$mrd_2yr,
                                mrd_df_clin$mrd_1yr))

code_mrd <- function(x) {
  out <- rep(NA_character_, length(x))
  out[!is.na(x) & x==1] <- "MRD-"
  out[!is.na(x) & x==0] <- "MRD+"
  out
}

mrd_long_clin <- data.frame(
  Study.ID     = rep(mrd_df_clin$Study.ID, times=5),
  timepoint    = factor(rep(tp_levels_long, each=nrow(mrd_df_clin)), levels=tp_levels_long),
  mrd_status   = c(rep(NA_character_, nrow(mrd_df_clin)),
                   rep(NA_character_, nrow(mrd_df_clin)),
                   code_mrd(mrd_df_clin$mrd_c8),
                   code_mrd(mrd_df_clin$mrd_c18),
                   code_mrd(mrd_df_clin$mrd_3yrFU)),
  best_response = rep(mrd_df_clin$best_response, times=5),
  stringsAsFactors=FALSE
)

patients_in_study <- sort(unique(pi_info$Study.ID))
grid_clin <- expand.grid(Study.ID=patients_in_study,
                         timepoint=factor(tp_levels_long, levels=tp_levels_long),
                         stringsAsFactors=FALSE)
grid_clin$timepoint <- factor(as.character(grid_clin$timepoint), levels=tp_levels_long)
heat_clin <- merge(grid_clin, mrd_long_clin, by=c("Study.ID","timepoint"), all.x=TRUE, sort=FALSE)
heat_clin$Study.ID     <- factor(heat_clin$Study.ID, levels=sort(patients_in_study))
heat_clin$cell_label   <- ifelse(is.na(heat_clin$mrd_status), "", heat_clin$mrd_status)
heat_clin$response_val <- ifelse(heat_clin$mrd_status=="MRD-", 1,
                          ifelse(heat_clin$mrd_status=="MRD+", 0, NA))

best_resp_df <- unique(mrd_df_clin[, c("Study.ID","best_response")])
heat_clin <- merge(heat_clin, best_resp_df, by="Study.ID", all.x=TRUE, suffixes=c("","_bo"))

resp_palette <- c("VGPR"="#9F6BA0","sCR"="#2C7A4B")
strip_df2 <- unique(heat_clin[, c("Study.ID","best_response")])
strip_df2$best_response <- factor(strip_df2$best_response, levels=c("VGPR","sCR"))
strip_df2$xcol <- "Best\nresponse"
strip_df2$Study.ID <- factor(strip_df2$Study.ID, levels=sort(patients_in_study))

theme_nat <- theme_bw(base_size=10) +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        panel.grid=element_blank())

p_strip6 <- ggplot(strip_df2, aes(x=xcol, y=Study.ID, fill=best_response)) +
  geom_tile(colour="white", linewidth=0.5) +
  scale_fill_manual(values=resp_palette, name="Best overall\nresponse",
                    na.value="grey90", drop=FALSE) +
  labs(x=NULL, y="Patient") +
  theme_nat +
  theme(plot.margin=margin(5.5,0,5.5,5.5), axis.title.x=element_blank())

p_main6 <- ggplot(heat_clin, aes(x=timepoint, y=Study.ID, fill=response_val)) +
  geom_tile(colour="white", linewidth=0.5) +
  geom_text(aes(label=cell_label), size=2.3, colour="black") +
  scale_fill_gradient2(low="#E63946", mid="#FFFFCC", high="#457B9D",
                       midpoint=0.5, limits=c(0,1),
                       breaks=c(0,1), labels=c("MRD+","MRD-"),
                       name="MRD by NGS",
                       na.value="grey90") +
  labs(x="Treatment Timepoint", y=NULL,
       title="Fig 5: Longitudinal MRD response and best overall response") +
  theme_nat +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        plot.margin=margin(5.5,5.5,5.5,0))

p_fig6 <- p_strip6 + p_main6 + plot_layout(widths=c(0.16, 1))
ggsave(file.path(fig_b, "krd_clinical_outcomes_heatmap.pdf"), p_fig6, width=8, height=8)
ggsave(file.path(fig_n, "fig6_clinical_outcomes_heatmap.pdf"), p_fig6, width=8, height=8)
ggsave(file.path(fig_n, "fig3b_response_heatmap.pdf"), p_fig6, width=8, height=8)
cat("  Saved fig6_clinical_outcomes_heatmap.pdf\n")


###############################################################################
## FIGURE 7 — Static theta posterior interval: gDNA vs cfDNA per subclone
##
## Static theta model: theta is [p x K], no time dimension.
## Shows gene-mean theta posterior distribution (median + 80% CI across MCMC draws)
## for gDNA and cfDNA side-by-side for each occupied subclone.
## theta_g_trace / theta_cf_trace are [Q x p x K] (3D, no T).
###############################################################################
cat("\n[Fig 7] Static theta posterior interval: gDNA vs cfDNA per subclone...\n")

## --- Build per-draw gene-mean theta from thinned theta traces (static, 3D) ---
tt_trace <- readRDS(file.path(res_dir, "krd_M1_fit_thetatrace.rds"))
trace_chains <- tt_trace$chains[clean1]

## For each chain, compute gene-mean per draw per subclone: [Q x K].
## There is one shared signature, so the old gDNA/cfDNA pair is replaced by
## (population signature, cross-subject SD of the subject-level signatures) --
## the pair that is actually informative under this model.
pool_draw_genemean_static <- function(trace) {
  tb <- trace$theta_bar_trace     # [Q x p x K]
  th <- trace$theta_trace         # [Q x p x K x n]
  Q <- dim(tb)[1]; K_tr <- dim(tb)[3]
  gm_g  <- matrix(NA_real_, Q, K_tr)
  gm_cf <- matrix(NA_real_, Q, K_tr)
  for(k in seq_len(K_tr)) {
    gm_g[, k] <- rowMeans(tb[, , k])                   # population gene-mean
    if (!is.null(th)) {
      ## between-subject SD of the per-draw subject gene-means
      sm <- apply(th[, , k, , drop=FALSE], c(1, 4), mean)   # [Q x n]
      gm_cf[, k] <- apply(sm, 1, sd)
    }
  }
  list(g=gm_g, cf=gm_cf)
}

draws_list <- lapply(trace_chains, pool_draw_genemean_static)

## Concatenate along draw dimension: [Q_total x K]
gm_g_all  <- do.call(rbind, lapply(draws_list, `[[`, "g"))
gm_cf_all <- do.call(rbind, lapply(draws_list, `[[`, "cf"))

## Summarise per (k, source): median, 10th, 90th percentile
theta_band_k <- do.call(rbind, lapply(k_occ, function(k) {
  vg  <- gm_g_all[,  k]
  vcf <- gm_cf_all[, k]
  rbind(
    data.frame(k=k, pi_k=final_pi[k], source="population",
               med=median(vg), lo=quantile(vg, 0.10), hi=quantile(vg, 0.90)),
    data.frame(k=k, pi_k=final_pi[k], source="between-subject SD",
               med=median(vcf, na.rm=TRUE),
               lo=quantile(vcf, 0.10, na.rm=TRUE), hi=quantile(vcf, 0.90, na.rm=TRUE))
  )
}))
theta_band_k$k_label <- sprintf("k=%d (pi=%.3f)", theta_band_k$k, theta_band_k$pi_k)
theta_band_k$k_label <- factor(theta_band_k$k_label,
  levels=sprintf("k=%d (pi=%.3f)", k_occ, final_pi[k_occ]))
theta_band_k$source  <- factor(theta_band_k$source, levels=c("population","between-subject SD"))

p_fig7 <- ggplot(theta_band_k, aes(x=source, y=med, colour=source)) +
  geom_point(size=2.5, position=position_dodge(0.4)) +
  geom_errorbar(aes(ymin=lo, ymax=hi), width=0.2,
                position=position_dodge(0.4), linewidth=0.8) +
  facet_wrap(~k_label, scales="free_y") +
  scale_colour_manual(values=c("gDNA"="#1D3557","cfDNA"="#E63946"), name=NULL) +
  labs(x=NULL, y="Gene-mean theta (posterior draw, median +/- 80% CI)",
       title="Fig 7: Subclonal epigenetic signature theta — gDNA vs cfDNA (static across cycles)",
       subtitle="Each point: median gene-mean theta across posterior draws; bar = 10th-90th pctile") +
  theme_bw(base_size=10) +
  theme(legend.position="bottom", panel.grid.minor=element_blank(),
        strip.text=element_text(face="bold"),
        plot.subtitle=element_text(size=8, colour="grey40"))

for(d in c(fig_n, fig_b)) {
  pfx <- if(d == fig_b) "krd_" else ""
  ggsave(file.path(d, sprintf("%sfig7_theta_posterior_interval.pdf", pfx)),
         p_fig7, width=8, height=6, limitsize=FALSE)
}
cat(sprintf("  Saved fig7_theta_posterior_interval.pdf (%d subclones, static theta)\n",
            length(k_occ)))


###############################################################################
## FIGURE 8 — PCA analyses
##
## Panel A: per-cycle observed 5-hMC PCA (cfDNA + gDNA, VST, faceted by cycle).
##   Data: Y_cf and Y_g [p x N_obs].  Combined [p x 2*N_obs] -> prcomp on
##   transposed matrix [2*N_obs x p]. Points = individual observations,
##   coloured by source, connected within pair by grey segment, faceted by cycle.
##
## Panels B+: per occupied subclone k, patient-level reconstruction PCA.
##   For each observation obs (patient i, time t), the p-dimensional RECONSTRUCTED
##   signature for source s is:
##     Yhat_s[j, obs] = omega_s_pm[k, obs] * theta_s[j, k, t]
##   (single-component contribution from subclone k).
##   This gives 84 cfDNA vectors and 84 gDNA vectors of length p = 11,837, one
##   per observation — the same structure as Panel A.  PCA on combined
##   [p x 2*N_obs] matrix, faceted by cycle, paired segments as in Panel A.
###############################################################################
cat("\n[Fig 7] PCA analyses (observed + per-subclone theta)...\n")

## ---- helper: build PCA plot from a combined [p x 2*N_obs] matrix ----
## mat_cf = [p x N_obs], mat_g = [p x N_obs]
make_pca_plot <- function(mat_cf, mat_g, title_str, base_size=9) {
  Y_comb <- cbind(mat_cf, mat_g)          # [p x 2*N_obs]
  pc     <- prcomp(t(Y_comb), center=TRUE, scale.=FALSE)
  ve     <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)

  pts <- data.frame(
    PC1       = pc$x[, 1],
    PC2       = pc$x[, 2],
    Source    = c(rep("cfDNA", n_obs), rep("gDNA", n_obs)),
    Timepoint = factor(rep(tp_map_long[as.character(tvec)], 2), levels=tp_levels_long),
    Study.ID  = c(pi_info$Study.ID, pi_info$Study.ID)
  )
  seg <- data.frame(
    PC1_cf = pc$x[seq_len(n_obs),              1],
    PC2_cf = pc$x[seq_len(n_obs),              2],
    PC1_g  = pc$x[seq_len(n_obs) + n_obs,      1],
    PC2_g  = pc$x[seq_len(n_obs) + n_obs,      2],
    Timepoint = factor(tp_map_long[as.character(tvec)], levels=tp_levels_long)
  )

  ggplot() +
    geom_segment(data=seg,
                 aes(x=PC1_cf, y=PC2_cf, xend=PC1_g, yend=PC2_g),
                 colour="grey70", linewidth=0.3) +
    geom_point(data=pts, aes(PC1, PC2, colour=Source), size=1.3) +
    facet_wrap(~ Timepoint, nrow=1) +
    scale_colour_manual(values=c(cfDNA="#E63946", gDNA="#457B9D")) +
    labs(x=sprintf("PC1 (%.1f%%)", ve[1]),
         y=sprintf("PC2 (%.1f%%)", ve[2]),
         title=title_str) +
    theme_bw(base_size=base_size) +
    theme(legend.position="bottom", panel.grid.minor=element_blank(),
          strip.text=element_text(size=8))
}

## ---- Panel A: observed VST data ----
p_pca_obs <- make_pca_plot(
  Ycf, Yg,
  "Fig 7a  Observed 5-hMC PCA: paired cfDNA (red) and gDNA (blue) per cycle"
)

## ---- Panels B+: per-subclone theta-only PCA ----
## For each subclone k, the data are the T unique signature vectors:
##   theta_g[,  k, t]  (gDNA,  one p-vector per timepoint)
##   theta_cf[, k, t]  (cfDNA, one p-vector per timepoint)
## Combined matrix: [2T x p] — exactly T blue and T red points.
## Segments connect cfDNA and gDNA at the same timepoint.
## Static theta: one signature vector per source per component.
## PCA across genes with gDNA vs cfDNA as the two points — replaced by a simple
## gene-level scatter of theta_g vs theta_cf for component k.
make_theta_pca_plot <- function(k, title_str, base_size=9) {
  ## One theta serves both sources, so a gDNA-vs-cfDNA signature PCA would be
  ## degenerate. Compare the population signature against the per-subject
  ## signatures instead.
  tg  <- tt_thbar[, k]
  tcf <- rowMeans(tt_th_i[, k, , drop=FALSE][, 1, , drop=TRUE])
  r_k <- if(sd(tg)>1e-12 && sd(tcf)>1e-12) cor(tg, tcf) else NA_real_
  df  <- data.frame(theta_g=tg, theta_cf=tcf)
  ggplot(df, aes(theta_g, theta_cf)) +
    geom_point(alpha=0.15, size=0.4, colour="steelblue") +
    geom_abline(slope=1, intercept=0, colour="red", linetype="dashed") +
    annotate("text", x=-Inf, y=Inf, hjust=-0.1, vjust=1.4,
             label=sprintf("r = %.4f", r_k), size=3) +
    labs(x=expression(theta[jk]^{g}), y=expression(theta[jk]^{cf}),
         title=title_str,
         subtitle="One point per gene; static (time-invariant) signatures") +
    theme_bw(base_size=base_size) +
    theme(panel.grid.minor=element_blank())
}

pca_theta_plots <- lapply(k_occ, function(k) {
  make_theta_pca_plot(
    k,
    sprintf("Fig 7  Subclone k=%d (pi=%.3f): theta PCA — 1 point per source per cycle",
            k, final_pi[k])
  )
})

## ---- Save ----
ggsave(file.path(fig_b, "krd_pca_observed_bycycle.pdf"),  p_pca_obs, width=12, height=3.5)
ggsave(file.path(fig_n, "fig8a_pca_observed_bycycle.pdf"), p_pca_obs, width=12, height=3.5)

for(idx in seq_along(k_occ)) {
  k <- k_occ[idx]
  ggsave(file.path(fig_b, sprintf("krd_pca_theta_k%02d.pdf", k)),
         pca_theta_plots[[idx]], width=12, height=3.5)
  ggsave(file.path(fig_n, sprintf("fig8b_pca_theta_k%02d.pdf", k)),
         pca_theta_plots[[idx]], width=12, height=3.5)
}

## Combined multi-page PDF
combined8_path <- file.path(fig_n, "fig8_pca_combined.pdf")
pdf(combined8_path, width=12, height=3.8)
print(p_pca_obs)
for(pl in pca_theta_plots) print(pl)
dev.off()
ggsave(file.path(fig_b, "krd_pca_combined.pdf"), p_pca_obs, width=12, height=3.5)
cat(sprintf("  Saved fig8_pca_combined.pdf (1 observed + %d subclone theta PCA panels)\n",
            length(k_occ)))

## ---- helper for segment start: prcomp already done above; expose for legends ----
## Also save per-subclone as separate files for Fig 8b1, 8b2 etc.
for(idx in seq_along(k_occ)) {
  k <- k_occ[idx]
  ggsave(file.path(fig_n, sprintf("fig8b%d_pca_theta_k%02d.pdf", idx, k)),
         pca_theta_plots[[idx]], width=12, height=3.5)
}

## Back-compat: preserve old filenames referenced in fig legends
if(length(k_occ) >= 1)
  ggsave(file.path(fig_n, sprintf("fig8b_pca_theta_k%02d.pdf", k_occ[1])),
         pca_theta_plots[[1]], width=12, height=3.5)
if(length(k_occ) >= 2)
  ggsave(file.path(fig_n, sprintf("fig8b_pca_theta_k%02d.pdf", k_occ[2])),
         pca_theta_plots[[2]], width=12, height=3.5)



###############################################################################
## SUPPLEMENTARY: THETA HEATMAP (top DE genes, dominant + 2nd component)
###############################################################################
cat("\n[S4] Theta heatmap (top DE genes per component)...\n")
n_top <- 50
for(k_show in unique(c(k_dom, k_2nd))) {
  scores  <- abs(tt_thbar[,k_show])
  top_idx <- order(scores, decreasing=TRUE)[1:n_top]
  top_genes <- genes[top_idx]
  ## Population signature, and the between-subject SD at each gene: the pair
  ## that replaces the old gDNA/cfDNA pair under a shared theta.
  heat_g   <- tt_thbar[top_idx, k_show]
  heat_cf  <- apply(tt_th_i[top_idx, k_show, , drop=FALSE][, 1, , drop=TRUE], 1, sd)
  heat_df <- rbind(
    data.frame(gene=top_genes, theta=heat_g,  source="gDNA"),
    data.frame(gene=top_genes, theta=heat_cf, source="cfDNA")
  )
  heat_df$gene <- factor(heat_df$gene, levels=top_genes)
  clamp <- quantile(abs(heat_df$theta), 0.97)
  heat_df$theta_c <- pmin(pmax(heat_df$theta, -clamp), clamp)
  p_heat <- ggplot(heat_df, aes(x=source, y=gene, fill=theta_c)) +
    geom_tile() +
    scale_fill_gradient2(low="#313695", mid="white", high="#A50026", midpoint=0,
                         name=expression(theta)) +
    labs(x="Source", y="Gene (top 50 by |theta_g| + |theta_cf|)",
         title=sprintf("Subclonal signature heatmap (component k=%d, static theta)", k_show)) +
    theme_bw(base_size=8) +
    theme(axis.text.y=element_text(size=4), panel.grid=element_blank())
  fname <- sprintf("krd_theta_heatmap_k%02d.pdf", k_show)
  ggsave(file.path(fig_b, fname), p_heat, width=8, height=10)
  ggsave(file.path(fig_n, fname), p_heat, width=8, height=10)
}
cat("  Saved krd_theta_heatmap_kXX.pdf\n")


###############################################################################
## SUPPLEMENTARY: THETA SCATTER GRID (population vs subject-mean signature)
###############################################################################
cat("\n[S5] Theta scatter (static, pooled across all genes)...\n")
## Static theta: no timepoint dimension — single scatter per dominant subclone
## Scatter the POPULATION signature against each subject's signature for the
## dominant subclone. (theta_g vs theta_cf is no longer meaningful: one theta
## serves both sources, so such a plot would be the identity line.)
tg_dom  <- tt_thbar[, k_dom]
tcf_dom <- rowMeans(tt_th_i[, k_dom, , drop=FALSE][, 1, , drop=TRUE])
r_pool_dom <- cor(tg_dom, tcf_dom)
qlo <- quantile(c(tg_dom, tcf_dom), 0.005)
qhi <- quantile(c(tg_dom, tcf_dom), 0.995)
scatter_all <- data.frame(tg=tg_dom, tcf=tcf_dom) %>%
  filter(tg > qlo, tg < qhi, tcf > qlo, tcf < qhi)

p_sgrid <- ggplot(scatter_all, aes(x=tg, y=tcf)) +
  geom_point(alpha=0.04, size=0.3, colour="steelblue") +
  geom_abline(slope=1, intercept=0, colour="red", linetype="dashed") +
  annotate("text", x=-Inf, y=Inf, hjust=-0.1, vjust=1.3, size=3.5,
           label=sprintf("r=%.3f", r_pool_dom)) +
  labs(x=expression(bar(theta)[jk]), y=expression(bar(theta)[jk]^{~subject-mean}),
       title=sprintf("Population vs subject-mean signature (k=%d, dominant)", k_dom)) +
  theme_bw(base_size=9)
ggsave(file.path(fig_b, "krd_theta_scatter_grid.pdf"),   p_sgrid, width=5, height=5)
ggsave(file.path(fig_n, "fig1d_theta_scatter_grid.pdf"), p_sgrid, width=5, height=5)
cat("  Saved krd_theta_scatter_grid.pdf (static theta, dominant subclone)\n")

## Save pooled concordance table (one r per subclone, static theta)
r_by_comp <- data.frame(k=k_occ, pi=final_pi[k_occ],
  ## cross-SOURCE r is identically 1 under shared theta; report the
  ## cross-SUBJECT agreement instead (see the Fig 3 note above).
  r_subject_vs_pop=sapply(k_occ, function(k)
    mean(vapply(seq_len(n_pat), function(i)
      if (sd(tt_th_i[,k,i]) > 1e-12 && sd(tt_thbar[,k]) > 1e-12)
        cor(tt_th_i[,k,i], tt_thbar[,k]) else NA_real_, numeric(1)), na.rm=TRUE)))
write.csv(r_by_comp, file.path(tab_n, "theta_concordance_by_component.csv"), row.names=FALSE)


###############################################################################
## SUPPLEMENTARY: KAPPA SENSITIVITY
###############################################################################
cat("\n[S6] Kappa sensitivity...\n")
kappa_vals <- c(0.5, 1.0, 2.0, 5.0)
kappa_res  <- lapply(kappa_vals, function(kv) {
  fname <- sprintf("fit_M1_kappa%.2f.rds", kv)
  fpath <- file.path(res_dir, "kappa", fname)
  if(!file.exists(fpath)) return(NULL)
  fit_k <- readRDS(fpath)
  ll_k  <- sapply(seq_along(fit_k$chains), function(i) mean(fit_k$chains[[i]]$samples$loglik))
  cl_k  <- which(ll_k >= max(ll_k) - 10000)
  mk    <- pool_clean(fit_k, cl_k)
  kp_k  <- smry(mk$K_plus)
  sn_k  <- smry(mk$sigma_0)
  K_k   <- fit_k$chains[[1]]$K
  pc_k  <- sapply(1:n_obs, function(obs) {
    cols <- ((seq_len(K_k)-1)*n_obs)+obs
    if(ncol(mk$omega_g_trace) < max(cols)) return(NA_real_)
    vg <- as.numeric(mk$omega_g_trace[,cols])
    vc <- as.numeric(mk$omega_cf_trace[,cols])
    if(sd(vg)>1e-12 && sd(vc)>1e-12) cor(vg,vc) else NA_real_
  })
  list(kappa=kv, K_plus_mean=kp_k["mean"], K_plus_lo=kp_k["lo"], K_plus_hi=kp_k["hi"],
       sigma_0_mean=sn_k["mean"], sigma_0_lo=sn_k["lo"], sigma_0_hi=sn_k["hi"],
       omega_post_corr_mean=mean(pc_k,na.rm=TRUE))
})
kappa_df <- do.call(rbind, lapply(Filter(Negate(is.null), kappa_res), as.data.frame))
if(nrow(kappa_df) > 0) {
  p_kap_r <- ggplot(kappa_df, aes(x=kappa, y=omega_post_corr_mean)) +
    geom_line(colour="#1D3557") + geom_point(size=3, colour="#1D3557") +
    geom_text(aes(label=sprintf("%.3f",omega_post_corr_mean)), vjust=-0.7, size=3) +
    labs(x=expression(kappa), y=expression("Mean posterior "*cor*"("*omega^{g}*","*omega^{cf}*")"),
         title=expression("Sensitivity to "*kappa)) +
    theme_bw(base_size=11)
  ggsave(file.path(fig_n, "figS1_kappa_sensitivity.pdf"), p_kap_r, width=5, height=4)
  p_kap_k <- ggplot(kappa_df, aes(x=kappa, y=K_plus_mean, ymin=K_plus_lo, ymax=K_plus_hi)) +
    geom_line(colour="#2A9D8F") + geom_point(size=3, colour="#2A9D8F") +
    geom_errorbar(width=0.1, colour="#2A9D8F") +
    labs(x=expression(kappa), y=expression("Mean "*K^"+"),
         title=expression("Sensitivity of "*K^"+"*" to "*kappa)) +
    theme_bw(base_size=11)
  ggsave(file.path(fig_n, "figS1b_kappa_sensitivity_kplus.pdf"), p_kap_k, width=5, height=4)
  write.csv(kappa_df, file.path(tab_n, "kappa_sensitivity.csv"), row.names=FALSE)
  cat("  Saved kappa sensitivity figures.\n")
} else cat("  No kappa sensitivity fits found — skipped.\n")


###############################################################################
## SUPPLEMENTARY: LOPO FOREST PLOT
###############################################################################
cat("\n[S7] LOPO cross-validation forest plot...\n")
lopo_cv_path <- file.path(res_dir, "lopo_cv.csv")
lopo_pp_path <- file.path(res_dir, "lopo_per_patient.csv")
if(!file.exists(lopo_cv_path)) {
  cat("  WARNING: Missing results/lopo_cv.csv — run slurm/run_lopo_aggregate.R first. Skipping.\n")
} else {
  cv_summary <- read.csv(lopo_cv_path)
  lopo_df    <- read.csv(lopo_pp_path)
  flagged_str <- as.character(cv_summary$folds_flagged)
  n_flagged   <- if (is.na(flagged_str) || nchar(trimws(flagged_str)) == 0) 0L else
                   length(strsplit(flagged_str, ";")[[1]])
  n_folds    <- cv_summary$folds_used + n_flagged
  cat(sprintf("  LOPO: %d folds total, %d used  delta_elpd=%.1f (SE %.1f)  %d/%d obs favour M1\n",
              n_folds, cv_summary$folds_used, cv_summary$delta_elpd, cv_summary$se_delta,
              cv_summary$n_obs_favouring_M1, cv_summary$n_obs))
  lopo_df$pat_label <- factor(paste0("Fold ", lopo_df$patient),
                              levels=paste0("Fold ", lopo_df$patient[order(lopo_df$d)]))
  p_lopo <- ggplot(lopo_df, aes(x=d, y=pat_label, colour=d>0)) +
    geom_vline(xintercept=0, colour="grey60", linetype="dashed") +
    geom_point(size=2.5) +
    scale_colour_manual(values=c(`TRUE`="#1D3557",`FALSE`="#E63946"),
                        labels=c(`TRUE`="favours M1",`FALSE`="favours M0"), name=NULL) +
    annotate("text", x=Inf, y=Inf, hjust=1.1, vjust=1.3, size=3.5,
             label=sprintf("total delta elpd = %.1f (SE %.1f)\n%d/%d obs favour M1",
                           cv_summary$delta_elpd, cv_summary$se_delta,
                           cv_summary$n_obs_favouring_M1, cv_summary$n_obs)) +
    labs(x=expression("per-patient held-out "*Delta*"elpd ("*M[1]-M[0]*")"), y="Patient fold",
         title="Leave-one-patient-out predictive comparison") +
    theme_bw(base_size=11) + theme(legend.position="bottom")
  ggsave(file.path(fig_b, "krd_lopo_forest.pdf"),  p_lopo, width=6, height=5)
  ggsave(file.path(fig_n, "figS2_lopo_forest.pdf"), p_lopo, width=6, height=5)
  lopo_tex <- c("\\begin{tabular}{lr}","\\toprule",
    "Summary & Value \\\\","\\midrule",
    sprintf("Folds used (not flagged) & %d of %d \\\\", cv_summary$folds_used, n_folds),
    sprintf("Total $\\Delta$elpd ($M_1 - M_0$) & %.1f \\\\", cv_summary$delta_elpd),
    sprintf("SE($\\Delta$elpd) & %.1f \\\\", cv_summary$se_delta),
    sprintf("$\\Delta$elpd / SE & %.2f \\\\", cv_summary$delta_elpd / cv_summary$se_delta),
    sprintf("Observations favouring $\\mathcal{M}_1$ & %d / %d \\\\",
            cv_summary$n_obs_favouring_M1, cv_summary$n_obs),
    sprintf("$\\Delta$elpd (cfDNA, pairing-specific) & %.1f (SE %.1f) \\\\",
            cv_summary$delta_elpd_cf, cv_summary$se_delta_cf),
    "\\bottomrule","\\end{tabular}")
  writeLines(lopo_tex, file.path(tab_b, "krd_lopo_comparison.tex"))
  write.csv(data.frame(
    folds_used=cv_summary$folds_used, folds_total=n_folds,
    delta_elpd=cv_summary$delta_elpd, se_delta=cv_summary$se_delta,
    ratio=cv_summary$delta_elpd/cv_summary$se_delta,
    n_obs_favouring_M1=cv_summary$n_obs_favouring_M1, n_obs=cv_summary$n_obs,
    delta_elpd_cf=cv_summary$delta_elpd_cf, se_delta_cf=cv_summary$se_delta_cf),
    file.path(tab_n, "lopo_comparison.csv"), row.names=FALSE)
  cat("  Saved figS2_lopo_forest.pdf\n")
}


###############################################################################
## SUPPLEMENTARY: PER-PATIENT CONTAMINATION TABLE
###############################################################################
cat("\n[S8] Per-patient contamination table...\n")
## Per-patient averages of the per-observation omega_0 (see the note at the
## omega_pm.rds block above).
wN_pm  <- as.numeric(tapply(colMeans(m1c$omega_0), pi_info$patient_id, mean))
wN_lo  <- as.numeric(tapply(apply(m1c$omega_0, 2, quantile, 0.025), pi_info$patient_id, mean))
wN_hi  <- as.numeric(tapply(apply(m1c$omega_0, 2, quantile, 0.975), pi_info$patient_id, mean))
pat_df <- data.frame(
  patient_id = seq_len(n_pat),
  Study.ID   = unique(pi_info$Study.ID[order(pi_info$patient_id)])[seq_len(n_pat)],
  wN_mean=wN_pm, wN_lo=wN_lo, wN_hi=wN_hi
)
wN_tex <- c("\\begin{tabular}{lrrr}","\\toprule",
  "Patient & Posterior mean $\\omega_{0}$ & 2.5\\% & 97.5\\% \\\\","\\midrule",
  sprintf("%s & %.2e & %.2e & %.2e \\\\", pat_df$Study.ID, pat_df$wN_mean, pat_df$wN_lo, pat_df$wN_hi),
  "\\bottomrule","\\end{tabular}")
writeLines(wN_tex, file.path(tab_b, "krd_omega0_by_patient.tex"))
write.csv(pat_df, file.path(tab_n, "omega0_by_patient.csv"), row.names=FALSE)
cat("  Saved krd_omega0_by_patient.tex\n")


###############################################################################
## SUPPLEMENTARY: DATA SNAPSHOT TABLE
###############################################################################
cat("\n[S9] Data snapshot table...\n")
data_snap <- pi_info %>%
  count(Timepoint) %>%
  mutate(Timepoint=factor(Timepoint, levels=c("Screening","C4","C8","C18","3 YR F/U")))
write.csv(data_snap, file.path(tab_n, "data_snapshot_by_timepoint.csv"), row.names=FALSE)
pat_tp <- pi_info %>% group_by(patient_id) %>%
  summarise(n_timepoints=n(), timepoints=paste(sort(unique(Timepoint)),collapse=", "), .groups="drop")
write.csv(pat_tp, file.path(tab_n, "data_snapshot_per_patient.csv"), row.names=FALSE)


###############################################################################
## SUPPLEMENTARY: MODEL-BASED PREDICTION CONCORDANCE HEATMAP
###############################################################################
cat("\n[S10] Model-based prediction concordance...\n")
pred_corr <- numeric(n_obs)
for(obs in seq_len(n_obs)) {
  ## Static theta: [p x K], no time dimension
  ## Both sources share theta, so the reconstructed profiles differ ONLY
  ## through omega -- which is precisely the concordance estimand. Use
  ## subject i's own signature slice.
  i_s <- pi_info$patient_id[obs]
  yhat_g  <- as.vector(tt_th_i[, , i_s] %*% omega_g_pm[,obs])
  yhat_cf <- as.vector(tt_th_i[, , i_s] %*% omega_cf_pm[,obs])
  if(sd(yhat_g)>1e-10 && sd(yhat_cf)>1e-10) {
    pred_corr[obs] <- cor(yhat_g, yhat_cf)
  } else pred_corr[obs] <- NA_real_
}
cat(sprintf("  Prediction concordance: mean=%.4f  median=%.4f\n",
            mean(pred_corr,na.rm=T), median(pred_corr,na.rm=T)))
pred_df <- data.frame(
  Study.ID  = pi_info$Study.ID,
  Timepoint = factor(tp_map[as.character(tvec)], levels=tp_levels),
  pred_corr = pred_corr
) %>% group_by(Study.ID, Timepoint) %>%
  summarise(pred_corr=mean(pred_corr,na.rm=TRUE), .groups="drop")
p_pred <- ggplot(pred_df, aes(x=Timepoint, y=Study.ID, fill=pred_corr)) +
  geom_tile(colour="white", linewidth=0.4) +
  geom_text(aes(label=sprintf("%.2f", pred_corr)), size=2.2, colour="black") +
  scale_fill_gradient2(low="#FEE090", mid="#4575B4", high="#313695",
                       midpoint=median(pred_corr,na.rm=T), name="Model r", na.value="grey90") +
  labs(x="Timepoint", y="Patient", title="Model-based prediction concordance") +
  theme_bw(base_size=9) +
  theme(axis.text.x=element_text(angle=45,hjust=1), panel.grid=element_blank())
ggsave(file.path(fig_b, "krd_prediction_concordance_heatmap.pdf"), p_pred, width=7, height=8)
ggsave(file.path(fig_n, "fig3b_prediction_concordance_heatmap.pdf"), p_pred, width=7, height=8)
write.csv(pred_df, file.path(tab_n, "prediction_concordance.csv"), row.names=FALSE)
cat("  Saved krd_prediction_concordance_heatmap.pdf\n")


###############################################################################
## FINAL SUMMARY
###############################################################################
cat("\n=======================================================\n")
cat("All outputs written:\n")
cat("  figures/      :", paste(sort(list.files(fig_b,pattern="\\.pdf$")),collapse="\n               "),"\n")
cat("  tables/       :", paste(sort(list.files(tab_b,pattern="\\.tex$")),collapse="\n               "),"\n")
cat("  nature/figs/  :", paste(sort(list.files(fig_n,pattern="\\.pdf$")),collapse="\n               "),"\n")
cat("  nature/tables/:", paste(sort(list.files(tab_n,pattern="\\.csv$")),collapse="\n               "),"\n")
cat("=======================================================\n")
kap_s2 <- smry(m1c$kappa)
cat(sprintf("\nKey results:\n"))
cat(sprintf("  K+ mode = %d\n", K_plus_mode))
cat(sprintf("  kappa (fixed) = %.1f  sigma_0 = %.3f (%.3f, %.3f)\n",
            kap_s2["mean"], sn_s["mean"], sn_s["lo"], sn_s["hi"]))
cat(sprintf("  Dominant k=%d  pi=%.4f\n", k_dom, final_pi[k_dom]))
cat(sprintf("  Theta concordance (k=%d, pooled j,t): r=%.4f\n", k_dom, r_pool_dom))
cat(sprintf("  Significant subclones (pi >= 0.01): k = %s\n",
            paste(k_occ, collapse=", ")))
cat(sprintf("  WAIC delta (M0-M1) = %.1f (%s)\n", delta_waic,
            if(delta_waic>0)"favours M1" else "favours M0"))
