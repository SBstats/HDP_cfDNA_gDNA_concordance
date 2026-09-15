###############################################################################
## R/12_subclonal_weights_vs_mrd.R
##
## GOAL (reframed predictor): the pre-specified concordance<->MRD association
## (R/11) was null. Here we test whether the clinical OUTCOME (MRD by NGS)
## tracks the SUBCLONAL MIXTURE WEIGHTS themselves and their LONGITUDINAL
## DYNAMICS, rather than cfDNA-gDNA concordance.
##
## Predictors (per patient-visit), from posterior-mean source-specific weights
## under M1 (canonical relabeled/aligned components; K+ = 2 occupied):
##   w1g, w1cf        weight on dominant subclone 1 (gDNA / cfDNA)
##   dom_g, dom_cf    subclonal dominance  = max_k w_k  (0.5 balanced .. 1 pure)
##   ent_g, ent_cf    subclonal entropy over occupied components (diversity)
## Longitudinal / patient-level:
##   d_w1g            w1g(t) - w1g(screening)          (shift from baseline)
##   sd_w1g           within-patient SD of w1g over t  (subclonal instability)
##   switch           dominant subclone changes over t (0/1)
##
## Outcome: MRD-positive (1) vs negative (0); windows C8/C18/3Y. Same coding
## and inference as R/11 (LMM + patient-level cluster permutation; logistic
## GLMM for OR; patient-level Wilcoxon; prognostic screening -> later MRD).
##
## gDNA (bone marrow) is the clinical ground truth, so gDNA weights are primary.
##
## Outputs (results/):
##   subclonal_weights_dataset.csv          per-visit weights + MRD
##   subclonal_weights_mrd_association.csv   tidy tests
## Outputs (nature_manuscript/figures/):
##   fig_weights_by_mrd.pdf
###############################################################################

set.seed(20260718)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readxl); library(lme4); library(ggplot2)
})

base_dir <- getwd()
res_dir  <- file.path(base_dir, "results")
fig_dir  <- file.path(base_dir, "nature_manuscript", "figures")
clin_file <- file.path(base_dir, "KRd trial", "Clinical data",
                       "KRd 12_1725 data_2_14_2023 from Ben without PHI.xlsx")
source(file.path(base_dir, "R", "01_lib_core.R"))   # relabel/align/merge helpers

## ---------------------------------------------------------------------------
## 1. Canonical posterior-mean weights per observation (match figure pipeline)
## ---------------------------------------------------------------------------
md  <- readRDS(file.path(res_dir, "model_data.rds"))
n_obs <- md$N_obs
krd <- readRDS(file.path(res_dir, "krd_M1_fit.rds"))
m1  <- krd
for (i in seq_along(m1$chains)) m1$chains[[i]] <- relabel_by_weight(m1$chains[[i]])
m1$chains <- align_chains_by_signature(m1$chains)
merged <- merge_chain_samples(m1$chains)

K_trace <- ncol(merged$omega_g_trace) / n_obs
occ <- 1:2                                  # occupied components (K+ = 2)
col_of <- function(k, obs) (k - 1) * n_obs + obs

wmean <- function(trace) {
  W <- matrix(NA_real_, n_obs, K_trace)
  for (k in 1:K_trace) for (o in 1:n_obs) W[o, k] <- mean(trace[, col_of(k, o)])
  W
}
Wg  <- wmean(merged$omega_g_trace)
Wcf <- wmean(merged$omega_cf_trace)

# renormalize over occupied components (tail comps are ~0 but keep exact)
ent <- function(w) { w <- w[w > 1e-12]; w <- w / sum(w); -sum(w * log(w)) }
obs_df <- data.frame(
  patient_id = md$patient, t_num = md$time,
  w1g = Wg[, 1], w2g = Wg[, 2], w1cf = Wcf[, 1], w2cf = Wcf[, 2]
) %>%
  mutate(
    dom_g  = pmax(w1g, w2g) / (w1g + w2g),
    dom_cf = pmax(w1cf, w2cf) / (w1cf + w2cf),
    ent_g  = apply(cbind(w1g, w2g), 1, ent),
    ent_cf = apply(cbind(w1cf, w2cf), 1, ent),
    dom_sub_g = ifelse(w1g >= w2g, 1L, 2L)
  )

cat("w1g  range:", paste(round(range(obs_df$w1g), 3), collapse = "-"),
    " mean:", round(mean(obs_df$w1g), 3), "\n")

# average replicate patient x t_num rows
obs_df <- obs_df %>% group_by(patient_id, t_num) %>%
  summarise(across(c(w1g, w1cf, dom_g, dom_cf, ent_g, ent_cf),
                   ~ mean(.x)),
            dom_sub_g = dom_sub_g[which.max(abs(w1g - 0.5))], .groups = "drop")

## id map + longitudinal features
id_map <- unique(md$paired_info[, c("patient_id", "Study.ID")])
obs_df <- left_join(obs_df, id_map, by = "patient_id")

scr <- obs_df %>% filter(t_num == 1) %>% select(patient_id, w1g_scr = w1g)
obs_df <- obs_df %>% left_join(scr, by = "patient_id") %>%
  mutate(d_w1g = w1g - w1g_scr)

pat_dyn <- obs_df %>% group_by(patient_id, Study.ID) %>%
  summarise(sd_w1g = ifelse(n() > 1, sd(w1g), NA_real_),
            switch  = as.integer(n_distinct(dom_sub_g) > 1),
            mean_w1g = mean(w1g), mean_dom_g = mean(dom_g),
            mean_ent_g = mean(ent_g), .groups = "drop")

## ---------------------------------------------------------------------------
## 2. MRD outcome (validated Fig-3b mapping; MRD_pos = 1 - "negative=1")
## ---------------------------------------------------------------------------
suppressWarnings({ clin <- read_excel(clin_file, sheet = "12_1725 cfDNA Dataset") })
mrd_df <- data.frame(
  Study.ID = clin$`Study number`, best_response = clin$`best overall response`,
  neg_c8 = clin$`MRD post 8 cycles by NGS (Adaptive) negative =1`,
  neg_c18 = clin$`MRD by NGS EoT Negative = 1`,
  neg_1yr = clin$`MRD by NGS 1yr f/u`, neg_2yr = clin$`MRD by NGS 2yr f/u`,
  neg_3yr = clin$`MRD by NGS 3yr f/u`, stringsAsFactors = FALSE)
mrd_df$neg_3yrFU <- ifelse(!is.na(mrd_df$neg_3yr), mrd_df$neg_3yr,
                    ifelse(!is.na(mrd_df$neg_2yr), mrd_df$neg_2yr, mrd_df$neg_1yr))
to_pos <- function(neg) ifelse(is.na(neg), NA_integer_, as.integer(1 - neg))
fin <- function(x) ifelse(is.finite(x), x, NA_real_)
mrd_long <- bind_rows(
  data.frame(Study.ID = mrd_df$Study.ID, t_num = 3L, MRD_pos = to_pos(mrd_df$neg_c8)),
  data.frame(Study.ID = mrd_df$Study.ID, t_num = 4L, MRD_pos = to_pos(mrd_df$neg_c18)),
  data.frame(Study.ID = mrd_df$Study.ID, t_num = 5L, MRD_pos = to_pos(mrd_df$neg_3yrFU)))
mrd_pat <- data.frame(
  Study.ID = mrd_df$Study.ID, best_response = mrd_df$best_response,
  MRD_ever_pos = to_pos(fin(pmin(mrd_df$neg_c8, mrd_df$neg_c18, mrd_df$neg_3yrFU, na.rm = TRUE))),
  MRD_persistent = to_pos(fin(pmax(mrd_df$neg_c8, mrd_df$neg_c18, mrd_df$neg_3yrFU, na.rm = TRUE))))

## ---------------------------------------------------------------------------
## 3. Contemporaneous dataset
## ---------------------------------------------------------------------------
dat <- obs_df %>%
  inner_join(mrd_long, by = c("Study.ID", "t_num")) %>%
  filter(!is.na(MRD_pos)) %>%
  mutate(patient = factor(patient_id),
         MRD = factor(ifelse(MRD_pos == 1, "MRD+", "MRD-"), levels = c("MRD-", "MRD+")))
write.csv(dat, file.path(res_dir, "subclonal_weights_dataset.csv"), row.names = FALSE)
cat("\nContemporaneous visits:", nrow(dat), "| patients:", n_distinct(dat$Study.ID),
    "| MRD+:", sum(dat$MRD_pos), " MRD-:", sum(dat$MRD_pos == 0), "\n")

## ---------------------------------------------------------------------------
## 4. Inference engine (LMM + patient-level cluster permutation)
## ---------------------------------------------------------------------------
n_perm <- 2000
perm_patient_labels <- function(df, seed_i) {
  set.seed(seed_i)
  pats <- unique(df$Study.ID); donor <- sample(pats)
  key <- df %>% distinct(Study.ID, t_num, MRD_pos)
  lut <- setNames(donor, pats)
  out <- df %>% mutate(donorID = lut[as.character(Study.ID)]) %>% select(-MRD_pos) %>%
    left_join(key %>% rename(donorID = Study.ID), by = c("donorID", "t_num"))
  out$MRD_pos
}
lmm_effect <- function(df, yvar) {
  d <- df[!is.na(df[[yvar]]) & !is.na(df$MRD_pos), ]; d$y <- d[[yvar]]
  if (n_distinct(d$MRD_pos) < 2 || nrow(d) < 6) return(NULL)
  ctl <- lmerControl(check.conv.singular = "ignore")
  fit  <- suppressMessages(lmer(y ~ MRD_pos + (1|patient), d, REML = FALSE, control = ctl))
  fit0 <- suppressMessages(lmer(y ~ 1 + (1|patient), d, REML = FALSE, control = ctl))
  beta <- fixef(fit)[["MRD_pos"]]; se <- sqrt(vcov(fit)["MRD_pos","MRD_pos"])
  p_lrt <- anova(fit0, fit)$`Pr(>Chisq)`[2]
  obs_stat <- abs(beta); ps <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    d$MRDp <- perm_patient_labels(d, i); dd <- d[!is.na(d$MRDp), ]
    ps[i] <- if (n_distinct(dd$MRDp) < 2) NA else
      tryCatch(abs(fixef(suppressMessages(lmer(y ~ MRDp + (1|patient), dd, REML = FALSE,
        control = ctl)))[["MRDp"]]), error = function(e) NA)
  }
  p_perm <- (1 + sum(ps >= obs_stat, na.rm = TRUE)) / (1 + sum(!is.na(ps)))
  data.frame(metric = yvar, n = nrow(d),
             mean_MRDneg = mean(d$y[d$MRD_pos == 0]), mean_MRDpos = mean(d$y[d$MRD_pos == 1]),
             beta = beta, ci_lo = beta - 1.96*se, ci_hi = beta + 1.96*se,
             p_LRT = p_lrt, p_perm = p_perm)
}
glmm_or <- function(df, yvar) {
  d <- df[!is.na(df[[yvar]]), ]; d$y <- d[[yvar]]
  fit <- tryCatch(suppressMessages(glmer(MRD_pos ~ scale(y) + (1|patient), d,
           family = binomial, control = glmerControl(optimizer = "bobyqa"))),
           error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  co <- summary(fit)$coefficients["scale(y)", ]
  data.frame(metric = yvar, OR_per_SD = exp(co[1]),
             OR_lo = exp(co[1] - 1.96*co[2]), OR_hi = exp(co[1] + 1.96*co[2]), p_wald = co[4])
}

vis_vars <- c("w1g","dom_g","ent_g","w1cf","dom_cf","d_w1g")
lab <- c(w1g="Dominant-subclone wt (gDNA)", dom_g="Subclonal dominance (gDNA)",
         ent_g="Subclonal entropy (gDNA)", w1cf="Dominant-subclone wt (cfDNA)",
         dom_cf="Subclonal dominance (cfDNA)", d_w1g="Shift in w1g from baseline")
cat("\n==== Contemporaneous LMM: weight-feature ~ MRD + (1|patient) ====\n")
lmm_res <- bind_rows(lapply(vis_vars, function(v) lmm_effect(dat, v)))
lmm_res$metric_label <- lab[lmm_res$metric]
print(lmm_res, digits = 3)
cat("\n==== Logistic GLMM: MRD ~ weight-feature(SD) + (1|patient) ====\n")
glmm_res <- bind_rows(lapply(vis_vars, function(v) glmm_or(dat, v)))
if (nrow(glmm_res)) { glmm_res$metric_label <- lab[glmm_res$metric]; print(glmm_res, digits = 3) }

## ---------------------------------------------------------------------------
## 5. Patient-level: static composition AND dynamics vs MRD / response
## ---------------------------------------------------------------------------
pat <- pat_dyn %>% left_join(mrd_pat, by = "Study.ID")
wilcox_by <- function(df, yvar, grp) {
  d <- df[!is.na(df[[yvar]]) & !is.na(df[[grp]]), ]; g <- factor(d[[grp]])
  if (nlevels(g) != 2) return(NULL)
  w <- suppressWarnings(wilcox.test(d[[yvar]] ~ g)); lv <- levels(g)
  data.frame(metric = yvar, group = grp,
             grp0 = lv[1], median0 = median(d[[yvar]][g==lv[1]]), n0 = sum(g==lv[1]),
             grp1 = lv[2], median1 = median(d[[yvar]][g==lv[2]]), n1 = sum(g==lv[2]),
             p_wilcox = w$p.value)
}
pat_vars <- c("mean_w1g","mean_dom_g","mean_ent_g","sd_w1g","switch")
cat("\n==== Patient-level Wilcoxon (composition + DYNAMICS) vs MRD / response ====\n")
pat_res <- bind_rows(
  lapply(pat_vars, wilcox_by, df = pat, grp = "MRD_persistent"),
  lapply(pat_vars, wilcox_by, df = pat, grp = "MRD_ever_pos"),
  lapply(pat_vars, function(v) wilcox_by(
     pat %>% mutate(resp2 = ifelse(best_response=="sCR","sCR","VGPR")), v, "resp2")))
print(pat_res, digits = 3)

## prognostic: screening composition -> later MRD
prog <- obs_df %>% filter(t_num == 1) %>% select(Study.ID, w1g, dom_g, ent_g) %>%
  left_join(mrd_pat, by = "Study.ID")
cat("\n==== Prognostic: screening composition vs MRD_ever_pos ====\n")
prog_res <- bind_rows(lapply(c("w1g","dom_g","ent_g"),
                             function(v) wilcox_by(prog, v, "MRD_ever_pos")))
if (nrow(prog_res)) print(prog_res, digits = 3)

## ---------------------------------------------------------------------------
## 6. Save tidy results + figure
## ---------------------------------------------------------------------------
out <- bind_rows(
  lmm_res %>% transmute(analysis = "contemporaneous_LMM", metric = metric_label, n,
    est = beta, ci_lo, ci_hi, p_primary = p_perm, p_secondary = p_LRT,
    detail = sprintf("mean MRD- %.3f vs MRD+ %.3f", mean_MRDneg, mean_MRDpos)),
  pat_res %>% transmute(analysis = paste0("patient_", group),
    metric = metric, n = n0 + n1, est = median1 - median0, ci_lo = NA, ci_hi = NA,
    p_primary = p_wilcox, p_secondary = NA,
    detail = sprintf("%s med %.3f (n=%d) vs %s med %.3f (n=%d)", grp0, median0, n0, grp1, median1, n1)))
write.csv(out, file.path(res_dir, "subclonal_weights_mrd_association.csv"), row.names = FALSE)
cat("\nSaved: results/subclonal_weights_mrd_association.csv\n")

pA <- ggplot(dat, aes(MRD, w1g, fill = MRD)) +
  geom_boxplot(width=.55, outlier.shape=NA, alpha=.7) +
  geom_jitter(width=.12, size=1.4, alpha=.6) +
  scale_fill_manual(values=c("MRD-"="#457B9D","MRD+"="#E63946"), guide="none") +
  labs(x=NULL, y="Dominant-subclone weight (gDNA)", title="Contemporaneous (visits)") +
  theme_bw(base_size=11)
pB <- pat %>% mutate(MRD=factor(ifelse(MRD_persistent==1,"MRD+","MRD-"),levels=c("MRD-","MRD+"))) %>%
  filter(!is.na(MRD), !is.na(sd_w1g)) %>%
  ggplot(aes(MRD, sd_w1g, fill=MRD)) +
  geom_boxplot(width=.55, outlier.shape=NA, alpha=.7) +
  geom_jitter(width=.12, size=1.8, alpha=.7) +
  scale_fill_manual(values=c("MRD-"="#457B9D","MRD+"="#E63946"), guide="none") +
  labs(x=NULL, y="Within-patient SD of w1g (instability)", title="Patient-level dynamics") +
  theme_bw(base_size=11)
fig <- gridExtra::grid.arrange(pA, pB, nrow = 1)
ggsave(file.path(fig_dir, "fig_weights_by_mrd.pdf"), fig, width = 8, height = 4)
cat("Saved: nature_manuscript/figures/fig_weights_by_mrd.pdf\n")

sink(file.path(res_dir, "sessionInfo_weights_assoc.txt")); print(sessionInfo()); sink()
cat("\nDone.\n")
