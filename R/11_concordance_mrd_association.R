###############################################################################
## R/11_concordance_mrd_association.R
##
## GOAL: formally establish the longitudinal association between model-based
## cfDNA-gDNA subclonal CONCORDANCE and the clinical OUTCOME (MRD by NGS).
## Figures 3a (concordance) and 3b (MRD) are currently read side-by-side only;
## this script turns that visual juxtaposition into a statistical test.
##
## Design (pre-specified):
##   Primary metric   : posterior correlation of omega_cf and omega_g (C_postcorr)
##                      -- the single subclonal-concordance measure in this paper
##                      (the L1 metric was removed by design)
##   Primary outcome  : MRD-positive (1) vs MRD-negative (0), contemporaneous
##   Primary design   : concordance@{C8,C18,3Y} matched to MRD@same window
##   Primary test     : LMM C ~ MRD + (1|patient) with CLUSTER-PERMUTATION p-value
##   Secondary        : prediction corr (theta-reconstructed, C_pred);
##                      patient-level Wilcoxon vs MRD-persistence / best response;
##                      prognostic screening concordance -> later MRD.
##
## Clustering (26 patients, 84 obs) is handled by (i) random patient intercepts
## and (ii) permuting MRD at the PATIENT level, which is exact under the null of
## no concordance-MRD association and robust to small n / separation.
##
## MRD coding: clinical file uses "negative = 1"; we recode MRD_pos = 1 - value
## so that 1 = residual disease present (MRD-positive).
##
## Outputs (results/):
##   concordance_mrd_association.csv        tidy table of all tests
##   concordance_mrd_analysis_dataset.csv   the merged obs-level analysis data
##   sessionInfo_assoc.txt
## Outputs (nature_manuscript/figures/):
##   fig_concordance_by_mrd.pdf/.png        concordance by MRD status (obs + patient)
###############################################################################

set.seed(20260718)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(lme4)
  library(ggplot2)
})

base_dir <- getwd()
res_dir  <- file.path(base_dir, "results")
fig_dir  <- file.path(base_dir, "nature_manuscript", "figures")
clin_file <- file.path(base_dir, "KRd trial", "Clinical data",
                       "KRd 12_1725 data_2_14_2023 from Ben without PHI.xlsx")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

tp_levels <- c("Screening", "C4", "C8", "C18", "3 YR F/U")   # t_num 1..5

## ---------------------------------------------------------------------------
## 1. id/timepoint maps from model_data
## ---------------------------------------------------------------------------
md <- readRDS(file.path(res_dir, "model_data.rds"))
paired <- md$paired_info
id_map <- unique(paired[, c("patient_id", "Study.ID")])          # int <-> Study.ID
id_map <- id_map[order(id_map$patient_id), ]

## normalize any timepoint label -> t_num (1..5)
norm_t <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    x %in% c("Screening", "SCR", "1")            ~ 1L,
    x %in% c("C4", "2")                           ~ 2L,
    x %in% c("C8", "3")                           ~ 3L,
    x %in% c("C18", "4")                          ~ 4L,
    x %in% c("3 YR F/U", "3Y", "3YR", "5")        ~ 5L,
    TRUE                                          ~ NA_integer_
  )
}

## ---------------------------------------------------------------------------
## 2. MRD outcome (validated Fig-3b mapping), long by patient x window
##    clinical: negative = 1  -> MRD_pos = 1 - value  (1 = residual disease)
## ---------------------------------------------------------------------------
suppressWarnings({ clin <- read_excel(clin_file, sheet = "12_1725 cfDNA Dataset") })
mrd_df <- data.frame(
  Study.ID      = clin$`Study number`,
  best_response = clin$`best overall response`,
  neg_c8        = clin$`MRD post 8 cycles by NGS (Adaptive) negative =1`,
  neg_c18       = clin$`MRD by NGS EoT Negative = 1`,
  neg_1yr       = clin$`MRD by NGS 1yr f/u`,
  neg_2yr       = clin$`MRD by NGS 2yr f/u`,
  neg_3yr       = clin$`MRD by NGS 3yr f/u`,
  stringsAsFactors = FALSE
)
mrd_df$neg_3yrFU <- ifelse(!is.na(mrd_df$neg_3yr), mrd_df$neg_3yr,
                    ifelse(!is.na(mrd_df$neg_2yr), mrd_df$neg_2yr, mrd_df$neg_1yr))
to_pos <- function(neg) ifelse(is.na(neg), NA_integer_, as.integer(1 - neg))

mrd_long <- bind_rows(
  data.frame(Study.ID = mrd_df$Study.ID, t_num = 3L, MRD_pos = to_pos(mrd_df$neg_c8)),
  data.frame(Study.ID = mrd_df$Study.ID, t_num = 4L, MRD_pos = to_pos(mrd_df$neg_c18)),
  data.frame(Study.ID = mrd_df$Study.ID, t_num = 5L, MRD_pos = to_pos(mrd_df$neg_3yrFU))
)

## patient-level MRD summaries
# MRD_ever_pos: positive at ANY assessed window (min neg = 0 -> pos). NA if all NA.
# MRD_persistent: positive at ALL assessed windows (max neg = 0 -> all pos).
fin <- function(x) ifelse(is.finite(x), x, NA_real_)
mrd_pat <- data.frame(
  Study.ID       = mrd_df$Study.ID,
  best_response  = mrd_df$best_response,
  MRD_ever_pos   = to_pos(fin(pmin(mrd_df$neg_c8, mrd_df$neg_c18, mrd_df$neg_3yrFU, na.rm = TRUE))),
  MRD_persistent = to_pos(fin(pmax(mrd_df$neg_c8, mrd_df$neg_c18, mrd_df$neg_3yrFU, na.rm = TRUE)))
)

## ---------------------------------------------------------------------------
## 3. Concordance metrics -> obs level (patient_id, t_num, metric)
## ---------------------------------------------------------------------------
# The L1 concordance metric was removed by design. The PRIMARY concordance
# predictor is the posterior correlation of omega_g and omega_cf (C_postcorr);
# the model-based prediction correlation (C_pred) is retained as a secondary,
# theta-reconstructed measure.
# prediction correlation (theta-reconstructed: Corr(Yhat_g, Yhat_cf))
pc <- read.csv(file.path(res_dir, "krd_prediction_correlation.csv"))
pc <- pc %>% transmute(patient_id = patient, t_num = time, C_pred = pred_correlation)

# posterior correlation of omega_g, omega_cf (PRIMARY; keyed by Study.ID + label)
poc <- read.csv(file.path(res_dir, "krd_posterior_correlation.csv"),
                stringsAsFactors = FALSE)
poc <- poc %>%
  rename(Study.ID = patient) %>%
  mutate(t_num = norm_t(timepoint)) %>%
  left_join(id_map, by = "Study.ID") %>%
  transmute(patient_id, t_num, C_postcorr = post_corr)

## merge metrics; average duplicate patient x t_num rows (replicate samples)
metrics <- poc %>%
  full_join(pc,    by = c("patient_id", "t_num")) %>%
  group_by(patient_id, t_num) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(across(everything(), ~ ifelse(is.nan(.x), NA, .x))) %>%
  left_join(id_map, by = "patient_id")

## ---------------------------------------------------------------------------
## 4. Analysis dataset: contemporaneous (MRD windows only)
## ---------------------------------------------------------------------------
dat <- metrics %>%
  inner_join(mrd_long, by = c("Study.ID", "t_num")) %>%   # only t in {3,4,5}
  left_join(mrd_pat %>% select(Study.ID, best_response), by = "Study.ID") %>%
  filter(!is.na(MRD_pos)) %>%
  mutate(window = c(`3` = "C8", `4` = "C18", `5` = "3Y")[as.character(t_num)],
         patient = factor(patient_id),
         MRD = factor(ifelse(MRD_pos == 1, "MRD+", "MRD-"), levels = c("MRD-", "MRD+")))

write.csv(dat, file.path(res_dir, "concordance_mrd_analysis_dataset.csv"),
          row.names = FALSE)

cat("\n==== Contemporaneous analysis dataset ====\n")
cat("Matched patient-visits:", nrow(dat),
    "| patients:", n_distinct(dat$Study.ID),
    "| MRD+:", sum(dat$MRD_pos), " MRD-:", sum(dat$MRD_pos == 0), "\n")
print(table(window = dat$window, MRD = dat$MRD))

## ---------------------------------------------------------------------------
## 5. Inference engine
##    - LMM  concordance ~ MRD_pos + (1|patient)   (lme4)
##    - cluster-permutation p-value: permute each patient's MRD vector
##    - GLMM logistic MRD ~ concordance + (1|patient) for OR
## ---------------------------------------------------------------------------
n_perm <- 2000

## permute MRD_pos WITHIN the design by shuffling patient-level status labels.
## Each patient keeps its (constant-or-varying) MRD vector, but which concordance
## profile it is attached to is randomized => exact under H0 of no association,
## preserving within-patient clustering.
perm_patient_labels <- function(df, seed_i) {
  set.seed(seed_i)
  pats <- unique(df$Study.ID)
  # map each patient to a permuted donor patient's MRD-by-window vector
  donor <- sample(pats)
  key <- df %>% distinct(Study.ID, t_num, MRD_pos)
  lut <- setNames(donor, pats)
  # attach donor's MRD for the same window
  don_mrd <- key %>% rename(donorID = Study.ID)
  out <- df %>%
    mutate(donorID = lut[as.character(Study.ID)]) %>%
    select(-MRD_pos) %>%
    left_join(key %>% rename(donorID = Study.ID), by = c("donorID", "t_num"))
  out$MRD_pos
}

lmm_effect <- function(df, yvar) {
  d <- df[!is.na(df[[yvar]]) & !is.na(df$MRD_pos), ]
  if (n_distinct(d$MRD_pos) < 2 || nrow(d) < 6) return(NULL)
  d$y <- d[[yvar]]
  fit  <- suppressMessages(lmer(y ~ MRD_pos + (1 | patient), data = d,
                                REML = FALSE,
                                control = lmerControl(check.conv.singular = "ignore")))
  fit0 <- suppressMessages(lmer(y ~ 1 + (1 | patient), data = d, REML = FALSE,
                                control = lmerControl(check.conv.singular = "ignore")))
  beta <- fixef(fit)[["MRD_pos"]]
  lrt  <- anova(fit0, fit)
  p_lrt <- lrt$`Pr(>Chisq)`[2]
  # cluster permutation on beta
  obs_stat <- abs(beta)
  perm_stat <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    d$MRDp <- perm_patient_labels(d, i)
    dd <- d[!is.na(d$MRDp), ]
    if (n_distinct(dd$MRDp) < 2) { perm_stat[i] <- NA; next }
    bi <- tryCatch(fixef(suppressMessages(lmer(y ~ MRDp + (1|patient), data = dd,
              REML = FALSE, control = lmerControl(check.conv.singular = "ignore"))))[["MRDp"]],
              error = function(e) NA)
    perm_stat[i] <- abs(bi)
  }
  p_perm <- (1 + sum(perm_stat >= obs_stat, na.rm = TRUE)) /
            (1 + sum(!is.na(perm_stat)))
  # 95% CI via profile-free Wald from lmer vcov
  se <- sqrt(vcov(fit)["MRD_pos", "MRD_pos"])
  data.frame(metric = yvar, n = nrow(d),
             mean_MRDneg = mean(d$y[d$MRD_pos == 0]),
             mean_MRDpos = mean(d$y[d$MRD_pos == 1]),
             beta = beta, ci_lo = beta - 1.96*se, ci_hi = beta + 1.96*se,
             p_LRT = p_lrt, p_perm = p_perm)
}

metric_vars <- c("C_postcorr", "C_pred")
metric_lab  <- c(C_postcorr = "Posterior correlation of omega (primary)",
                 C_pred = "Prediction correlation (theta-reconstructed)")

cat("\n==== PRIMARY + secondary: LMM concordance ~ MRD + (1|patient) ====\n")
lmm_res <- bind_rows(lapply(metric_vars, function(v) lmm_effect(dat, v)))
lmm_res$metric_label <- metric_lab[lmm_res$metric]
lmm_res$analysis <- "contemporaneous_LMM"
print(lmm_res, digits = 3)

## logistic GLMM for OR (primary metric only), permutation p
glmm_or <- function(df, yvar) {
  d <- df[!is.na(df[[yvar]]), ]; d$y <- d[[yvar]]
  fit <- tryCatch(suppressMessages(glmer(MRD_pos ~ scale(y) + (1|patient),
           data = d, family = binomial,
           control = glmerControl(optimizer = "bobyqa"))), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  co <- summary(fit)$coefficients["scale(y)", ]
  data.frame(metric = yvar, OR_per_SD = exp(co["Estimate"]),
             OR_lo = exp(co["Estimate"] - 1.96*co["Std. Error"]),
             OR_hi = exp(co["Estimate"] + 1.96*co["Std. Error"]),
             p_wald = co["Pr(>|z|)"])
}
cat("\n==== Logistic GLMM: MRD ~ concordance(SD) + (1|patient) ====\n")
glmm_res <- bind_rows(lapply(metric_vars, function(v) glmm_or(dat, v)))
if (nrow(glmm_res)) print(glmm_res, digits = 3)

## ---------------------------------------------------------------------------
## 6. Patient-level analyses (no pseudo-replication)
## ---------------------------------------------------------------------------
pat_lvl <- metrics %>%
  group_by(Study.ID) %>%
  summarise(across(all_of(metric_vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  left_join(mrd_pat, by = "Study.ID")

wilcox_by <- function(df, yvar, grp) {
  d <- df[!is.na(df[[yvar]]) & !is.na(df[[grp]]), ]
  g <- factor(d[[grp]])
  if (nlevels(g) != 2) return(NULL)
  w <- suppressWarnings(wilcox.test(d[[yvar]] ~ g))
  lv <- levels(g)
  data.frame(metric = yvar, group = grp,
             grp0 = lv[1], median0 = median(d[[yvar]][g == lv[1]]), n0 = sum(g == lv[1]),
             grp1 = lv[2], median1 = median(d[[yvar]][g == lv[2]]), n1 = sum(g == lv[2]),
             p_wilcox = w$p.value)
}
cat("\n==== Patient-level Wilcoxon: mean concordance by MRD group / response ====\n")
pat_res <- bind_rows(
  lapply(metric_vars, wilcox_by, df = pat_lvl, grp = "MRD_persistent"),
  lapply(metric_vars, wilcox_by, df = pat_lvl, grp = "MRD_ever_pos"),
  lapply(metric_vars, function(v) wilcox_by(
     pat_lvl %>% mutate(resp2 = ifelse(best_response == "sCR", "sCR", "VGPR")),
     v, "resp2"))
)
print(pat_res, digits = 3)

## ---------------------------------------------------------------------------
## 7. Prognostic: SCREENING concordance -> later MRD (patient level)
## ---------------------------------------------------------------------------
scr <- metrics %>% filter(t_num == 1) %>%
  select(Study.ID, all_of(metric_vars)) %>%
  rename_with(~ paste0("scr_", .x), all_of(metric_vars)) %>%
  left_join(mrd_pat, by = "Study.ID")
cat("\n==== Prognostic: screening concordance vs MRD_ever_pos ====\n")
prog_res <- bind_rows(lapply(metric_vars, function(v) {
  yv <- paste0("scr_", v)
  r <- wilcox_by(scr, yv, "MRD_ever_pos"); if (!is.null(r)) r$metric <- v; r
}))
if (nrow(prog_res)) print(prog_res, digits = 3)

## ---------------------------------------------------------------------------
## 8. Save tidy results
## ---------------------------------------------------------------------------
out <- bind_rows(
  lmm_res %>% transmute(analysis, metric = metric_label, n,
                        est = beta, ci_lo, ci_hi, p_primary = p_perm, p_secondary = p_LRT,
                        detail = sprintf("mean MRD- %.3f vs MRD+ %.3f", mean_MRDneg, mean_MRDpos)),
  pat_res %>% transmute(analysis = paste0("patient_", group),
                        metric = metric_lab[metric], n = n0 + n1,
                        est = median1 - median0, ci_lo = NA, ci_hi = NA,
                        p_primary = p_wilcox, p_secondary = NA,
                        detail = sprintf("%s med %.3f (n=%d) vs %s med %.3f (n=%d)",
                                         grp0, median0, n0, grp1, median1, n1))
)
write.csv(out, file.path(res_dir, "concordance_mrd_association.csv"), row.names = FALSE)
cat("\nSaved: results/concordance_mrd_association.csv\n")

## ---------------------------------------------------------------------------
## 9. Figure: concordance by MRD status (obs-level + patient-level)
## ---------------------------------------------------------------------------
pA <- ggplot(dat, aes(MRD, C_postcorr, fill = MRD)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, height = 0, size = 1.4, alpha = 0.6) +
  scale_fill_manual(values = c("MRD-" = "#457B9D", "MRD+" = "#E63946"), guide = "none") +
  labs(x = NULL, y = expression("Posterior correlation "*rho(omega^g, omega^cf)),
       title = "Contemporaneous (patient-visits)") +
  theme_bw(base_size = 11)

pB <- pat_lvl %>%
  mutate(MRD = factor(ifelse(MRD_persistent == 1, "MRD+", "MRD-"),
                      levels = c("MRD-", "MRD+"))) %>%
  filter(!is.na(MRD)) %>%
  ggplot(aes(MRD, C_postcorr, fill = MRD)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, height = 0, size = 1.8, alpha = 0.7) +
  scale_fill_manual(values = c("MRD-" = "#457B9D", "MRD+" = "#E63946"), guide = "none") +
  labs(x = NULL, y = "Mean posterior correlation of omega",
       title = "Patient-level (persistent MRD)") +
  theme_bw(base_size = 11)

fig <- gridExtra::grid.arrange(pA, pB, nrow = 1)
ggsave(file.path(fig_dir, "fig_concordance_by_mrd.pdf"), fig, width = 8, height = 4)
tryCatch(
  ggsave(file.path(fig_dir, "fig_concordance_by_mrd.png"), fig,
         width = 8, height = 4, dpi = 300, type = "cairo"),
  error = function(e) message("PNG export skipped (", conditionMessage(e), "); PDF written."))
cat("Saved: nature_manuscript/figures/fig_concordance_by_mrd.pdf\n")

sink(file.path(res_dir, "sessionInfo_assoc.txt")); print(sessionInfo()); sink()
cat("\nDone.\n")
