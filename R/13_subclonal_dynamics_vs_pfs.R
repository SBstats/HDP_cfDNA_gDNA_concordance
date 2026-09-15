###############################################################################
## R/13_subclonal_dynamics_vs_pfs.R
##
## GOAL (new endpoint = PFS): the MRD binary endpoint gave no robust signal
## (R/11 concordance, R/12 weights/dynamics). PFS is a patient-level
## time-to-event endpoint with more information (12 events / 25 patients) and
## may have more power for a subclonal-instability predictor.
##
## Predictors (patient-level, from M1 canonical weights; gDNA = marrow truth):
##   Composition : mean_w1g, mean_dom_g (dominance), mean_ent_g (entropy)
##   Dynamics    : sd_w1g (within-patient SD of w1g = instability),
##                 switch (dominant subclone changes over t, argmax),
##                 sw10   (robust switch: w1g crosses 0.4<->0.6)
##   Baseline    : scr_w1g, scr_dom_g (screening composition, prognostic)
##   Concordance : mean_postcorr (patient-mean posterior correlation of omega,
##                 from krd_posterior_correlation.csv) -- closes the loop on the
##                 original concordance question against PFS. (L1 removed by design.)
##
## Endpoint: PFS = Surv(PFS_mo, PFS_event) from clinical file (25 pts; 101-08
## absent from clinical sheet -> dropped).
##
## Tests: univariable Cox PH per predictor (HR per SD for continuous; HR for
## binary), Wald + likelihood-ratio p; log-rank + Kaplan-Meier for switch.
## Small n / 12 events => EXPLORATORY (>=1 covariate per 10 events rule).
##
## Outputs (results/):  subclonal_dynamics_pfs_dataset.csv,
##                      subclonal_dynamics_pfs_cox.csv
## Outputs (figures/):  fig_pfs_by_switch_km.pdf, fig_pfs_by_instability_km.pdf
###############################################################################

set.seed(20260718)
suppressPackageStartupMessages({
  library(dplyr); library(readxl); library(survival); library(ggplot2)
})

base_dir <- getwd()
res_dir  <- file.path(base_dir, "results")
fig_dir  <- file.path(base_dir, "nature_manuscript", "figures")
clin_file <- file.path(base_dir, "KRd trial", "Clinical data",
                       "KRd 12_1725 data_2_14_2023 from Ben without PHI.xlsx")
source(file.path(base_dir, "R", "01_lib_core.R"))

## ---------------------------------------------------------------------------
## 1. Canonical posterior-mean weights per observation -> patient features
## ---------------------------------------------------------------------------
md  <- readRDS(file.path(res_dir, "model_data.rds")); n_obs <- md$N_obs
krd <- readRDS(file.path(res_dir, "krd_M1_fit.rds")); m1 <- krd
for (i in seq_along(m1$chains)) m1$chains[[i]] <- relabel_by_weight(m1$chains[[i]])
m1$chains <- align_chains_by_signature(m1$chains)
merged <- merge_chain_samples(m1$chains)
col_of <- function(k, o) (k - 1) * n_obs + o
w1g <- sapply(1:n_obs, function(o) mean(merged$omega_g_trace[, col_of(1, o)]))
w2g <- sapply(1:n_obs, function(o) mean(merged$omega_g_trace[, col_of(2, o)]))
ent <- function(a, b) { w <- c(a, b); w <- w[w > 1e-12]; w <- w / sum(w); -sum(w * log(w)) }

obs <- data.frame(patient_id = md$patient, t_num = md$time, w1g = w1g, w2g = w2g) %>%
  mutate(dom_g = pmax(w1g, w2g) / (w1g + w2g),
         ent_g = mapply(ent, w1g, w2g),
         dom_sub = ifelse(w1g >= w2g, 1L, 2L)) %>%
  group_by(patient_id, t_num) %>%              # average replicate rows
  summarise(w1g = mean(w1g), dom_g = mean(dom_g), ent_g = mean(ent_g),
            dom_sub = dom_sub[which.max(abs(w1g - 0.5))], .groups = "drop")

id_map <- unique(md$paired_info[, c("patient_id", "Study.ID")])
scr <- obs %>% filter(t_num == 1) %>% transmute(patient_id, scr_w1g = w1g, scr_dom_g = dom_g)

feat <- obs %>% group_by(patient_id) %>%
  summarise(nvis = n(),
            mean_w1g = mean(w1g), mean_dom_g = mean(dom_g), mean_ent_g = mean(ent_g),
            sd_w1g = ifelse(n() > 1, sd(w1g), NA_real_),
            switch = as.integer(n_distinct(dom_sub) > 1),
            sw10   = as.integer(any(w1g > 0.60) & any(w1g < 0.40)),
            .groups = "drop") %>%
  left_join(scr, by = "patient_id") %>%
  left_join(id_map, by = "patient_id")

## patient-mean concordance (posterior correlation of omega) -> loop back on
## the original question. (The L1 concordance metric was removed by design.)
cpc <- read.csv(file.path(res_dir, "krd_posterior_correlation.csv"), stringsAsFactors = FALSE)
# posterior_correlation.csv is keyed by Study.ID + timepoint label; average per patient
cpc <- cpc %>% rename(Study.ID = patient) %>%
  group_by(Study.ID) %>% summarise(mean_postcorr = mean(post_corr, na.rm = TRUE), .groups = "drop")
feat <- left_join(feat, cpc, by = "Study.ID")

## ---------------------------------------------------------------------------
## 2. PFS endpoint
## ---------------------------------------------------------------------------
clin <- suppressWarnings(read_excel(clin_file, sheet = "12_1725 cfDNA Dataset"))
clin <- as.data.frame(clin)
pfs <- data.frame(Study.ID = clin[[grep("Study number", names(clin))[1]]],
                  PFS_mo = clin[[41]], PFS_event = clin[[34]],
                  best_response = clin$`best overall response`)
pfs <- pfs[!is.na(pfs$PFS_mo) & !is.na(pfs$PFS_event), ]

dat <- inner_join(feat, pfs, by = "Study.ID")
write.csv(dat, file.path(res_dir, "subclonal_dynamics_pfs_dataset.csv"), row.names = FALSE)
cat("\nPFS analysis set: n =", nrow(dat), " events =", sum(dat$PFS_event),
    " censored =", sum(dat$PFS_event == 0), "\n")
cat("switch: ", sum(dat$switch == 1), "switchers /", nrow(dat),
    " | sw10:", sum(dat$sw10 == 1), "\n")

## ---------------------------------------------------------------------------
## 3. Univariable Cox per predictor
## ---------------------------------------------------------------------------
cont_vars <- c("mean_w1g","mean_dom_g","mean_ent_g","sd_w1g",
               "scr_w1g","scr_dom_g","mean_postcorr")
bin_vars  <- c("switch","sw10")

cox_one <- function(v, binary = FALSE) {
  d <- dat[!is.na(dat[[v]]), ]
  x <- if (binary) d[[v]] else as.numeric(scale(d[[v]]))   # per-SD for continuous
  fit <- tryCatch(coxph(Surv(PFS_mo, PFS_event) ~ x, data = d), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s <- summary(fit)
  data.frame(predictor = v, type = ifelse(binary, "binary(1 vs 0)", "per_SD"),
             n = nrow(d), events = sum(d$PFS_event),
             HR = s$coef[1, "exp(coef)"],
             ci_lo = s$conf.int[1, "lower .95"], ci_hi = s$conf.int[1, "upper .95"],
             p_wald = s$coef[1, "Pr(>|z|)"], p_LRT = s$logtest["pvalue"])
}
cox_res <- bind_rows(
  bind_rows(lapply(cont_vars, cox_one, binary = FALSE)),
  bind_rows(lapply(bin_vars,  cox_one, binary = TRUE))
)
cat("\n==== Univariable Cox PH: PFS ~ predictor ====\n")
print(cox_res, digits = 3, row.names = FALSE)

## log-rank for the binary switch predictors + median PFS by group
cat("\n==== Log-rank (PFS by dominant-subclone switching) ====\n")
for (v in bin_vars) {
  d <- dat[!is.na(dat[[v]]), ]; d$grp <- d[[v]]
  sd_ <- survdiff(Surv(PFS_mo, PFS_event) ~ grp, data = d)
  p <- 1 - pchisq(sd_$chisq, df = 1)
  km <- survfit(Surv(PFS_mo, PFS_event) ~ grp, data = d)
  med <- summary(km)$table[, "median"]
  cat(sprintf("  %-7s log-rank p = %.4f | median PFS(mo): %s\n",
              v, p, paste(sprintf("%s=%.1f", names(med), med), collapse = "  ")))
}

write.csv(cox_res, file.path(res_dir, "subclonal_dynamics_pfs_cox.csv"), row.names = FALSE)
cat("\nSaved: results/subclonal_dynamics_pfs_cox.csv\n")

## ---------------------------------------------------------------------------
## 4. Kaplan-Meier figures (switch; instability median split)
## ---------------------------------------------------------------------------
km_plot <- function(group, glabels, title, file) {
  d <- dat[!is.na(group), ]; d$g <- factor(group[!is.na(group)], labels = glabels)
  km <- survfit(Surv(PFS_mo, PFS_event) ~ g, data = d)
  sd_ <- survdiff(Surv(PFS_mo, PFS_event) ~ g, data = d); p <- 1 - pchisq(sd_$chisq, df = 1)
  # build step data
  sm <- data.frame(time = km$time, surv = km$surv,
                   grp = rep(glabels, km$strata))
  ggplot(sm, aes(time, surv, color = grp)) +
    geom_step(linewidth = 0.9) +
    scale_color_manual(values = c("#457B9D", "#E63946"), name = NULL) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "PFS (months)", y = "Progression-free survival",
         title = title, subtitle = sprintf("log-rank p = %.3f", p)) +
    theme_bw(base_size = 11) + theme(legend.position = c(0.8, 0.85))
}
p_sw <- km_plot(dat$switch, c("no switch", "switch"),
                "PFS by dominant-subclone switching", NULL)
ggsave(file.path(fig_dir, "fig_pfs_by_switch_km.pdf"), p_sw, width = 5.5, height = 4)

inst_hi <- as.integer(dat$sd_w1g > median(dat$sd_w1g, na.rm = TRUE))
p_in <- km_plot(inst_hi, c("low instability", "high instability"),
                "PFS by subclonal instability (SD of w1g, median split)", NULL)
ggsave(file.path(fig_dir, "fig_pfs_by_instability_km.pdf"), p_in, width = 5.5, height = 4)
cat("Saved: fig_pfs_by_switch_km.pdf, fig_pfs_by_instability_km.pdf\n")

sink(file.path(res_dir, "sessionInfo_pfs.txt")); print(sessionInfo()); sink()
cat("\nDone.\n")
