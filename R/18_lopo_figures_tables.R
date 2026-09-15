###############################################################################
## R/18_lopo_figures_tables.R
## Tables + figures for the Tier-4 leave-one-patient-out (LOPO) predictive
## comparison and the theta split-Rhat audit. Reads the fresh results and emits
## outputs for BOTH manuscript tracks. Presents whatever the numbers show
## (favours is computed from the sign of delta_elpd, not assumed).
##
## Inputs : results/lopo_cv.csv, lopo_per_patient.csv, lopo_convergence_gate.csv,
##          results/tier2/theta_rhat_M1.csv
## Outputs: tables/krd_lopo_comparison.tex, tables/krd_theta_rhat.tex          (methods paper)
##          figures/krd_lopo_per_patient.pdf                                    (methods paper)
##          nature_manuscript/table_lopo_comparison.csv, table_theta_rhat.csv   (applied paper)
##          nature_manuscript/figures/figS4_lopo_per_patient.pdf                (applied paper)
###############################################################################

suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
base_dir <- getwd()  # run from project root
setwd(base_dir)
res <- file.path(base_dir, "results")
fig_b <- file.path(base_dir, "figures");            dir.create(fig_b, showWarnings = FALSE)
tab_b <- file.path(base_dir, "tables");             dir.create(tab_b, showWarnings = FALSE)
nm    <- file.path(base_dir, "nature_manuscript")
fig_n <- file.path(nm, "figures")

cv  <- read.csv(file.path(res, "lopo_cv.csv"))
pp  <- read.csv(file.path(res, "lopo_per_patient.csv"))
gate<- read.csv(file.path(res, "lopo_convergence_gate.csv"))
th  <- read.csv(file.path(res, "tier2", "theta_rhat_M1.csv"))
fav <- function(x) ifelse(x > 0, "$\\mathcal{M}_1$", "$\\mathcal{M}_0$")
favc<- function(x) ifelse(x > 0, "M1", "M0")

cat("=== LOPO summary (from results/lopo_cv.csv) ===\n"); print(cv)
cat(sprintf("\nfolds used=%s  flagged={%s}  n_obs=%d  favouring M1=%d/%d\n",
            cv$folds_used, cv$folds_flagged, cv$n_obs, cv$n_obs_favouring_M1, cv$n_obs))

## ---------- Model-comparison table: total + source-specific delta elpd ----------
comp <- data.frame(
  Quantity = c("$\\Delta$elpd (total)", "$\\Delta$elpd (gDNA)", "$\\Delta$elpd (cfDNA)"),
  Delta = c(cv$delta_elpd, cv$delta_elpd_g, cv$delta_elpd_cf),
  SE    = c(cv$se_delta,   cv$se_delta_g,   cv$se_delta_cf))
comp$ratio <- comp$Delta / comp$SE
comp$Favours <- fav(comp$Delta)

## Methods paper LaTeX
tex <- c("\\begin{tabular}{lrrrl}", "\\toprule",
  "Predictive quantity & $\\Delta$elpd & SE & $\\Delta$/SE & Favours \\\\", "\\midrule",
  sprintf("%s & %s & %s & %.2f & %s \\\\",
          comp$Quantity, formatC(comp$Delta, format="f", big.mark=",", digits=1),
          formatC(comp$SE, format="f", big.mark=",", digits=1), comp$ratio, comp$Favours),
  "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(tab_b, "krd_lopo_comparison.tex"))

## Applied paper CSV
nat <- data.frame(Quantity = c("Delta elpd (total)","Delta elpd (gDNA)","Delta elpd (cfDNA)"),
                  Delta_elpd = round(comp$Delta,1), SE = round(comp$SE,1),
                  ratio = round(comp$ratio,2), Favours = favc(comp$Delta))
write.csv(nat, file.path(nm, "table_lopo_comparison.csv"), row.names = FALSE)

## ---------- theta Rhat summary table (occupied components) ----------
tho <- th[th$comp %in% c(1,2), ]
tho_sum <- tho %>% group_by(comp) %>%
  summarise(rhat_median = round(mean(rhat_median),3), rhat_p95 = round(max(rhat_p95),3),
            rhat_max = round(max(rhat_max),2), pct_gt_1.1 = round(max(pct_gt_1.1),2), .groups="drop")
write.csv(as.data.frame(tho_sum), file.path(nm, "table_theta_rhat.csv"), row.names = FALSE)
tex2 <- c("\\begin{tabular}{lrrrr}", "\\toprule",
  "Occupied component & median $\\hat R$ & 95th pctile $\\hat R$ & max $\\hat R$ & \\% $\\hat R>1.1$ \\\\",
  "\\midrule",
  sprintf("Component %d & %.3f & %.3f & %.2f & %.2f \\\\", tho_sum$comp, tho_sum$rhat_median,
          tho_sum$rhat_p95, tho_sum$rhat_max, tho_sum$pct_gt_1.1),
  "\\bottomrule", "\\end{tabular}")
writeLines(tex2, file.path(tab_b, "krd_theta_rhat.tex"))

## ---------- per-patient forest plot of delta elpd (M1 - M0) ----------
pp$patient <- factor(pp$patient, levels = pp$patient[order(pp$d)])
pfor <- ggplot(pp, aes(d, patient, color = d > 0)) +
  geom_vline(xintercept = 0, color = "grey60", linetype = "dashed") +
  geom_point(size = 2) +
  scale_color_manual(values = c(`TRUE` = "#1D3557", `FALSE` = "#E63946"),
                     labels = c(`TRUE` = "favours M1", `FALSE` = "favours M0"), name = NULL) +
  labs(x = expression("per-patient held-out "*Delta*"elpd ("*M[1]-M[0]*")"), y = "Patient",
       title = sprintf("Leave-one-patient-out: total delta-elpd = %.0f (favours %s)",
                       cv$delta_elpd, favc(cv$delta_elpd))) +
  theme_bw(base_size = 9) + theme(axis.text.y = element_text(size = 6))
ggsave(file.path(fig_b, "krd_lopo_per_patient.pdf"), pfor, width = 6, height = 6.5)
ggsave(file.path(fig_n, "figS4_lopo_per_patient.pdf"), pfor, width = 6, height = 6.5)

cat("\nWrote:\n  tables/krd_lopo_comparison.tex, tables/krd_theta_rhat.tex\n",
    "  figures/krd_lopo_per_patient.pdf\n",
    "  nature_manuscript/table_lopo_comparison.csv, table_theta_rhat.csv\n",
    "  nature_manuscript/figures/figS4_lopo_per_patient.pdf\n", sep="")
cat(sprintf("\nHEADLINE: total delta_elpd = %.1f (SE %.1f, ratio %.2f) -> favours %s | cfDNA-specific %.1f (SE %.1f)\n",
            cv$delta_elpd, cv$se_delta, cv$delta_elpd/cv$se_delta, favc(cv$delta_elpd),
            cv$delta_elpd_cf, cv$se_delta_cf))
