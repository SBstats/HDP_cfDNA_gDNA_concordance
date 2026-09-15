
###############################################################################
## ---- 17_manuscript_tables.R
###############################################################################

###############################################################################
## 17_manuscript_tables.R
## Generate all LaTeX tables for the manuscript
##
## Reads aggregated simulation results and KRd analysis output,
## produces .tex files that can be \input{} directly in the manuscript.
##
## Usage: Rscript R/17_manuscript_tables.R [base_dir]
###############################################################################

suppressPackageStartupMessages({
  library(xtable)
})

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[1] else getwd()

tab_dir <- file.path(base_dir, "tables")
dir.create(tab_dir, showWarnings = FALSE)

sim_dir <- file.path(base_dir, "results", "sim_aggregated")

# ===========================================================================
# SCENARIO 1: Signature estimation table (formerly Scenario 2)
# ===========================================================================
s2_file <- file.path(sim_dir, "scenario_2_summary.csv")
if (file.exists(s2_file)) {
  s2 <- read.csv(s2_file, stringsAsFactors = FALSE)

  tab2 <- data.frame(
    `$p$` = s2$p,
    `$s_k/p$` = s2$s_k_frac,
    MSE = sprintf("%.2f (%.2f)", s2$mse_mean, s2$mse_sd),
    `Norm. MSE` = sprintf("%.3f", s2$norm_mse_mean),
    TPR = sprintf("%.2f", s2$tpr_mean),
    FPR = sprintf("%.3f", s2$fpr_mean),
    check.names = FALSE
  )

  xt2 <- xtable(tab2, caption = "Scenario 2: Signature estimation.",
                 label = "tab:sim_s2_detail")
  print(xt2, file = file.path(tab_dir, "sim_scenario2_theta.tex"),
        include.rownames = FALSE, sanitize.text.function = identity,
        floating = FALSE)
  cat("Saved: sim_scenario2_theta.tex\n")
}

# ===========================================================================
# SCENARIO 3: Kappa/composition recovery table
# ===========================================================================
s3_file <- file.path(sim_dir, "scenario_3_summary.csv")
if (file.exists(s3_file)) {
  s3 <- read.csv(s3_file, stringsAsFactors = FALSE)
  s3a <- s3[s3$sub == "3a", ]

  if (nrow(s3a) > 0) {
    tab3 <- data.frame(
      `$\\kappa^*$` = s3a$kappa_true,
      Bias = sprintf("%.2f", s3a$kappa_bias),
      RMSE = sprintf("%.2f", s3a$kappa_rmse),
      `Coverage` = sprintf("%.2f", s3a$kappa_cover),
      `$\\hat{K}^+$` = sprintf("%.1f", s3a$Kplus_mean),
      `$\\|\\hat{\\pi}-\\pi^*\\|_1$` = sprintf("%.3f", s3a$pi_l1_mean),
      `Corr. error` = sprintf("%.3f", s3a$corr_error_mean),
      check.names = FALSE
    )

    xt3 <- xtable(tab3,
                   caption = "Scenario 3a: Recovery of $\\kappa$, $K^+$, and tracking correlation.",
                   label = "tab:sim_s3_detail")
    print(xt3, file = file.path(tab_dir, "sim_scenario3_kappa.tex"),
          include.rownames = FALSE, sanitize.text.function = identity,
          floating = FALSE)
    cat("Saved: sim_scenario3_kappa.tex\n")
  }
}

# ===========================================================================
# KRD: Posterior summary table
# ===========================================================================
krd_csv <- file.path(base_dir, "results", "krd_posterior_summary.csv")
if (file.exists(krd_csv)) {
  krd_summary <- read.csv(krd_csv, stringsAsFactors = FALSE)

  tab_krd <- data.frame(
    Parameter = krd_summary$parameter,
    Mean = sprintf("%.2f", krd_summary$mean),
    Median = sprintf("%.2f", krd_summary$median),
    `95\\% CI` = sprintf("(%.2f, %.2f)", krd_summary$ci_lo, krd_summary$ci_hi),
    check.names = FALSE
  )

  xt_krd <- xtable(tab_krd,
                    caption = "Posterior summaries under $\\mathcal{M}_1$ for the KRd trial.",
                    label = "tab:krd_summary_detail")
  print(xt_krd, file = file.path(tab_dir, "krd_posterior_summary.tex"),
        include.rownames = FALSE, sanitize.text.function = identity,
        floating = FALSE)
  cat("Saved: krd_posterior_summary.tex\n")
}

cat("\nTable generation complete.\n")


###############################################################################
## ---- 19_nature_tables.R
###############################################################################

###############################################################################
## 19_nature_tables.R
## Generate tables for the applied paper version
##
## Table 1 (study characteristics) uses observed data only.
## Tables 2-4 require model fitting results and will be generated later.
##
## Output: CSV files in nature_manuscript/ for integration into tables.md
##
## Usage: Rscript R/19_nature_tables.R [base_dir]
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
})

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[1] else
  getwd()

out_dir <- file.path(base_dir, "nature_manuscript")
res_dir <- file.path(base_dir, "results")
data_dir <- file.path(base_dir, "KRd trial", "5hmC data")
clin_dir <- file.path(base_dir, "KRd trial", "Clinical data")

cat("=== Applied Paper Table Generation ===\n\n")

# ===========================================================================
# TABLE 1: Patient Demographics & Study Characteristics
# ===========================================================================
cat("--- Table 1: Study Characteristics ---\n")

# Load sample keys
cfDNA_key <- read.csv(file.path(data_dir, "kRd-cfDNA_sample_key.csv"),
                      stringsAsFactors = FALSE)
gDNA_key  <- read.csv(file.path(data_dir, "kRd-gDNA_sample_key.csv"),
                      stringsAsFactors = FALSE)

# Load preprocessed model data for final dimensions
md <- readRDS(file.path(res_dir, "model_data.rds"))

# Identify paired observations
paired <- inner_join(
  cfDNA_key %>% select(Study.ID, Timepoint),
  gDNA_key  %>% select(Study.ID, Timepoint),
  by = c("Study.ID", "Timepoint")
)

# Timepoint distribution
tp_counts <- table(paired$Timepoint)
tp_order <- c("Screening", "C4", "C8", "C18", "3 YR F/U")
tp_counts_ordered <- tp_counts[tp_order]

# Per-patient sample counts
patient_counts <- paired %>%
  group_by(Study.ID) %>%
  summarize(n_timepoints = n(), .groups = "drop")

# Build Table 1
table1 <- data.frame(
  Characteristic = c(
    "Study design",
    "Treatment regimen",
    "Total patients enrolled",
    "Patients with paired samples",
    "Total cfDNA samples",
    "Total gDNA samples",
    "Total paired observations",
    "Timepoints per patient, median (range)",
    "",
    "Samples per timepoint:",
    paste0("  ", tp_order[1]),
    paste0("  ", tp_order[2]),
    paste0("  ", tp_order[3]),
    paste0("  ", tp_order[4]),
    paste0("  ", tp_order[5]),
    "",
    "Genomic features:",
    "  Genes profiled (total)",
    "  Genes after QC filtering",
    "  Filtering criterion",
    "",
    "Normalization",
    "5-hMC profiling method"
  ),
  Value = c(
    "Phase II clinical trial (KRd)",
    "Carfilzomib, lenalidomide, dexamethasone",
    as.character(nrow(cfDNA_key %>% distinct(Study.ID))),
    as.character(md$n),
    as.character(nrow(cfDNA_key)),
    as.character(nrow(gDNA_key)),
    as.character(md$N_obs),
    paste0(median(patient_counts$n_timepoints), " (",
           min(patient_counts$n_timepoints), "-",
           max(patient_counts$n_timepoints), ")"),
    "",
    "",
    as.character(tp_counts_ordered[1]),
    as.character(tp_counts_ordered[2]),
    as.character(tp_counts_ordered[3]),
    as.character(tp_counts_ordered[4]),
    as.character(tp_counts_ordered[5]),
    "",
    "",
    "19,100",
    as.character(md$p),
    "<10 counts in >5% of samples",
    "",
    "DESeq2 variance-stabilizing transformation (VST)",
    "hMe-Seal gene body 5-hMC"
  ),
  stringsAsFactors = FALSE
)

write.csv(table1, file.path(out_dir, "table1_study_characteristics.csv"),
          row.names = FALSE)

cat("  Patients with paired data:", md$n, "\n")
cat("  Total paired observations:", md$N_obs, "\n")
cat("  cfDNA samples:", nrow(cfDNA_key), "\n")
cat("  gDNA samples:", nrow(gDNA_key), "\n")
cat("  Genes after filtering:", md$p, "\n")
cat("  Timepoints per patient: median =", median(patient_counts$n_timepoints),
    ", range =", min(patient_counts$n_timepoints), "-",
    max(patient_counts$n_timepoints), "\n")
cat("  Timepoint distribution:\n")
print(tp_counts_ordered)
cat("  Saved: table1_study_characteristics.csv\n")

# ===========================================================================
# TABLE 2: Posterior Parameter Summaries (REQUIRES MODEL RESULTS)
# ===========================================================================
cat("\n--- Table 2: Posterior Parameter Summaries ---\n")

m1_file <- file.path(res_dir, "krd_M1_fit.rds")
if (file.exists(m1_file)) {
  m1 <- readRDS(m1_file)
  merged <- m1$merged

  # Both alpha_dp and kappa are fixed — report as constants
  param_summary <- data.frame(
    Parameter = c("K_plus", "sigma_0",
                  "omega_0 (cohort mean, per-observation)",
                  "kappa", "alpha", "rho (tracking correlation)"),
    Description = c("Occupied components",
                    "Background noise SD",
                    "Normal contamination fraction",
                    "Tracking precision",
                    "DP concentration",
                    "Cross-source correlation"),
    `Posterior Mean` = c(
      round(mean(merged$K_plus), 1),
      round(mean(merged$sigma_0), 3),
      ## Posterior mean of the cohort-average contamination: average over
      ## observations within each draw, then over draws.
      round(mean(rowMeans(merged$omega_0)), 3),
      1.0,
      1.0,
      round(2/3, 3)
    ),
    `95% CI` = c(
      paste0("(", round(quantile(merged$K_plus, 0.025), 0), ", ",
             round(quantile(merged$K_plus, 0.975), 0), ")"),
      paste0("(", round(quantile(merged$sigma_0, 0.025), 3), ", ",
             round(quantile(merged$sigma_0, 0.975), 3), ")"),
      ## A genuine POSTERIOR interval for the cohort mean, consistent with the
      ## K_plus and sigma_0 rows above. The previous version took quantiles of
      ## the per-observation posterior MEANS, which measures between-observation
      ## heterogeneity, not posterior uncertainty -- and silently mixed the two
      ## kinds of interval in one column.
      paste0("(", round(quantile(rowMeans(merged$omega_0), 0.025), 3), ", ",
             round(quantile(rowMeans(merged$omega_0), 0.975), 3), ")"),
      "Fixed",
      "Fixed",
      "Exact (kappa=alpha=1)"
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  write.csv(param_summary, file.path(out_dir, "table2_posterior_summary.csv"),
            row.names = FALSE)
  cat("  Saved: table2_posterior_summary.csv\n")
} else {
  cat("  KRd production results not found; Table 2 will be generated after analysis.\n")
}

# ===========================================================================
# TABLE 3: Model Comparison (REQUIRES MODEL RESULTS)
# ===========================================================================
cat("\n--- Table 3: Model Comparison ---\n")

mc_file <- file.path(res_dir, "krd_model_comparison.csv")
if (file.exists(mc_file)) {
  mc <- read.csv(mc_file, stringsAsFactors = FALSE)
  # Expected rows: delta_WAIC_all, delta_LOO_all, delta_WAIC_clean, delta_LOO_clean
  label_map <- c(delta_WAIC_all   = "$\\Delta$WAIC (all chains)",
                 delta_LOO_all    = "$\\Delta$PSIS-LOO (all chains)",
                 delta_WAIC_clean = "$\\Delta$WAIC (clean chains)",
                 delta_LOO_clean  = "$\\Delta$PSIS-LOO (clean chains)")
  mc$label <- ifelse(mc$metric %in% names(label_map),
                     label_map[mc$metric], mc$metric)
  model_comp <- data.frame(
    Criterion = mc$label,
    Value  = round(mc$value, 1),
    Favors = mc$favors,
    stringsAsFactors = FALSE
  )
  write.csv(model_comp, file.path(out_dir, "table3_model_comparison.csv"), row.names = FALSE)
  cat("  Saved: table3_model_comparison.csv\n")
} else {
  cat("  Model comparison CSV not found; Table 3 will be generated after analysis.\n")
}

# ===========================================================================
# TABLE 4: Per-Patient Concordance Metrics (REQUIRES MODEL RESULTS)
# ===========================================================================
cat("\n--- Table 4: Per-Patient Concordance ---\n")

wN_file <- file.path(res_dir, "krd_omega0_by_patient.csv")
if (file.exists(wN_file)) {
  wN_df <- read.csv(wN_file, stringsAsFactors = FALSE)
  ## Fail loudly on schema drift. sprintf() on a NULL column returns
  ## character(0), which only surfaces later as an opaque
  ## "arguments imply differing number of rows" from data.frame().
  if (!all(c("patient", "post_mean", "ci_lo", "ci_hi") %in% names(wN_df)))
    stop("krd_omega0_by_patient.csv is missing ci_lo/ci_hi; regenerate it with ",
         "the current slurm/run_krd_combine.R.")
  wN_df$post_mean_pct <- sprintf("%.4f%%", wN_df$post_mean * 100)
  wN_df$ci_pct <- sprintf("(%.5f%%--%.4f%%)", wN_df$ci_lo * 100, wN_df$ci_hi * 100)
  table4 <- data.frame(
    Patient = wN_df$patient,
    `Posterior Mean` = wN_df$post_mean_pct,
    `95% CI` = wN_df$ci_pct,
    check.names = FALSE, stringsAsFactors = FALSE
  )
  write.csv(table4, file.path(out_dir, "table4_contamination.csv"), row.names = FALSE)
  cat("  Saved: table4_contamination.csv\n")
} else {
  cat("  Per-patient wN CSV not found; Table 4 will be generated after analysis.\n")
}

cat("\n=== Applied paper table generation complete ===\n")

