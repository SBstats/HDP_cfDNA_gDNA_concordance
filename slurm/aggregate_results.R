## aggregate_results.R
## Reads all rep_*.rds files from results/sim/scenario_*/config_*/ and writes
## one CSV per scenario to results/sim_aggregated/.
## Usage: Rscript slurm/aggregate_results.R [base_dir]

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[1] else getwd()

raw_dir <- file.path(base_dir, "results", "sim")
out_dir <- file.path(base_dir, "results", "sim_aggregated")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mcse_mean <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
mcse_prop <- function(x) {
  p <- mean(x, na.rm = TRUE); n <- sum(!is.na(x))
  sqrt(p * (1 - p) / n)
}

load_reps <- function(scen) {
  sdir <- file.path(raw_dir, sprintf("scenario_%d", scen))
  if (!dir.exists(sdir)) return(list())
  cfg_dirs <- sort(list.dirs(sdir, recursive = FALSE))
  reps <- list()
  for (cd in cfg_dirs) {
    fls <- list.files(cd, pattern = "^rep_.*\\.rds$", full.names = TRUE)
    for (f in fls) {
      x <- tryCatch(readRDS(f), error = function(e) NULL)
      if (!is.null(x)) reps[[length(reps) + 1]] <- x
    }
  }
  reps
}

## ---- Scenario 1: signature estimation -----------------------------------
agg_s1 <- function(reps) {
  rows <- lapply(reps, function(x) {
    th <- x$theta; comp <- x$competitors
    data.frame(
      config_id    = x$config_id,
      p            = x$config$p,
      s_k_frac     = x$config$s_k_frac,
      mse          = if (!is.null(th)) th$mse else NA,
      normalized_mse = if (!is.null(th)) th$normalized_mse else NA,
      tpr          = if (!is.null(th)) th$tpr else NA,
      fpr          = if (!is.null(th)) th$fpr else NA,
      recovery_rate = if (!is.null(th)) th$recovery_rate else NA,
      coverage_active   = if (!is.null(th)) th$coverage_active else NA,
      coverage_inactive = if (!is.null(th)) th$coverage_inactive else NA,
      lasso_mse    = if (!is.null(comp) && !is.null(comp$lasso))         comp$lasso$mse          else NA,
      spike_slab_mse = if (!is.null(comp) && !is.null(comp$spike_slab)) comp$spike_slab$mse     else NA,
      unregularized_mse = if (!is.null(comp) && !is.null(comp$unregularized)) comp$unregularized$mse else NA,
      lasso_norm_mse = if (!is.null(comp) && !is.null(comp$lasso))      comp$lasso$normalized_mse else NA,
      spike_slab_norm_mse = if (!is.null(comp) && !is.null(comp$spike_slab)) comp$spike_slab$normalized_mse else NA,
      unregularized_norm_mse = if (!is.null(comp) && !is.null(comp$unregularized)) comp$unregularized$normalized_mse else NA,
      correct_selection = as.integer(x$correct_selection),
      delta_waic   = x$delta_waic,
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)
  keys <- c("config_id", "p", "s_k_frac")
  agg <- do.call(rbind, lapply(split(df, df$config_id), function(g) {
    n <- nrow(g)
    data.frame(
      config_id = g$config_id[1], p = g$p[1], s_k_frac = g$s_k_frac[1],
      n_reps = n,
      mse_mean = mean(g$mse, na.rm = TRUE), mse_mcse = mcse_mean(g$mse),
      norm_mse_mean = mean(g$normalized_mse, na.rm = TRUE),
      norm_mse_mcse = mcse_mean(g$normalized_mse),
      tpr_mean = mean(g$tpr, na.rm = TRUE), fpr_mean = mean(g$fpr, na.rm = TRUE),
      recovery_rate_mean = mean(g$recovery_rate, na.rm = TRUE),
      coverage_active_mean = mean(g$coverage_active, na.rm = TRUE),
      coverage_inactive_mean = mean(g$coverage_inactive, na.rm = TRUE),
      lasso_mse_mean = mean(g$lasso_mse, na.rm = TRUE),
      spike_slab_mse_mean = mean(g$spike_slab_mse, na.rm = TRUE),
      unregularized_mse_mean = mean(g$unregularized_mse, na.rm = TRUE),
      lasso_norm_mse_mean = mean(g$lasso_norm_mse, na.rm = TRUE),
      spike_slab_norm_mse_mean = mean(g$spike_slab_norm_mse, na.rm = TRUE),
      unregularized_norm_mse_mean = mean(g$unregularized_norm_mse, na.rm = TRUE),
      model_sel_rate = mean(g$correct_selection, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  agg[order(agg$s_k_frac, agg$p), ]
}

## ---- Scenario 2: signature recovery vs n (REDESIGNED) ------------------
agg_s2 <- function(reps) {
  rows <- lapply(reps, function(x) {
    th <- x$theta
    rt <- x$ident   # identifiability metrics (rho_theta retired)
    data.frame(
      config_id    = x$config_id,
      n            = x$config$n,
      p            = x$config$p,
      recovery_rate = if (!is.null(th)) th$recovery_rate else NA,
      tpr          = if (!is.null(th)) th$tpr else NA,
      fpr          = if (!is.null(th)) th$fpr else NA,
      mse          = if (!is.null(th)) th$mse else NA,
      coverage_active = if (!is.null(th)) th$coverage_active else NA,
      thetabar_cor   = if (!is.null(rt)) rt$thetabar_cor_mean else NA,
      thetabar_rmse  = if (!is.null(rt)) rt$thetabar_rmse else NA,
      varsigma_k_bias= if (!is.null(rt)) rt$varsigma_k_bias else NA,
      omega0_rmse    = if (!is.null(rt)) rt$omega0_rmse else NA,
      pi_l1_ident    = if (!is.null(rt)) rt$pi_l1 else NA,
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)
  agg <- do.call(rbind, lapply(split(df, df$config_id), function(g) {
    data.frame(
      config_id    = g$config_id[1], n = g$n[1], p = g$p[1],
      n_reps       = nrow(g),
      recovery_rate_mean = mean(g$recovery_rate, na.rm = TRUE),
      recovery_rate_mcse = mcse_mean(g$recovery_rate),
      tpr_mean     = mean(g$tpr, na.rm = TRUE),
      fpr_mean     = mean(g$fpr, na.rm = TRUE),
      mse_mean     = mean(g$mse, na.rm = TRUE),
      coverage_active_mean = mean(g$coverage_active, na.rm = TRUE),
      thetabar_cor_mean   = mean(g$thetabar_cor, na.rm = TRUE),
      thetabar_rmse_mean  = mean(g$thetabar_rmse, na.rm = TRUE),
      varsigma_k_bias_mean= mean(g$varsigma_k_bias, na.rm = TRUE),
      omega0_rmse_mean    = mean(g$omega0_rmse, na.rm = TRUE),
      pi_l1_ident_mean    = mean(g$pi_l1_ident, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  agg[order(agg$n), ]
}

## ---- Scenario 3: sample-size power curve (REDESIGNED) ------------------
agg_s3 <- function(reps) {
  rows <- lapply(reps, function(x) {
    sel <- x$selection
    data.frame(
      config_id  = x$config_id,
      sub        = x$config$sub,
      true_model = x$config$true_model,
      n          = x$config$n,
      kappa_true = x$config$kappa_true,
      favors_M1  = as.integer(x$favors_M1),
      correct_selection = as.integer(x$correct_selection),
      false_positive = if (!is.null(sel)) as.integer(sel$false_positive) else
                       as.integer(x$config$true_model == "M0" && !is.na(x$favors_M1) && x$favors_M1),
      incr_excess = if (!is.null(sel)) sel$incr_excess else NA,
      incr_p      = if (!is.null(sel)) sel$incr_p      else NA,
      cor_excess  = if (!is.null(sel)) sel$cor_excess   else NA,
      conc_excess = if (!is.null(sel)) sel$conc_excess  else NA,
      tracking_score = if (!is.null(sel)) sel$tracking_score else NA,
      delta_waic  = x$delta_waic,
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)
  agg <- do.call(rbind, lapply(split(df, df$config_id), function(g) {
    data.frame(
      config_id   = g$config_id[1],
      sub         = g$sub[1],
      true_model  = g$true_model[1],
      n           = g$n[1],
      kappa_true  = g$kappa_true[1],
      n_reps      = nrow(g),
      waic_detection_rate = mean(g$favors_M1, na.rm = TRUE),
      waic_detection_mcse = mcse_prop(g$favors_M1),
      correct_sel_rate = mean(g$correct_selection, na.rm = TRUE),
      false_positive_rate = mean(g$false_positive, na.rm = TRUE),
      incr_excess_mean = mean(g$incr_excess, na.rm = TRUE),
      incr_excess_mcse = mcse_mean(g$incr_excess),
      incr_detection_rate = mean(g$incr_p < 0.05, na.rm = TRUE),
      cor_excess_mean = mean(g$cor_excess, na.rm = TRUE),
      conc_excess_mean = mean(g$conc_excess, na.rm = TRUE),
      tracking_score_mean = mean(g$tracking_score, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  agg[order(agg$true_model == "M0", agg$n), ]
}

## ---- Scenario 4: misspecification robustness ----------------------------
agg_s4 <- function(reps) {
  rows <- lapply(reps, function(x) {
    conc <- x$concordance; incr <- x$increment
    rt <- x$ident   # identifiability metrics (rho_theta retired)
    data.frame(
      config_id  = x$config_id,
      sub        = x$config$sub,
      cor_excess = if (!is.null(conc)) conc$cor_excess else NA,
      conc_excess = if (!is.null(conc)) conc$excess     else NA,
      incr_excess = if (!is.null(incr)) incr$excess     else NA,
      Kplus_median = if (!is.null(x$kplus_median)) x$kplus_median else x$basic$K_plus_median,
      thetabar_cor  = if (!is.null(rt)) rt$thetabar_cor_mean else NA,
      omega0_rmse   = if (!is.null(rt)) rt$omega0_rmse else NA,
      correct_selection = as.integer(x$correct_selection),
      delta_waic  = x$delta_waic,
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)
  agg <- do.call(rbind, lapply(split(df, df$config_id), function(g) {
    data.frame(
      config_id    = g$config_id[1],
      sub          = g$sub[1],
      n_reps       = nrow(g),
      cor_excess_mean  = mean(g$cor_excess,  na.rm = TRUE),
      cor_excess_mcse  = mcse_mean(g$cor_excess),
      conc_excess_mean = mean(g$conc_excess, na.rm = TRUE),
      incr_excess_mean = mean(g$incr_excess, na.rm = TRUE),
      Kplus_median_mean = mean(g$Kplus_median, na.rm = TRUE),
      thetabar_cor_mean = mean(g$thetabar_cor, na.rm = TRUE),
      omega0_rmse_mean  = mean(g$omega0_rmse,  na.rm = TRUE),
      correct_sel_rate  = mean(g$correct_selection, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  agg[order(agg$config_id), ]
}

## ---- Scenario 5: identifiability of the hierarchical signature and the
##      time-varying background (varsigma_theta sweep). Replaces the retired
##      compartment-specific-signature (delta_theta / rho^theta) scenario. ---
agg_s5 <- function(reps) {
  rows <- lapply(reps, function(x) {
    rt <- x$ident   # identifiability metrics (rho_theta retired)
    data.frame(
      config_id    = x$config_id,
      sub          = x$config$sub,
      varsigma_theta_config = if (!is.null(x$config$varsigma_theta) && length(x$config$varsigma_theta)) x$config$varsigma_theta else if (length(x$varsigma_theta)) x$varsigma_theta else NA_real_,
      varsigma_theta_rep    = if (length(x$varsigma_theta)) x$varsigma_theta else NA_real_,
      kplus_median = if (length(x$kplus_median)) x$kplus_median else if (length(x$basic$K_plus_median)) x$basic$K_plus_median else NA_real_,
      thetabar_cor      = if (!is.null(rt)) rt$thetabar_cor_mean else NA,
      thetabar_rmse     = if (!is.null(rt)) rt$thetabar_rmse else NA,
      theta_subj_rmse   = if (!is.null(rt)) rt$theta_subj_rmse else NA,
      subj_vs_pop_ratio = if (!is.null(rt)) rt$theta_subj_vs_pop_ratio else NA,
      varsigma_k_est    = if (!is.null(rt)) rt$varsigma_k_est else NA,
      varsigma_k_bias   = if (!is.null(rt)) rt$varsigma_k_bias else NA,
      omega0_rmse       = if (!is.null(rt)) rt$omega0_rmse else NA,
      omega0_cor        = if (!is.null(rt)) rt$omega0_cor else NA,
      omega0_frac_bdry  = if (!is.null(rt)) rt$omega0_frac_boundary else NA,
      theta0_rmse       = if (!is.null(rt)) rt$theta0_rmse else NA,
      pi_l1_ident       = if (!is.null(rt)) rt$pi_l1 else NA,
      label_switch_rate = if (!is.null(rt)) rt$label_switch_rate else NA,
      A1_satisfied      = if (!is.null(rt)) rt$A1_satisfied else NA,
      favors_M1    = as.integer(x$favors_M1),
      correct_selection = as.integer(x$correct_selection),
      delta_waic   = x$delta_waic,
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)
  agg <- do.call(rbind, lapply(split(df, df$config_id), function(g) {
    data.frame(
      config_id    = g$config_id[1],
      sub          = g$sub[1],
      varsigma_theta = mean(g$varsigma_theta_rep, na.rm = TRUE),
      n_reps       = nrow(g),
      thetabar_cor_mean    = mean(g$thetabar_cor,    na.rm = TRUE),
      thetabar_cor_mcse    = mcse_mean(g$thetabar_cor),
      thetabar_rmse_mean   = mean(g$thetabar_rmse,   na.rm = TRUE),
      theta_subj_rmse_mean = mean(g$theta_subj_rmse, na.rm = TRUE),
      subj_vs_pop_ratio    = mean(g$subj_vs_pop_ratio, na.rm = TRUE),
      varsigma_k_est_mean  = mean(g$varsigma_k_est,  na.rm = TRUE),
      varsigma_k_bias_mean = mean(g$varsigma_k_bias, na.rm = TRUE),
      varsigma_k_bias_mcse = mcse_mean(g$varsigma_k_bias),
      omega0_rmse_mean     = mean(g$omega0_rmse,     na.rm = TRUE),
      omega0_cor_mean      = mean(g$omega0_cor,      na.rm = TRUE),
      omega0_frac_bdry_mean= mean(g$omega0_frac_bdry, na.rm = TRUE),
      theta0_rmse_mean     = mean(g$theta0_rmse,     na.rm = TRUE),
      ## THE critical readout: weight-level accuracy should hold even where the
      ## signatures are not identified (large varsigma_theta).
      pi_l1_ident_mean     = mean(g$pi_l1_ident,     na.rm = TRUE),
      pi_l1_ident_mcse     = mcse_mean(g$pi_l1_ident),
      label_switch_rate    = mean(g$label_switch_rate, na.rm = TRUE),
      A1_satisfied_frac    = mean(as.logical(g$A1_satisfied), na.rm = TRUE),
      kplus_median_mean    = mean(g$kplus_median,  na.rm = TRUE),
      detection_rate       = mean(g$favors_M1,     na.rm = TRUE),
      correct_sel_rate     = mean(g$correct_selection, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  agg[order(agg$varsigma_theta), ]
}

## ---- Run aggregation for all 5 scenarios --------------------------------
cat("Aggregating scenarios...\n")
for (s in 1:5) {
  cat(sprintf("  Scenario %d ... ", s))
  reps <- load_reps(s)
  cat(sprintf("%d reps loaded\n", length(reps)))
  if (length(reps) == 0) next
  agg <- switch(as.character(s),
    "1" = agg_s1(reps),
    "2" = agg_s2(reps),
    "3" = agg_s3(reps),
    "4" = agg_s4(reps),
    "5" = agg_s5(reps)
  )
  out_f <- file.path(out_dir, sprintf("scenario_%d_summary.csv", s))
  write.csv(agg, out_f, row.names = FALSE)
  cat(sprintf("    -> wrote %s (%d rows)\n", out_f, nrow(agg)))
}
cat("Done.\n")
