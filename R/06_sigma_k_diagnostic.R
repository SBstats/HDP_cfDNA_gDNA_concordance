###############################################################################
## 25_per_component_variance_diagnostic.R
##
## Round 3 diagnostic (Priority 1b): test whether the residual 5-9k log-unit
## chain-to-chain log-likelihood spread is explained by per-component sigma_k^2
## "locking" -- i.e. different chains settling on different component variances
## (especially for small/infrequently-allocated components) because sigma_k^2
## is updated independently per k with no pooling across components.
##
## Mechanism: the data log-likelihood contribution from a gene allocated to
## component k is proportional to  -0.5*log(sigma_k^2) - resid^2/(2*sigma_k^2).
## Summed over ~2*p*N_obs ~ 2e6 gene allocations, chain-level differences in
## sigma_k^2 accumulate into very large log-likelihood differences.
##
## Estimand: for each component k, compute
##   (A) split-R-hat and bulk-ESS of sigma_k across chains
##   (B) per-chain mean sigma_k[k] (and log sigma_k^2)
##   (C) approximate allocation mass N_k per chain (from omega traces)
##   (D) implied per-chain LL offset: -0.5 * N_k * log(sigma_k^2)
##       (this is the -0.5 log sigma_k^2 piece; the residual piece is
##        self-cancelling in Gibbs at stationarity.)
## Compare max(D) - min(D) across chains to the observed LL spread.
##
## Inputs:
##   results/krd_M1_fit.rds, results/krd_M0_fit.rds, results/krd_chain_health.csv
##
## Outputs:
##   results/krd_sigma_k_by_chain_M1.csv, krd_sigma_k_by_chain_M0.csv
##   results/krd_sigma_k_locking_verdict.csv
##   figures/krd_sigma_k_boxplot_by_chain.pdf
###############################################################################

suppressPackageStartupMessages({
  library(posterior)
  library(ggplot2)
})

base_dir <- getwd()
results_dir <- file.path(base_dir, "results")
figures_dir <- file.path(base_dir, "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== sigma_k locking diagnostic ===\n\n")

chain_health <- read.csv(file.path(results_dir, "krd_chain_health.csv"),
                         stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# Helper: per-model analysis
# ---------------------------------------------------------------------------
analyze_model <- function(multi_fit, model_label) {
  chains <- multi_fit$chains
  n_ch <- length(chains)
  K <- chains[[1]]$K
  n_save <- chains[[1]]$n_save
  p <- chains[[1]]$p
  N_obs <- chains[[1]]$N_obs
  K_trace <- chains[[1]]$K_trace

  cat(sprintf("\n--- %s: %d chains, K=%d, %d saved iters, p=%d, N_obs=%d ---\n",
              model_label, n_ch, K, n_save, p, N_obs))

  # Per-component sigma_k trace: array [n_save, n_ch, K]
  sigma_k_arr <- array(NA, dim = c(n_save, n_ch, K))
  for (ci in seq_len(n_ch)) {
    sigma_k_arr[, ci, ] <- chains[[ci]]$samples$sigma_k
  }

  # Per-chain posterior-mean allocation mass per component.
  # omega_g_trace is stored as [n_save x (N_obs * K_trace)]: column
  # (k-1)*N_obs + obs holds omega_g[k, obs].
  N_k_per_chain <- matrix(0, n_ch, K)
  for (ci in seq_len(n_ch)) {
    og <- chains[[ci]]$samples$omega_g_trace
    oc <- chains[[ci]]$samples$omega_cf_trace
    for (k in seq_len(min(K, K_trace))) {
      col_start <- (k - 1) * N_obs + 1
      col_end   <- k * N_obs
      # Expected allocations (gDNA + cfDNA) to component k, aggregated
      # across iterations and observations. gDNA contributes p genes per
      # observation; cfDNA contributes p*(1 - wN_i) genes per observation.
      # We approximate by ignoring wN (negligible in this dataset, < 0.1%).
      N_k_per_chain[ci, k] <-
        p * (mean(rowMeans(og[, col_start:col_end, drop = FALSE])) +
             mean(rowMeans(oc[, col_start:col_end, drop = FALSE]))) * N_obs
    }
  }

  # Per-component split-Rhat and ESS across chains
  rhat_vec <- numeric(K)
  ess_vec  <- numeric(K)
  for (k in seq_len(K)) {
    draws <- sigma_k_arr[, , k, drop = TRUE]  # [n_save x n_ch]
    rhat_vec[k] <- posterior::rhat(draws)
    ess_vec[k]  <- posterior::ess_bulk(draws)
  }

  # Per-chain mean sigma_k[k] and implied LL offset
  # sigma_k stored by sampler as sqrt(sigma_k2), so sigma_k2 = sigma_k^2.
  mean_sigma_k  <- apply(sigma_k_arr, c(2, 3), mean)  # [n_ch x K]
  log_sigma_k2  <- 2 * log(mean_sigma_k)              # [n_ch x K]
  # Implied LL offset per chain, per k: -0.5 * N_k * log(sigma_k^2)
  ll_offset_k   <- -0.5 * N_k_per_chain * log_sigma_k2

  # Per-chain total LL offset (summed over k)
  ll_offset_total <- rowSums(ll_offset_k)
  ll_offset_gap   <- diff(range(ll_offset_total))

  # Build per-(chain, k) table
  out_rows <- data.frame(
    model     = model_label,
    chain     = rep(seq_len(n_ch), each = K),
    component = rep(seq_len(K), times = n_ch),
    mean_sigma_k = as.vector(t(mean_sigma_k)),
    N_k_alloc    = as.vector(t(N_k_per_chain)),
    ll_offset    = as.vector(t(ll_offset_k)),
    stringsAsFactors = FALSE
  )

  # Per-component summary
  summary_k <- data.frame(
    model     = model_label,
    component = seq_len(K),
    rhat      = rhat_vec,
    bulk_ess  = ess_vec,
    min_mean  = apply(mean_sigma_k, 2, min),
    max_mean  = apply(mean_sigma_k, 2, max),
    ratio_max_min = apply(mean_sigma_k, 2, function(x) max(x)/min(x)),
    N_k_avg   = colMeans(N_k_per_chain),
    stringsAsFactors = FALSE
  )

  cat("\nPer-component diagnostics:\n")
  print(summary_k, digits = 4, row.names = FALSE)

  cat(sprintf("\nPer-chain implied LL offset (sum over k):\n"))
  for (ci in seq_len(n_ch)) {
    cat(sprintf("  Chain %d: %.1f\n", ci, ll_offset_total[ci]))
  }
  cat(sprintf("Implied LL offset gap (max - min): %.1f\n", ll_offset_gap))

  # Observed LL spread from chain_health
  obs_gap <- diff(range(
    chain_health$mean_loglik[chain_health$model == model_label]))
  cat(sprintf("Observed LL spread (from chain_health): %.1f\n", obs_gap))
  frac_explained <- ll_offset_gap / obs_gap
  cat(sprintf("Fraction of observed LL gap explained by sigma_k offset: %.2f\n",
              frac_explained))

  list(
    per_k   = out_rows,
    summary_k = summary_k,
    ll_offset_total = ll_offset_total,
    ll_offset_gap   = ll_offset_gap,
    observed_gap    = obs_gap,
    fraction_explained = frac_explained,
    sigma_k_arr = sigma_k_arr
  )
}

# ---------------------------------------------------------------------------
# Run for both models
# ---------------------------------------------------------------------------
multi_M1 <- readRDS(file.path(results_dir, "krd_M1_fit.rds"))
multi_M0 <- readRDS(file.path(results_dir, "krd_M0_fit.rds"))

res_M1 <- analyze_model(multi_M1, "M1")
res_M0 <- analyze_model(multi_M0, "M0")

write.csv(res_M1$per_k, file.path(results_dir, "krd_sigma_k_by_chain_M1.csv"),
          row.names = FALSE)
write.csv(res_M0$per_k, file.path(results_dir, "krd_sigma_k_by_chain_M0.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# Verdict CSV
# ---------------------------------------------------------------------------
verdict_row <- function(res, label) {
  any_rhat_high <- any(res$summary_k$rhat > 1.10, na.rm = TRUE)
  # Locking confirmed if: any per-component Rhat > 1.10 AND
  # the sigma_k offset explains >= 50% of the observed LL gap
  locking_confirmed <- any_rhat_high && res$fraction_explained >= 0.5
  data.frame(
    model = label,
    max_rhat_sigma_k = max(res$summary_k$rhat, na.rm = TRUE),
    n_components_rhat_gt_1p10 = sum(res$summary_k$rhat > 1.10, na.rm = TRUE),
    max_ratio_chain_sigma_k = max(res$summary_k$ratio_max_min, na.rm = TRUE),
    observed_ll_gap = res$observed_gap,
    implied_ll_gap_from_sigma_k = res$ll_offset_gap,
    fraction_explained = res$fraction_explained,
    locking_confirmed = locking_confirmed,
    stringsAsFactors = FALSE
  )
}

verdict_df <- rbind(verdict_row(res_M1, "M1"),
                    verdict_row(res_M0, "M0"))
write.csv(verdict_df, file.path(results_dir, "krd_sigma_k_locking_verdict.csv"),
          row.names = FALSE)

cat("\n=== VERDICT ===\n")
print(verdict_df, digits = 4, row.names = FALSE)

# ---------------------------------------------------------------------------
# Boxplot: sigma_k per chain per component
# ---------------------------------------------------------------------------
make_boxplot_df <- function(res, label) {
  K <- dim(res$sigma_k_arr)[3]
  n_ch <- dim(res$sigma_k_arr)[2]
  n_save <- dim(res$sigma_k_arr)[1]
  df <- data.frame(
    model = label,
    chain = factor(rep(seq_len(n_ch), each = n_save * K)),
    component = factor(rep(rep(seq_len(K), each = n_save), n_ch)),
    sigma_k = as.vector(res$sigma_k_arr)
  )
  df
}

box_df <- rbind(make_boxplot_df(res_M1, "M1"),
                make_boxplot_df(res_M0, "M0"))

p <- ggplot(box_df, aes(x = component, y = sigma_k, fill = chain)) +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.25) +
  facet_wrap(~ model, ncol = 1, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Component k", y = expression(sigma[k]),
       title = "Per-component sigma_k posterior by chain",
       fill = "Chain") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")

ggsave(file.path(figures_dir, "krd_sigma_k_boxplot_by_chain.pdf"),
       p, width = 8, height = 7)

cat("\nFiles written:\n")
cat("  results/krd_sigma_k_by_chain_M1.csv\n")
cat("  results/krd_sigma_k_by_chain_M0.csv\n")
cat("  results/krd_sigma_k_locking_verdict.csv\n")
cat("  figures/krd_sigma_k_boxplot_by_chain.pdf\n")
