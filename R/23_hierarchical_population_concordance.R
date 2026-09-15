###############################################################################
## R/23_hierarchical_population_concordance.R
##
## The "one number" endpoint of the concordance hierarchy (storyline item 7):
## pool the PER-PATIENT longitudinal concordance correlations into a single
## POPULATION-LEVEL estimate with a confidence interval, via a Fisher-z random-
## effects meta-analysis (DerSimonian-Laird), instead of treating each patient's
## correlation as an independent replicate or reporting a single pooled number
## that ignores between-patient heterogeneity.
##
## For each patient i we form a change-based concordance rho_i:
##   (A) MODEL-FREE gene-level:   rho_i = mean over i's within-patient transitions
##       of cor(DeltaYcf, DeltaYg) across the 11,837 genes.
##   (B) LATENT M1 / (C) LATENT M0: rho_i = cor of the dominant-subclone weight
##       increments Dwcf vs Dwg across i's transitions (posterior-mean weights;
##       gDNA-dominant reference so M0 is a valid cross-source comparison).
##
## Then z_i = atanh(rho_i), pooled as z_i ~ N(mu, tau^2 + v_i) with DL tau^2 and
## v_i approximated from each patient's information (n_i transitions; genes for A).
## Report tanh(mu) as the population longitudinal concordance + 95% CI, plus a
## per-patient forest.
##
## Runs on existing fits. Outputs:
##   results/population_concordance.csv
##   results/population_concordance_per_patient.csv
##   figures/fig_population_concordance.pdf  (+ nature_manuscript/figures/)
###############################################################################

set.seed(20260808)
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
base_dir <- getwd()  # run from project root
setwd(base_dir)
suppressWarnings(suppressMessages(source(file.path(base_dir, "R", "01_lib_core.R"))))
res_dir <- file.path(base_dir, "results")
fig_dirs <- c(file.path(base_dir, "figures"), file.path(base_dir, "nature_manuscript", "figures"))
for (d in fig_dirs) dir.create(d, showWarnings = FALSE, recursive = TRUE)

md  <- readRDS(file.path(res_dir, "model_data.rds"))
Yg  <- md$Y_g; Ycf <- md$Y_cf; N <- md$N_obs; tvec <- md$time; pvec <- md$patient

## within-patient transitions (strictly increasing time), grouped by patient
trans_by_patient <- lapply(sort(unique(pvec)), function(pt) {
  idx <- which(pvec == pt); idx <- idx[order(tvec[idx])]
  if (length(idx) < 2) return(NULL)
  do.call(rbind, lapply(1:(length(idx) - 1), function(j) {
    a <- idx[j]; b <- idx[j + 1]
    if (tvec[b] == tvec[a]) return(NULL)
    data.frame(a = a, b = b)
  }))
})
names(trans_by_patient) <- sort(unique(pvec))
trans_by_patient <- trans_by_patient[!sapply(trans_by_patient, is.null)]
patients <- as.integer(names(trans_by_patient))
n_i <- sapply(trans_by_patient, nrow)          # transitions per patient

## ------- (A) per-patient model-free gene-level change correlation -----------
rho_gene <- sapply(trans_by_patient, function(tb) {
  vals <- sapply(1:nrow(tb), function(j) cor(Yg[, tb$b[j]] - Yg[, tb$a[j]],
                                             Ycf[, tb$b[j]] - Ycf[, tb$a[j]]))
  mean(vals)
})
## variance of Fisher-z: gene-level r per transition has ~p obs, but transitions
## within a patient are correlated; use n_eff = n_i (independent transitions) so
## v_i reflects how many transitions the patient contributes (conservative).
v_gene <- 1 / pmax(n_i, 1)                       # heterogeneity dominated by between-patient

## ------- (B)/(C) per-patient latent dominant-weight change correlation ------
## VALID cross-source reference (as in R/15): within EACH posterior draw take
## gDNA's dominant signature and read cfDNA's loading on the SAME raw component
## (omega_g[k] and omega_cf[k] share theta_k in the stored traces), so M0's
## independent per-source relabeling cannot corrupt the comparison. Average each
## patient's Delta-weight correlation across draws.
latent_rho <- function(fit_path) {
  fit <- readRDS(fit_path)
  G  <- do.call(rbind, lapply(fit$chains, function(ch) ch$samples$omega_g_trace))   # S x (K*N)
  CF <- do.call(rbind, lapply(fit$chains, function(ch) ch$samples$omega_cf_trace))
  K  <- fit$chains[[1]]$K_trace; S <- nrow(G)
  # accumulate per-patient correlation over draws
  acc <- matrix(NA, S, length(trans_by_patient))
  for (s in 1:S) {
    Gm  <- matrix(G[s, ],  nrow = N, ncol = K)
    CFm <- matrix(CF[s, ], nrow = N, ncol = K)
    ref <- which.max(colMeans(Gm))
    wg <- Gm[, ref]; wcf <- CFm[, ref]
    for (pi in seq_along(trans_by_patient)) {
      tb <- trans_by_patient[[pi]]
      dwg  <- wg[tb$b]  - wg[tb$a]
      dwcf <- wcf[tb$b] - wcf[tb$a]
      if (length(dwg) >= 2 && sd(dwg) > 1e-9 && sd(dwcf) > 1e-9)
        acc[s, pi] <- cor(dwg, dwcf)
    }
  }
  colMeans(acc, na.rm = TRUE)   # posterior-mean per-patient Delta-weight correlation
}
rho_M1 <- latent_rho(file.path(res_dir, "krd_M1_fit.rds"))
rho_M0 <- latent_rho(file.path(res_dir, "krd_M0_fit.rds"))

## ------- DerSimonian-Laird random-effects pooling on Fisher-z --------------
pool_z <- function(rho, vi) {
  ok <- is.finite(rho) & rho > -0.999 & rho < 0.999 & is.finite(vi)
  rho <- rho[ok]; vi <- vi[ok]
  z <- atanh(rho); w <- 1 / vi
  mu_fe <- sum(w * z) / sum(w)
  Q <- sum(w * (z - mu_fe)^2); dfree <- length(z) - 1
  C <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0, (Q - dfree) / C)
  wr <- 1 / (vi + tau2)
  mu <- sum(wr * z) / sum(wr); se <- sqrt(1 / sum(wr))
  ci <- mu + c(-1, 1) * 1.96 * se
  list(k = length(z), rho_pool = tanh(mu), lo = tanh(ci[1]), hi = tanh(ci[2]),
       mu = mu, se = se, tau2 = tau2, Q = Q, p = 2 * pnorm(-abs(mu / se)),
       I2 = max(0, (Q - dfree) / Q) * 100)
}
res <- rbind(
  data.frame(estimator = "Model-free gene-level change", pool_z(rho_gene, v_gene)),
  data.frame(estimator = "Latent M1 (coupled) dom-weight change", pool_z(rho_M1, 1 / pmax(n_i, 1))),
  data.frame(estimator = "Latent M0 (independent) dom-weight change", pool_z(rho_M0, 1 / pmax(n_i, 1)))
)
res[ , c("rho_pool","lo","hi","mu","se","tau2","Q","p","I2")] <-
  round(res[ , c("rho_pool","lo","hi","mu","se","tau2","Q","p","I2")], 4)
write.csv(res, file.path(res_dir, "population_concordance.csv"), row.names = FALSE)
cat("\n=== Population-level longitudinal concordance (Fisher-z random effects) ===\n")
print(res[ , c("estimator","k","rho_pool","lo","hi","p","I2")], row.names = FALSE)

per_pt <- data.frame(patient = patients, n_transitions = n_i,
                     rho_gene = round(rho_gene, 4),
                     rho_M1 = round(rho_M1, 4), rho_M0 = round(rho_M0, 4))
write.csv(per_pt, file.path(res_dir, "population_concordance_per_patient.csv"), row.names = FALSE)

## ------- forest plot: per-patient model-free rho_i + pooled diamond --------
pl <- data.frame(patient = factor(patients, levels = patients[order(rho_gene)]),
                 rho = rho_gene, n = n_i)
pooled <- pool_z(rho_gene, v_gene)
p <- ggplot(pl, aes(rho, patient)) +
  geom_vline(xintercept = 0, color = "grey60", linetype = "dashed") +
  geom_vline(xintercept = pooled$rho_pool, color = "#E63946", linewidth = 0.5) +
  geom_point(aes(size = n), color = "#457B9D") +
  annotate("rect", xmin = pooled$lo, xmax = pooled$hi, ymin = -Inf, ymax = Inf,
           alpha = 0.10, fill = "#E63946") +
  scale_size_continuous(range = c(1, 3.5), name = "transitions") +
  labs(x = "Per-patient gene-level change concordance (rho_i)", y = "Patient",
       title = "Population-level longitudinal concordance",
       subtitle = sprintf("Pooled rho = %.3f (95%% CI %.3f-%.3f), random-effects Fisher-z; red band = 95%% CI",
                          pooled$rho_pool, pooled$lo, pooled$hi)) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 8),
        axis.text.y = element_text(size = 6))
for (d in fig_dirs) ggsave(file.path(d, "fig_population_concordance.pdf"), p, width = 7, height = 5.5)

cat(sprintf("\nModel-free pooled rho = %.3f (95%% CI %.3f-%.3f), I2 = %.0f%%\n",
            pooled$rho_pool, pooled$lo, pooled$hi, pooled$I2))
cat("Saved: results/population_concordance.csv, population_concordance_per_patient.csv, fig_population_concordance.pdf\n")
