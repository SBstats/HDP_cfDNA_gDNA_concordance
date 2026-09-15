###############################################################################
## R/17_normalization_robustness.R
## Robustness of the model-free longitudinal Δ-concordance finding to the
## normalization scheme. The paper normalizes cfDNA and gDNA JOINTLY (one pooled
## DESeq2 VST), which a reviewer may argue manufactures cross-source similarity.
## Here we re-derive the 5-hMC matrices with cfDNA and gDNA normalized
## SEPARATELY (two independent VST fits) on the SAME 11,837 genes and SAME
## pairing, and re-run the model-free within-patient increment test. No mixture
## model, no MCMC. If the pairing-specific longitudinal signal survives, the
## joint-normalization objection is refuted.
##
## Outputs: results/robustness/normalization_delta_concordance.csv
###############################################################################

set.seed(20260724)
suppressPackageStartupMessages({ library(DESeq2); library(dplyr) })
base_dir <- getwd()  # run from project root
setwd(base_dir)
res_dir <- file.path(base_dir, "results")
out_dir <- file.path(res_dir, "robustness"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
data_dir <- file.path(base_dir, "KRd trial", "5hmC data")

md    <- readRDS(file.path(res_dir, "model_data.rds"))
genes <- md$genes; pinfo <- md$paired_info; tvec <- md$time; pvec <- md$patient; N <- md$N_obs

## ---- raw counts, subset to the SAME genes and paired barcodes ----
raw_cf <- readRDS(file.path(data_dir, "kRd-cfDNA_genebody_count.RDS"))
raw_g  <- readRDS(file.path(data_dir, "kRd-gDNA_genebody_count.RDS"))
cf_counts <- round(as.matrix(raw_cf[genes, pinfo$barcode_cf]))
g_counts  <- round(as.matrix(raw_g [genes, pinfo$barcode_g ]))
stopifnot(ncol(cf_counts) == N, ncol(g_counts) == N, nrow(cf_counts) == length(genes))

## ---- SEPARATE per-source VST (each source normalized on its own) ----
vst_sep <- function(counts) {
  cd <- data.frame(row.names = colnames(counts), grp = factor(rep("a", ncol(counts))))
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = cd, design = ~ 1)
  assay(vst(dds, blind = FALSE))
}
Ycf_sep <- vst_sep(cf_counts); Yg_sep <- vst_sep(g_counts)
## joint-normalized matrices as used in the paper (from model_data.rds)
Yg_joint <- md$Y_g; Ycf_joint <- md$Y_cf

## ---- within-patient transitions (consecutive, strictly increasing time) ----
TR <- do.call(rbind, lapply(unique(pvec), function(pt){
  idx <- which(pvec == pt); idx <- idx[order(tvec[idx])]
  if (length(idx) < 2) return(NULL)
  do.call(rbind, lapply(1:(length(idx)-1), function(j){
    a <- idx[j]; b <- idx[j+1]; if (tvec[b] == tvec[a]) return(NULL)
    data.frame(a=a, b=b, st=paste0(tvec[a],"_",tvec[b])) }))
}))
M <- nrow(TR)

## ---- model-free level + increment concordance with stratified permutation null ----
delta_test <- function(Yg, Ycf, n_perm = 5000) {
  level <- mean(sapply(1:N, function(i) cor(Yg[, i], Ycf[, i])))         # cross-sectional r
  dG  <- sapply(1:M, function(m) Yg [, TR$b[m]] - Yg [, TR$a[m]])
  dCF <- sapply(1:M, function(m) Ycf[, TR$b[m]] - Ycf[, TR$a[m]])
  C   <- cor(dG, dCF); obs <- mean(diag(C))
  null <- replicate(n_perm, { q <- 1:M
    for (s in unique(TR$st)) { ii <- which(TR$st == s); if (length(ii) > 1) q[ii] <- sample(ii) }
    mean(C[cbind(1:M, q)]) })
  data.frame(level_r = round(level,4), delta_true = round(obs,4),
             delta_null = round(mean(null),4), frac_pos = round(mean(diag(C) > 0),3),
             p = (1 + sum(null >= obs))/(1 + n_perm),
             z = round((obs - mean(null))/sd(null), 2))
}

joint <- cbind(normalization = "joint (pooled VST, as in paper)", delta_test(Yg_joint, Ycf_joint))
sep   <- cbind(normalization = "separate (per-source VST)",       delta_test(Yg_sep,   Ycf_sep))
tab <- rbind(joint, sep)

cat("\n===== Normalization robustness: model-free within-patient Δ concordance =====\n")
cat(sprintf("genes=%d  paired obs=%d  within-patient intervals=%d\n\n", length(genes), N, M))
print(tab, row.names = FALSE)
cat("\nInterpretation: 'level_r' is the cross-sectional (non-specific) correlation; ",
    "'delta_true' vs 'delta_null' is the pairing-specific longitudinal signal.\n", sep = "")
write.csv(tab, file.path(out_dir, "normalization_delta_concordance.csv"), row.names = FALSE)
cat("Saved results/robustness/normalization_delta_concordance.csv\n")
