###############################################################################
## R/02_preprocess.R
## Preprocess KRd trial 5-hMC data for the HDP concordance model
##
## Gene selection: ~318 candidate genes identified by pooled limma-voom DE
## analysis (design: ~ timepoint + source, blocking on patient_id via
## duplicateCorrelation) across all paired observations jointly.
## Threshold: adj.P.Val < 0.05 & |logFC| > 0.5 (gDNA vs cfDNA effect).
## This pooled approach maximises statistical power by sharing information
## across cycles and yields a stable, reproducible gene set.
##
## Inputs (all within this project root):
##   KRd trial/5hmC data/kRd-cfDNA_genebody_count.RDS
##   KRd trial/5hmC data/kRd-gDNA_genebody_count.RDS
##   KRd trial/5hmC data/kRd-cfDNA_sample_key.csv
##   KRd trial/5hmC data/kRd-gDNA_sample_key.csv
##
## Outputs:
##   results/candidate_genes_union.rds   DE gene names (character vector)
##   results/candidate_genes_union.csv   DE gene table with logFC / adj.P.Val
##   results/model_data.rds              model input (see structure below)
##
## model_data.rds fields:
##   Y_g     [p x N_obs]  gDNA  joint VST
##   Y_cf    [p x N_obs]  cfDNA joint VST
##   patient [N_obs]      integer patient ID (1..n)
##   time    [N_obs]      integer timepoint  (1..T_max)
##   n                    unique patients
##   p                    number of DE genes
##   T_max                number of timepoints
##   N_obs                total paired observations (after deduplication)
##   genes   character    gene names
##   paired_info          data.frame with Study.ID, Timepoint, barcodes, etc.
##
## Run from project root:
##   Rscript R/02_preprocess.R
###############################################################################

set.seed(20260328)

suppressPackageStartupMessages({
  library(dplyr)
  library(DESeq2)
  library(limma)
})

# ---------------------------------------------------------------------------
# Paths — all relative to project root
# ---------------------------------------------------------------------------
base_dir <- getwd()
data_dir <- file.path(base_dir, "KRd trial", "5hmC data")
out_dir  <- file.path(base_dir, "results")
dir.create(out_dir, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. Load sample keys
# ---------------------------------------------------------------------------
cfDNA_key <- read.csv(file.path(data_dir, "kRd-cfDNA_sample_key.csv"),
                      stringsAsFactors = FALSE) %>%
  mutate(barcode = Assigned.ID)

gDNA_key <- read.csv(file.path(data_dir, "kRd-gDNA_sample_key.csv"),
                     stringsAsFactors = FALSE) %>%
  mutate(barcode = Assigned.ID)

cat("cfDNA samples:", nrow(cfDNA_key), "\n")
cat("gDNA samples: ", nrow(gDNA_key),  "\n")

# ---------------------------------------------------------------------------
# 2. Identify paired observations (same patient + timepoint)
# ---------------------------------------------------------------------------
paired <- inner_join(
  cfDNA_key %>% select(Study.ID, Timepoint, barcode_cf = barcode),
  gDNA_key  %>% select(Study.ID, Timepoint, barcode_g  = barcode),
  by = c("Study.ID", "Timepoint")
)

cat("\nRaw paired observations:", nrow(paired), "\n")
cat("Unique patients:", n_distinct(paired$Study.ID), "\n")
print(table(paired$Timepoint))

# ---------------------------------------------------------------------------
# 3. Load count matrices and resolve duplicate gDNA entries at 3YR F/U
#    Patients 101-17 and 101-73 each have 2 gDNA barcodes at 3YR F/U.
#    Average their raw count vectors into a synthetic barcode, keep one row.
# ---------------------------------------------------------------------------
cat("\nResolving duplicate gDNA entries...\n")
cfDNA_counts_full <- as.matrix(readRDS(
  file.path(data_dir, "kRd-cfDNA_genebody_count.RDS")))
gDNA_counts_full  <- as.matrix(readRDS(
  file.path(data_dir, "kRd-gDNA_genebody_count.RDS")))

stopifnot(all(rownames(cfDNA_counts_full) == rownames(gDNA_counts_full)))
cat("Raw count matrix:", nrow(cfDNA_counts_full), "genes\n")

# Identify duplicated patient x timepoint rows
dup_key <- paste(paired$Study.ID, paired$Timepoint)
dup_rows <- which(duplicated(dup_key) | duplicated(dup_key, fromLast = TRUE))

if (length(dup_rows) > 0) {
  dup_groups <- unique(dup_key[dup_rows])
  rows_to_drop <- integer(0)
  for (grp in dup_groups) {
    rows <- which(dup_key == grp)
    g_bcs <- paired$barcode_g[rows]
    cat("  Duplicate:", grp, "| gDNA barcodes:", paste(g_bcs, collapse = ", "), "\n")
    avg_g <- rowMeans(gDNA_counts_full[, g_bcs, drop = FALSE])
    synth_bc <- paste0(g_bcs[1], "_avg")
    gDNA_counts_full <- cbind(gDNA_counts_full, avg_g)
    colnames(gDNA_counts_full)[ncol(gDNA_counts_full)] <- synth_bc
    paired$barcode_g[rows[1]] <- synth_bc
    rows_to_drop <- c(rows_to_drop, rows[-1])
  }
  paired <- paired[-rows_to_drop, ]
  cat("After deduplication:", nrow(paired), "paired observations\n")
} else {
  cat("No duplicates found.\n")
}

cat("Timepoint breakdown after deduplication:\n")
print(table(paired$Timepoint))

# ---------------------------------------------------------------------------
# 5. Select DE genes: pooled limma-voom (design: ~ timepoint + source,
#    blocking on patient_id). Threshold: adj.P.Val < 0.05 & |logFC| > 0.5.
# ---------------------------------------------------------------------------
cat("\nSelecting DE genes via pooled limma-voom...\n")

# Build pooled raw count matrix (all paired cfDNA + gDNA, all genes)
raw_cf_all <- cfDNA_counts_full[, paired$barcode_cf, drop = FALSE]
raw_g_all  <- gDNA_counts_full[,  paired$barcode_g,  drop = FALSE]
colnames(raw_cf_all) <- paste0("CF_", paired$barcode_cf)
colnames(raw_g_all)  <- paste0("G_",  paired$barcode_g)
pooled_all <- cbind(raw_cf_all, raw_g_all)

# Gene-level QC defining the universe for the DE test: keep genes with >= 10
# counts in >= 95% of the pooled cfDNA + gDNA samples.
#
# KRd DATA ONLY. This is deliberate and is a change from the historical
# pipeline, which is worth recording because it changes the gene count.
#
# The original selection (parent folder R/01_match_and_normalize.R ->
# R/04_DE_genes_and_concordance.R) tested an 11,430-gene universe and selected
# 310 genes. That 11,430 was NOT reachable from KRd data alone: it was the
# intersection of the KRd QC set with a separate breast-cancer case-control
# cohort, after per-source VST had dropped all-zero genes in each. Depending on
# an unrelated cohort to fix this analysis's feature set is both scientifically
# unmotivated and fatal to this folder being self-contained, so that dependency
# has been removed.
#
# Consequence: the universe is 11,816 genes and the DE set is 318, not 310.
# The two sets overlap in 308 genes (10 gained, 2 lost) -- 318 is NOT a
# superset of 310, because the universe size feeds the BH correction and hence
# the adjusted p-values. Any manuscript text or downstream artifact quoting
# p = 310 must be updated to p = 318.
n_all_samples <- ncol(pooled_all)
keep_genes    <- rownames(pooled_all)[rowSums(pooled_all >= 10) >= 0.95 * n_all_samples]
pooled_filt   <- pooled_all[keep_genes, , drop = FALSE]
cat("Genes passing QC filter (KRd only, >=10 counts in >=95% of samples):",
    length(keep_genes), "\n")

# Joint DESeq2 VST on the QC-filtered pooled matrix (needed for voom weights)
mat_qc <- round(pooled_filt)
storage.mode(mat_qc) <- "integer"
mat_qc <- mat_qc[rowSums(mat_qc) > 0, , drop = FALSE]
colnames(mat_qc) <- make.unique(colnames(mat_qc))
coldata_qc <- data.frame(row.names = colnames(mat_qc),
                         condition = factor(rep("x", ncol(mat_qc))))
dds_qc   <- DESeqDataSetFromMatrix(mat_qc, colData = coldata_qc, design = ~1)
vst_qc   <- assay(varianceStabilizingTransformation(dds_qc, blind = TRUE))
genes_qc <- rownames(vst_qc)
cat("Genes after joint VST:", length(genes_qc), "\n")

# Sample metadata for limma design
n_obs      <- nrow(paired)
timepoint_order_de <- c("Screening" = 1, "C4" = 2, "C8" = 3, "C18" = 4, "3 YR F/U" = 5)
sample_df <- data.frame(
  col_name   = colnames(pooled_all),
  source     = factor(c(rep("cfDNA", n_obs), rep("gDNA", n_obs)),
                      levels = c("cfDNA", "gDNA")),
  patient_id = c(paired$Study.ID, paired$Study.ID),
  timepoint  = factor(c(as.character(paired$Timepoint),
                         as.character(paired$Timepoint)),
                      levels = names(timepoint_order_de)),
  stringsAsFactors = FALSE
)

# Pooled limma-voom: ~ timepoint + source, blocking on patient_id
pooled_de   <- pooled_all[genes_qc, , drop = FALSE]
design_all  <- model.matrix(~ timepoint + source, data = sample_df)
v_all   <- voom(pooled_de, design_all, plot = FALSE)
corf1   <- duplicateCorrelation(v_all, design_all, block = sample_df$patient_id)
v_all2  <- voom(pooled_de, design_all, plot = FALSE,
                block = sample_df$patient_id, correlation = corf1$consensus)
corf2   <- duplicateCorrelation(v_all2, design_all, block = sample_df$patient_id)
fit_all <- lmFit(v_all2, design_all,
                 block = sample_df$patient_id, correlation = corf2$consensus)
fit_all <- eBayes(fit_all)
tt_all  <- topTable(fit_all, coef = "sourcegDNA", number = Inf,
                    sort.by = "none", adjust.method = "BH")
tt_all$gene <- rownames(tt_all)
tt_all$DE   <- tt_all$adj.P.Val < 0.05 & abs(tt_all$logFC) > 0.5

n_de <- sum(tt_all$DE, na.rm = TRUE)
cat("Pooled DE genes (adj.P<0.05 & |logFC|>0.5):", n_de, "\n")

de_genes <- tt_all$gene[tt_all$DE & !is.na(tt_all$DE)]
cat("Candidate gene set:", length(de_genes), "genes\n")

# Save gene list and table for reference
saveRDS(de_genes, file.path(out_dir, "candidate_genes_union.rds"))
de_table <- data.frame(
  gene      = de_genes,
  logFC     = tt_all$logFC[tt_all$DE & !is.na(tt_all$DE)],
  adj.P.Val = tt_all$adj.P.Val[tt_all$DE & !is.na(tt_all$DE)],
  stringsAsFactors = FALSE
) %>% arrange(desc(abs(logFC)))
write.csv(de_table, file.path(out_dir, "candidate_genes_union.csv"), row.names = FALSE)
cat("Saved: results/candidate_genes_union.rds and .csv\n")

# ---------------------------------------------------------------------------
# 6. Subset count matrices to paired samples and DE genes
# ---------------------------------------------------------------------------
cat("\nSubsetting to DE genes and paired samples...\n")

missing_de <- setdiff(de_genes, rownames(cfDNA_counts_full))
if (length(missing_de) > 0) {
  warning("DE genes missing from count matrix: ", paste(missing_de, collapse = ", "))
  de_genes <- intersect(de_genes, rownames(cfDNA_counts_full))
}

cfDNA_paired <- cfDNA_counts_full[de_genes, paired$barcode_cf, drop = FALSE]
gDNA_paired  <- gDNA_counts_full[ de_genes, paired$barcode_g,  drop = FALSE]
cat("cfDNA subset:", nrow(cfDNA_paired), "genes x", ncol(cfDNA_paired), "samples\n")
cat("gDNA  subset:", nrow(gDNA_paired),  "genes x", ncol(gDNA_paired),  "samples\n")

# ---------------------------------------------------------------------------
# 7. Pool cfDNA + gDNA for joint DESeq2 VST (common scale for both sources)
# ---------------------------------------------------------------------------
cat("\nBuilding pooled count matrix for joint VST...\n")
pooled <- cbind(
  cfDNA_paired,
  gDNA_paired
)
colnames(pooled) <- c(
  paste0("cf_", paired$barcode_cf),
  paste0("gd_", paired$barcode_g)
)
cat("Pooled matrix:", nrow(pooled), "genes x", ncol(pooled), "samples\n")

# Confirm all genes pass the count filter at this reduced gene set
# (DE genes were already identified from QC-filtered data, so this is expected)
n_total_samples <- ncol(pooled)
low_count_per_gene <- rowSums(pooled < 10)
threshold_samples <- ceiling(n_total_samples * 0.05)
genes_pass_filter <- rownames(pooled)[low_count_per_gene <= threshold_samples]
genes_drop <- setdiff(rownames(pooled), genes_pass_filter)
if (length(genes_drop) > 0) {
  cat("Dropping", length(genes_drop), "DE genes that fail count filter:",
      paste(genes_drop, collapse = ", "), "\n")
  pooled <- pooled[genes_pass_filter, , drop = FALSE]
} else {
  cat("All", nrow(pooled), "DE genes pass count filter (>=10 counts in >=95% of samples).\n")
}

# ---------------------------------------------------------------------------
# 8. DESeq2 joint VST
# ---------------------------------------------------------------------------
cat("\nApplying DESeq2 joint VST...\n")
mat_int <- round(pooled)
storage.mode(mat_int) <- "integer"
mat_int <- mat_int[rowSums(mat_int) > 0, , drop = FALSE]
colnames(mat_int) <- make.unique(colnames(mat_int))

coldata <- data.frame(
  row.names = colnames(mat_int),
  condition = factor(rep("x", ncol(mat_int)))
)
dds <- DESeqDataSetFromMatrix(mat_int, colData = coldata, design = ~1)
vst_mat <- assay(varianceStabilizingTransformation(dds, blind = TRUE))

cat("VST complete:", nrow(vst_mat), "genes x", ncol(vst_mat), "samples\n")
cat("Value range: [", round(min(vst_mat), 2), ",", round(max(vst_mat), 2), "]\n")

# ---------------------------------------------------------------------------
# 9. Split back into cfDNA and gDNA
# ---------------------------------------------------------------------------
n_paired <- nrow(paired)
# make.unique in the VST step may append suffixes; use positional split
cfDNA_final <- vst_mat[, seq_len(n_paired), drop = FALSE]
gDNA_final  <- vst_mat[, seq_len(n_paired) + n_paired, drop = FALSE]
colnames(cfDNA_final) <- paired$barcode_cf
colnames(gDNA_final)  <- paired$barcode_g

p <- nrow(cfDNA_final)
cat("\nFinal cfDNA:", p, "genes x", ncol(cfDNA_final), "samples\n")
cat("Final gDNA: ", p, "genes x", ncol(gDNA_final),  "samples\n")

# ---------------------------------------------------------------------------
# 10. Build model data structure
# ---------------------------------------------------------------------------
timepoint_order <- c("Screening" = 1, "C4" = 2, "C8" = 3,
                     "C18" = 4, "3 YR F/U" = 5)

paired <- paired %>%
  mutate(
    t_num      = timepoint_order[Timepoint],
    patient_id = as.integer(factor(Study.ID))
  )

n     <- n_distinct(paired$patient_id)
T_max <- max(paired$t_num)

model_data <- list(
  Y_g         = gDNA_final,
  Y_cf        = cfDNA_final,
  patient     = paired$patient_id,
  time        = paired$t_num,
  n           = n,
  p           = p,
  T_max       = T_max,
  N_obs       = ncol(cfDNA_final),
  genes       = rownames(cfDNA_final),
  paired_info = paired
)

cat("\n=== Final Model Data ===\n")
cat("Patients (n):              ", n, "\n")
cat("Genes (p):                 ", p, "(DE genes only)\n")
cat("Max timepoints (T_max):    ", T_max, "\n")
cat("Total paired observations: ", model_data$N_obs, "\n")

# ---------------------------------------------------------------------------
# 11. Save
# ---------------------------------------------------------------------------
saveRDS(model_data, file.path(out_dir, "model_data.rds"))
cat("\nSaved: results/model_data.rds\n")

sink(file.path(out_dir, "sessionInfo_preprocess.txt"))
cat("Preprocessing completed:", format(Sys.time()), "\n\n")
cat("Gene set: pooled limma-voom DE genes (p =", p, ")\n")
cat("Selection: adj.P.Val < 0.05 & |logFC| > 0.5, design ~ timepoint + source\n\n")
sessionInfo()
sink()

cat("Preprocessing complete.\n")
