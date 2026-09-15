###############################################################################
## R/29d_cfdna_vs_gdna_separate_norm.R
##
## Comparison scatter: joint VST (existing) vs. source-separate VST (new)
##
## Generates a 6-page PDF with 3 panels per page:
##   Left:   Joint DESeq2 VST on pooled 168 samples (existing pipeline)
##   Middle: Separate DESeq2 VST — cfDNA normalized alone, gDNA normalized alone
##   Right:  Raw counts (no normalization, 19,100 genes)
##
## Output: KRD_cfdna_vs_gdna_scatter_separate_norm.pdf  (project root)
## Existing KRD_cfdna_vs_gdna_scatter.pdf is NOT touched.
## Run from project root: Rscript R/29d_cfdna_vs_gdna_separate_norm.R
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(DESeq2)
})

base_dir <- getwd()
res_dir  <- file.path(base_dir, "results")
data_dir <- file.path(base_dir, "KRd trial", "5hmC data")
out_pdf  <- file.path(base_dir, "KRD_cfdna_vs_gdna_scatter_separate_norm.pdf")

tp_levels <- c("Screening", "C4", "C8", "C18", "3 YR F/U")
tp_map    <- c("Screening" = 1L, "C4" = 2L, "C8" = 3L, "C18" = 4L, "3 YR F/U" = 5L)

## ---------------------------------------------------------------------------
## 1. Load joint-VST data (existing pipeline output)
## ---------------------------------------------------------------------------
cat("Loading joint-VST model_data.rds...\n")
md <- readRDS(file.path(res_dir, "model_data.rds"))

## ---------------------------------------------------------------------------
## 2. Load raw counts and build paired sample index
## ---------------------------------------------------------------------------
cat("Loading raw counts...\n")
cfDNA_key <- read.csv(file.path(data_dir, "kRd-cfDNA_sample_key.csv"), stringsAsFactors = FALSE)
gDNA_key  <- read.csv(file.path(data_dir, "kRd-gDNA_sample_key.csv"),  stringsAsFactors = FALSE)
cfDNA_key$barcode <- cfDNA_key$Assigned.ID
gDNA_key$barcode  <- gDNA_key$Assigned.ID

paired <- inner_join(
  cfDNA_key |> select(Study.ID, Timepoint, barcode_cf = barcode),
  gDNA_key  |> select(Study.ID, Timepoint, barcode_g  = barcode),
  by = c("Study.ID", "Timepoint")
)
paired$t_num <- tp_map[paired$Timepoint]

cfDNA_counts <- readRDS(file.path(data_dir, "kRd-cfDNA_genebody_count.RDS"))
gDNA_counts  <- readRDS(file.path(data_dir, "kRd-gDNA_genebody_count.RDS"))

cf_raw_paired <- as.matrix(cfDNA_counts[, paired$barcode_cf])  # [19100 x 84]
gd_raw_paired <- as.matrix(gDNA_counts[,  paired$barcode_g])

## ---------------------------------------------------------------------------
## 3. Apply the same gene filter as 02_preprocess.R on the pooled matrix,
##    then run SEPARATE VST for cfDNA and gDNA independently
## ---------------------------------------------------------------------------
cat("Applying joint gene filter (same 11,837 genes as existing pipeline)...\n")
pooled_raw <- cbind(cf_raw_paired, gd_raw_paired)
n_total    <- ncol(pooled_raw)
thresh     <- ceiling(n_total * 0.05)
keep       <- rowSums(pooled_raw < 10) <= thresh
cat(sprintf("  Genes retained: %d / %d\n", sum(keep), nrow(pooled_raw)))

cf_filt <- cf_raw_paired[keep, ]
gd_filt <- gd_raw_paired[keep, ]

cat("Running SEPARATE DESeq2 VST for cfDNA (84 samples)...\n")
dds_cf <- DESeqDataSetFromMatrix(
  countData = round(cf_filt),
  colData   = data.frame(sample = colnames(cf_filt), row.names = colnames(cf_filt)),
  design    = ~ 1
)
vst_cf <- assay(vst(dds_cf, blind = FALSE))  # [11837 x 84]

cat("Running SEPARATE DESeq2 VST for gDNA (84 samples)...\n")
dds_gd <- DESeqDataSetFromMatrix(
  countData = round(gd_filt),
  colData   = data.frame(sample = colnames(gd_filt), row.names = colnames(gd_filt)),
  design    = ~ 1
)
vst_gd <- assay(vst(dds_gd, blind = FALSE))  # [11837 x 84]

cat(sprintf("  Separate-VST cfDNA range: [%.2f, %.2f]\n", min(vst_cf), max(vst_cf)))
cat(sprintf("  Separate-VST gDNA range:  [%.2f, %.2f]\n", min(vst_gd), max(vst_gd)))
cat(sprintf("  Joint-VST   cfDNA range:  [%.2f, %.2f]\n", min(md$Y_cf), max(md$Y_cf)))
cat(sprintf("  Joint-VST   gDNA range:   [%.2f, %.2f]\n", min(md$Y_g),  max(md$Y_g)))

## ---------------------------------------------------------------------------
## 4. Shared axis limits per normalisation (computed once for comparability)
## ---------------------------------------------------------------------------
joint_lim <- range(c(md$Y_g, md$Y_cf))
sep_lim   <- range(c(vst_gd, vst_cf))
raw_lim   <- range(c(gd_raw_paired, cf_raw_paired))

## ---------------------------------------------------------------------------
## 5. Helper: build one 3-panel page
## ---------------------------------------------------------------------------
make_page <- function(obs_idx,        # integer vector: columns to use (paired index)
                      page_title, n_label) {

  ## --- joint VST ---
  g_j  <- as.vector(md$Y_g[,  obs_idx])
  cf_j <- as.vector(md$Y_cf[, obs_idx])
  r_j  <- cor(g_j, cf_j)

  ## --- separate VST ---
  g_s  <- as.vector(vst_gd[, obs_idx])
  cf_s <- as.vector(vst_cf[, obs_idx])
  r_s  <- cor(g_s, cf_s)

  ## --- raw counts (all 19,100 genes) ---
  g_r  <- as.vector(gd_raw_paired[, obs_idx])
  cf_r <- as.vector(cf_raw_paired[, obs_idx])
  r_r  <- cor(g_r, cf_r)

  n_pts_filt <- length(g_j)
  n_pts_raw  <- length(g_r)

  panel_joint <- ggplot(data.frame(gDNA = g_j, cfDNA = cf_j), aes(gDNA, cfDNA)) +
    geom_bin2d(bins = 80) +
    geom_abline(slope = 1, intercept = 0, colour = "firebrick",
                linewidth = 0.5, linetype = "dashed") +
    scale_fill_viridis_c(option = "magma", trans = "log10", name = "Count\n(log10)") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
             label = sprintf("r = %.4f\n%s gene×obs", r_j,
                             format(n_pts_filt, big.mark = ",")),
             size = 3) +
    coord_cartesian(xlim = joint_lim, ylim = joint_lim) +
    labs(title = "Joint VST (existing pipeline)",
         subtitle = sprintf("11,837 genes, %d obs, pooled 168-sample normalisation", n_label),
         x = "gDNA 5-hMC (joint VST)", y = "cfDNA 5-hMC (joint VST)") +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank())

  panel_sep <- ggplot(data.frame(gDNA = g_s, cfDNA = cf_s), aes(gDNA, cfDNA)) +
    geom_bin2d(bins = 80) +
    geom_abline(slope = 1, intercept = 0, colour = "firebrick",
                linewidth = 0.5, linetype = "dashed") +
    scale_fill_viridis_c(option = "magma", trans = "log10", name = "Count\n(log10)") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
             label = sprintf("r = %.4f\n%s gene×obs", r_s,
                             format(n_pts_filt, big.mark = ",")),
             size = 3) +
    coord_cartesian(xlim = sep_lim, ylim = sep_lim) +
    labs(title = "Separate VST (cfDNA alone / gDNA alone)",
         subtitle = sprintf("11,837 genes, %d obs, each source normalised independently", n_label),
         x = "gDNA 5-hMC (separate VST)", y = "cfDNA 5-hMC (separate VST)") +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank())

  panel_raw <- ggplot(data.frame(gDNA = g_r, cfDNA = cf_r), aes(gDNA, cfDNA)) +
    geom_bin2d(bins = 80) +
    geom_abline(slope = 1, intercept = 0, colour = "firebrick",
                linewidth = 0.5, linetype = "dashed") +
    scale_fill_viridis_c(option = "magma", trans = "log10", name = "Count\n(log10)") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
             label = sprintf("r = %.4f\n%s gene×obs", r_r,
                             format(n_pts_raw, big.mark = ",")),
             size = 3) +
    coord_cartesian(xlim = raw_lim, ylim = raw_lim) +
    labs(title = "Raw counts (no normalisation)",
         subtitle = sprintf("19,100 genes, %d obs", n_label),
         x = "gDNA 5-hMC (raw count)", y = "cfDNA 5-hMC (raw count)") +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank())

  panel_joint + panel_sep + panel_raw +
    plot_annotation(
      title    = page_title,
      subtitle = "Each point = one gene in one paired observation  |  dashed line = identity",
      theme    = theme(plot.title    = element_text(size = 11, face = "bold"),
                       plot.subtitle = element_text(size = 8, colour = "grey40"))
    )
}

## ---------------------------------------------------------------------------
## 6. Build pages and write PDF
## ---------------------------------------------------------------------------
cat("\nBuilding PDF...\n")
pdf(out_pdf, width = 15, height = 5.5)

## Page 1: all 84 observations pooled
cat("  Page 1: all cycles pooled...\n")
print(make_page(seq_len(ncol(md$Y_g)),
                "KRd Trial | cfDNA vs gDNA 5-hMC — all cycles pooled (n = 84 paired observations)",
                ncol(md$Y_g)))

## Pages 2-6: one per cycle
for (tt in seq_along(tp_levels)) {
  tp_label <- tp_levels[tt]
  obs_idx  <- which(paired$t_num == tt)
  if (length(obs_idx) == 0) next
  cat(sprintf("  Page %d: %s (n = %d obs)...\n", tt + 1, tp_label, length(obs_idx)))
  print(make_page(obs_idx,
                  sprintf("KRd Trial | cfDNA vs gDNA 5-hMC — %s (n = %d paired observations)",
                          tp_label, length(obs_idx)),
                  length(obs_idx)))
}

dev.off()
cat(sprintf("\nSaved: %s\n", out_pdf))
cat("Existing KRD_cfdna_vs_gdna_scatter.pdf was not modified.\n")
