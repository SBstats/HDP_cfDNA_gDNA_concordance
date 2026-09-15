###############################################################################
## R/29_raw_value_scatter_grid.R
##
## GOAL: show WHY cfDNA and gDNA correlate so highly at the global gene level.
## For every paired observation (patient x cycle) we overlay the raw 5-hMC
## values of the two sources on a shared axis: each gene is plotted twice ---
## once for cfDNA (red) and once for gDNA (blue) --- with genes ordered on the x
## axis by their within-observation mean value. If the two sources share the
## same value distribution (the driver of the near-1 gene-level correlation),
## the two colored clouds sit essentially on top of one another and trace the
## same monotone sweep from low- to high-signal genes.
##
## Two versions are produced (a grid: one row per patient, one column per cycle):
##   (1) VST version   -- the analysis-ready DESeq2 variance-stabilized values,
##                        11,837 filtered genes (results/model_data.rds).
##   (2) RAW version    -- the raw gene-body counts BEFORE DESeq2 normalization
##                        and BEFORE low-count filtering, all 19,100 genes,
##                        plotted on log10(count + 1). This shows the same
##                        distributional overlap is present in the raw data and
##                        is not an artifact of the joint normalization.
##
## Outputs: nature_manuscript/figures/fig_raw_value_scatter_grid.pdf   (VST)
##          nature_manuscript/figures/fig_raw_value_scatter_grid_rawcounts.pdf
###############################################################################

suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })

base_dir <- getwd()  # run from project root
setwd(base_dir)
res_dir  <- file.path(base_dir, "results")
fig_dir  <- file.path(base_dir, "nature_manuscript", "figures")
data_dir <- file.path(base_dir, "KRd trial", "5hmC data")

tp_levels <- c("Screening", "C4", "C8", "C18", "3 YR F/U")
set.seed(20260718)
n_show <- 2500                      # genes sub-sampled per cell for legibility

## ---------------------------------------------------------------------------
## Shared plotting routine: given aligned [p x N_obs] cfDNA and gDNA matrices,
## per-observation patient labels and cycle indices, build the overlaid grid.
## ---------------------------------------------------------------------------
make_grid <- function(Ycf, Yg, plabel, tnum, pat_levels, ylab, title, subtitle, out_pdf) {
  p <- nrow(Yg); N <- ncol(Yg)
  gene_idx <- sort(sample.int(p, min(n_show, p)))

  rows <- vector("list", N)
  for (o in 1:N) {
    g  <- Yg[gene_idx,  o]
    cf <- Ycf[gene_idx, o]
    ord <- order((g + cf) / 2)          # order genes by within-obs mean value
    rows[[o]] <- data.frame(
      rank    = rep(seq_along(ord), 2),
      value   = c(g[ord], cf[ord]),
      Source  = rep(c("gDNA", "cfDNA"), each = length(ord)),
      Patient = factor(plabel[o], levels = pat_levels),
      Cycle   = factor(tp_levels[tnum[o]], levels = tp_levels)
    )
  }
  df <- bind_rows(rows)

  p_fig <- ggplot(df, aes(rank, value, color = Source)) +
    geom_point(size = 0.15, alpha = 0.35) +
    scale_color_manual(values = c(cfDNA = "#E63946", gDNA = "#457B9D")) +
    facet_grid(Patient ~ Cycle, switch = "y") +
    labs(x = "Genes (ordered by within-observation mean value)",
         y = ylab, title = title, subtitle = subtitle) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    theme_bw(base_size = 6) +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 4),
          strip.text.y.left = element_text(angle = 0, size = 5),
          strip.text.x = element_text(size = 5))

  n_pat <- nlevels(df$Patient)
  ## Landscape aspect (~13 x 8.5, ratio ~1.5) so the grid scales down to fill a
  ## landscape Letter printable area (~9.5in wide x ~6.5in tall) without becoming
  ## a thin strip. Width driven; height held to keep 26 rows legible.
  ggsave(out_pdf, p_fig, width = 13, height = 8.5, limitsize = FALSE)
  cat(sprintf("Saved: %s  (%d patients x %d cycles; %d genes/cell)\n",
              out_pdf, n_pat, length(tp_levels), length(gene_idx)))
}

## ---------------------------------------------------------------------------
## Same as make_grid but genes are shown in their natural (sub-sampled) index
## order rather than sorted by within-observation mean value.
## ---------------------------------------------------------------------------
make_grid_unordered <- function(Ycf, Yg, plabel, tnum, pat_levels, ylab, title, subtitle, out_pdf) {
  p <- nrow(Yg); N <- ncol(Yg)
  gene_idx <- sort(sample.int(p, min(n_show, p)))

  rows <- vector("list", N)
  for (o in 1:N) {
    g  <- Yg[gene_idx,  o]
    cf <- Ycf[gene_idx, o]
    # no sorting — natural gene index order
    rows[[o]] <- data.frame(
      rank    = rep(seq_along(gene_idx), 2),
      value   = c(g, cf),
      Source  = rep(c("gDNA", "cfDNA"), each = length(gene_idx)),
      Patient = factor(plabel[o], levels = pat_levels),
      Cycle   = factor(tp_levels[tnum[o]], levels = tp_levels)
    )
  }
  df <- bind_rows(rows)

  p_fig <- ggplot(df, aes(rank, value, color = Source)) +
    geom_point(size = 0.15, alpha = 0.35) +
    scale_color_manual(values = c(cfDNA = "#E63946", gDNA = "#457B9D")) +
    facet_grid(Patient ~ Cycle, switch = "y") +
    labs(x = "Genes (natural index order, no sorting)",
         y = ylab, title = title, subtitle = subtitle) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    theme_bw(base_size = 6) +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 4),
          strip.text.y.left = element_text(angle = 0, size = 5),
          strip.text.x = element_text(size = 5))

  n_pat <- nlevels(df$Patient)
  ggsave(out_pdf, p_fig, width = 13, height = 8.5, limitsize = FALSE)
  cat(sprintf("Saved: %s  (%d patients x %d cycles; %d genes/cell)\n",
              out_pdf, n_pat, length(tp_levels), length(gene_idx)))
}

## ---------------------------------------------------------------------------
## Per-cycle Pearson correlation heatmap: for each paired observation, the
## cfDNA-vs-gDNA gene-level Pearson r (over ALL genes, not the sub-sample) laid
## out as a patient x cycle grid. This is the correlation the paired scatterplot
## in make_grid() visualizes, made quantitative per cell.
## ---------------------------------------------------------------------------
make_corr_heatmap <- function(Ycf, Yg, plabel, tnum, pat_levels, title, out_pdf) {
  N <- ncol(Yg)
  r_obs <- vapply(1:N, function(o) {
    a <- Yg[, o]; b <- Ycf[, o]
    if (sd(a) < 1e-12 || sd(b) < 1e-12) NA_real_ else cor(a, b)
  }, numeric(1))

  hd <- data.frame(
    Patient = factor(plabel, levels = rev(pat_levels)),   # top-to-bottom = grid order
    Cycle   = factor(tp_levels[tnum], levels = tp_levels),
    r       = r_obs
  )

  p_fig <- ggplot(hd, aes(Cycle, Patient, fill = r)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.3f", r)), size = 2.0, color = "black") +
    scale_fill_gradient2(low = "#E63946", mid = "#FFFFCC", high = "#457B9D",
                         midpoint = 0.9, limits = c(min(0.5, min(r_obs, na.rm = TRUE)), 1),
                         name = "Pearson r", na.value = "grey90") +
    labs(x = "Treatment timepoint", y = "Patient", title = title) +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank())

  n_pat <- length(pat_levels)
  ggsave(out_pdf, p_fig, width = 6.5, height = 0.32 * n_pat + 1.2, limitsize = FALSE)
  cat(sprintf("Saved: %s  (per-cycle Pearson r; mean=%.3f, range %.3f-%.3f)\n",
              out_pdf, mean(r_obs, na.rm = TRUE),
              min(r_obs, na.rm = TRUE), max(r_obs, na.rm = TRUE)))
}

## ---------------------------------------------------------------------------
## (1) VST version -- analysis-ready values (11,837 filtered genes)
## ---------------------------------------------------------------------------
md  <- readRDS(file.path(res_dir, "model_data.rds"))
id_map <- unique(md$paired_info[, c("patient_id", "Study.ID")])
plab_vst <- id_map$Study.ID[match(md$patient, id_map$patient_id)]
plab_vst[is.na(plab_vst)] <- as.character(md$patient[is.na(plab_vst)])

## ---------------------------------------------------------------------------
## CANONICAL patient row order, shared by ALL FOUR figures so rows line up for
## head-to-head comparison. Ordered by the internal patient_id (analysis order),
## labelled by Study.ID. Any label not in this list (should be none) is appended.
## ---------------------------------------------------------------------------
PAT_LEVELS <- unique(plab_vst[order(md$patient)])

make_grid(
  Ycf = md$Y_cf, Yg = md$Y_g, plabel = plab_vst, tnum = md$time, pat_levels = PAT_LEVELS,
  ylab = "VST 5-hMC value",
  title = "Raw cfDNA (red) and gDNA (blue) 5-hMC values overlaid per patient and cycle (DESeq2 VST, filtered genes)",
  subtitle = sprintf("Each gene plotted twice; %d of 11,837 filtered genes sub-sampled per cell. Overlapping clouds indicate the two sources share the same value distribution.", n_show),
  out_pdf = file.path(fig_dir, "fig_raw_value_scatter_grid.pdf")
)

make_corr_heatmap(
  Ycf = md$Y_cf, Yg = md$Y_g, plabel = plab_vst, tnum = md$time, pat_levels = PAT_LEVELS,
  title = "Per-cycle cfDNA-gDNA gene-level Pearson correlation (DESeq2 VST, 11,837 filtered genes)",
  out_pdf = file.path(fig_dir, "fig_raw_value_corr_heatmap.pdf")
)

## ---------------------------------------------------------------------------
## (2) RAW version -- gene-body counts BEFORE DESeq2 and BEFORE filtering
##     (all 19,100 genes). Rebuild the identical patient x cycle pairing from
##     the sample keys (same logic as R/02_preprocess.R steps 1-3).
## ---------------------------------------------------------------------------
cfDNA_key <- read.csv(file.path(data_dir, "kRd-cfDNA_sample_key.csv"), stringsAsFactors = FALSE)
gDNA_key  <- read.csv(file.path(data_dir, "kRd-gDNA_sample_key.csv"),  stringsAsFactors = FALSE)
cfDNA_key$barcode <- cfDNA_key$Assigned.ID
gDNA_key$barcode  <- gDNA_key$Assigned.ID

paired <- inner_join(
  cfDNA_key %>% select(Study.ID, Timepoint, barcode_cf = barcode),
  gDNA_key  %>% select(Study.ID, Timepoint, barcode_g  = barcode),
  by = c("Study.ID", "Timepoint")
)

cfDNA_counts <- readRDS(file.path(data_dir, "kRd-cfDNA_genebody_count.RDS"))
gDNA_counts  <- readRDS(file.path(data_dir, "kRd-gDNA_genebody_count.RDS"))
stopifnot(all(rownames(cfDNA_counts) == rownames(gDNA_counts)))

cf_raw <- as.matrix(cfDNA_counts[, paired$barcode_cf])   # [19100 x N_paired]
g_raw  <- as.matrix(gDNA_counts[,  paired$barcode_g])

## map Timepoint label -> cycle index; patient rows use the SAME canonical order
## as the VST figures (PAT_LEVELS). Any raw-only Study.ID not in the analysis set
## is appended at the end so no patient is dropped, but the shared patients keep
## identical row positions across all four figures.
tnum_raw <- match(paired$Timepoint, tp_levels)
plab_raw <- paired$Study.ID
PAT_LEVELS_RAW <- c(PAT_LEVELS, setdiff(unique(plab_raw), PAT_LEVELS))

cat(sprintf("\nRaw-count grid: %d paired observations, %d genes (pre-filter, pre-VST)\n",
            ncol(cf_raw), nrow(cf_raw)))
extra <- setdiff(unique(plab_raw), PAT_LEVELS)
if (length(extra)) cat("  raw-only patients appended (not in analysis set):",
                       paste(extra, collapse = ", "), "\n")

make_grid(
  Ycf = log10(cf_raw + 1), Yg = log10(g_raw + 1),
  plabel = plab_raw, tnum = tnum_raw, pat_levels = PAT_LEVELS_RAW,
  ylab = expression(log[10](raw~count + 1)),
  title = "Raw cfDNA (red) and gDNA (blue) 5-hMC gene-body counts overlaid per patient and cycle (BEFORE DESeq2 and low-count filtering)",
  subtitle = sprintf("Each gene plotted twice; %d of 19,100 unfiltered genes sub-sampled per cell, on log10(count+1). The distributional overlap is present in the raw counts, not created by normalization.", n_show),
  out_pdf = file.path(fig_dir, "fig_raw_value_scatter_grid_rawcounts.pdf")
)

make_corr_heatmap(
  Ycf = log10(cf_raw + 1), Yg = log10(g_raw + 1),
  plabel = plab_raw, tnum = tnum_raw, pat_levels = PAT_LEVELS_RAW,
  title = "Per-cycle cfDNA-gDNA gene-level Pearson correlation (raw log10 counts, 19,100 unfiltered genes)",
  out_pdf = file.path(fig_dir, "fig_raw_value_corr_heatmap_rawcounts.pdf")
)

## ---------------------------------------------------------------------------
## (3) Unordered variants — same data, genes in natural index order (no x-axis sort)
## ---------------------------------------------------------------------------
make_grid_unordered(
  Ycf = md$Y_cf, Yg = md$Y_g, plabel = plab_vst, tnum = md$time, pat_levels = PAT_LEVELS,
  ylab = "VST 5-hMC value",
  title = "Raw cfDNA (red) and gDNA (blue) 5-hMC values overlaid per patient and cycle (DESeq2 VST, filtered genes)",
  subtitle = sprintf("Each gene plotted twice; %d of 11,837 filtered genes sub-sampled per cell. Genes in natural index order (no sorting by mean).", n_show),
  out_pdf = file.path(fig_dir, "fig_raw_value_scatter_grid_unordered.pdf")
)

make_grid_unordered(
  Ycf = log10(cf_raw + 1), Yg = log10(g_raw + 1),
  plabel = plab_raw, tnum = tnum_raw, pat_levels = PAT_LEVELS_RAW,
  ylab = expression(log[10](raw~count + 1)),
  title = "Raw cfDNA (red) and gDNA (blue) 5-hMC gene-body counts overlaid per patient and cycle (BEFORE DESeq2 and low-count filtering)",
  subtitle = sprintf("Each gene plotted twice; %d of 19,100 unfiltered genes sub-sampled per cell, on log10(count+1). Genes in natural index order (no sorting by mean).", n_show),
  out_pdf = file.path(fig_dir, "fig_raw_value_scatter_grid_rawcounts_unordered.pdf")
)

cat("\nDone.\n")
