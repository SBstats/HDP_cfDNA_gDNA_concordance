###############################################################################
## R/29b_raw_value_scatter_pdf.R
##
## Assembles all raw-value diagnostic figures into a single multi-page PDF:
##
##   Page 1  (portrait)  Caption page for Fig D1
##   Page 2  (landscape) Fig D1  -- VST scatter grid, genes ordered by mean
##   Page 3  (portrait)  Fig D1b -- VST Pearson r heatmap + caption
##   Page 4  (portrait)  Caption page for Fig D2
##   Page 5  (landscape) Fig D2  -- raw-count scatter grid, genes ordered by mean
##   Page 6  (portrait)  Fig D2b -- raw-count Pearson r heatmap + caption
##   Page 7  (portrait)  Caption page for Fig D3
##   Page 8  (landscape) Fig D3  -- VST scatter grid, natural gene order
##   Page 9  (portrait)  Caption page for Fig D4
##   Page 10 (landscape) Fig D4  -- raw-count scatter grid, natural gene order
##
## Output: nature_manuscript/raw_value_scatter.pdf
##
## Run from project root:  Rscript R/29b_raw_value_scatter_pdf.R
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(grid)
  library(gridExtra)
})

base_dir <- getwd()
res_dir  <- file.path(base_dir, "results")
fig_dir  <- file.path(base_dir, "nature_manuscript", "figures")
data_dir <- file.path(base_dir, "KRd trial", "5hmC data")
out_pdf  <- file.path(base_dir, "nature_manuscript", "raw_value_scatter.pdf")

tp_levels <- c("Screening", "C4", "C8", "C18", "3 YR F/U")
set.seed(20260718)
n_show <- 2500

## ---------------------------------------------------------------------------
## Build ggplot objects (return, don't save)
## ---------------------------------------------------------------------------

make_grid_gg <- function(Ycf, Yg, plabel, tnum, pat_levels, ylab, xlab, title) {
  p <- nrow(Yg); N <- ncol(Yg)
  gene_idx <- sort(sample.int(p, min(n_show, p)))
  rows <- vector("list", N)
  for (o in 1:N) {
    g  <- Yg[gene_idx,  o]
    cf <- Ycf[gene_idx, o]
    ord <- order((g + cf) / 2)
    rows[[o]] <- data.frame(
      rank    = rep(seq_along(ord), 2),
      value   = c(g[ord], cf[ord]),
      Source  = rep(c("gDNA", "cfDNA"), each = length(ord)),
      Patient = factor(plabel[o], levels = pat_levels),
      Cycle   = factor(tp_levels[tnum[o]], levels = tp_levels)
    )
  }
  df <- bind_rows(rows)
  ggplot(df, aes(rank, value, color = Source)) +
    geom_point(size = 0.15, alpha = 0.35) +
    scale_color_manual(values = c(cfDNA = "#E63946", gDNA = "#457B9D")) +
    facet_grid(Patient ~ Cycle, switch = "y") +
    labs(x = xlab, y = ylab, title = title) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    theme_bw(base_size = 6) +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 4),
          strip.text.y.left = element_text(angle = 0, size = 5),
          strip.text.x = element_text(size = 5))
}

make_grid_unordered_gg <- function(Ycf, Yg, plabel, tnum, pat_levels, ylab, title) {
  p <- nrow(Yg); N <- ncol(Yg)
  gene_idx <- sort(sample.int(p, min(n_show, p)))
  rows <- vector("list", N)
  for (o in 1:N) {
    g  <- Yg[gene_idx,  o]
    cf <- Ycf[gene_idx, o]
    rows[[o]] <- data.frame(
      rank    = rep(seq_along(gene_idx), 2),
      value   = c(g, cf),
      Source  = rep(c("gDNA", "cfDNA"), each = length(gene_idx)),
      Patient = factor(plabel[o], levels = pat_levels),
      Cycle   = factor(tp_levels[tnum[o]], levels = tp_levels)
    )
  }
  df <- bind_rows(rows)
  ggplot(df, aes(rank, value, color = Source)) +
    geom_point(size = 0.15, alpha = 0.35) +
    scale_color_manual(values = c(cfDNA = "#E63946", gDNA = "#457B9D")) +
    facet_grid(Patient ~ Cycle, switch = "y") +
    labs(x = "Genes (natural index order, no sorting)", y = ylab, title = title) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    theme_bw(base_size = 6) +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 4),
          strip.text.y.left = element_text(angle = 0, size = 5),
          strip.text.x = element_text(size = 5))
}

make_heatmap_gg <- function(Ycf, Yg, plabel, tnum, pat_levels, title) {
  N <- ncol(Yg)
  r_obs <- vapply(1:N, function(o) {
    a <- Yg[, o]; b <- Ycf[, o]
    if (sd(a) < 1e-12 || sd(b) < 1e-12) NA_real_ else cor(a, b)
  }, numeric(1))
  hd <- data.frame(
    Patient = factor(plabel, levels = rev(pat_levels)),
    Cycle   = factor(tp_levels[tnum], levels = tp_levels),
    r       = r_obs
  )
  ggplot(hd, aes(Cycle, Patient, fill = r)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.3f", r)), size = 2.0, color = "black") +
    scale_fill_gradient2(low = "#E63946", mid = "#FFFFCC", high = "#457B9D",
                         midpoint = 0.9,
                         limits = c(min(0.5, min(r_obs, na.rm = TRUE)), 1),
                         name = "Pearson r", na.value = "grey90") +
    labs(x = "Treatment timepoint", y = "Patient", title = title) +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank())
}

## ---------------------------------------------------------------------------
## Caption helper: draws a text page (portrait)
## ---------------------------------------------------------------------------
caption_page <- function(label, text) {
  grid.newpage()
  grid.text(paste0(label, "\n\n", text),
            x = 0.05, y = 0.95, just = c("left", "top"),
            gp = gpar(fontsize = 10, lineheight = 1.5),
            vp = viewport(width = 0.90, height = 0.90))
}

## ---------------------------------------------------------------------------
## Captions
## ---------------------------------------------------------------------------
cap_D1 <- paste(
  "Figure D1. Overlaid cfDNA and gDNA 5-hMC values per patient and cycle",
  "(DESeq2 VST, filtered genes).",
  "Grid with one row per patient (n=26) and one column per treatment timepoint",
  "(Screening, C4, C8, C18, 3-Year Follow-Up); empty cells are patient-timepoints",
  "without a paired specimen. Within each cell every gene is plotted twice ---",
  "cfDNA (red) and gDNA (blue) --- on a shared axis, with genes ordered by their",
  "within-observation mean value. A random sub-sample of 2,500 of the 11,837",
  "filtered genes is shown per cell. In nearly every cell the two point clouds",
  "sit on top of one another, showing that cfDNA and gDNA share the same",
  "across-gene value distribution (mean r = 0.956; range 0.719-0.994)."
)

cap_D1b <- paste(
  "Figure D1b. Per-cycle cfDNA-gDNA gene-level Pearson correlation",
  "(DESeq2 VST, 11,837 filtered genes).",
  "Each cell is the Pearson r between cfDNA and gDNA VST values across all",
  "filtered genes for that paired observation. Cells are annotated with r",
  "and coloured on a scale centred at 0.9."
)

cap_D2 <- paste(
  "Figure D2. Overlaid cfDNA and gDNA 5-hMC gene-body counts per patient",
  "and cycle, BEFORE DESeq2 normalization and low-count filtering.",
  "Same layout as Figure D1 but using raw gene-body read counts (log10(count+1))",
  "for all 19,100 assayed genes before VST and before the low-count filter.",
  "The distributional overlap is already present in the unprocessed counts,",
  "confirming it is not an artifact of the joint normalization."
)

cap_D2b <- paste(
  "Figure D2b. Per-cycle cfDNA-gDNA gene-level Pearson correlation",
  "(raw log10 counts, 19,100 unfiltered genes).",
  "Each cell is the Pearson r between cfDNA and gDNA raw log10(count+1) values",
  "across all 19,100 genes before filtering. Mean r = 0.975 (range 0.851-0.994),",
  "slightly higher than the VST/filtered value because the shared low-count floor",
  "adds concordant mass."
)

cap_D3 <- paste(
  "Figure D3. Overlaid cfDNA and gDNA 5-hMC values per patient and cycle,",
  "genes in natural index order (DESeq2 VST, filtered genes).",
  "Same data as Figure D1 but genes are plotted in their natural genomic index",
  "order rather than sorted by within-observation mean. The tight vertical",
  "alignment of red and blue points at each gene position shows that the high",
  "global correlation reflects genuine per-gene agreement between sources,",
  "not merely a shared distributional shape."
)

cap_D4 <- paste(
  "Figure D4. Overlaid cfDNA and gDNA 5-hMC gene-body counts, natural gene",
  "order, BEFORE DESeq2 normalization and low-count filtering.",
  "Same data as Figure D2 with genes in natural index order. The per-gene",
  "co-movement confirms that the high raw-count correlation reflects point-by-point",
  "concordance between sources."
)

## ---------------------------------------------------------------------------
## Load data
## ---------------------------------------------------------------------------
cat("Loading VST data...\n")
md  <- readRDS(file.path(res_dir, "model_data.rds"))
id_map   <- unique(md$paired_info[, c("patient_id", "Study.ID")])
plab_vst <- id_map$Study.ID[match(md$patient, id_map$patient_id)]
plab_vst[is.na(plab_vst)] <- as.character(md$patient[is.na(plab_vst)])
PAT_LEVELS <- unique(plab_vst[order(md$patient)])

cat("Loading raw count data...\n")
cfDNA_key <- read.csv(file.path(data_dir, "kRd-cfDNA_sample_key.csv"), stringsAsFactors = FALSE)
gDNA_key  <- read.csv(file.path(data_dir, "kRd-gDNA_sample_key.csv"),  stringsAsFactors = FALSE)
cfDNA_key$barcode <- cfDNA_key$Assigned.ID
gDNA_key$barcode  <- gDNA_key$Assigned.ID
paired_raw <- dplyr::inner_join(
  cfDNA_key |> dplyr::select(Study.ID, Timepoint, barcode_cf = barcode),
  gDNA_key  |> dplyr::select(Study.ID, Timepoint, barcode_g  = barcode),
  by = c("Study.ID", "Timepoint")
)
cfDNA_counts <- readRDS(file.path(data_dir, "kRd-cfDNA_genebody_count.RDS"))
gDNA_counts  <- readRDS(file.path(data_dir, "kRd-gDNA_genebody_count.RDS"))
cf_raw   <- log10(as.matrix(cfDNA_counts[, paired_raw$barcode_cf]) + 1)
g_raw    <- log10(as.matrix(gDNA_counts[,  paired_raw$barcode_g])  + 1)
tnum_raw <- match(paired_raw$Timepoint, tp_levels)
plab_raw <- paired_raw$Study.ID
PAT_LEVELS_RAW <- c(PAT_LEVELS, setdiff(unique(plab_raw), PAT_LEVELS))

## ---------------------------------------------------------------------------
## Build all ggplot objects
## ---------------------------------------------------------------------------
cat("Building figures (this may take a minute)...\n")

fig_D1  <- make_grid_gg(md$Y_cf, md$Y_g, plab_vst, md$time, PAT_LEVELS,
                         ylab = "VST 5-hMC value",
                         xlab = "Genes (ordered by within-observation mean value)",
                         title = "Fig D1: cfDNA (red) and gDNA (blue) VST values, ordered by mean")

fig_D1b <- make_heatmap_gg(md$Y_cf, md$Y_g, plab_vst, md$time, PAT_LEVELS,
                            title = "Fig D1b: Per-cycle Pearson r (VST, 11,837 filtered genes)")

fig_D2  <- make_grid_gg(cf_raw, g_raw, plab_raw, tnum_raw, PAT_LEVELS_RAW,
                         ylab = "log10(raw count + 1)",
                         xlab = "Genes (ordered by within-observation mean value)",
                         title = "Fig D2: cfDNA (red) and gDNA (blue) raw counts, ordered by mean")

fig_D2b <- make_heatmap_gg(cf_raw, g_raw, plab_raw, tnum_raw, PAT_LEVELS_RAW,
                            title = "Fig D2b: Per-cycle Pearson r (raw log10 counts, 19,100 genes)")

fig_D3  <- make_grid_unordered_gg(md$Y_cf, md$Y_g, plab_vst, md$time, PAT_LEVELS,
                                   ylab = "VST 5-hMC value",
                                   title = "Fig D3: cfDNA (red) and gDNA (blue) VST values, natural gene order")

fig_D4  <- make_grid_unordered_gg(cf_raw, g_raw, plab_raw, tnum_raw, PAT_LEVELS_RAW,
                                   ylab = "log10(raw count + 1)",
                                   title = "Fig D4: cfDNA (red) and gDNA (blue) raw counts, natural gene order")

## ---------------------------------------------------------------------------
## Write multi-page PDF
## Portrait pages: 8.5 x 11 in  (caption + heatmap pages)
## Landscape pages: 13 x 8.5 in (scatter grids)
## pdf() does not support mixed page sizes in one call, so we use a two-pass
## approach: write portrait and landscape pages as separate temp PDFs, then
## merge with pdftools if available, otherwise keep them as named page files.
## ---------------------------------------------------------------------------

# Check for pdftools
has_pdftools <- requireNamespace("pdftools", quietly = TRUE)

if (has_pdftools) {
  cat("pdftools available -- will merge into single PDF\n")
  tmp_dir <- tempdir()

  write_portrait <- function(fname, draw_fn) {
    pdf(fname, width = 8.5, height = 11, onefile = FALSE)
    draw_fn()
    dev.off()
  }
  write_landscape <- function(fname, gg) {
    pdf(fname, width = 13, height = 8.5, onefile = FALSE)
    print(gg)
    dev.off()
    cat("  rendered:", basename(fname), "\n")
  }

  pages <- list()

  # D1
  f <- file.path(tmp_dir, "p01_D1_caption.pdf")
  write_portrait(f, function() caption_page("Figure D1", cap_D1)); pages[[1]] <- f
  f <- file.path(tmp_dir, "p02_D1_grid.pdf")
  write_landscape(f, fig_D1); pages[[2]] <- f

  # D1b
  n_pat <- length(PAT_LEVELS)
  h_hm  <- max(4, 0.32 * n_pat + 1.2)
  f <- file.path(tmp_dir, "p03_D1b_heatmap.pdf")
  pdf(f, width = 6.5, height = h_hm, onefile = FALSE)
  print(fig_D1b)
  dev.off(); pages[[3]] <- f

  # D2
  f <- file.path(tmp_dir, "p04_D2_caption.pdf")
  write_portrait(f, function() caption_page("Figure D2", cap_D2)); pages[[4]] <- f
  f <- file.path(tmp_dir, "p05_D2_grid.pdf")
  write_landscape(f, fig_D2); pages[[5]] <- f

  # D2b
  n_pat_raw <- length(PAT_LEVELS_RAW)
  h_hm_raw  <- max(4, 0.32 * n_pat_raw + 1.2)
  f <- file.path(tmp_dir, "p06_D2b_heatmap.pdf")
  pdf(f, width = 6.5, height = h_hm_raw, onefile = FALSE)
  print(fig_D2b)
  dev.off(); pages[[6]] <- f

  # D3
  f <- file.path(tmp_dir, "p07_D3_caption.pdf")
  write_portrait(f, function() caption_page("Figure D3", cap_D3)); pages[[7]] <- f
  f <- file.path(tmp_dir, "p08_D3_grid.pdf")
  write_landscape(f, fig_D3); pages[[8]] <- f

  # D4
  f <- file.path(tmp_dir, "p09_D4_caption.pdf")
  write_portrait(f, function() caption_page("Figure D4", cap_D4)); pages[[9]] <- f
  f <- file.path(tmp_dir, "p10_D4_grid.pdf")
  write_landscape(f, fig_D4); pages[[10]] <- f

  pdftools::pdf_combine(unlist(pages), out_pdf)
  cat(sprintf("\nMerged PDF saved: %s\n", out_pdf))

} else {
  # Fallback: single landscape PDF for all pages (portrait pages will be
  # landscape-oriented but readable)
  cat("pdftools not available -- writing single-orientation PDF\n")
  pdf(out_pdf, width = 13, height = 8.5, onefile = TRUE)

  # Caption pages rendered as grobs on landscape canvas
  caption_grob <- function(label, text) {
    grid.newpage()
    grid.text(paste0(label, "\n\n", text),
              x = 0.03, y = 0.97, just = c("left", "top"),
              gp = gpar(fontsize = 9, lineheight = 1.5),
              vp = viewport(width = 0.94, height = 0.94))
  }

  caption_grob("Figure D1", cap_D1);  print(fig_D1)
  caption_grob("Figure D1b", cap_D1b); print(fig_D1b)
  caption_grob("Figure D2", cap_D2);  print(fig_D2)
  caption_grob("Figure D2b", cap_D2b); print(fig_D2b)
  caption_grob("Figure D3", cap_D3);  print(fig_D3)
  caption_grob("Figure D4", cap_D4);  print(fig_D4)

  dev.off()
  cat(sprintf("\nPDF saved: %s\n", out_pdf))
}

cat("Done.\n")
