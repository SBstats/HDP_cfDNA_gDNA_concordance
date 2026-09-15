###############################################################################
## R/29c_cfdna_vs_gdna_scatter.R
##
## cfDNA vs gDNA 5-hMC scatter, 6-page PDF:
##   Page 1: all 84 observations pooled
##   Pages 2-6: one page per cycle (Screening, C4, C8, C18, 3 YR F/U)
##
## Each page = two panels:
##   Left:  after DESeq2 VST + low-count filtering (11,837 genes)
##   Right: before DESeq2 and filtering, raw read counts (19,100 genes)
##
## Output: nature_manuscript/figures/fig_cfdna_vs_gdna_scatter.pdf
## Run from project root: Rscript R/29c_cfdna_vs_gdna_scatter.R
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(pdftools)
})

base_dir <- getwd()
res_dir  <- file.path(base_dir, "results")
data_dir <- file.path(base_dir, "KRd trial", "5hmC data")
fig_dir  <- file.path(base_dir, "nature_manuscript", "figures")
out_pdf  <- file.path(fig_dir, "fig_cfdna_vs_gdna_scatter.pdf")
tmp_dir  <- tempdir()

tp_levels <- c("Screening", "C4", "C8", "C18", "3 YR F/U")

## ---------------------------------------------------------------------------
## Load data
## ---------------------------------------------------------------------------
cat("Loading VST data...\n")
md <- readRDS(file.path(res_dir, "model_data.rds"))

cat("Loading raw count data...\n")
cfDNA_key <- read.csv(file.path(data_dir, "kRd-cfDNA_sample_key.csv"), stringsAsFactors = FALSE)
gDNA_key  <- read.csv(file.path(data_dir, "kRd-gDNA_sample_key.csv"),  stringsAsFactors = FALSE)
cfDNA_key$barcode <- cfDNA_key$Assigned.ID
gDNA_key$barcode  <- gDNA_key$Assigned.ID

paired_raw <- inner_join(
  cfDNA_key |> select(Study.ID, Timepoint, barcode_cf = barcode),
  gDNA_key  |> select(Study.ID, Timepoint, barcode_g  = barcode),
  by = c("Study.ID", "Timepoint")
)
cfDNA_counts <- readRDS(file.path(data_dir, "kRd-cfDNA_genebody_count.RDS"))
gDNA_counts  <- readRDS(file.path(data_dir, "kRd-gDNA_genebody_count.RDS"))

## Full raw matrices [genes x N_obs], aligned to paired_raw order
g_raw_mat  <- as.matrix(gDNA_counts[,  paired_raw$barcode_g])
cf_raw_mat <- as.matrix(cfDNA_counts[, paired_raw$barcode_cf])

## ---------------------------------------------------------------------------
## Shared axis limits (computed once so all pages are comparable)
## ---------------------------------------------------------------------------
vst_xlim <- range(md$Y_g)
vst_ylim <- range(md$Y_cf)
raw_lim  <- range(c(g_raw_mat, cf_raw_mat))

## ---------------------------------------------------------------------------
## Helper: build one two-panel page
## ---------------------------------------------------------------------------
make_page <- function(g_vst, cf_vst, g_raw, cf_raw,
                      page_title, n_obs_vst, n_obs_raw) {
  r_vst <- cor(g_vst, cf_vst)
  r_raw <- cor(g_raw,  cf_raw)

  panel_A <- ggplot(data.frame(gDNA = g_vst, cfDNA = cf_vst), aes(gDNA, cfDNA)) +
    geom_bin2d(bins = 80) +
    geom_abline(slope = 1, intercept = 0, color = "firebrick",
                linewidth = 0.5, linetype = "dashed") +
    scale_fill_viridis_c(option = "magma", trans = "log10", name = "Count\n(log10)") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
             label = sprintf("r = %.4f\nn genes x obs = %s",
                             r_vst, format(length(g_vst), big.mark = ",")),
             size = 3, color = "black") +
    coord_cartesian(xlim = vst_xlim, ylim = vst_ylim) +
    labs(title = "After DESeq2 VST + filtering",
         subtitle = sprintf("11,837 genes x %d observations", n_obs_vst),
         x = "gDNA 5-hMC (VST)", y = "cfDNA 5-hMC (VST)") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), legend.position = "right")

  panel_B <- ggplot(data.frame(gDNA = g_raw, cfDNA = cf_raw), aes(gDNA, cfDNA)) +
    geom_bin2d(bins = 80) +
    geom_abline(slope = 1, intercept = 0, color = "firebrick",
                linewidth = 0.5, linetype = "dashed") +
    scale_fill_viridis_c(option = "magma", trans = "log10", name = "Count\n(log10)") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
             label = sprintf("r = %.4f\nn genes x obs = %s",
                             r_raw, format(length(g_raw), big.mark = ",")),
             size = 3, color = "black") +
    coord_cartesian(xlim = raw_lim, ylim = raw_lim) +
    labs(title = "Before DESeq2 and filtering",
         subtitle = sprintf("19,100 genes x %d observations, raw counts", n_obs_raw),
         x = "gDNA 5-hMC (raw count)", y = "cfDNA 5-hMC (raw count)") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), legend.position = "right")

  panel_A + panel_B +
    plot_annotation(
      title = page_title,
      subtitle = "Each point = one gene in one paired observation; dashed line = identity",
      theme = theme(plot.title    = element_text(size = 12, face = "bold"),
                    plot.subtitle = element_text(size = 9))
    )
}

## ---------------------------------------------------------------------------
## Page 1: all observations pooled
## ---------------------------------------------------------------------------
cat("Building page 1 (all cycles pooled)...\n")
pages <- list()

fig1 <- make_page(
  g_vst    = as.vector(md$Y_g),
  cf_vst   = as.vector(md$Y_cf),
  g_raw    = as.vector(g_raw_mat),
  cf_raw   = as.vector(cf_raw_mat),
  page_title = "KRd Trial | cfDNA vs gDNA 5-hMC -- all cycles pooled (n = 84 paired observations)",
  n_obs_vst = ncol(md$Y_g),
  n_obs_raw = ncol(g_raw_mat)
)
f <- file.path(tmp_dir, "page_01_all.pdf")
ggsave(f, fig1, width = 10, height = 5)
pages[[1]] <- f

## ---------------------------------------------------------------------------
## Pages 2-6: one per cycle
## ---------------------------------------------------------------------------
for (tt in seq_along(tp_levels)) {
  tp_label <- tp_levels[tt]
  cat(sprintf("Building page %d (%s)...\n", tt + 1, tp_label))

  # VST: subset observations at this timepoint
  obs_vst <- which(md$time == tt)
  g_vst_t  <- as.vector(md$Y_g[,  obs_vst])
  cf_vst_t <- as.vector(md$Y_cf[, obs_vst])
  n_vst_t  <- length(obs_vst)

  # Raw: subset paired_raw to this timepoint
  obs_raw  <- which(paired_raw$Timepoint == tp_label)
  g_raw_t  <- as.vector(g_raw_mat[,  obs_raw])
  cf_raw_t <- as.vector(cf_raw_mat[, obs_raw])
  n_raw_t  <- length(obs_raw)

  if (n_vst_t == 0 && n_raw_t == 0) next

  fig_t <- make_page(
    g_vst    = g_vst_t,
    cf_vst   = cf_vst_t,
    g_raw    = g_raw_t,
    cf_raw   = cf_raw_t,
    page_title = sprintf("KRd Trial | cfDNA vs gDNA 5-hMC -- %s (n = %d paired observations)",
                         tp_label, n_vst_t),
    n_obs_vst = n_vst_t,
    n_obs_raw = n_raw_t
  )
  f <- file.path(tmp_dir, sprintf("page_%02d_%s.pdf", tt + 1,
                                  gsub("[^A-Za-z0-9]", "_", tp_label)))
  ggsave(f, fig_t, width = 10, height = 5)
  pages[[tt + 1]] <- f
}

## ---------------------------------------------------------------------------
## Merge into single PDF
## ---------------------------------------------------------------------------
pdf_combine(unlist(pages), out_pdf)
cat(sprintf("\nSaved: %s  (%d pages)\n", out_pdf, length(pages)))
