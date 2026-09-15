#!/usr/bin/env Rscript
###############################################################################
## R/32_omega_trajectory_plots.R
##
## Grid plots of the longitudinal subclonal-weight trajectories.
##
## Reads:  results/krd_omega_trajectories.csv        (written by slurm/run_krd_combine.R)
##         results/krd_omega_trajectories_wide.csv
## Writes: figures/omega_traj_grid_by_patient.pdf    patient grid, all series
##         figures/omega_traj_grid_gdna_vs_cfdna.pdf patient grid, source-paired
##         figures/omega_traj_background.pdf         cfDNA non-tumour fraction
##         figures/omega_traj_cohort.pdf             cohort-level summary
##         figures/omega_cf_vs_g_scatter.pdf         matched-observation scatter
##
## Usage (from project root, after run_krd_combine.R):
##   Rscript R/32_omega_trajectory_plots.R
##
## This script deliberately does NOT refit anything; it is cheap to rerun while
## iterating on plot aesthetics.
##
## RELATIONSHIP TO R/30. R/30_generate_all_manuscript_outputs.R already draws a
## manuscript Figure 4 trajectory grid (krd_omega_trajectory_grid.pdf, faceted
## Study.ID x subclone, gDNA vs cfDNA). This script is not a replacement for it.
## The differences that motivated adding it:
##   (a) R/30 builds its trajectory frame in memory and never writes it out, so
##       the numbers behind the figure cannot be inspected, re-plotted, or
##       shared without rerunning the whole manuscript pipeline. The CSVs read
##       here are a durable artifact of run_krd_combine.R.
##   (b) R/30's trajectory omits the cfDNA non-tumour background entirely, so
##       the contamination trajectory -- the quantity whose earlier "negligible"
##       estimate was retracted -- had no per-patient longitudinal view at all.
##   (c) It adds credible-interval ribbons, a cohort-level summary, and the
##       matched-observation cfDNA-vs-gDNA scatter.
## Use R/30's version for the manuscript figure; use these for exploration and
## for anything involving the background.
##
## Design notes:
##  - The panel is UNBALANCED (T_i between 1 and 5). Gaps are missing data, not
##    absent subclones, so lines are simply not drawn across them and never
##    imputed to zero.
##  - The x-axis uses the real trial labels (Screening / C4 / C8 / C18 /
##    3 YR F/U) ordered by t_num, not integers.
##  - The cfDNA background is plotted as its own series, because it sums with
##    the tumour block to 1 rather than living inside it.
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

base_dir <- getwd()
res_dir  <- file.path(base_dir, "results")
fig_dir  <- file.path(base_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

traj_file <- file.path(res_dir, "krd_omega_trajectories.csv")
if (!file.exists(traj_file))
  stop("results/krd_omega_trajectories.csv not found. Run slurm/run_krd_combine.R first.")

traj <- read.csv(traj_file, stringsAsFactors = FALSE)
stopifnot(all(c("patient", "timepoint", "timepoint_label", "source",
                "subclone", "post_mean", "ci_lo", "ci_hi") %in% names(traj)))

## Chronological factor levels for the x-axis, taken from the data itself.
tp_levels <- traj %>%
  distinct(timepoint, timepoint_label) %>%
  arrange(timepoint) %>%
  pull(timepoint_label)
traj$timepoint_label <- factor(traj$timepoint_label, levels = tp_levels)

## Label panels by the trial Study.ID when available; fall back to the integer.
traj$panel <- if ("study_id" %in% names(traj)) traj$study_id else
              paste("Patient", traj$patient)

## A stable series key so colour never depends on row order.
traj$series <- ifelse(traj$source == "cfDNA background",
                      "cfDNA background",
                      paste0(traj$source, " · ", traj$subclone))

subclones <- sort(unique(traj$subclone[traj$source != "cfDNA background"]))
cat(sprintf("Loaded %d rows | %d patients | %d timepoints | subclones: %s\n",
            nrow(traj), length(unique(traj$patient)), length(tp_levels),
            paste(subclones, collapse = ", ")))

## ---------------------------------------------------------------------------
## Palette. Colour encodes the SUBCLONE; linetype/shape encodes the SOURCE.
## This is the key readability decision: a reader comparing compartments wants
## to see the same subclone in the same colour in both, with the compartment
## distinguished by a non-colour channel.
## ---------------------------------------------------------------------------
sub_cols <- setNames(
  c("#2C6E9B", "#C1553B", "#4C9A6A", "#9B7BB8", "#D89C3F")[seq_along(subclones)],
  subclones)
bg_col <- "#7A7A7A"

theme_traj <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        strip.background  = element_rect(fill = "grey93", colour = NA),
        strip.text        = element_text(size = 7.5, face = "bold"),
        axis.text.x       = element_text(angle = 45, hjust = 1, size = 6),
        legend.position   = "bottom",
        legend.key.width  = unit(1.1, "lines"))

n_pat  <- length(unique(traj$panel))
n_cols <- 6
n_rows <- ceiling(n_pat / n_cols)

## ---------------------------------------------------------------------------
## FIGURE 1. Patient grid, every series on one panel per patient.
## One panel per patient; colour = subclone, linetype = source, grey = background.
## ---------------------------------------------------------------------------
tum <- traj %>% filter(source != "cfDNA background")
bg  <- traj %>% filter(source == "cfDNA background")

p1 <- ggplot() +
  geom_line(data = tum,
            aes(timepoint_label, post_mean, colour = subclone,
                linetype = source, group = interaction(subclone, source)),
            linewidth = 0.45) +
  geom_point(data = tum,
             aes(timepoint_label, post_mean, colour = subclone, shape = source),
             size = 0.9) +
  geom_line(data = bg,
            aes(timepoint_label, post_mean, group = 1),
            colour = bg_col, linewidth = 0.45, linetype = "dotted") +
  geom_point(data = bg, aes(timepoint_label, post_mean),
             colour = bg_col, size = 0.9, shape = 4) +
  facet_wrap(~ panel, ncol = n_cols) +
  scale_colour_manual(values = sub_cols, name = "Subclone") +
  scale_linetype_manual(values = c(gDNA = "solid", cfDNA = "dashed"),
                        name = "Source") +
  scale_shape_manual(values = c(gDNA = 16, cfDNA = 17), name = "Source") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = NULL, y = expression("Posterior mean " * omega),
       title = "Longitudinal subclonal weight trajectories, by patient",
       subtitle = paste("Colour = subclone; solid/dashed = gDNA/cfDNA;",
                        "grey dotted (×) = cfDNA non-tumour background")) +
  theme_traj
ggsave(file.path(fig_dir, "omega_traj_grid_by_patient.pdf"), p1,
       width = 13, height = 2.1 * n_rows + 1.6, limitsize = FALSE)
cat("  Saved: figures/omega_traj_grid_by_patient.pdf\n")

## ---------------------------------------------------------------------------
## FIGURE 2. Source-paired grid: rows = subclone, columns = patient, with a
## credible-interval ribbon. This is the figure for reading concordance --
## whether the cfDNA curve moves with the gDNA curve within each patient.
## ---------------------------------------------------------------------------
p2 <- ggplot(tum, aes(timepoint_label, post_mean,
                      colour = source, fill = source, group = source)) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.16, colour = NA) +
  geom_line(linewidth = 0.45) +
  geom_point(size = 0.8) +
  facet_grid(subclone ~ panel, switch = "y") +
  scale_colour_manual(values = c(gDNA = "#2C6E9B", cfDNA = "#C1553B"),
                      name = "Source") +
  scale_fill_manual(values = c(gDNA = "#2C6E9B", cfDNA = "#C1553B"),
                    name = "Source") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = NULL, y = expression("Posterior mean " * omega * "  (95% CrI)"),
       title = "gDNA vs cfDNA subclonal weight trajectories",
       subtitle = "Ribbon = 95% credible interval. Gaps are unobserved timepoints, not zero weight.") +
  theme_traj +
  theme(strip.text.x = element_text(size = 5.5),
        axis.text.x  = element_text(size = 4.5))
ggsave(file.path(fig_dir, "omega_traj_grid_gdna_vs_cfdna.pdf"), p2,
       width = 1.05 * n_pat + 2, height = 2.3 * length(subclones) + 1.6,
       limitsize = FALSE)
cat("  Saved: figures/omega_traj_grid_gdna_vs_cfdna.pdf\n")

## ---------------------------------------------------------------------------
## FIGURE 3. cfDNA non-tumour background trajectory, per patient.
## Reported on its own because it is the quantity whose earlier "negligible"
## estimate was retracted; it deserves to be inspected directly.
## ---------------------------------------------------------------------------
p3 <- ggplot(bg, aes(timepoint_label, post_mean, group = 1)) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.2,
              fill = bg_col, colour = NA) +
  geom_line(colour = bg_col, linewidth = 0.5) +
  geom_point(colour = bg_col, size = 1) +
  facet_wrap(~ panel, ncol = n_cols) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = NULL, y = expression("Posterior mean " * omega[0] * "  (95% CrI)"),
       title = "cfDNA non-tumour (normal tissue) fraction over time",
       subtitle = "Per subject-timepoint, with 95% credible intervals") +
  theme_traj
ggsave(file.path(fig_dir, "omega_traj_background.pdf"), p3,
       width = 13, height = 2.1 * n_rows + 1.4, limitsize = FALSE)
cat("  Saved: figures/omega_traj_background.pdf\n")

## ---------------------------------------------------------------------------
## FIGURE 4. Cohort-level summary: median across patients at each timepoint,
## with the IQR band. Patient counts differ by timepoint (unbalanced panel), so
## n is printed on the axis to keep that visible.
## ---------------------------------------------------------------------------
coh <- traj %>%
  group_by(series, subclone, source, timepoint, timepoint_label) %>%
  summarise(med = median(post_mean), lo = quantile(post_mean, 0.25),
            hi = quantile(post_mean, 0.75), n = dplyr::n(), .groups = "drop")
n_by_tp <- traj %>% group_by(timepoint_label) %>%
  summarise(n = n_distinct(patient), .groups = "drop") %>% arrange(timepoint_label)
lab_tp <- setNames(sprintf("%s\n(n=%d)", n_by_tp$timepoint_label, n_by_tp$n),
                   as.character(n_by_tp$timepoint_label))

p4 <- ggplot(coh, aes(timepoint_label, med, colour = series, fill = series,
                      group = series)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.8) +
  scale_x_discrete(labels = lab_tp) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = expression("Median " * omega * " across patients (IQR)"),
       title = "Cohort-level subclonal weight trajectories",
       subtitle = "Median over patients observed at each timepoint; band = IQR") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "right",
        legend.title = element_blank())
ggsave(file.path(fig_dir, "omega_traj_cohort.pdf"), p4, width = 9, height = 5)
cat("  Saved: figures/omega_traj_cohort.pdf\n")

## ---------------------------------------------------------------------------
## FIGURE 5. cfDNA against gDNA at MATCHED observations, faceted by subclone.
## The concordance question in its most direct form: points on the diagonal
## mean the two compartments agree at that subject-timepoint.
## ---------------------------------------------------------------------------
wide_file <- file.path(res_dir, "krd_omega_trajectories_wide.csv")
if (file.exists(wide_file)) {
  w <- read.csv(wide_file, stringsAsFactors = FALSE)
  ks <- as.integer(sub("^omega_g_k", "", grep("^omega_g_k", names(w), value = TRUE)))
  sc <- do.call(rbind, lapply(ks, function(k) data.frame(
    patient   = w$patient,
    timepoint = w$timepoint,
    subclone  = paste0("subclone ", k),
    g         = w[[sprintf("omega_g_k%d",  k)]],
    cf        = w[[sprintf("omega_cf_k%d", k)]],
    stringsAsFactors = FALSE)))

  cors <- sc %>% group_by(subclone) %>%
    summarise(r = suppressWarnings(cor(g, cf, use = "complete.obs")),
              .groups = "drop") %>%
    mutate(lab = sprintf("r = %.3f", r))

  p5 <- ggplot(sc, aes(g, cf)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = "dashed") +
    geom_point(aes(colour = factor(timepoint)), size = 1.6, alpha = 0.85) +
    geom_text(data = cors, aes(x = 0.04, y = 0.96, label = lab),
              hjust = 0, size = 3.2, inherit.aes = FALSE) +
    facet_wrap(~ subclone) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    scale_colour_viridis_d(name = "Timepoint", option = "D", end = 0.9) +
    labs(x = expression("gDNA " * omega), y = expression("cfDNA " * omega),
         title = "cfDNA vs gDNA subclonal weights at matched observations",
         subtitle = "One point per (patient, timepoint); dashed line is identity") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank())
  ggsave(file.path(fig_dir, "omega_cf_vs_g_scatter.pdf"), p5,
         width = 4.2 * min(length(ks), 3) + 1.2, height = 4.6)
  cat("  Saved: figures/omega_cf_vs_g_scatter.pdf\n")
  print(as.data.frame(cors[, c("subclone", "r")]), row.names = FALSE)
} else {
  cat("  (wide file absent; skipping the matched-observation scatter)\n")
}

cat("\nDone. 5 trajectory figures written to figures/.\n")
