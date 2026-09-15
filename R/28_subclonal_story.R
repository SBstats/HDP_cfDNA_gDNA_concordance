###############################################################################
## R/28_subclonal_story.R
##
## GOAL (#5, descriptive / hypothesis-generating): tell the per-patient "story"
## linking longitudinal SUBCLONAL CHANGE to the clinical outcome (MRD by NGS and
## best overall response). With n = 26 (few with MRD conversions) this is NOT a
## powered test -- we surface individual trajectories where the dominant-subclone
## weight shift aligns with MRD status / response, as illustrative narrative.
##
## For each patient we plot the dominant-subclone weight w1(t) across treatment
## timepoints for cfDNA (solid) and gDNA (dashed), one small panel per patient,
## with the panel annotated by best overall response and MRD status points
## overlaid at the MRD-measured cycles (C8/C18/3Y). A companion table flags a
## handful of illustrative patients (largest |w1 shift|, or a shift co-occurring
## with an MRD conversion).
##
## Reuses the R/12 weight-extraction + clinical-loader recipe verbatim.
##
## Outputs:
##   nature_manuscript/figures/fig_subclonal_story.pdf
##   results/subclonal_story_illustrative.csv
###############################################################################

set.seed(20260718)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readxl); library(ggplot2)
})

base_dir  <- getwd()
res_dir   <- file.path(base_dir, "results")
fig_dir   <- file.path(base_dir, "nature_manuscript", "figures")
clin_file <- file.path(base_dir, "KRd trial", "Clinical data",
                       "KRd 12_1725 data_2_14_2023 from Ben without PHI.xlsx")
source(file.path(base_dir, "R", "01_lib_core.R"))

## ---------------------------------------------------------------------------
## 1. Canonical posterior-mean weights per observation (same as R/12)
## ---------------------------------------------------------------------------
md    <- readRDS(file.path(res_dir, "model_data.rds"))
n_obs <- md$N_obs
m1    <- readRDS(file.path(res_dir, "krd_M1_fit.rds"))
for (i in seq_along(m1$chains)) m1$chains[[i]] <- relabel_by_weight(m1$chains[[i]])
m1$chains <- align_chains_by_signature(m1$chains)
merged <- merge_chain_samples(m1$chains)

K_trace <- ncol(merged$omega_g_trace) / n_obs
col_of  <- function(k, obs) (k - 1) * n_obs + obs
wmean <- function(trace) {
  W <- matrix(NA_real_, n_obs, K_trace)
  for (k in 1:K_trace) for (o in 1:n_obs) W[o, k] <- mean(trace[, col_of(k, o)])
  W
}
Wg  <- wmean(merged$omega_g_trace)
Wcf <- wmean(merged$omega_cf_trace)

id_map <- unique(md$paired_info[, c("patient_id", "Study.ID")])
obs_df <- data.frame(
  patient_id = md$patient, t_num = md$time,
  w1g = Wg[, 1], w1cf = Wcf[, 1]
) %>%
  group_by(patient_id, t_num) %>%
  summarise(w1g = mean(w1g), w1cf = mean(w1cf), .groups = "drop") %>%
  left_join(id_map, by = "patient_id")

tp_levels <- c("Screening", "C4", "C8", "C18", "3 YR F/U")
obs_df$Timepoint <- factor(tp_levels[obs_df$t_num], levels = tp_levels)

## ---------------------------------------------------------------------------
## 2. Clinical: MRD (per window) + best overall response (same coding as R/12)
## ---------------------------------------------------------------------------
suppressWarnings({ clin <- read_excel(clin_file, sheet = "12_1725 cfDNA Dataset") })
mrd_df <- data.frame(
  Study.ID = clin$`Study number`, best_response = clin$`best overall response`,
  neg_c8 = clin$`MRD post 8 cycles by NGS (Adaptive) negative =1`,
  neg_c18 = clin$`MRD by NGS EoT Negative = 1`,
  neg_1yr = clin$`MRD by NGS 1yr f/u`, neg_2yr = clin$`MRD by NGS 2yr f/u`,
  neg_3yr = clin$`MRD by NGS 3yr f/u`, stringsAsFactors = FALSE)
mrd_df$neg_3yrFU <- ifelse(!is.na(mrd_df$neg_3yr), mrd_df$neg_3yr,
                    ifelse(!is.na(mrd_df$neg_2yr), mrd_df$neg_2yr, mrd_df$neg_1yr))
to_pos <- function(neg) ifelse(is.na(neg), NA_integer_, as.integer(1 - neg))

mrd_long <- bind_rows(
  data.frame(Study.ID = mrd_df$Study.ID, t_num = 3L, MRD_pos = to_pos(mrd_df$neg_c8)),
  data.frame(Study.ID = mrd_df$Study.ID, t_num = 4L, MRD_pos = to_pos(mrd_df$neg_c18)),
  data.frame(Study.ID = mrd_df$Study.ID, t_num = 5L, MRD_pos = to_pos(mrd_df$neg_3yrFU)))
mrd_long$Timepoint <- factor(tp_levels[mrd_long$t_num], levels = tp_levels)

resp_df <- mrd_df %>% distinct(Study.ID, best_response)

## MRD points to overlay on each patient panel (placed at w1cf height for that cycle)
mrd_pts <- obs_df %>%
  inner_join(mrd_long %>% select(Study.ID, t_num, MRD_pos), by = c("Study.ID", "t_num")) %>%
  filter(!is.na(MRD_pos)) %>%
  mutate(MRD = factor(ifelse(MRD_pos == 1, "MRD+", "MRD-"), levels = c("MRD-", "MRD+")))

## Facet label: "<ID>  (<best response>)"
panel_lab <- obs_df %>% distinct(Study.ID) %>%
  left_join(resp_df, by = "Study.ID") %>%
  mutate(panel = ifelse(is.na(best_response),
                        as.character(Study.ID),
                        sprintf("%s (%s)", Study.ID, best_response)))
obs_df  <- left_join(obs_df,  panel_lab, by = "Study.ID")
mrd_pts <- left_join(mrd_pts, panel_lab %>% select(Study.ID, panel, best_response),
                     by = "Study.ID")

## ---------------------------------------------------------------------------
## 3. Small-multiples: dominant-subclone weight trajectory per patient
## ---------------------------------------------------------------------------
long_w <- obs_df %>%
  select(panel, Study.ID, Timepoint, w1g, w1cf) %>%
  pivot_longer(c(w1g, w1cf), names_to = "Source", values_to = "w1") %>%
  mutate(Source = recode(Source, w1g = "gDNA", w1cf = "cfDNA"))

p <- ggplot(long_w, aes(Timepoint, w1, group = Source,
                        color = Source, linetype = Source)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.1) +
  geom_point(data = mrd_pts, aes(Timepoint, w1cf, shape = MRD),
             inherit.aes = FALSE, size = 2.2, stroke = 0.7) +
  scale_color_manual(values = c(cfDNA = "#E63946", gDNA = "#457B9D")) +
  scale_linetype_manual(values = c(cfDNA = "solid", gDNA = "dashed")) +
  scale_shape_manual(values = c("MRD-" = 1, "MRD+" = 4), name = "MRD by NGS") +
  facet_wrap(~ panel, ncol = 5) +
  labs(x = NULL, y = expression("Dominant-subclone weight " * omega[1]),
       title = "Per-patient dominant-subclone trajectory (cfDNA solid, gDNA dashed), annotated by response and MRD") +
  ylim(0, 1) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.text = element_text(size = 6.5))

n_pat <- n_distinct(obs_df$Study.ID)
ggsave(file.path(fig_dir, "fig_subclonal_story.pdf"), p,
       width = 11, height = 2.2 * ceiling(n_pat / 5), limitsize = FALSE)
cat("Saved: nature_manuscript/figures/fig_subclonal_story.pdf\n")

## ---------------------------------------------------------------------------
## 4. Narrative table: illustrative patients
##    - per-patient dominant-weight shift (max - min over cycles) in each source
##    - MRD trajectory summary + best response
##    - flag "conversion" patients (MRD status changes across measured cycles)
## ---------------------------------------------------------------------------
shift_tbl <- obs_df %>% group_by(Study.ID, panel, best_response) %>%
  summarise(nvis = n(),
            w1g_shift  = ifelse(n() > 1, max(w1g)  - min(w1g),  NA_real_),
            w1cf_shift = ifelse(n() > 1, max(w1cf) - min(w1cf), NA_real_),
            w1g_scr  = w1g[which.min(t_num)],  w1g_last  = w1g[which.max(t_num)],
            w1cf_scr = w1cf[which.min(t_num)], w1cf_last = w1cf[which.max(t_num)],
            .groups = "drop")

mrd_traj <- mrd_long %>% filter(!is.na(MRD_pos)) %>%
  arrange(Study.ID, t_num) %>%
  group_by(Study.ID) %>%
  summarise(mrd_cycles = paste(tp_levels[t_num], collapse = ","),
            mrd_seq    = paste(ifelse(MRD_pos == 1, "+", "-"), collapse = ""),
            mrd_conversion = as.integer(n_distinct(MRD_pos) > 1),
            .groups = "drop")

story <- shift_tbl %>%
  left_join(mrd_traj, by = "Study.ID") %>%
  mutate(max_shift = pmax(w1g_shift, w1cf_shift, na.rm = TRUE)) %>%
  arrange(desc(mrd_conversion), desc(max_shift))

write.csv(story, file.path(res_dir, "subclonal_story_illustrative.csv"),
          row.names = FALSE)
cat("Saved: results/subclonal_story_illustrative.csv\n\n")

cat("==== Illustrative patients (MRD conversions and/or largest subclonal shift) ====\n")
ill <- story %>%
  filter(mrd_conversion == 1 | max_shift >= quantile(max_shift, 0.75, na.rm = TRUE)) %>%
  select(Study.ID, best_response, nvis, w1g_shift, w1cf_shift,
         mrd_cycles, mrd_seq, mrd_conversion)
print(as.data.frame(ill), digits = 3, row.names = FALSE)

sink(file.path(res_dir, "sessionInfo_subclonal_story.txt")); print(sessionInfo()); sink()
cat("\nDone.\n")
