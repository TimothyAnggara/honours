# =============================================================================
# figures_results.R
#
# Generates all Results chapter figures for the thesis.
# Assumes you have already run the benchmarking pipeline and saved:
#   results/simulation/sensitivity_single_metrics.rds
#   results/simulation/sensitivity_repeated_metrics.rds
#
# For the Keren figure, paste in the kerenBench$raw object from your keren.Rmd.
#
# Output: results/figures/  (PDF + PNG for each figure)
# =============================================================================

library(tidyverse)
library(patchwork)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# ── Consistent method ordering and colours ────────────────────────────────────
methodOrder = c("TreeBH", "DART", "pfilter", "Ours", "ClusterTrim_BY", "Flat BH", "Flat BY")

methodColours = c(
  "TreeBH"        = "#534AB7",
  "DART"          = "#D85A30",
  "pfilter"       = "#1D9E75",
  "Ours"          = "#BA7517",
  "ClusterTrim_BY"= "#D4537E",
  "Flat BH"       = "#888780",
  "Flat BY"       = "#B4B2A9"
)

# Shared theme
plotTheme = theme_minimal(base_size = 11) +
  theme(
    legend.position  = "right",
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold"),
    axis.title       = element_text(face = "bold"),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, colour = "grey40"),
    legend.title     = element_text(face = "bold", size = 9)
  )

# Helper to save both PDF and PNG
saveFig = function(plot, name, width = 8, height = 5) {
  ggsave(sprintf("results/figures/%s.pdf", name), plot, width = width, height = height)
  ggsave(sprintf("results/figures/%s.png", name), plot, width = width, height = height, dpi = 300)
  message("Saved: ", name)
}

# ── Load data ─────────────────────────────────────────────────────────────────
sensitivitySingleMetrics   = readRDS("results/simulation/sensitivity_single_metrics.rds")
sensitivityRepeatedMetrics = readRDS("results/simulation/sensitivity_repeated_metrics.rds")

# ── Summarise ─────────────────────────────────────────────────────────────────
summarise_metrics = function(df) {
  df %>%
    group_by(method, signal_weight) %>%
    summarise(
      mean_TPR = mean(TPR, na.rm = TRUE),
      sd_TPR   = sd(TPR,   na.rm = TRUE),
      mean_FDR = mean(FDR, na.rm = TRUE),
      sd_FDR   = sd(FDR,   na.rm = TRUE),
      n_reps   = n(),
      .groups  = "drop"
    ) %>%
    mutate(
      se_TPR = sd_TPR / sqrt(n_reps),
      se_FDR = sd_FDR / sqrt(n_reps),
      method = factor(method, levels = methodOrder),
      signal_weight = as.numeric(signal_weight)
    )
}

singleSummary   = summarise_metrics(sensitivitySingleMetrics)
repeatedSummary = summarise_metrics(sensitivityRepeatedMetrics)

# =============================================================================
# Figure 1: Power vs signal weight — single measure
# =============================================================================
figPowerSingle = ggplot(singleSummary,
                        aes(x = signal_weight, y = mean_TPR, colour = method, group = method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  geom_errorbar(
    aes(ymin = pmax(mean_TPR - se_TPR, 0), ymax = pmin(mean_TPR + se_TPR, 1)),
    width = 0.025, alpha = 0.5
  ) +
  scale_x_continuous(breaks = sort(unique(singleSummary$signal_weight))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_colour_manual(values = methodColours, name = "Method") +
  labs(
    title    = "Power across increasing signal strength — single measure",
    subtitle = "Mean true positive rate across 20 simulation replicates; error bars show ±1 SE",
    x        = "Signal weight",
    y        = "Mean power"
  ) +
  plotTheme

saveFig(figPowerSingle, "fig_power_single")

# =============================================================================
# Figure 2: FDR vs signal weight — single measure
# =============================================================================
figFDRSingle = ggplot(singleSummary,
                      aes(x = signal_weight, y = mean_FDR, colour = method, group = method)) +
  geom_hline(yintercept = 0.25, linetype = "dashed", linewidth = 0.6, colour = "black") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  geom_errorbar(
    aes(ymin = pmax(mean_FDR - se_FDR, 0), ymax = pmin(mean_FDR + se_FDR, 1)),
    width = 0.025, alpha = 0.5
  ) +
  scale_x_continuous(breaks = sort(unique(singleSummary$signal_weight))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_colour_manual(values = methodColours, name = "Method") +
  labs(
    title    = "False discovery rate across increasing signal strength — single measure",
    subtitle = "Dashed line shows nominal FDR level q = 0.05",
    x        = "Signal weight",
    y        = "Mean FDR"
  ) +
  plotTheme

saveFig(figFDRSingle, "fig_fdr_single")

# =============================================================================
# Figure 3: Power vs signal weight — repeated measure
# =============================================================================
figPowerRepeated = ggplot(repeatedSummary,
                          aes(x = signal_weight, y = mean_TPR, colour = method, group = method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  geom_errorbar(
    aes(ymin = pmax(mean_TPR - se_TPR, 0), ymax = pmin(mean_TPR + se_TPR, 1)),
    width = 0.025, alpha = 0.5
  ) +
  scale_x_continuous(breaks = sort(unique(repeatedSummary$signal_weight))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_colour_manual(values = methodColours, name = "Method") +
  labs(
    title    = "Power across increasing signal strength — repeated measure",
    subtitle = "Mean true positive rate across 20 simulation replicates; error bars show ±1 SE",
    x        = "Signal weight",
    y        = "Mean power"
  ) +
  plotTheme

saveFig(figPowerRepeated, "fig_power_repeated")

# =============================================================================
# Figure 4: FDR vs signal weight — repeated measure
# =============================================================================
figFDRRepeated = ggplot(repeatedSummary,
                        aes(x = signal_weight, y = mean_FDR, colour = method, group = method)) +
  geom_hline(yintercept = 0.25, linetype = "dashed", linewidth = 0.6, colour = "black") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  geom_errorbar(
    aes(ymin = pmax(mean_FDR - se_FDR, 0), ymax = pmin(mean_FDR + se_FDR, 1)),
    width = 0.025, alpha = 0.5
  ) +
  scale_x_continuous(breaks = sort(unique(repeatedSummary$signal_weight))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_colour_manual(values = methodColours, name = "Method") +
  labs(
    title    = "False discovery rate across increasing signal strength — repeated measure",
    subtitle = "Dashed line shows nominal FDR level q = 0.05",
    x        = "Signal weight",
    y        = "Mean FDR"
  ) +
  plotTheme

saveFig(figFDRRepeated, "fig_fdr_repeated")

# =============================================================================
# Figure 5: Power–FDR trade-off — single measure
# =============================================================================
tradeoffSingle = singleSummary %>%
  mutate(swLabel = as.character(signal_weight))

figTradeoffSingle = ggplot(tradeoffSingle,
                           aes(x = mean_FDR, y = mean_TPR, colour = method, label = swLabel)) +
  geom_vline(xintercept = 0.25, linetype = "dashed", linewidth = 0.6, colour = "black") +
  geom_path(aes(group = method), linewidth = 0.5, alpha = 0.4) +
  geom_point(size = 2.8, alpha = 0.9) +
  geom_text(nudge_y = 0.04, size = 2.8, show.legend = FALSE) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_colour_manual(values = methodColours, name = "Method") +
  labs(
    title    = "Power–FDR trade-off — single measure",
    subtitle = "Each point is one method at one signal weight; labels show signal weight",
    x        = "Mean FDR",
    y        = "Mean power"
  ) +
  plotTheme

saveFig(figTradeoffSingle, "fig_tradeoff_single")

# =============================================================================
# Figure 6: Power–FDR trade-off — repeated measure
# =============================================================================
tradeoffRepeated = repeatedSummary %>%
  mutate(swLabel = as.character(signal_weight))

figTradeoffRepeated = ggplot(tradeoffRepeated,
                             aes(x = mean_FDR, y = mean_TPR, colour = method, label = swLabel)) +
  geom_vline(xintercept = 0.25, linetype = "dashed", linewidth = 0.6, colour = "black") +
  geom_path(aes(group = method), linewidth = 0.5, alpha = 0.4) +
  geom_point(size = 2.8, alpha = 0.9) +
  geom_text(nudge_y = 0.04, size = 2.8, show.legend = FALSE) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_colour_manual(values = methodColours, name = "Method") +
  labs(
    title    = "Power–FDR trade-off — repeated measure",
    subtitle = "Each point is one method at one signal weight; labels show signal weight",
    x        = "Mean FDR",
    y        = "Mean power"
  ) +
  plotTheme

saveFig(figTradeoffRepeated, "fig_tradeoff_repeated")

# =============================================================================
# Figure 7: Combined single vs repeated comparison (power + FDR, 2x2)
# =============================================================================
combinedDf = bind_rows(
  singleSummary   %>% mutate(design = "Single measure"),
  repeatedSummary %>% mutate(design = "Repeated measure")
) %>%
  mutate(design = factor(design, levels = c("Single measure", "Repeated measure")))

figCombinedPower = ggplot(combinedDf,
                          aes(x = signal_weight, y = mean_TPR, colour = method, group = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  facet_wrap(~design) +
  scale_x_continuous(breaks = sort(unique(combinedDf$signal_weight))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_colour_manual(values = methodColours, name = "Method") +
  labs(x = "Signal weight", y = "Mean power") +
  plotTheme

figCombinedFDR = ggplot(combinedDf,
                        aes(x = signal_weight, y = mean_FDR, colour = method, group = method)) +
  geom_hline(yintercept = 0.25, linetype = "dashed", linewidth = 0.6, colour = "black") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  facet_wrap(~design) +
  scale_x_continuous(breaks = sort(unique(combinedDf$signal_weight))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_colour_manual(values = methodColours, name = "Method") +
  labs(x = "Signal weight", y = "Mean FDR") +
  plotTheme

figCombined = (figCombinedPower / figCombinedFDR) +
  plot_annotation(
    title    = "Power and FDR across designs and signal strengths",
    subtitle = "Dashed line in FDR panels shows nominal level q = 0.05",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, colour = "grey40")
    )
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

saveFig(figCombined, "fig_combined_designs", width = 10, height = 8)

# =============================================================================
# Figure 8: Bar chart at signal weight 0.99 — single and repeated side by side
# =============================================================================

make_bar = function(df, metric_col, se_col, metric_label, show_dashed = FALSE) {
  p = ggplot(df, aes(x = method, y = .data[[metric_col]], fill = method)) +
    geom_col(width = 0.7) +
    geom_errorbar(
      aes(
        ymin = pmax(.data[[metric_col]] - .data[[se_col]], 0),
        ymax = pmin(.data[[metric_col]] + .data[[se_col]], 1)
      ),
      width = 0.3,
      linewidth = 0.6,
      colour = "grey30"
    ) +
    scale_fill_manual(values = methodColours, guide = "none") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    coord_flip() +
    labs(x = NULL, y = metric_label) +
    plotTheme +
    theme(panel.grid.major.y = element_blank())
  
  if (show_dashed) {
    p = p + geom_hline(yintercept = 0.25, linetype = "dashed",
                       linewidth = 0.6, colour = "black")
  }
  p
}

make_tradeoff = function(df) {
  ggplot(df, aes(x = mean_FDR, y = mean_TPR, colour = method, label = method)) +
    geom_vline(xintercept = 0.25, linetype = "dashed",
               linewidth = 0.6, colour = "black") +
    geom_point(size = 3.5) +
    geom_text(nudge_y = 0.06, size = 2.8, show.legend = FALSE) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_colour_manual(values = methodColours, name = "Method") +
    labs(x = "FDR", y = "Power") +
    plotTheme
}
make_bar = function(df, metric_col, se_col, metric_label, show_dashed = FALSE) {
  p = ggplot(df, aes(x = method, y = .data[[metric_col]], fill = method)) +
    geom_col(width = 0.7) +
    geom_errorbar(
      aes(
        ymin = pmax(.data[[metric_col]] - .data[[se_col]], 0),
        ymax = pmin(.data[[metric_col]] + .data[[se_col]], 1)
      ),
      width = 0.3,
      linewidth = 0.6,
      colour = "grey30"
    ) +
    scale_fill_manual(values = methodColours, guide = "none") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    coord_flip() +
    labs(x = NULL, y = metric_label) +
    plotTheme +
    theme(panel.grid.major.y = element_blank())
  
  if (show_dashed) {
    p = p + geom_hline(yintercept = 0.25, linetype = "dashed",
                       linewidth = 0.6, colour = "black")
  }
  p
}

singleW99   = singleSummary   %>% filter(signal_weight == 0.99)
repeatedW99 = repeatedSummary %>% filter(signal_weight == 0.99)

# Figure 8a: single measure
figBar8a = (
  (make_bar(singleW99, "mean_TPR", "se_TPR", "Power") |
     make_bar(singleW99, "mean_FDR", "se_FDR", "FDR", show_dashed = TRUE)) /
    (plot_spacer() | make_tradeoff(singleW99) | plot_spacer())
) +
  plot_annotation(
    title    = "Method performance at signal weight 0.99 - single measure",
    subtitle = "Dashed lines show nominal FDR level q = 0.25",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, colour = "grey40")
    )
  ) +
  plot_layout(widths = c(0.5, 5, 0.5)) &
  theme(legend.position = "right")

figBar8a

saveFig(figBar8a, "fig_bar_w099_single", width = 10, height = 10)

# Figure 8b: repeated measure
figBar8b = (
  (make_bar(repeatedW99, "mean_TPR", "se_TPR", "Power") |
     make_bar(repeatedW99, "mean_FDR", "se_FDR", "FDR", show_dashed = TRUE)) /
    (plot_spacer() | make_tradeoff(repeatedW99) | plot_spacer())
) +
  plot_annotation(
    title    = "Method performance at signal weight 0.99 - repeated measure",
    subtitle = "Dashed lines show nominal FDR level q = 0.25",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, colour = "grey40")
    )
  ) +
  plot_layout(widths = c(0.5, 5, 0.5)) &
  theme(legend.position = "right")

figBar8b

saveFig(figBar8b, "fig_bar_w099_repeated", width = 10, height = 10)

# =============================================================================
# Figure A: Single measure sensitivity — power | FDR | trade-off
# =============================================================================
figSensitivitySingle = (
  (figSinglePower + guides(colour = "none") |
     figSingleFDR  + guides(colour = "none")) /
    (plot_spacer() | figSingleTradeoff | plot_spacer()) +
    plot_layout(widths = c(0.5, 3, 0.5))
) +
  plot_annotation(
    title    = "Sensitivity analysis: single measure",
    subtitle = "Rows: power (top-left), FDR (top-right), power–FDR trade-off (bottom); dashed lines show q = 0.05; error bars ±1 SE",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 9, colour = "grey40")
    )
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

saveFig(figSensitivitySingle, "fig_sensitivity_single", width = 12, height = 10)

# =============================================================================
# Figure B: Repeated measure sensitivity — power | FDR | trade-off
# =============================================================================
figSensitivityRepeated = (
  (figRepeatedPower + guides(colour = "none") |
     figRepeatedFDR  + guides(colour = "none")) /
    (plot_spacer() | figRepeatedTradeoff | plot_spacer()) +
    plot_layout(widths = c(0.5, 3, 0.5))
) +
  plot_annotation(
    title    = "Sensitivity analysis: repeated measure",
    subtitle = "Rows: power (top-left), FDR (top-right), power–FDR trade-off (bottom); dashed lines show q = 0.05; error bars ±1 SE",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 9, colour = "grey40")
    )
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

saveFig(figSensitivityRepeated, "fig_sensitivity_repeated", width = 12, height = 10)

# =============================================================================
# Figure 9: Keren et al. rejection heatmap
# Paste in kerenBench$raw from your keren.Rmd before running this section
# =============================================================================

# ── Replace this line with your actual object ─────────────────────────────────
kerenRaw = kerenBench$raw

kerenHeatmap = kerenRaw %>%
  mutate(method = factor(method, levels = methodOrder)) %>%
  group_by(method) %>%
  mutate(
    family = sub("__.*", "", leaf),
    reject = as.integer(reject)
  ) %>%
  ungroup() %>%
  group_by(leaf) %>%
  filter(any(reject == 1)) %>%       # keep only pairs rejected by at least one method
  ungroup() %>%
  ggplot(aes(x = leaf, y = method, fill = factor(reject))) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_manual(
    values = c("0" = "#F1EFE8", "1" = "#534AB7"),
    labels = c("Not rejected", "Rejected"),
    name   = NULL
  ) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  labs(
    title    = "Rejected interactions — Keren et al. (2018) TNBC dataset",
    subtitle = "Cold vs compartmentalised tumour-immune archetypes; q = 0.05",
    x        = "Cell-type pair",
    y        = NULL
  ) +
  plotTheme +
  theme(
    axis.text.x  = element_text(size = 6),
    panel.grid   = element_blank()
  )

saveFig(kerenHeatmap, "fig_keren_heatmap", width = 14, height = 5)

message("All figures saved to results/figures/")

# =============================================================================
# Figure A: Single measure sensitivity — power | FDR | trade-off
# =============================================================================

figSinglePower = ggplot(
  singleSummary,
  aes(x = signal_weight, y = mean_TPR, colour = method, group = method)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  geom_errorbar(
    aes(ymin = pmax(mean_TPR - se_TPR, 0), ymax = pmin(mean_TPR + se_TPR, 1)),
    width = 0.025, alpha = 0.5
  ) +
  scale_x_continuous(breaks = sort(unique(singleSummary$signal_weight))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_colour_manual(values = methodColours, name = "Method") +
  labs(x = "Signal weight", y = "Mean power") +
  plotTheme

figSingleFDR = ggplot(
  singleSummary,
  aes(x = signal_weight, y = mean_FDR, colour = method, group = method)
) +
  geom_hline(yintercept = 0.25, linetype = "dashed", linewidth = 0.6, colour = "black") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  geom_errorbar(
    aes(ymin = pmax(mean_FDR - se_FDR, 0), ymax = pmin(mean_FDR + se_FDR, 1)),
    width = 0.025, alpha = 0.5
  ) +
  scale_x_continuous(breaks = sort(unique(singleSummary$signal_weight))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_colour_manual(values = methodColours, name = "Method") +
  labs(x = "Signal weight", y = "Mean FDR") +
  plotTheme

figSingleTradeoff = ggplot(
  singleSummary %>% mutate(swLabel = as.character(signal_weight)),
  aes(x = mean_FDR, y = mean_TPR, colour = method, label = swLabel)
) +
  geom_vline(xintercept = 0.25, linetype = "dashed", linewidth = 0.6, colour = "black") +
  geom_path(aes(group = method), linewidth = 0.5, alpha = 0.4) +
  geom_point(size = 2.2, alpha = 0.9) +
  geom_text(nudge_y = 0.045, size = 2.5, show.legend = FALSE) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_colour_manual(values = methodColours, name = "Method") +
  labs(x = "Mean FDR", y = "Mean power") +
  plotTheme

figSensitivitySingle = (
  (figSinglePower + guides(colour = "none") |
     figSingleFDR  + guides(colour = "none")) /
    (plot_spacer() | figSingleTradeoff | plot_spacer()) +
    plot_layout(widths = c(0.5, 3, 0.5))
) +
  plot_annotation(
    title    = "Sensitivity analysis: single measure",
    subtitle = "Rows: power (top-left), FDR (top-right), power–FDR trade-off (bottom); dashed lines show theoretical minimum FDR of 0.25; error bars ±1 SE",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 9, colour = "grey40")
    )
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

saveFig(figSensitivitySingle, "fig_sensitivity_single", width = 16, height = 13)

# =============================================================================
# Figure B: Repeated measure sensitivity — power | FDR | trade-off
# =============================================================================

figRepeatedPower = ggplot(
  repeatedSummary,
  aes(x = signal_weight, y = mean_TPR, colour = method, group = method)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  geom_errorbar(
    aes(ymin = pmax(mean_TPR - se_TPR, 0), ymax = pmin(mean_TPR + se_TPR, 1)),
    width = 0.025, alpha = 0.5
  ) +
  scale_x_continuous(breaks = sort(unique(repeatedSummary$signal_weight))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_colour_manual(values = methodColours, name = "Method") +
  labs(x = "Signal weight", y = "Mean power") +
  plotTheme

figRepeatedFDR = ggplot(
  repeatedSummary,
  aes(x = signal_weight, y = mean_FDR, colour = method, group = method)
) +
  geom_hline(yintercept = 0.25, linetype = "dashed", linewidth = 0.6, colour = "black") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  geom_errorbar(
    aes(ymin = pmax(mean_FDR - se_FDR, 0), ymax = pmin(mean_FDR + se_FDR, 1)),
    width = 0.025, alpha = 0.5
  ) +
  scale_x_continuous(breaks = sort(unique(repeatedSummary$signal_weight))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_colour_manual(values = methodColours, name = "Method") +
  labs(x = "Signal weight", y = "Mean FDR") +
  plotTheme

figRepeatedTradeoff = ggplot(
  repeatedSummary %>% mutate(swLabel = as.character(signal_weight)),
  aes(x = mean_FDR, y = mean_TPR, colour = method, label = swLabel)
) +
  geom_vline(xintercept = 0.25, linetype = "dashed", linewidth = 0.6, colour = "black") +
  geom_path(aes(group = method), linewidth = 0.5, alpha = 0.4) +
  geom_point(size = 2.2, alpha = 0.9) +
  geom_text(nudge_y = 0.045, size = 2.5, show.legend = FALSE) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_colour_manual(values = methodColours, name = "Method") +
  labs(x = "Mean FDR", y = "Mean power") +
  plotTheme

figSensitivityRepeated = (
  (figRepeatedPower + guides(colour = "none") |
     figRepeatedFDR  + guides(colour = "none")) /
    (plot_spacer() | figRepeatedTradeoff | plot_spacer()) +
    plot_layout(widths = c(0.5, 3, 0.5))
) +
  plot_annotation(
    title    = "Sensitivity analysis: repeated measure",
    subtitle = "Rows: power (top-left), FDR (top-right), power–FDR trade-off (bottom); dashed lines show theoretical minimum FDR of 0.25; error bars ±1 SE",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 9, colour = "grey40")
    )
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

saveFig(figSensitivityRepeated, "fig_sensitivity_repeated", width = 16, height = 13)
