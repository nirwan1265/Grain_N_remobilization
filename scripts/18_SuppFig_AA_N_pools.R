################################################################################
### SUPPLEMENTARY FIGURE: AMINO-ACID-DERIVED N POOLS
################################################################################

library(tidyverse)
library(ggplot2)
library(patchwork)

################################################################################
### CONFIGURATION
################################################################################

samples_file <- "tables/supplementary/SuppTable_IC_J_AA_N_biological_means.csv"
summary_file <- "tables/supplementary/SuppTable_IC_J_AA_N_summary.csv"
output_dir <- "Figs/Supplementary"
output_file <- file.path(output_dir, "SuppFig_IC_J_AA_N_pools.png")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

group_colors <- c(
  "Indian Chief G0" = "#92C5DE",
  "Indian Chief G14" = "#2166AC",
  "Jarvis G0" = "#FDB863",
  "Jarvis G14" = "#B35806"
)

metric_labels <- c(
  FAA_N_per_mg = "FAA-derived N per mg",
  PBAA_N_per_mg_est = "Protein-bound AA-derived N per mg",
  Total_AA_N_per_mg_est = "Total AA-derived N per mg"
)

group_levels <- c(
  "Indian Chief G0",
  "Indian Chief G14",
  "Jarvis G0",
  "Jarvis G14"
)

plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.tag = element_text(size = 24, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 15, face = "bold", color = "black", angle = 18, hjust = 1),
    axis.text.y = element_text(size = 16, face = "bold", color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.7),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.35),
    legend.position = "none",
    plot.margin = margin(15, 15, 15, 15)
  )

################################################################################
### LOAD
################################################################################

samples_df <- readr::read_csv(samples_file, show_col_types = FALSE) %>%
  mutate(group = factor(paste(population, generation), levels = group_levels))

summary_df <- readr::read_csv(summary_file, show_col_types = FALSE)

plot_df <- samples_df %>%
  select(population, generation, biological_sample, group, all_of(names(metric_labels))) %>%
  pivot_longer(
    cols = all_of(names(metric_labels)),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(metric_label = factor(metric_labels[metric], levels = metric_labels))

annot_df <- summary_df %>%
  filter(metric %in% names(metric_labels)) %>%
  mutate(
    metric_label = factor(metric_labels[metric], levels = metric_labels),
    label = if_else(raw_sig == "" | is.na(raw_sig), "ns", raw_sig),
    x_start = if_else(population == "Indian Chief", 1, 3),
    x_end = if_else(population == "Indian Chief", 2, 4)
  ) %>%
  left_join(
    plot_df %>%
      group_by(metric_label) %>%
      summarise(y_max = max(value, na.rm = TRUE), .groups = "drop"),
    by = "metric_label"
  ) %>%
  group_by(metric_label) %>%
  mutate(
    y = y_max * c(1.08, 1.18)
  ) %>%
  ungroup()

################################################################################
### PLOT
################################################################################

make_metric_plot <- function(metric_name, letter_tag) {
  metric_dat <- plot_df %>%
    filter(metric == metric_name)

  metric_annot <- annot_df %>%
    filter(metric == metric_name)

  ggplot(metric_dat, aes(x = group, y = value, fill = group, color = group)) +
    geom_boxplot(
      width = 0.62,
      alpha = 0.28,
      linewidth = 0.7,
      outlier.shape = NA
    ) +
    geom_jitter(
      width = 0.11,
      size = 3,
      alpha = 0.9
    ) +
    geom_segment(
      data = metric_annot,
      aes(x = x_start, xend = x_end, y = y, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.7,
      color = "black"
    ) +
    geom_segment(
      data = metric_annot,
      aes(x = x_start, xend = x_start, y = y * 0.985, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.7,
      color = "black"
    ) +
    geom_segment(
      data = metric_annot,
      aes(x = x_end, xend = x_end, y = y * 0.985, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.7,
      color = "black"
    ) +
    geom_text(
      data = metric_annot,
      aes(x = (x_start + x_end) / 2, y = y * 1.01, label = label),
      inherit.aes = FALSE,
      size = 5.2,
      fontface = "bold"
    ) +
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.18))) +
    labs(
      title = metric_labels[[metric_name]],
      y = "Estimated N units per mg",
      tag = letter_tag
    ) +
    plot_theme
}

p1 <- make_metric_plot("FAA_N_per_mg", "A")
p2 <- make_metric_plot("PBAA_N_per_mg_est", "B")
p3 <- make_metric_plot("Total_AA_N_per_mg_est", "C")

supp_fig <- p1 + p2 + p3 + plot_layout(ncol = 3)

ggsave(output_file, supp_fig, width = 20, height = 7.5, dpi = 300, bg = "white")

cat("\nSaved supplementary figure to:\n")
cat("  ", output_file, "\n", sep = "")
