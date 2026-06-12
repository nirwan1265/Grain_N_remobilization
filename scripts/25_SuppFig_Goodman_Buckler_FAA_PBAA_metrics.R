################################################################################
### SUPPLEMENTARY FIGURE: GOODMAN-BUCKLER FAA / PBAA PANEL METRICS
################################################################################

library(tidyverse)
library(ggplot2)
library(patchwork)

################################################################################
### CONFIGURATION
################################################################################

taxa_means_file <- "tables/supplementary/SuppTable_12_Goodman_Buckler_FAA_PBAA_taxa_means.csv"
paired_tests_file <- "tables/supplementary/SuppTable_14_Goodman_Buckler_FAA_PBAA_paired_tests.csv"
output_dir <- "Figs/Supplementary"
output_file <- file.path(output_dir, "SuppFig_Goodman_Buckler_FAA_PBAA_metrics.png")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

group_colors <- c(
  "FAA" = "#2166AC",
  "PBAA" = "#B35806"
)

metric_specs <- tribble(
  ~metric, ~source_group, ~panel_label, ~y_label, ~plot_tag, ~use_log10,
  "Total_FAA", "FAA", "Broad amino acid pool", "Assay-scale abundance", "A", TRUE,
  "Total_PBAA", "PBAA", "Broad amino acid pool", "Assay-scale abundance", "A", TRUE,
  "FAA_N_proxy", "FAA", "Estimated amino-acid-derived N proxy", "Estimated N units (assay scale)", "B", TRUE,
  "PBAA_N_proxy", "PBAA", "Estimated amino-acid-derived N proxy", "Estimated N units (assay scale)", "B", TRUE,
  "Proline_fraction_FAA", "FAA", "Proline fraction", "Fraction of total pool", "C", FALSE,
  "Proline_fraction_PBAA", "PBAA", "Proline fraction", "Fraction of total pool", "C", FALSE,
  "N_rich_FAA_fraction", "FAA", "Nitrogen-rich fraction", "Fraction of total pool", "D", FALSE,
  "N_rich_PBAA_fraction", "PBAA", "Nitrogen-rich fraction", "Fraction of total pool", "D", FALSE
)

plot_theme <- theme_minimal(base_size = 22) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.tag = element_text(size = 24, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 17, face = "bold"),
    axis.text.x = element_text(size = 15, face = "bold", color = "black"),
    axis.text.y = element_text(size = 15, face = "bold", color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.7),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.35),
    legend.position = "none",
    plot.margin = margin(14, 14, 14, 14)
  )

################################################################################
### LOAD
################################################################################

taxa_df <- readr::read_csv(taxa_means_file, show_col_types = FALSE) %>%
  filter(matched_FAA_PBAA)

paired_df <- readr::read_csv(paired_tests_file, show_col_types = FALSE)

plot_df <- taxa_df %>%
  select(
    taxa,
    Total_FAA,
    Total_PBAA,
    FAA_N_proxy,
    PBAA_N_proxy,
    Proline_fraction_FAA,
    Proline_fraction_PBAA,
    N_rich_FAA_fraction,
    N_rich_PBAA_fraction
  ) %>%
  pivot_longer(
    cols = -taxa,
    names_to = "metric",
    values_to = "value"
  ) %>%
  left_join(metric_specs, by = "metric") %>%
  mutate(source_group = factor(source_group, levels = c("FAA", "PBAA")))

annotation_df <- paired_df %>%
  transmute(
    panel_label = case_when(
      comparison_label == "Broad amino acid pool: FAA vs PBAA" ~ "Broad amino acid pool",
      comparison_label == "Estimated amino-acid-derived nitrogen proxy: FAA vs PBAA" ~ "Estimated amino-acid-derived N proxy",
      comparison_label == "Proline fraction: FAA vs PBAA" ~ "Proline fraction",
      comparison_label == "Nitrogen-rich fraction: FAA vs PBAA" ~ "Nitrogen-rich fraction"
    ),
    label = sig
  ) %>%
  left_join(
    plot_df %>%
      group_by(panel_label) %>%
      summarise(y_max = max(value, na.rm = TRUE), .groups = "drop"),
    by = "panel_label"
  ) %>%
  mutate(
    x_start = 1,
    x_end = 2,
    y = y_max * 1.10
  )

################################################################################
### PLOT
################################################################################

make_panel_plot <- function(panel_name, tag_name, y_axis_label, use_log) {
  metric_dat <- plot_df %>%
    filter(panel_label == panel_name)

  metric_annot <- annotation_df %>%
    filter(panel_label == panel_name)

  p <- ggplot(metric_dat, aes(x = source_group, y = value, fill = source_group, color = source_group)) +
    geom_violin(
      width = 0.9,
      alpha = 0.20,
      linewidth = 0.9,
      trim = FALSE
    ) +
    geom_boxplot(
      width = 0.28,
      alpha = 0.25,
      linewidth = 0.9,
      outlier.shape = NA
    ) +
    geom_jitter(
      width = 0.11,
      size = 1.5,
      alpha = 0.35
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
      aes(x = 1.5, y = y * 1.01, label = label),
      inherit.aes = FALSE,
      size = 5.0,
      fontface = "bold"
    ) +
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors) +
    labs(
      title = panel_name,
      y = y_axis_label,
      tag = tag_name
    ) +
    plot_theme

  if (use_log) {
    p <- p +
      scale_y_log10(expand = expansion(mult = c(0.03, 0.18)))
  } else {
    p <- p +
      scale_y_continuous(expand = expansion(mult = c(0.03, 0.18)))
  }

  p
}

p1 <- make_panel_plot("Broad amino acid pool", "A", "Assay-scale abundance", TRUE)
p2 <- make_panel_plot("Estimated amino-acid-derived N proxy", "B", "Estimated N units (assay scale)", TRUE)
p3 <- make_panel_plot("Proline fraction", "C", "Fraction of total pool", FALSE)
p4 <- make_panel_plot("Nitrogen-rich fraction", "D", "Fraction of total pool", FALSE)

supp_fig <- (p1 + p2) / (p3 + p4)

ggsave(output_file, supp_fig, width = 16, height = 11, dpi = 300, bg = "white")

cat("\nSaved supplementary figure to:\n")
cat("  ", output_file, "\n", sep = "")
