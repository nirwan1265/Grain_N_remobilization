################################################################################
### INDIAN JARVIS LOW-COVERAGE SEQUENCING QC FIGURES
################################################################################

library(tidyverse)
library(scales)
library(patchwork)

################################################################################
### CONFIGURATION
################################################################################

qc_dir <- "/Users/nirwantandukar/Documents/Research/data/Indian_Jarvis/qc_stats"
output_dir <- "Figs/Supplementary"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

variant_file <- file.path(output_dir, "SuppFig_seq_statistics_variant_level.png")
sample_file <- file.path(output_dir, "SuppFig_seq_statistics_sample_level.png")

plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 10)
    ),
    plot.tag = element_text(size = 18, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    legend.background = element_rect(fill = "white", color = "grey70", linewidth = 0.4),
    legend.title = element_blank(),
    legend.text = element_text(size = 13),
    plot.margin = margin(15, 15, 15, 15)
  )

save_png <- function(plot_obj, file_path, width, height) {
  ggsave(
    filename = file_path,
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 300,
    bg = "white",
    device = "png"
  )
}

################################################################################
### VARIANT-LEVEL QC
################################################################################

dp <- readr::read_table(
  file.path(qc_dir, "all_dp_clean.txt"),
  col_names = "depth",
  show_col_types = FALSE
) %>%
  dplyr::filter(depth <= 10)

p_dp <- ggplot(dp, aes(x = depth)) +
  geom_histogram(binwidth = 1, fill = "#2C7FB8", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#D55E00", linewidth = 0.5) +
  scale_x_continuous(breaks = 0:10) +
  scale_y_continuous(labels = comma, trans = "log10") +
  labs(
    title = "Variant Depth Distribution",
    x = "Depth (x)",
    y = "Count (log10)"
  ) +
  plot_theme

af_df <- readr::read_table(
  file.path(qc_dir, "af_bins_final.txt"),
  col_names = c("bin_max", "site_count"),
  show_col_types = FALSE
) %>%
  dplyr::filter(bin_max >= 0.01) %>%
  dplyr::mutate(
    af_bin = factor(
      percent(bin_max, accuracy = 0.1),
      levels = percent(bin_max, accuracy = 0.1)
    )
  )

p_af <- ggplot(af_df, aes(x = af_bin, y = site_count)) +
  geom_col(fill = "#E69F00", width = 0.8) +
  scale_y_continuous(labels = label_scientific()) +
  labs(
    title = "Minor Allele Frequency",
    x = "Allele Frequency Bin",
    y = "Number of Variants"
  ) +
  plot_theme

miss_df <- readr::read_table(
  file.path(qc_dir, "missing_counts.tsv"),
  col_names = c("sample", "missing_gt"),
  show_col_types = FALSE
) %>%
  dplyr::arrange(dplyr::desc(missing_gt)) %>%
  dplyr::mutate(
    sample_rank = dplyr::row_number(),
    outlier = missing_gt >= quantile(missing_gt, 0.95)
  )

p_miss <- ggplot(miss_df, aes(x = sample_rank, y = missing_gt, fill = outlier)) +
  geom_col(width = 1) +
  scale_fill_manual(values = c("FALSE" = "#009E73", "TRUE" = "#D55E00"), guide = "none") +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Missing Genotypes per Sample",
    x = "Sample Rank",
    y = "Missing Genotypes (./.)"
  ) +
  plot_theme

titv_chr <- readr::read_table(
  file.path(qc_dir, "titv_per_chr.tsv"),
  col_names = c("raw_chr", "titv"),
  show_col_types = FALSE
) %>%
  dplyr::mutate(chr = factor(seq_len(dplyr::n())))

p_titv <- ggplot(titv_chr, aes(x = chr, y = titv)) +
  geom_col(fill = "#8E63B0", width = 0.8) +
  geom_hline(
    yintercept = mean(titv_chr$titv, na.rm = TRUE),
    linetype = "dashed",
    color = "grey35",
    linewidth = 0.45
  ) +
  labs(
    title = "Ti/Tv Ratio by Chromosome",
    x = "Chromosome",
    y = "Ti/Tv Ratio"
  ) +
  plot_theme

variant_plot <- (p_dp + p_af) / (p_miss + p_titv) +
  plot_annotation(tag_levels = "A")

save_png(variant_plot, variant_file, width = 12, height = 9)

################################################################################
### SAMPLE-LEVEL QC
################################################################################

depth_df <- readr::read_tsv(
  file.path(qc_dir, "per_sample_mean_depth.tsv"),
  col_names = c("sample", "mean_depth", "n_genotypes"),
  show_col_types = FALSE
) %>%
  dplyr::mutate(
    coverage_group = cut(
      mean_depth,
      breaks = c(0, 0.5, 1, 1.5, Inf),
      labels = c("<0.5x", "0.5-1x", "1-1.5x", ">1.5x"),
      include.lowest = TRUE
    )
  )

coverage_colors <- c(
  "<0.5x" = "#56B4E9",
  "0.5-1x" = "#E69F00",
  "1-1.5x" = "#CC79A7",
  ">1.5x" = "#009E73"
)

global_avg_depth <- mean(depth_df$mean_depth, na.rm = TRUE)

p_density <- ggplot(depth_df, aes(x = mean_depth)) +
  geom_histogram(binwidth = 0.05, fill = "#4C78A8", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "#D55E00", linewidth = 0.5) +
  geom_vline(xintercept = global_avg_depth, linetype = "solid", color = "grey35", linewidth = 0.45) +
  labs(
    title = "Sample Mean Depth",
    x = "Mean Depth (x)",
    y = "Number of Samples"
  ) +
  plot_theme

ranked_df <- depth_df %>%
  dplyr::arrange(dplyr::desc(mean_depth)) %>%
  dplyr::mutate(sample_rank = dplyr::row_number())

p_ranked <- ggplot(ranked_df, aes(x = sample_rank, y = mean_depth, fill = coverage_group)) +
  geom_col(width = 1) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "#D55E00", linewidth = 0.5) +
  scale_fill_manual(values = coverage_colors, name = "Coverage Group") +
  labs(
    title = "Ranked Sample Depth",
    x = "Sample Rank",
    y = "Mean Depth (x)"
  ) +
  plot_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

sample_plot <- (p_density / p_ranked) +
  plot_layout(guides = "collect", heights = c(1, 1.15)) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "top")

save_png(sample_plot, sample_file, width = 10, height = 9)

cat("\nSaved variant-level QC figure to:\n")
cat("  ", variant_file, "\n", sep = "")
cat("Saved sample-level QC figure to:\n")
cat("  ", sample_file, "\n", sep = "")
