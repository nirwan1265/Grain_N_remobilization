################################################################################
### SUPPLEMENTARY FIGURE: ASPARAGINE / PROLINE FAA RANK COUNTS
################################################################################

library(tidyverse)
library(ggplot2)
library(patchwork)
library(scales)

################################################################################
### CONFIGURATION
################################################################################

output_dir <- "Figs/Supplementary"
output_file <- file.path(output_dir, "SuppFig_FAA_Asn_Pro_panel.png")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 8)),
    plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 10)),
    plot.tag = element_text(size = 22, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 15, face = "bold", color = "black"),
    axis.text.y = element_text(size = 15, face = "bold", color = "black"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.3),
    legend.position = "none",
    plot.margin = margin(15, 15, 15, 15)
  )

rank1_colors <- c(
  "Asparagine" = "#2166AC",
  "Proline" = "#B35806",
  "Other amino acids" = "#9E9E9E"
)

pair_colors <- c(
  "Asn + Pro top two" = "#6A3D9A",
  "Other top-two pairings" = "#BDBDBD"
)

################################################################################
### SOURCE COUNTS FROM MANUSCRIPT SUMMARY
################################################################################

n_genotypes <- 280

rank1_df <- tibble(
  category = factor(
    c("Asparagine", "Proline", "Other amino acids"),
    levels = c("Asparagine", "Proline", "Other amino acids")
  ),
  n_count = c(144, 116, 20)
) %>%
  mutate(prop = n_count / n_genotypes)

top2_df <- tibble(
  category = factor(
    c("Asn + Pro top two", "Other top-two pairings"),
    levels = c("Asn + Pro top two", "Other top-two pairings")
  ),
  n_count = c(180, 100)
) %>%
  mutate(prop = n_count / n_genotypes)

################################################################################
### BUILD FIGURE
################################################################################

p_rank1 <- ggplot(rank1_df, aes(x = category, y = n_count, fill = category)) +
  geom_col(width = 0.68, alpha = 0.92) +
  geom_text(
    aes(
      label = paste0(n_count, " (", percent(prop, accuracy = 1), ")")
    ),
    vjust = -0.6,
    size = 5
  ) +
  scale_fill_manual(values = rank1_colors) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.08)),
    breaks = pretty_breaks(n = 5)
  ) +
  labs(
    tag = "A",
    title = "Most Abundant FAA by Genotype",
    subtitle = "Counts of genotypes in which each amino acid ranked first",
    x = NULL,
    y = "Number of Genotypes"
  ) +
  plot_theme

p_top2 <- ggplot(top2_df, aes(x = category, y = n_count, fill = category)) +
  geom_col(width = 0.62, alpha = 0.92) +
  geom_text(
    aes(
      label = paste0(n_count, " (", percent(prop, accuracy = 1), ")")
    ),
    vjust = -0.6,
    size = 5
  ) +
  scale_fill_manual(values = pair_colors) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.08)),
    breaks = pretty_breaks(n = 5)
  ) +
  labs(
    tag = "B",
    title = "Asparagine and Proline as the Top-Two FAA",
    subtitle = "Counts of genotypes where Asn and Pro were the top two, regardless of order",
    x = NULL,
    y = "Number of Genotypes"
  ) +
  plot_theme

combined_plot <- p_rank1 + p_top2 + plot_layout(widths = c(1.2, 1))

################################################################################
### SAVE
################################################################################

ggsave(
  filename = output_file,
  plot = combined_plot,
  width = 13,
  height = 6.8,
  dpi = 300,
  bg = "white"
)

cat("Saved figure:\n  ", output_file, "\n", sep = "")
