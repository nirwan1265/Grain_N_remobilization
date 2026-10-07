################################################################################
### SUPPLEMENTARY FIGURE: GOODMAN-BUCKLER FAA DIVERSITY
################################################################################

library(tidyverse)
library(ggplot2)
library(patchwork)
library(scales)

################################################################################
### CONFIGURATION
################################################################################

faa_file <- "data/FAA_Goodman_Buckler.csv"
output_dir <- "Figs/Supplementary"
output_file <- file.path(output_dir, "SuppFig_Goodman_Buckler_FAA_diversity.png")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

faa_cols <- c(
  "A", "R", "N", "D", "Q", "E", "G", "H", "I", "L",
  "K", "M", "F", "P", "S", "W", "T", "Y", "V", "C"
)

aa_labels <- c(
  A = "Ala", R = "Arg", N = "Asn", D = "Asp", Q = "Gln", E = "Glu",
  G = "Gly", H = "His", I = "Ile", L = "Leu", K = "Lys", M = "Met",
  F = "Phe", P = "Pro", S = "Ser", W = "Trp", T = "Thr", Y = "Tyr",
  V = "Val", C = "Cys"
)

highlight_colors <- c(
  "Asparagine" = "#9C2F2F",
  "Proline" = "#C77C1A",
  "Other FAA" = "#4A78A6"
)

plot_theme <- theme_minimal(base_size = 22) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.tag = element_text(size = 24, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    axis.title.x = element_text(size = 17, face = "bold"),
    axis.title.y = element_text(size = 17, face = "bold"),
    axis.text.x = element_text(size = 13, face = "bold", color = "black"),
    axis.text.y = element_text(size = 15, face = "bold", color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.7),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.35),
    legend.background = element_rect(fill = "white", color = "grey70", linewidth = 0.4),
    legend.title = element_blank(),
    legend.text = element_text(size = 13),
    plot.margin = margin(14, 14, 14, 14)
  )

################################################################################
### LOAD AND SUMMARISE
################################################################################

faa_raw <- readr::read_csv(faa_file, show_col_types = FALSE)

faa_taxa <- faa_raw %>%
  group_by(taxa) %>%
  summarise(
    across(all_of(faa_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(Total_FAA = rowSums(across(all_of(faa_cols)), na.rm = TRUE))

mean_order <- faa_taxa %>%
  summarise(across(all_of(faa_cols), ~ mean(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "aa", values_to = "mean_value") %>%
  arrange(desc(mean_value)) %>%
  pull(aa)

abs_long <- faa_taxa %>%
  select(taxa, all_of(faa_cols)) %>%
  pivot_longer(
    cols = -taxa,
    names_to = "aa",
    values_to = "value"
  ) %>%
  mutate(
    aa = factor(aa, levels = mean_order),
    aa_label = factor(aa_labels[as.character(aa)], levels = aa_labels[mean_order]),
    highlight = case_when(
      as.character(aa) == "N" ~ "Asparagine",
      as.character(aa) == "P" ~ "Proline",
      TRUE ~ "Other FAA"
    )
  )

rel_long <- faa_taxa %>%
  mutate(across(all_of(faa_cols), ~ .x / Total_FAA)) %>%
  select(taxa, all_of(faa_cols)) %>%
  pivot_longer(
    cols = -taxa,
    names_to = "aa",
    values_to = "fraction"
  )

variability_df <- rel_long %>%
  group_by(aa) %>%
  summarise(
    mean_fraction = mean(fraction, na.rm = TRUE),
    sd_fraction = sd(fraction, na.rm = TRUE),
    iqr_fraction = IQR(fraction, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(sd_fraction)) %>%
  slice_head(n = 10) %>%
  mutate(
    aa = factor(aa, levels = aa),
    aa_label = factor(aa_labels[as.character(aa)], levels = aa_labels[as.character(aa)]),
    highlight = case_when(
      as.character(aa) == "N" ~ "Asparagine",
      as.character(aa) == "P" ~ "Proline",
      TRUE ~ "Other FAA"
    )
  )

################################################################################
### PLOT
################################################################################

p_dist <- ggplot(abs_long, aes(x = aa_label, y = value, fill = highlight, color = highlight)) +
  geom_violin(
    alpha = 0.20,
    linewidth = 0.7,
    trim = FALSE
  ) +
  geom_boxplot(
    width = 0.18,
    alpha = 0.28,
    linewidth = 0.8,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    size = 0.65,
    alpha = 0.18
  ) +
  scale_fill_manual(values = highlight_colors) +
  scale_color_manual(values = highlight_colors) +
  scale_y_log10(labels = label_number(accuracy = 0.1)) +
  labs(
    title = "FAA abundance distributions across 280 genotypes",
    x = NULL,
    y = "Genotype mean FAA abundance",
    tag = "A"
  ) +
  plot_theme +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_var <- ggplot(variability_df, aes(x = aa_label, y = sd_fraction, fill = highlight, color = highlight)) +
  geom_col(width = 0.72, linewidth = 0.7, alpha = 0.92) +
  geom_text(
    aes(label = percent(sd_fraction, accuracy = 0.1)),
    vjust = -0.35,
    size = 4.2,
    fontface = "bold",
    color = "black"
  ) +
  scale_fill_manual(values = highlight_colors) +
  scale_color_manual(values = highlight_colors) +
  scale_y_continuous(
    labels = label_percent(accuracy = 1),
    expand = expansion(mult = c(0, 0.14))
  ) +
  labs(
    title = "Top 10 most variable FAA by relative composition",
    x = NULL,
    y = "SD of genotype-level FAA fraction",
    tag = "B"
  ) +
  plot_theme +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

final_plot <- p_dist / p_var +
  plot_layout(heights = c(1, 1.2))

ggsave(output_file, final_plot, width = 16, height = 12, dpi = 300, bg = "white")

cat("\nSaved Goodman-Buckler FAA diversity figure to:\n")
cat("  ", output_file, "\n", sep = "")
