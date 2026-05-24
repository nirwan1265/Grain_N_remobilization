################################################################################
### FIGURE 4: RELATIVE FREE AMINO ACID COMPOSITION
################################################################################

library(tidyverse)
library(ggplot2)
library(patchwork)

################################################################################
### CONFIGURATION
################################################################################

main_dir <- "Figs/main"
main_file <- file.path(main_dir, "Fig4.png")

dir.create(main_dir, showWarnings = FALSE, recursive = TRUE)

plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 10)
    ),
    plot.tag = element_text(size = 24, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, face = "bold", color = "black"),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    legend.background = element_rect(fill = "white", color = "grey70", linewidth = 0.4),
    legend.title = element_blank(),
    legend.text = element_text(size = 13),
    plot.margin = margin(15, 15, 15, 15)
  )

source("scripts/13_IC_J_aa_plots.R")

################################################################################
### BUILD FIGURE
################################################################################

cat("\n=== Building Fig4 relative FAA composition panels ===\n")

amino_outputs <- build_amino_profile_outputs(
  base_theme = plot_theme,
  faa_amino_label_size = 7.8,
  faa_sig_label_size = 9.4,
  faa_grid_label_size = 4.8,
  faa_trace_linewidth = 2.4,
  faa_amino_label_face = "bold",
  faa_sig_label_face = "bold",
  faa_plot_title_size = 17
)

fig4 <- (amino_outputs$faa_ic_plot + amino_outputs$faa_j_plot) +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A")

ggsave(main_file, fig4, width = 18, height = 9.5, dpi = 300, bg = "white")

cat("\nSaved main figure to:\n")
cat("  ", main_file, "\n", sep = "")
