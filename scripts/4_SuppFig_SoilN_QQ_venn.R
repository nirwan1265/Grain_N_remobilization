################################################################################
##### QQ PLOT WITH GGPLOT2
################################################################################

if (!exists("plot_theme")) {
  plot_theme <- theme_minimal(base_size = 24) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 16, color = "black"),
      axis.text.y = element_text(size = 16, color = "black"),
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.95, 0.95),
      legend.justification = c("right", "top"),
      legend.background = element_rect(fill = "white", color = "grey70", linewidth = 0.4),
      legend.direction = "vertical",
      legend.spacing.y = unit(0.2, "cm"),
      legend.title = element_blank(),
      legend.text = element_text(size = 16),
      plot.margin = margin(15, 15, 15, 15)
    )
}

# Prepare QQ data
n_snps <- nrow(best_per_snp)
expected <- -log10(ppoints(n_snps))
observed <- sort(best_per_snp$logp, decreasing = TRUE)

qq_df <- data.frame(
  expected = expected,
  observed = observed
)

# Calculate 95% confidence interval
conf_int <- data.frame(
  expected = expected,
  lower = -log10(qbeta(0.975, 1:n_snps, n_snps:1 + 1)),
  upper = -log10(qbeta(0.025, 1:n_snps, n_snps:1 + 1))
)

# Add best_model info to QQ data for coloring
# Sort best_per_snp by logp descending to match qq_df order
sorted_models <- best_per_snp %>%
  arrange(desc(logp)) %>%
  pull(best_model)

qq_df <- qq_df %>%
  mutate(best_model = sorted_models)

# QQ plot with points colored by model
p_qq <- ggplot(qq_df, aes(x = expected, y = observed)) +
  # Confidence interval ribbon
  geom_ribbon(data = conf_int, aes(x = expected, ymin = lower, ymax = upper),
              fill = "grey85", alpha = 0.5, inherit.aes = FALSE) +
  # Reference line
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", linewidth = 0.8) +
  # Points colored by best model
  geom_point(aes(color = best_model), alpha = 0.6, size = 1.5) +
  scale_color_manual(
    values = model_colors,
    name = "Best Model",
    limits = names(model_colors)
  ) +
  # Labels
  labs(
    x = expression(Expected~~-log[10](italic(p))),
    y = expression(Observed~~-log[10](italic(p)))
  ) +
  plot_theme +
  theme(
    legend.position = "top"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))





################################################################################
#### VENN DIAGRAM
################################################################################

# Define threshold
# Suggestive threshold: 1e-5 (-log10 p = 5)
# For Bonferroni-style threshold, use 1e-7 (-log10 p = 7)
p_thresh <- 1e-5

# Filter significant SNPs
mlm_sig     <- na.omit(mlm_raw$SNP     [mlm_raw$P.value     < p_thresh])
mlmm_sig    <- na.omit(MLMM_raw$SNP    [MLMM_raw$P.value    < p_thresh])
blink_sig   <- na.omit(BLINK_raw$SNP   [BLINK_raw$P.value   < p_thresh])
farmcpu_sig <- na.omit(farmcpu_raw$SNP[farmcpu_raw$P.value < p_thresh])


# Create named list for input
venn_input <- list(
  MLM     = mlm_sig,
  MLMM    = mlmm_sig,
  BLINK   = blink_sig,
  FarmCPU = farmcpu_sig
)



venn.plot <- venn.diagram(
  x = venn_input,
  filename = NULL,  # Draw to R object instead of file
  output = TRUE,    # Return grob object
  
  # Color settings - use same colors as Manhattan plot
  fill = model_colors[names(venn_input)],
  alpha = 0.65,     # Slightly more transparency for better overlap visibility
  
  # Circle borders
  lwd = 2,          # Thicker border lines
  lty = "solid",    # Solid line style
  col = "white",    # White borders for clean look
  
  # Numeric labels
  cex = 1.6,        # Larger size for numbers
  fontface = "bold",
  
  # Category labels
  cat.cex = 1.5,
  cat.fontface = "bold",
  
  # Margins and scaling
  margin = 0.08,
  height = 6,       # Inches
  width = 6,
  
  # Rotation (if needed)
  rotation.degree = 0,
  
  # No internal title; use panel tags in the combined figure
  main = NULL
)

# Combine QQ and Venn side by side
p_venn <- patchwork::wrap_elements(full = venn.plot)

suppfig_combined <- (p_qq | p_venn) +
  plot_layout(widths = c(1.35, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 18, face = "bold"),
    plot.tag.position = c(0.02, 0.98)
  )

output_file <- "Figs/Supplementary/SuppFig_QQ_venn.png"
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
ggsave(output_file, plot = suppfig_combined, width = 14, height = 7, dpi = 300, bg = "white")

if (interactive()) {
  quartz(width = 14, height = 7)
  print(suppfig_combined)
}

# ggsave("results/figures/SuppFig_QQ_venn.png",
#        plot = suppfig_combined,
#        width = 14, height = 7, dpi = 300, bg = "white")
