#!/usr/bin/env Rscript
################################################################################
## SUPP FIGURE -- sensitivity of the pathway recurrence test to candidate thresholds
##
## L-homoserine biosynthesis across 12 combinations of GWAS threshold (FDR 5%, 3%, 1%,
## Bonferroni; applied to both the amino acid and soil-N GWAS) and F_ST window threshold
## (top 5%, 3%, 1% in Indian Chief and Jarvis).
## A: one tile per combination. Shade = -log10 pathway-specific permutation P; label = rank of
##    L-homoserine biosynthesis among the 468 CornCyc pathways; * = corrected P < 0.05
##    (min-P correction across all pathways).
## B: number of the 12 combinations in which L-homoserine biosynthesis ranked 1st, 2nd, or lower.
## Permutations: 20,000 per combination; 200,000 for Bonferroni with top 5% and top 3%.
## Input : final/results/pathway_recurrence/sensitivity/pathway_threshold_sensitivity_grid.csv
## Output: final/supp_figs/SuppFig_pathway_threshold_sensitivity.{png,pdf}
################################################################################
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
ROOT <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "."
IN   <- file.path(ROOT, "final", "results", "pathway_recurrence", "sensitivity", "pathway_threshold_sensitivity_grid.csv")
OUT  <- file.path(ROOT, "final", "supp_figs")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.tag     = element_text(size = 24, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x  = element_text(size = 16, colour = "black"),
    axis.text.y  = element_text(size = 16, colour = "black"),
    axis.line    = element_line(colour = "black"),
    panel.grid   = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 15, face = "bold"),
    legend.text  = element_text(size = 14),
    plot.margin  = margin(15, 15, 15, 15)
  )

d <- fread(IN)
d[, GWAS := factor(GWAS_threshold, levels = c("FDR 5%", "FDR 3%", "FDR 1%", "Bonferroni"))]
d[, FST := factor(paste0("Top ", FST_top_pct, "%"), levels = c("Top 1%", "Top 3%", "Top 5%"))]
d[, lp := -log10(homoserine_P)]
d[, lab := paste0("#", homoserine_rank, ifelse(homoserine_P_corrected < 0.05, "*", ""))]
d[, txtcol := ifelse(lp > 1.8, "white", "black")]

## A: tile grid
pA <- ggplot(d, aes(GWAS, FST)) +
  geom_tile(aes(fill = lp), colour = "black", linewidth = 0.8) +
  geom_text(aes(label = lab, colour = txtcol), size = 8, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "black", limits = c(0, 3.6),
                      name = expression(bold(-log[10](italic(P))))) +
  scale_colour_identity() +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = "GWAS threshold", y = expression(bold(italic(F)[ST]~windows))) +
  plot_theme + theme(axis.line = element_blank(), axis.ticks = element_blank())

## B: rank frequency
d[, rank_class := fifelse(homoserine_rank == 1, "1st", fifelse(homoserine_rank == 2, "2nd", "Lower"))]
rc <- d[, .N, by = rank_class][, rank_class := factor(rank_class, levels = c("1st", "2nd", "Lower"))]
rc <- merge(data.table(rank_class = factor(c("1st", "2nd", "Lower"), levels = c("1st", "2nd", "Lower"))), rc, all.x = TRUE)
rc[is.na(N), N := 0]
pB <- ggplot(rc, aes(rank_class, N)) +
  geom_col(fill = "black", width = 0.6) +
  geom_text(aes(label = N), vjust = -0.5, size = 7) +
  scale_y_continuous(limits = c(0, 7), breaks = 0:7, expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Rank of L-homoserine biosynthesis", y = "Threshold combinations") + plot_theme

fig <- (pA | pB) + plot_layout(widths = c(1.6, 1)) + plot_annotation(tag_levels = "A")
for (ext in c("png", "pdf"))
  ggsave(file.path(OUT, paste0("SuppFig_pathway_threshold_sensitivity.", ext)), fig, width = 17, height = 6.5,
         dpi = 300, bg = "white", device = if (ext == "png") ragg::agg_png else "pdf")
print(d[, .(GWAS, FST, homoserine_P, homoserine_P_corrected, homoserine_rank)])
print(rc)
