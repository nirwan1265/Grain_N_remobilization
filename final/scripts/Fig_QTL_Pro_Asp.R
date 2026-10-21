suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
g <- fread("final/results/qtl/lod.csv"); res <- jsonlite::fromJSON("final/results/qtl/res.json")
plot_theme <- theme_minimal(base_size = 24) + theme(
  plot.tag = element_text(size = 24, face = "bold"), plot.tag.position = c(0.01, 0.99),
  axis.title.x = element_text(size = 18, face = "bold"), axis.title.y = element_text(size = 18, face = "bold"),
  axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
  axis.line = element_line(colour = "black"), panel.grid = element_blank(), legend.position = "none",
  plot.margin = margin(15, 15, 15, 15))
setorder(g, chr, cM)
g[, cM0 := cM - min(cM), by = chr]
off <- g[, .(len = max(cM0)), by = chr][, start := shift(cumsum(len + 20), fill = 0)]
g <- merge(g, off[, .(chr, start)], by = "chr"); g[, x := cM0 + start]
mid <- g[, .(mid = (min(x) + max(x)) / 2, lo = min(x), hi = max(x)), by = chr]
mk <- g[is_marker == TRUE]
pan <- function(tr, tag) {
  thr <- c(Proline = 3.63, Aspartate = 3.70)[[tr]]  # permutation thresholds reported in the paper (R/qtl, 1,000 permutations)
  ggplot(g, aes(x, get(tr), group = chr)) +
    annotate("rect", xmin = mid[chr == 3, lo] - 5, xmax = mid[chr == 3, hi] + 5, ymin = -Inf, ymax = Inf, fill = "grey88") +
    geom_line(colour = "black", linewidth = 0.6) +
    geom_hline(yintercept = thr, linetype = "dashed", colour = "black", linewidth = 0.6) +
    geom_rug(data = mk, aes(x = x), inherit.aes = FALSE, sides = "b", length = unit(0.015, "npc"), colour = "grey40") +
    scale_x_continuous(breaks = mid$mid, labels = mid$chr, expand = expansion(mult = 0.01)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Chromosome", y = "LOD", tag = tag) + plot_theme
}
fig <- pan("Proline", "A") / pan("Aspartate", "B")
for (e in c("png", "pdf")) ggsave(paste0("final/main_figs/Fig_QTL_Pro_Asp.", e), fig, width = 14, height = 9, dpi = 300, bg = "white")
