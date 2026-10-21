#!/usr/bin/env Rscript
## FIGURE -- free Asp and Pro: spread in the Goodman-Buckler panel and QTL in B64 x CML14
## A, C: genotype means (black) and plot values (grey) of free Asp / Pro across the 280 GB lines,
##       ranked by mean; B64 and CML14 (parents of the mapping population) marked.
## B, D: genome-wide Haley-Knott LOD scans for Asp / Pro BLUPs in the F2:3 population; shaded band = overlap of the Asp and Pro 1.5-LOD support intervals on chr3;
##       dashed line = permutation threshold reported in the text (3.70 Asp, 3.63 Pro).
## Input : data/FAA_Goodman_Buckler.csv; final/results/qtl/lod.csv (QTL_HK_scan.py)
## Output: final/supp_figs/SuppFig_QTL_Asp_Pro_panel.{png,pdf}
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
args <- commandArgs(TRUE); ROOT <- if (length(args)) args[1] else "."
FAA <- file.path(ROOT, "data/FAA_Goodman_Buckler.csv"); LOD <- file.path(ROOT, "final/results/qtl/lod.csv")
OUT <- file.path(ROOT, "final/supp_figs")
plot_theme <- theme_minimal(base_size = 24) + theme(
  plot.tag = element_text(size = 24, face = "bold"), plot.tag.position = c(0.01, 0.99),
  axis.title.x = element_text(size = 18, face = "bold"), axis.title.y = element_text(size = 18, face = "bold"),
  axis.text.x = element_text(size = 15, colour = "black"), axis.text.y = element_text(size = 15, colour = "black"),
  axis.line = element_line(colour = "black"), panel.grid = element_blank(), legend.position = "none",
  plot.margin = margin(15, 15, 15, 15))
r <- fread(FAA, encoding = "UTF-8"); setnames(r, 1, "taxa")
spread <- function(aa, lab, tag) {
  d <- r[, .(taxa, v = as.numeric(get(aa)))][is.finite(v)]
  m <- d[, .(mean = mean(v)), by = taxa][order(mean)][, rank := .I]
  d <- merge(d, m[, .(taxa, rank)], by = "taxa")
  par <- m[taxa %in% c("B64", "CML14")]
  ggplot() +
    geom_point(data = d, aes(rank, v), colour = "grey70", size = 0.9) +
    geom_line(data = m, aes(rank, mean), colour = "black", linewidth = 0.8) +
    geom_point(data = par, aes(rank, mean), shape = 21, fill = "white", colour = "black", size = 4.5, stroke = 1.2) +
    geom_text(data = par, aes(rank, mean, label = taxa), vjust = -1.1, size = 6, fontface = "bold") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    labs(x = "Goodman-Buckler lines, ranked by mean", y = lab, tag = tag) + plot_theme +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}
g <- fread(LOD); setorder(g, chr, cM); g[, cM0 := cM - min(cM), by = chr]
off <- g[, .(len = max(cM0)), by = chr][, start := shift(cumsum(len + 20), fill = 0)]
g <- merge(g, off[, .(chr, start)], by = "chr"); g[, x := cM0 + start]
mid <- g[, .(mid = (min(x) + max(x)) / 2, lo = min(x), hi = max(x)), by = chr]
mk <- g[is_marker == TRUE]
## shared 1.5-LOD support interval on chr3 (overlap of the Asp and Pro intervals)
lodint <- function(tr) { c3 <- g[chr == 3][order(cM0)]; i <- which.max(c3[[tr]]); L <- c3[[tr]][i]
  lo <- i; hi <- i
  while (lo > 1 && c3[[tr]][lo - 1] >= L - 1.5) lo <- lo - 1
  while (hi < nrow(c3) && c3[[tr]][hi + 1] >= L - 1.5) hi <- hi + 1
  c(c3$x[lo], c3$x[hi]) }
ia <- lodint("Aspartate"); ip <- lodint("Proline")
shade <- c(max(ia[1], ip[1]), min(ia[2], ip[2]))
message(sprintf("Asp interval %.1f-%.1f, Pro %.1f-%.1f, shared %.1f-%.1f (plot x)", ia[1], ia[2], ip[1], ip[2], shade[1], shade[2]))
scan <- function(tr, thr, tag) {
  ggplot(g, aes(x, get(tr), group = chr)) +
    annotate("rect", xmin = shade[1], xmax = shade[2], ymin = -Inf, ymax = Inf, fill = "grey70") +
    geom_line(colour = "black", linewidth = 0.6) +
    geom_hline(yintercept = thr, linetype = "dashed", colour = "black", linewidth = 0.6) +
    geom_rug(data = mk, aes(x = x), inherit.aes = FALSE, sides = "b", length = unit(0.015, "npc"), colour = "grey40") +
    scale_x_continuous(breaks = mid$mid, labels = mid$chr, expand = expansion(mult = 0.01)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Chromosome", y = "LOD", tag = tag) + plot_theme
}
fig <- (spread("D", "Free Asp", "A") | scan("Aspartate", 3.70, "B")) /
       (spread("P", "Free Pro", "C") | scan("Proline", 3.63, "D")) + plot_layout(widths = c(1, 2))
for (e in c("png", "pdf")) ggsave(file.path(OUT, paste0("SuppFig_QTL_Asp_Pro_panel.", e)), fig, width = 20, height = 10, dpi = 300, bg = "white")
