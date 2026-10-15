#!/usr/bin/env Rscript
################################################################################
## FIGURE -- CornCyc pathways recurring across the four candidate sets
##
## Every pathway (distinct gene sets) with pathway-specific permutation P < 0.05, from
## pathway_recurrence_permutation.py (Bonferroni GWAS lists; 20,000 genome-wide
## circular shifts of each set's intervals, chromosome by chromosome).
## Left: which candidate sets contain a gene of the pathway (black = present, white = absent).
## Right: -log10 pathway-specific permutation P (uncorrected; dashed line P = 0.05).
## One figure, no panel letters.
## Pathways with identical gene sets (phosphatidylglycerol biosynthesis I and II) are
## shown once.
## Input : final/results/pathway_recurrence/pathway_permutation_all_Bonf.csv
## Output: final/main_figs/Fig_pathway_recurrence.{png,pdf}
################################################################################
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
ROOT <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "."
IN   <- file.path(ROOT, "final", "results", "pathway_recurrence", "pathway_permutation_all_Bonf.csv")
OUT  <- file.path(ROOT, "final", "main_figs")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.tag     = element_text(size = 24, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x  = element_text(size = 16, colour = "black"),
    axis.text.y  = element_text(size = 16, face = "bold", colour = "black"),
    axis.line    = element_line(colour = "black"),
    panel.grid   = element_blank(),
    legend.position = "none",
    plot.margin  = margin(15, 15, 15, 15)
  )

d <- fread(IN)
d[, pathway := gsub("<[^>]+>", "", pathway)]
## one row per distinct gene set (identical hit pattern and size, consecutive in rank)
d[, key := paste(n_genes, GB, soilN, IC, J, P_specific)]
d <- d[!duplicated(key)]
d[pathway %like% "^phosphatidylglycerol biosynthesis I \\(", pathway := "phosphatidylglycerol biosynthesis I and II"]
top <- d[P_specific < 0.05][order(P_specific)]          # every pathway with uncorrected P < 0.05
top[, label := sprintf("%s\n(%d genes)", pathway, n_genes)]
top[, label := factor(label, levels = rev(label))]

## A: hit matrix
SETS <- c("Amino acid GWAS", "Soil-N GWAS", "Indian Chief", "Jarvis")
m <- melt(top[, .(label, `Amino acid GWAS` = GB, `Soil-N GWAS` = soilN, `Indian Chief` = IC, Jarvis = J)],
          id.vars = "label", variable.name = "set", value.name = "hit")
m[, set := factor(set, levels = SETS)]
pA <- ggplot(m, aes(set, label)) +
  geom_point(aes(fill = factor(hit)), shape = 21, size = 11, colour = "black", stroke = 1) +
  scale_fill_manual(values = c("1" = "black", "0" = "white")) +
  scale_x_discrete(position = "top", labels = c("Amino acid\nGWAS", "Soil-N\nGWAS", "Indian\nChief", "Jarvis")) +
  labs(x = NULL, y = NULL) + plot_theme +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(size = 15, face = "bold", colour = "black"))

## B: -log10 P
top[, lp := -log10(P_specific)]
pB <- ggplot(top, aes(lp, label)) +
  geom_segment(aes(x = 0, xend = lp, yend = label), linewidth = 1.4, colour = "grey45") +
  geom_point(size = 7, colour = "#7F2704") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", colour = "grey40", linewidth = 0.7) +
  annotate("text", x = -log10(0.05), y = 0.45, label = "P = 0.05", hjust = -0.1, vjust = 0, size = 4.6, colour = "grey30") +
  geom_text(aes(label = paste0("P = ", formatC(signif(P_specific, 2), digits = 2, format = "fg", flag = "#"))),
            hjust = -0.18, size = 5.2) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.5)), limits = c(0, NA)) +
  labs(x = expression(bold(-log[10](italic(P)))), y = NULL) +
  plot_theme + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

fig <- (pA | pB) + plot_layout(widths = c(1, 0.8))
for (ext in c("png", "pdf"))
  ggsave(file.path(OUT, paste0("Fig_pathway_recurrence.", ext)), fig, width = 16, height = 1.35 * nrow(top) + 2,
         dpi = 300, bg = "white", device = if (ext == "png") ragg::agg_png else "pdf")
print(top[, .(pathway, n_genes, R_obs, P_specific)])
