#!/usr/bin/env Rscript
################################################################################
## SUPPLEMENTARY FIGURE 1, PANEL C -- abundance and variability of the 20 free
## amino acids in the Goodman-Buckler panel
##
## Panels A-B of this figure (SuppFig1_FAA_Asn_Pro_rank.R) establish that Asn
## and Pro dominate the RANK ordering. Panel C gives the quantity behind that
## ordering: how much of each amino acid is actually present, and how variable
## it is across the panel.
##
## Why this and not 20 density plots: concentrations span roughly four orders of
## magnitude between amino acids, so overlaid distributions on a common axis
## would be unreadable and separate axes would stop them being comparable. A
## log-scale interval plot shows the location and spread of every trait on one
## comparable axis, and the CV panel beside it separates "abundant" from
## "variable" -- which is the property that matters for mapping.
##
## C1  genotype means per amino acid, log10 scale: median, interquartile range
##     and 5th-95th percentile, ordered by median. Asn and Pro sit at the top,
##     which is the same fact panels A-B make by rank.
## C2  coefficient of variation across genotypes. Abundance and variability are
##     not the same ordering, and traits with the most genetic signal to map are
##     not necessarily the most abundant ones.
##
## Input : data/FAA_Goodman_Buckler.csv   (raw plot records)
## Output: final/supp_figs/SuppFig1C_FAA_abundance.png / .pdf
##         plus a CSV of the per-trait summary for the supplement
################################################################################

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

ROOT <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "."
IN   <- file.path(ROOT, "final", "results", "_stage", "FAA_Goodman_Buckler.csv")
if (!file.exists(IN)) IN <- file.path(ROOT, "data", "FAA_Goodman_Buckler.csv")
OUT  <- file.path(ROOT, "final", "supp_figs")
TAB  <- file.path(ROOT, "final", "supp_tables")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB, showWarnings = FALSE, recursive = TRUE)

plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title   = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.tag     = element_text(size = 24, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x  = element_text(size = 16, colour = "black"),
    axis.text.y  = element_text(size = 16, face = "bold", colour = "black"),
    axis.line    = element_line(colour = "black"),
    panel.grid   = element_blank(),
    legend.background = element_rect(fill = "white", colour = "grey70", linewidth = 0.4),
    legend.title = element_blank(),
    legend.text  = element_text(size = 13),
    plot.margin  = margin(15, 15, 15, 15)
  )

AA <- c("A","R","N","D","Q","E","G","H","I","L","K","M","F","P","S","W","T","Y","V","C")
FULL <- c(A="Ala", R="Arg", N="Asn", D="Asp", Q="Gln", E="Glu", G="Gly", H="His",
          I="Ile", L="Leu", K="Lys", M="Met", F="Phe", P="Pro", S="Ser", W="Trp",
          T="Thr", Y="Tyr", V="Val", C="Cys")
HILITE <- c("N", "P")      # Asn and Pro, the two panels A-B are about

raw <- fread(IN, check.names = FALSE)
setnames(raw, sub("^﻿", "", trimws(names(raw))))
stopifnot(all(c("taxa", AA) %in% names(raw)))

## genotype means across the ~4 replicates, then summarise across genotypes
gm <- raw[, lapply(.SD, mean, na.rm = TRUE), by = taxa, .SDcols = AA]
long <- melt(gm, id.vars = "taxa", variable.name = "aa", value.name = "conc")
long <- long[is.finite(conc) & conc > 0]

summ <- long[, .(n_taxa = .N,
                 median = median(conc),
                 q25 = quantile(conc, .25), q75 = quantile(conc, .75),
                 p05 = quantile(conc, .05), p95 = quantile(conc, .95),
                 mean = mean(conc), sd = sd(conc)), by = aa]
summ[, CV := 100 * sd / mean]
summ[, label := FULL[as.character(aa)]]
summ[, hilite := as.character(aa) %in% HILITE]
setorder(summ, median)
summ[, aa_f := factor(label, levels = label)]

fwrite(summ[order(-median), .(amino_acid = label, code = aa, n_taxa, median, q25, q75,
                              p05, p95, mean, sd, CV)],
       file.path(TAB, "SuppTable1b_FAA_abundance_summary.csv"))

cat("\n--- FAA abundance (genotype means) ---\n")
print(summ[order(-median), .(label, median = round(median, 1),
                             IQR = paste0(round(q25,1), "-", round(q75,1)),
                             CV = round(CV, 1))])
cat(sprintf("\nfold range between most and least abundant: %.0fx\n",
            summ[, max(median) / min(median)]))

pal <- c("TRUE" = "#B2182B", "FALSE" = "#2166AC")

c1 <- ggplot(summ, aes(y = aa_f, colour = hilite)) +
  geom_linerange(aes(xmin = p05, xmax = p95), linewidth = .7, alpha = .45) +
  geom_linerange(aes(xmin = q25, xmax = q75), linewidth = 2.2) +
  geom_point(aes(x = median), size = 3.1) +
  scale_colour_manual(values = pal, guide = "none") +
  scale_x_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  ## headroom inside the panel so the tag clears the top trait label
  scale_y_discrete(expand = expansion(add = c(0.6, 1.4))) +
  labs(x = "Concentration (genotype mean, log scale)", y = NULL, tag = "C") +
  plot_theme +
  ## extra headroom so the panel tag does not sit on the top trait label
  theme(axis.text.y = element_text(size = 14),
        plot.margin = margin(15, 15, 15, 15))

c2 <- ggplot(summ, aes(CV, aa_f, colour = hilite)) +
  geom_segment(aes(x = 0, xend = CV, yend = aa_f), linewidth = .5, alpha = .45) +
  geom_point(size = 3.1) +
  scale_colour_manual(values = pal, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, .06))) +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.4))) +
  labs(x = "CV across genotypes (%)", y = NULL) +
  plot_theme +
  theme(axis.text.y = element_blank(),
        plot.margin = margin(15, 15, 15, 15))

fig <- c1 + c2 + plot_layout(widths = c(1.25, 1))
ggsave(file.path(OUT, "SuppFig1C_FAA_abundance.png"), fig,
       width = 15, height = 9, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(OUT, "SuppFig1C_FAA_abundance.pdf"), fig,
       width = 15, height = 9, bg = "white", limitsize = FALSE)
message("Saved: ", file.path(OUT, "SuppFig1C_FAA_abundance.png"))
