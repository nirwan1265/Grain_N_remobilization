#!/usr/bin/env Rscript
################################################################################
## SUPPLEMENTARY FIGURE -- genotype PCA scree for the Romero-Navarro landrace
## panel, justifying the use of three principal components as GWAS covariates.
##
## Input : Grain_Nitrogen/results/PCA/PCA_variance_explained_maize_romerro.csv
##         SNPRelate varprop, as PERCENT of total variance. The file lists 3,539
##         rows but only the first 32 carry values -- SNPRelate returns the top
##         32 eigenvalues by default. Rows 33+ are NA and are dropped.
##
## WHERE THE ELBOW IS -- read before changing the annotation
##   The knee (maximum perpendicular distance from the chord joining the first
##   and last point) is at PC3 when computed over the first 10, 15 or 20 PCs,
##   and at PC5 over all 32. PC3 is defensible but not unique, so the figure
##   shows BOTH the per-PC drop and the cumulative curve and lets the reader
##   see the flattening rather than asserting a single elbow. The honest
##   statement is that PC2 -> PC3 is the last large decrement (0.194 points)
##   and every decrement from PC3 onward is <= 0.06.
################################################################################
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})

ROOT <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "."
IN   <- Sys.getenv("GN_PCA_VAR",
         file.path(ROOT, "PCA_variance_explained_maize_romerro.csv"))
OUT  <- Sys.getenv("GN_OUT", file.path(ROOT, "SuppFig_genotype_PCA_scree.png"))
NPC  <- 20L      # PCs to display
KEEP <- 3L       # PCs used as covariates

v <- fread(IN); setnames(v, c("PC","Variance"))
v <- v[!is.na(Variance)][, i := .I][, cum := cumsum(Variance)]
cat(sprintf("eigenvalues available: %d | first %d PCs explain %.2f%% of variance\n",
            nrow(v), KEEP, v$cum[KEEP]))

col_keep <- "#2166AC"; col_rest <- "#9E9E9E"; col_mark <- "#B35806"
th <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey88", linewidth = .3),
        axis.line  = element_line(colour = "black"),
        axis.text  = element_text(size = 10, colour = "black", face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.tag   = element_text(size = 16, face = "bold"),
        plot.margin = margin(10, 12, 10, 10))

d <- v[1:NPC][, kept := i <= KEEP]

pA <- ggplot(d, aes(i, Variance, fill = kept)) +
  geom_col(width = .68) +
  geom_vline(xintercept = KEEP + 0.5, linetype = "dashed",
             colour = col_mark, linewidth = .6) +
  annotate("text", x = KEEP + 0.9, y = max(d$Variance) * .85,
           label = sprintf("%d PCs used\nas covariates", KEEP),
           hjust = 0, size = 3.6, fontface = "bold", colour = col_mark) +
  scale_fill_manual(values = c(`TRUE` = col_keep, `FALSE` = col_rest), guide = "none") +
  scale_x_continuous(breaks = c(1, seq(5, NPC, 5))) +
  labs(x = "Genotype principal component", y = "Variance explained (%)", tag = "A") + th

pB <- ggplot(d, aes(i, cum)) +
  geom_line(colour = col_keep, linewidth = .8) +
  geom_point(aes(colour = kept), size = 2.2) +
  geom_vline(xintercept = KEEP + 0.5, linetype = "dashed",
             colour = col_mark, linewidth = .6) +
  annotate("text", x = KEEP + 0.9, y = v$cum[KEEP] - 0.18,
           label = sprintf("PC1-PC%d: %.2f%%", KEEP, v$cum[KEEP]),
           hjust = 0, size = 3.6, fontface = "bold", colour = col_mark) +
  scale_colour_manual(values = c(`TRUE` = col_keep, `FALSE` = col_rest), guide = "none") +
  scale_x_continuous(breaks = c(1, seq(5, NPC, 5))) +
  labs(x = "Genotype principal component",
       y = "Cumulative variance explained (%)", tag = "B") + th

ggsave(OUT, pA | pB, width = 11, height = 4.4, dpi = 300, bg = "white")
cat("Saved:", OUT, "\n")
