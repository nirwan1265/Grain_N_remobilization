#!/usr/bin/env Rscript
################################################################################
## SUPPLEMENTARY FIGURE -- genotype PCA scree for the Romero-Navarro landrace
## panel, justifying three principal components as GWAS covariates.
##
## Input : Grain_Nitrogen/results/PCA/PCA_variance_explained_maize_romerro.csv
##         SNPRelate varprop (% of total variance); only the first 32 rows carry
##         values (SNPRelate default), rows 33+ are NA and dropped.
## Elbow : PC2 -> PC3 is the last large decrement; every decrement from PC3 on
##         is <= 0.06 points (knee at PC3 over the first 10-20 PCs).
## Output: final/supp_figs/SuppFig_genotype_PCA_scree.{png,pdf}
################################################################################
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
gv   <- function(k, d) { v <- Sys.getenv(k, ""); if (nzchar(v)) v else d }
REPO <- gv("GN_REPO", "/Users/nirwantandukar/Documents/Github/Grain_N_remobilization")
IN   <- gv("GN_PCA_VAR", "/Users/nirwantandukar/Documents/Github/Grain_Nitrogen/results/PCA/PCA_variance_explained_maize_romerro.csv")
FIG  <- gv("GN_FIGDIR", file.path(REPO, "final", "supp_figs"))
dir.create(FIG, FALSE, TRUE)
NPC  <- 20L; KEEP <- 3L

plot_theme <- theme_minimal(base_size = 24) +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 10)),
        plot.tag = element_text(size = 24, face = "bold"), plot.tag.position = c(0.01, 0.99),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, face = "bold", colour = "black"),
        axis.line = element_line(colour = "black"), panel.grid = element_blank(),
        legend.background = element_rect(fill = "white", colour = "grey70", linewidth = 0.4),
        legend.title = element_blank(), legend.text = element_text(size = 13),
        plot.margin = margin(15, 15, 15, 15))

v <- fread(IN); setnames(v, c("PC", "Variance"))
v <- v[!is.na(Variance)][, i := .I][, cum := cumsum(Variance)]
cat(sprintf("eigenvalues: %d | PC1-PC%d explain %.2f%%\n", nrow(v), KEEP, v$cum[KEEP]))

col_keep <- "#2166AC"; col_rest <- "#9E9E9E"; col_mark <- "#B35806"
d <- v[1:NPC][, kept := i <= KEEP]

pA <- ggplot(d, aes(i, Variance, fill = kept)) +
  geom_col(width = .7) +
  geom_vline(xintercept = KEEP + 0.5, linetype = "dashed", colour = col_mark, linewidth = 1) +
  annotate("text", x = KEEP + 1, y = max(d$Variance) * .85,
           label = sprintf("%d PCs used\nas covariates", KEEP),
           hjust = 0, size = 6, fontface = "bold", colour = col_mark) +
  scale_fill_manual(values = c(`TRUE` = col_keep, `FALSE` = col_rest), guide = "none") +
  scale_x_continuous(breaks = c(1, seq(5, NPC, 5))) +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  labs(x = "Genotype principal component", y = "Variance explained (%)") + plot_theme

pB <- ggplot(d, aes(i, cum)) +
  geom_line(colour = col_keep, linewidth = 1.2) +
  geom_point(aes(colour = kept), size = 4) +
  geom_vline(xintercept = KEEP + 0.5, linetype = "dashed", colour = col_mark, linewidth = 1) +
  annotate("text", x = KEEP + 1, y = v$cum[KEEP] - 0.18,
           label = sprintf("PC1-PC%d: %.2f%%", KEEP, v$cum[KEEP]),
           hjust = 0, size = 6, fontface = "bold", colour = col_mark) +
  scale_colour_manual(values = c(`TRUE` = col_keep, `FALSE` = col_rest), guide = "none") +
  scale_x_continuous(breaks = c(1, seq(5, NPC, 5))) +
  labs(x = "Genotype principal component", y = "Cumulative variance explained (%)") + plot_theme

fig <- (pA | pB) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 28, face = "bold"))
for (ext in c("png", "pdf"))
  ggsave(file.path(FIG, paste0("SuppFig_genotype_PCA_scree.", ext)), fig,
         width = 18, height = 7.5, dpi = 300, bg = "white")
cat("written:", FIG, "\n")
