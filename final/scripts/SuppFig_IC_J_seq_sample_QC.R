#!/usr/bin/env Rscript
################################################################################
## SUPPLEMENTARY FIGURE -- per-sample sequencing depth, Indian Chief and Jarvis
## skim sequencing (379 plants). Shows the depth difference between cycles that
## motivates genotype-likelihood methods for Fig. 4.
## A  mean depth per sample by population-cycle group
## B  missing genotypes per sample (of 512,667 filtered SNPs) against mean depth
## Inputs: qc_stats/per_sample_mean_depth.tsv; per-sample missingness computed from the
##         filtered PLINK genotypes (final/results/ic_j_structure/IC_J_per_sample_missingness.csv).
##         NB qc_stats/missing_counts.tsv is NOT missing genotypes (it rises with depth) -- not used.
## Output: final/supp_figs/SuppFig_seq_statistics_sample_level.{png,pdf}
################################################################################
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
gv   <- function(k, d) { v <- Sys.getenv(k, ""); if (nzchar(v)) v else d }
REPO <- gv("GN_REPO", "/Users/nirwantandukar/Documents/Github/Grain_N_remobilization")
F_DP <- gv("GN_DEPTH", "/Users/nirwantandukar/Documents/Research/data/Indian_Jarvis/qc_stats/per_sample_mean_depth.tsv")
F_MS <- gv("GN_MISSING", file.path(REPO, "final/results/ic_j_structure/IC_J_per_sample_missingness.csv"))
OUT  <- gv("GN_SUPPDIR", file.path(REPO, "final/supp_figs"))
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
GRP <- c(I01 = "Indian Chief C0", I14 = "Indian Chief C14", J01 = "Jarvis C0", J14 = "Jarvis C14")
COL <- c("Indian Chief C0" = "#92C5DE", "Indian Chief C14" = "#2166AC", "Jarvis C0" = "#FDB863", "Jarvis C14" = "#B35806")
d <- fread(F_DP, header = FALSE, col.names = c("sample", "depth", "n_geno"))
m <- fread(F_MS)[, .(sample, missing = missing_genotypes, missing_fraction)]
d <- merge(d, m, by = "sample")
cat(sprintf("missing genotypes per sample: median %s of 512,667 (%.0f%%)\n", format(median(d$missing), big.mark = ","), 100 * median(d$missing_fraction)))
d[, group := factor(GRP[sub("_.*", "", sample)], levels = GRP)]
s <- d[, .(n = .N, mean = mean(depth), median = median(depth)), by = group][order(group)]
print(s); cat(sprintf("all: n = %d, mean %.2f, median %.2f, range %.2f-%.2f\n", nrow(d), mean(d$depth),
                      median(d$depth), min(d$depth), max(d$depth)))
pA <- ggplot(d, aes(group, depth, colour = group, fill = group)) +
  geom_violin(alpha = 0.15, linewidth = 0.7) +
  geom_boxplot(width = 0.18, alpha = 0.4, outlier.shape = NA, linewidth = 0.7) +
  geom_jitter(width = 0.1, size = 1.4, alpha = 0.5) +
  geom_text(data = s, aes(group, max(d$depth) * 1.06, label = sprintf("%.2fx", mean)),
            inherit.aes = FALSE, size = 5.5, fontface = "bold") +
  scale_colour_manual(values = COL, guide = "none") + scale_fill_manual(values = COL, guide = "none") +
  scale_x_discrete(labels = function(z) sub(" C", "\nC", z)) +
  labs(x = NULL, y = "Mean sequencing depth (x)", title = "Per-sample depth") + plot_theme +
  theme(axis.text.x = element_text(size = 15, face = "bold"), plot.title = element_text(size = 17, face = "bold"))
pB <- ggplot(d, aes(depth, missing / 1e3, colour = group)) +
  geom_point(size = 2.2, alpha = 0.75) +
  scale_colour_manual(values = COL) +
  labs(x = "Mean sequencing depth (x)", y = "Missing genotypes (thousands)", title = "Missing genotypes per sample") +
  plot_theme + theme(legend.position = "bottom", plot.title = element_text(size = 17, face = "bold")) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 5)))
fig <- (pA | pB) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 28, face = "bold"))
for (ext in c("png", "pdf")) ggsave(file.path(OUT, paste0("SuppFig_seq_statistics_sample_level.", ext)), fig,
                                    width = 18, height = 7.5, dpi = 300, bg = "white")
