#!/usr/bin/env Rscript
################################################################################
## SUPPLEMENTARY FIGURE -- variant-level QC for the 512,667 filtered SNPs used in
## every Indian Chief / Jarvis analysis (the old version summarised the raw
## 52-59 M variants and used a mislabelled missing-genotype file).
## A  mean depth per sample at each site (INFO DP / 379)
## B  minor allele frequency spectrum (equal 0.025 bins, called genotypes)
## C  missing genotypes per sample, ranked, coloured by group
## D  transition/transversion ratio per chromosome
## Inputs (final/results/ic_j_structure/): qc/site_depth.csv, qc/maf_spectrum.csv,
##   qc/titv_by_chr.csv (seq_variant_QC_summary.py), IC_J_per_sample_missingness.csv
## Output: final/supp_figs/SuppFig_seq_statistics_variant_level.{png,pdf}
################################################################################
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
gv   <- function(k, d) { v <- Sys.getenv(k, ""); if (nzchar(v)) v else d }
REPO <- gv("GN_REPO", "/Users/nirwantandukar/Documents/Github/Grain_N_remobilization")
RES  <- gv("GN_RESDIR", file.path(REPO, "final/results/ic_j_structure"))
OUT  <- gv("GN_SUPPDIR", file.path(REPO, "final/supp_figs"))
plot_theme <- theme_minimal(base_size = 24) +
  theme(plot.title = element_text(size = 17, face = "bold", hjust = 0.5, margin = margin(b = 10)),
        plot.tag = element_text(size = 24, face = "bold"), plot.tag.position = c(0.01, 0.99),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, face = "bold", colour = "black"),
        axis.line = element_line(colour = "black"), panel.grid = element_blank(),
        legend.background = element_rect(fill = "white", colour = "grey70", linewidth = 0.4),
        legend.title = element_blank(), legend.text = element_text(size = 13),
        plot.margin = margin(15, 15, 15, 15))
legend_inside <- function(x, y, just) {
  if (utils::packageVersion("ggplot2") >= "3.5.0")
    theme(legend.position = "inside", legend.position.inside = c(x, y), legend.justification = just)
  else theme(legend.position = c(x, y), legend.justification = just)
}
GRP <- c(I01 = "Indian Chief C0", I14 = "Indian Chief C14", J01 = "Jarvis C0", J14 = "Jarvis C14")
COL <- c("Indian Chief C0" = "#92C5DE", "Indian Chief C14" = "#2166AC", "Jarvis C0" = "#FDB863", "Jarvis C14" = "#B35806")
dp <- fread(file.path(RES, "qc/site_depth.csv")); mf <- fread(file.path(RES, "qc/maf_spectrum.csv"))
tv <- fread(file.path(RES, "qc/titv_by_chr.csv")); ms <- fread(file.path(RES, "IC_J_per_sample_missingness.csv"))
ms[, group := factor(GRP[group], levels = GRP)]; setorder(ms, -missing_genotypes); ms[, rank := .I]
tot <- tv[, sum(transitions) / sum(transversions)]
pA <- ggplot(dp, aes(mean_depth_bin_start + 0.05, n_sites / 1e3)) + geom_col(width = 0.09, fill = "#2C7FB8") +
  labs(x = "Mean depth per sample at site (x)", y = "SNPs (thousands)", title = "Site depth") + coord_cartesian(xlim = c(0.5, 1.5)) + plot_theme
pB <- ggplot(mf, aes(maf_bin_start + 0.0125, n_sites / 1e3)) + geom_col(width = 0.023, fill = "#E69F00") +
  labs(x = "Minor allele frequency", y = "SNPs (thousands)", title = "Allele frequency spectrum") + plot_theme
pC <- ggplot(ms, aes(rank, missing_genotypes / 1e3, fill = group)) + geom_col(width = 1) +
  scale_fill_manual(values = COL) +
  labs(x = "Samples ranked by missing genotypes", y = "Missing genotypes (thousands)", title = "Missing genotypes per sample") +
  plot_theme + legend_inside(0.98, 0.98, c(1, 1))
pD <- ggplot(tv, aes(factor(chromosome, levels = chromosome), titv)) + geom_col(width = 0.75, fill = "#8E63B0") +
  geom_hline(yintercept = tot, linetype = "dashed", linewidth = 0.7) +
  annotate("text", x = 10.4, y = tot, label = sprintf("genome %.2f", tot), hjust = 1, vjust = -0.6, size = 5.5, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Chromosome", y = "Ti/Tv", title = "Transition/transversion ratio") + plot_theme
fig <- ((pA | pB) / (pC | pD)) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 28, face = "bold"))
for (ext in c("png", "pdf")) ggsave(file.path(OUT, paste0("SuppFig_seq_statistics_variant_level.", ext)), fig,
                                    width = 18, height = 13, dpi = 300, bg = "white")
