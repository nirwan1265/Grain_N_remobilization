################################################################################
## Regional association plots with LD for NMX (TATC, Zm00001eb399680, chr9) and
## ZmFd9 (Zm00001eb421930, chr10), using the GBS landrace genotypes for LD because
## the PHG genotypes are not available. Inputs come from LD_GBS_regional_NMX_ZmFd9.py.
##
## Two versions are written:
##  _GBS : soil-N GWAS on the 495,755 GBS SNPs (lowest P of MLM, BLINK, FarmCPU per
##         SNP), every SNP coloured by r^2 with the GBS lead SNP.
##  _PHG : PHG soil-N GLM P values. Only PHG SNPs at GBS positions have an r^2
##         (with the GBS lead SNP), the rest are grey.
## Grey vertical line = gene; red dashed line = Bonferroni 0.05 / number of SNPs.
##
## Usage: Rscript SuppFig_LD_regional_NMX_ZmFd9_GBS.R <in_dir> <out_dir>
################################################################################
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
args <- commandArgs(trailingOnly = TRUE)
IN  <- if (length(args) >= 1) args[1] else "final/results/soil_nitrogen_gwas/LD_GBS_regional"
OUT <- if (length(args) >= 2) args[2] else "final/supp_figs"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
L <- fread(file.path(IN, "leads.csv"))
RAMP <- c("#FFFFCC", "#FED976", "#FD8D3C", "#E31A1C", "#800026")
THR <- c(GBS = -log10(0.05 / 495755), PHG = -log10(0.05 / 2381512))

panel <- function(i, set) {
  l <- L[i]
  if (set == "GBS") {
    d <- fread(file.path(IN, paste0(l$locus, "_gbs_r2.csv")))
    d[, p := pmin(P_MLM, P_BLINK, P_FarmCPU, na.rm = TRUE)]
  } else {
    d <- fread(file.path(IN, paste0(l$locus, "_phg_glm.csv")))
  }
  d <- d[!is.na(p)]
  d[, `:=`(mb = pos / 1e6, lp = -log10(p))]
  setorder(d, r2, na.last = FALSE)
  ld <- d[pos == l$gbs_lead]
  ymax <- max(d$lp, THR[[set]]) * 1.1
  g <- ggplot(d, aes(mb, lp)) +
    geom_vline(xintercept = (l$gene_start + l$gene_end) / 2e6, colour = "grey55", linewidth = 0.6) +
    geom_hline(yintercept = THR[[set]], colour = "#E31A1C", linetype = "dashed", linewidth = 0.6) +
    geom_point(data = d[is.na(r2)], colour = "grey80", size = 1.4) +
    geom_point(data = d[!is.na(r2) & pos != l$gbs_lead], aes(fill = r2), shape = 21,
               colour = "grey35", stroke = 0.25, size = 2.2) +
    annotate("text", x = (l$gbs_lead - 1e6) / 1e6, y = ymax * 0.96, label = paste0("Chr ", l$chr), hjust = 0, size = 4.2) +
    scale_fill_gradientn(colours = RAMP, limits = c(0, 1), name = expression(r^2)) +
    scale_x_continuous(limits = c(l$gbs_lead - 1e6, l$gbs_lead + 1e6) / 1e6, expand = expansion(mult = 0.01)) +
    scale_y_continuous(limits = c(0, ymax), expand = expansion(mult = c(0, 0.02))) +
    labs(title = l$locus, x = sprintf("Position on chromosome %d (Mb)", l$chr), y = expression(-log[10](italic(P)))) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(face = "italic", hjust = 0.5, size = 15),
          axis.text = element_text(colour = "black"),
          legend.key.height = unit(14, "pt"), legend.key.width = unit(8, "pt"))
  if (nrow(ld) == 1)
    g <- g + geom_point(data = ld, fill = "#800026", shape = 23, colour = "black", stroke = 0.6, size = 4) +
      geom_label(data = ld, aes(label = sprintf("chr%d:%s", l$chr, format(l$gbs_lead, big.mark = ","))),
                 hjust = -0.08, vjust = 0.2, size = 3.4, label.size = 0.3, label.r = unit(0, "pt"))
  g
}
for (set in c("GBS", "PHG")) {
  fig <- (panel(1, set) | panel(2, set)) + plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 16))
  for (ext in c("png", "pdf"))
    ggsave(file.path(OUT, sprintf("SuppFig_LD_regional_NMX_ZmFd9_%s.%s", set, ext)), fig,
           width = 11, height = 3.8, dpi = 300, bg = "white")
}
cat("done\n")
