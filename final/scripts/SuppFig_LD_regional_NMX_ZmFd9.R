################################################################################
## Regional association plots with LD for two soil-N loci
##   NMX   (TATC, Zm00001eb399680)  chr9
##   ZmFd9 (Zm00001eb421930)        chr10
## -log10(P) from the PHG soil-N GLM (GLM_SoilN_maize.txt, Pos column = B73 v5;
## the SNP names in that file are wrong and are not used). Points are coloured by
## r^2 with the lead SNP (PLINK --r2 output from LD_PHG_regional_NMX_ZmFd9.sh).
## Grey vertical line = gene; red dashed line = Bonferroni 0.05 / number of SNPs.
##
## Inputs : GLM file, final/results/soil_nitrogen_gwas/LD_PHG/{NMX,ZmFd9}_r2.ld
## Output : final/supp_figs/SuppFig_LD_regional_NMX_ZmFd9.{png,pdf}
## Usage  : Rscript SuppFig_LD_regional_NMX_ZmFd9.R   (paths can be overridden
##          with GN_REPO, GN_GLM, GN_LDDIR, GN_FIGDIR, GN_FLANK)
################################################################################
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
gv    <- function(k, d) { v <- Sys.getenv(k, ""); if (nzchar(v)) v else d }
REPO  <- gv("GN_REPO", "/Users/nirwantandukar/Documents/Github/Grain_N_remobilization")
GLM   <- gv("GN_GLM", "/Users/nirwantandukar/Documents/Research/results/GWAS/GLM/nitrogen/PHG/TN_PHG/GLM_SoilN_maize.txt")
LDD   <- gv("GN_LDDIR", file.path(REPO, "final/results/soil_nitrogen_gwas/LD_PHG"))
OUT   <- gv("GN_FIGDIR", file.path(REPO, "final/supp_figs"))
FLANK <- as.numeric(gv("GN_FLANK", "1e6"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

LOCI <- data.table(name   = c("NMX", "ZmFd9"),
                   gene   = c("Zm00001eb399680", "Zm00001eb421930"),
                   chr    = c(9L, 10L),
                   gstart = c(151731980, 116703178),
                   gend   = c(151737594, 116706315),
                   lead   = c(151733953, 116706616))

g <- fread(GLM, select = c("Chr", "Pos", "p"))
g[, chr := as.integer(sub("^chr", "", Chr))]
THR <- -log10(0.05 / nrow(g))
RAMP <- c("#FFFFCC", "#FED976", "#FD8D3C", "#E31A1C", "#800026")

panel <- function(i) {
  L <- LOCI[i]
  d <- g[chr == L$chr & Pos >= L$lead - FLANK & Pos <= L$lead + FLANK, .(Pos, p)]
  d <- d[, .(p = min(p)), by = Pos]
  r <- fread(file.path(LDD, paste0(L$name, "_r2.ld")))[, .(Pos = BP_B, R2)]
  r <- r[, .(R2 = max(R2)), by = Pos]
  d <- merge(d, r, by = "Pos", all.x = TRUE)
  d[Pos == L$lead, R2 := 1]
  d[, `:=`(mb = Pos / 1e6, lp = -log10(p))]
  setorder(d, R2, na.last = FALSE)            # strongest LD drawn on top
  ld <- d[Pos == L$lead]
  stopifnot(nrow(ld) == 1)
  cat(sprintf("%s: %d SNPs, %d with r2, lead -log10P %.1f\n", L$name, nrow(d), d[!is.na(R2), .N], ld$lp))
  ymax <- max(d$lp, THR) * 1.08
  ggplot(d, aes(mb, lp)) +
    geom_vline(xintercept = (L$gstart + L$gend) / 2e6, colour = "grey55", linewidth = 0.6) +
    geom_hline(yintercept = THR, colour = "#E31A1C", linetype = "dashed", linewidth = 0.6) +
    geom_point(data = d[is.na(R2)], colour = "grey80", size = 1.6) +
    geom_point(data = d[!is.na(R2) & Pos != L$lead], aes(fill = R2), shape = 21,
               colour = "grey35", stroke = 0.25, size = 2.2) +
    geom_point(data = ld, fill = "#800026", shape = 23, colour = "black", stroke = 0.6, size = 4) +
    geom_label(data = ld, aes(label = sprintf("chr%d:%s", L$chr, format(L$lead, big.mark = ","))),
               hjust = -0.08, vjust = 0.2, size = 3.4, label.size = 0.3, label.r = unit(0, "pt")) +
    annotate("text", x = (L$lead - FLANK) / 1e6, y = ymax * 0.96, label = paste0("Chr ", L$chr),
             hjust = 0, size = 4.2) +
    scale_fill_gradientn(colours = RAMP, limits = c(0, 1), name = expression(r^2)) +
    scale_x_continuous(limits = c(L$lead - FLANK, L$lead + FLANK) / 1e6,
                       expand = expansion(mult = 0.01)) +
    scale_y_continuous(limits = c(0, ymax), expand = expansion(mult = c(0, 0.02))) +
    labs(title = L$name, x = sprintf("Position on chromosome %d (Mb)", L$chr),
         y = expression(-log[10](italic(P)))) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(face = "italic", hjust = 0.5, size = 15),
          axis.text = element_text(colour = "black"),
          legend.key.height = unit(14, "pt"), legend.key.width = unit(8, "pt"))
}

fig <- (panel(1) | panel(2)) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 16))
for (ext in c("png", "pdf"))
  ggsave(file.path(OUT, paste0("SuppFig_LD_regional_NMX_ZmFd9.", ext)), fig,
         width = 11, height = 3.8, dpi = 300, bg = "white")
cat("wrote", file.path(OUT, "SuppFig_LD_regional_NMX_ZmFd9.{png,pdf}"), "\n")
