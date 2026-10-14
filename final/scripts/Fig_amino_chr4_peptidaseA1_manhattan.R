#!/usr/bin/env Rscript
################################################################################
## FIGURE 1 -- a locus on chromosome 4 shared by eight amino acid phenotypes,
## next to a peptidase A1 (aspartic protease) gene
##
## Locus chr4:174.55-174.57 Mb (B73 RefGen_v5). It is FDR-significant (5% BH,
## FarmCPU) for eight phenotypes: A, LAV, IVL, V, K, F, E.Total, T -- five
## absolute amino acids from four biosynthetic families (Ala, Val, Lys, Phe,
## Thr) and three composites/ratios. Two SNP groups:
##   * 174,551,750-174,555,207 (three SNPs, r2 0.92-0.94; E.Total, K, F, IVL,
##     LAV, T), inside Zm00001eb190790 (DUF2404) and 2.6-6.1 kb upstream of
##     Zm00001eb190800 (peptidase A1; + strand, 174,557,846-174,559,120);
##   * 174,565,066 (A, V; P = 6.0e-15), 5.9 kb downstream of Zm00001eb190800,
##     r2 0.25-0.27 with the first group.
## The arrow marks the gene; the gene label is printed once (first panel).
##
## Points: every SNP is plotted (4,177,796 per phenotype; nothing thinned).
## Point size 0.9 below the panel's FDR line and 1.8 above it.
##
## COORDINATES. The genotypes are HapMap 3.2.1 (hmp321_282_agpv4), so the scan
## positions are AGPv4; every SNP is lifted to v5 with the chain file
## (liftover_v4_v5.R) before plotting, so the x axis and the gene are both v5.
##
## Input : final/results/manhattan_full/<trait>__<MODEL>.csv.gz  (Chr,Pos,P; no header;
##         the model's rows of the full per-trait scan, AGPv4 positions)
##         final/results/goodman_buckler_amino_acid_gwas/by_model_fdr/by_model (FDR cutoffs)
## Output: final/main_figs/Fig_amino_chr4_peptidaseA1_manhattan.png
################################################################################

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

ROOT     <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "."
gv       <- function(k, d) { v <- Sys.getenv(k, ""); if (nzchar(v)) v else d }
IN_DIR   <- gv("GN_SCANS", file.path(ROOT, "final", "results", "manhattan_full"))
OUT_FILE <- gv("GN_OUT", file.path(ROOT, "final", "main_figs", "Fig_amino_chr4_peptidaseA1_manhattan.png"))
CHAIN    <- gv("GN_CHAIN", "/Users/nirwantandukar/Documents/Research/data/Sarah_amino_acid/chain/B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain.txt")
source(gv("GN_LIFT", file.path(ROOT, "final", "scripts", "liftover_v4_v5.R")))
BLOCKS   <- read_chain_blocks(CHAIN)

## THRESHOLD: Benjamini-Hochberg FDR 5%, the threshold used throughout the
## paper. BH gives each phenotype x model its own cutoff, so each panel draws
## its own line rather than one shared line. In practice the nine cutoffs span
## -log10 6.05 to 6.97, under one log unit, so the panels stay comparable by
## eye. The Bonferroni and suggestive lines of the previous version are gone:
## neither threshold is used anywhere in the manuscript any more.
FDR_DIR <- gv("GN_FDR", file.path(ROOT, "final", "results", "goodman_buckler_amino_acid_gwas",
                     "by_model_fdr", "by_model"))

fdr_cutoff <- function(trait, model) {
  f <- file.path(FDR_DIR, sprintf("%s_%s.RData", trait, model))
  if (!file.exists(f)) stop("no FDR file for ", trait, " / ", model, ": ", f)
  e <- new.env(parent = emptyenv()); nm <- load(f, envir = e)
  d <- as.data.table(get(nm[1], envir = e))
  if (!nrow(d[fdr_sig == TRUE])) stop(trait, "/", model, ": no FDR-significant SNP")
  max(d[fdr_sig == TRUE]$P.value)
}

GENE       <- "Zm00001eb190800"
GENE_LAB   <- "Peptidase A1 (aspartic protease)"
GENE_CHR   <- 4L
GENE_POS   <- 174558483       # gene midpoint, v5
PEAK_WINDOW <- 250e3          # anchor the arrow on the lead SNP within this

## two columns of four, ordered by the P value at this locus
## (A 6.0e-15, LAV 4.5e-14, IVL 2.3e-10, V 3.2e-10 | K 4.3e-10, F 3.1e-9,
##  E.Total 2.1e-8, T 1.5e-7); all FarmCPU, the model that detected each one
COLUMNS <- list(
  list(panels = list(
    c("A",   "FarmCPU"), c("LAV", "FarmCPU"),
    c("IVL", "FarmCPU"), c("V",   "FarmCPU"))),
  list(panels = list(
    c("K",       "FarmCPU"), c("F", "FarmCPU"),
    c("E.Total", "FarmCPU"), c("T", "FarmCPU")))
)
## draft mode: GN_ALLOW_MISSING=1 draws an empty placeholder for a missing scan
ALLOW_MISSING <- Sys.getenv("GN_ALLOW_MISSING", "") == "1"

band_dark  <- "#3B3B3B"
band_light <- "#9E9E9E"
fdr_col    <- "#C0392B"
marker_col <- "#1A5276"

theme_gwas <- theme_minimal(base_size = 11) +
  theme(axis.text          = element_text(size = 8, colour = "black"),
        axis.title         = element_text(size = 10, face = "bold"),
        axis.line          = element_line(colour = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        legend.position    = "none",
        plot.margin        = margin(4, 6, 4, 6))

read_scan <- function(trait, model) {
  f <- file.path(IN_DIR, sprintf("%s__%s.csv.gz", trait, model))
  if (!file.exists(f)) stop("missing: ", f)
  d <- fread(f, header = FALSE, col.names = c("Chr", "Pos", "P"), showProgress = FALSE)
  d <- d[is.finite(P) & P > 0 & Chr >= 1 & Chr <= 10]
  d[, Pos := lift_v4_v5(Chr, Pos, BLOCKS)$pos_v5]          # AGPv4 -> v5
  d
}

## genome coordinates from the first scan
first  <- read_scan(COLUMNS[[1]]$panels[[1]][1], COLUMNS[[1]]$panels[[1]][2])
chrlen <- first[, .(len = max(Pos)), by = Chr][order(Chr)]
chrlen[, offset := cumsum(as.numeric(len)) - as.numeric(len)]
axis_df <- chrlen[, .(Chr, centre = offset + len / 2)]
rm(first); invisible(gc())

marker_g   <- chrlen[Chr == GENE_CHR, offset] + GENE_POS
genome_end <- chrlen[, max(offset + len)]

label_hjust <- function(g) { f <- g / genome_end
  if (f > 0.85) 1.02 else if (f < 0.15) -0.02 else 0.5 }

panel <- function(trait, model, show_x, gene_lab = NA_character_) {
  if (ALLOW_MISSING && !file.exists(file.path(IN_DIR, sprintf("%s__%s.csv.gz", trait, model))))
    return(ggplot() + annotate("text", x = 0, y = 0, label = sprintf("%s (%s): scan not yet available", trait, model),
                               size = 3, colour = "grey40") + theme_void())
  d <- merge(read_scan(trait, model), chrlen[, .(Chr, offset)], by = "Chr")
  d[, gpos := Pos + offset][, logp := -log10(P)]
  d[, band := factor(Chr %% 2)]
  thr <- -log10(fdr_cutoff(trait, model))
  d[, sig := logp >= thr]

  near <- d[abs(gpos - marker_g) <= PEAK_WINDOW]
  px <- if (nrow(near)) near$gpos[which.max(near$logp)] else marker_g
  py <- if (nrow(near)) max(near$logp) else 0
  ytop <- max(d$logp)
  a_lo <- py + 0.05 * ytop
  a_hi <- py + 0.20 * ytop
  head_room <- if (is.na(gene_lab)) 0.26 else 0.36

  ggplot(d, aes(gpos, logp, colour = band)) +
    geom_point(data = d[sig == FALSE], size = .9, alpha = .85, shape = 16) +
    geom_point(data = d[sig == TRUE],  size = 1.8, alpha = 1, shape = 16) +   # above the FDR line: larger
    geom_hline(yintercept = thr, linetype = "dashed",
               colour = fdr_col, linewidth = .6) +
    annotate("segment", x = px, xend = px, y = a_hi, yend = a_lo,
             colour = marker_col, linewidth = .7,
             arrow = arrow(type = "closed", length = unit(0.07, "inches"))) +
    scale_colour_manual(values = c("0" = band_dark, "1" = band_light)) +
    scale_x_continuous(breaks = axis_df$centre, labels = axis_df$Chr,
                       expand = expansion(mult = .01)) +
    scale_y_continuous(expand = expansion(mult = c(0, head_room))) +
    annotate("text", x = -Inf, y = Inf,
             label = sprintf("%s  (%s)", trait, model),
             hjust = -0.06, vjust = 1.5, size = 2.9,
             fontface = "bold", colour = "grey15") +
    {if (!is.na(gene_lab))
       annotate("text", x = px, y = a_hi, label = gene_lab, vjust = -0.45,
                hjust = label_hjust(px), size = 2.8, fontface = "bold",
                colour = marker_col) else NULL} +
    coord_cartesian(clip = "off") +
    labs(x = if (show_x) "Chromosome" else NULL,
         y = expression(bold(-log[10](italic(p))))) +
    theme_gwas +
    theme(axis.text.x = if (show_x) element_text(size = 7.5) else element_blank())
}

## MEMORY: eight full scans (8 x 4.18M points) do not fit in memory as ggplot
## objects at once, so each panel is rendered to its own PNG at the final size
## and resolution, then the eight PNGs are tiled. Nothing is thinned; the tiles
## are placed 1:1 (11 x 11 in at 300 dpi -> 5.5 x 2.75 in per panel).
PANEL_W <- 11 / length(COLUMNS); PANEL_H <- 11 / 4; DPI <- 300
tmpdir <- tempfile("panels_"); dir.create(tmpdir)
tiles <- list()
for (ci in seq_along(COLUMNS)) {
  cl <- COLUMNS[[ci]]
  for (i in seq_along(cl$panels)) {
    p <- cl$panels[[i]]; message("  ", p[1], " / ", p[2])
    g <- panel(p[1], p[2],
               show_x   = (i == length(cl$panels)),
               ## label the gene once, on the very first panel
               gene_lab = if (i == 1 && ci == 1) sprintf("%s  %s", GENE_LAB, GENE) else NA_character_)
    f <- file.path(tmpdir, sprintf("c%d_p%d.png", ci, i))
    ggsave(f, g, width = PANEL_W, height = PANEL_H, dpi = DPI, bg = "white", device = ragg::agg_png)
    rm(g); invisible(gc())
    tiles[[length(tiles) + 1]] <- list(ci = ci, i = i, f = f)
  }
}
fig <- wrap_plots(lapply(seq_along(COLUMNS), function(ci)
  wrap_plots(lapply(Filter(function(t) t$ci == ci, tiles), function(t)
    wrap_elements(full = grid::rasterGrob(png::readPNG(t$f), interpolate = FALSE))), ncol = 1)), nrow = 1)
dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)
ggsave(OUT_FILE, fig & theme(plot.margin = margin(0, 0, 0, 0)), width = 11, height = 11, dpi = DPI, bg = "white", limitsize = FALSE, device = ragg::agg_png)
message("Saved: ", OUT_FILE)
