#!/usr/bin/env Rscript
################################################################################
## FIGURE -- one candidate locus recurring across eight amino acid phenotypes
##
## Zm00001eb393450 (lysine-rich arabinogalactan protein 17), chr9:130,359,239-
## 130,360,506 in B73 RefGen_v5, carries a Bonferroni-significant association in
## eight of the 76 amino acid phenotypes. Eight phenotypes is too many for one
## column, so they run four deep in two columns, ordered by significance down
## column A and then down column B, and the recurring peak reads straight down
## each column.
##
## There are no panel letters. All eight plots show the same quantity for the
## same locus, differing only in phenotype, so lettering them would imply a
## grouping that does not exist and would have to be referenced one by one in
## the text. Each plot names its own phenotype in the corner, and the model that
## produced the association, because the models were not all run for every
## phenotype -- E.Total is BLINK, the rest FarmCPU, and a reader comparing plots
## needs to know that. The two columns are a space-saving layout, nothing more.
##
## Input : final/results/manhattan_thin/<trait>__<MODEL>.csv  (Chr, Pos, P)
##         produced on the HPC by 32_thin_for_manhattan.R
## Output: final/main_figs/Fig_amino_Zm00001eb393450_manhattan.png
################################################################################

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

ROOT     <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "."
IN_DIR   <- file.path(ROOT, "final", "results", "manhattan_thin")
OUT_FILE <- file.path(ROOT, "final", "main_figs", "Fig_amino_Zm00001eb393450_manhattan.png")

## THRESHOLD: Benjamini-Hochberg FDR 5%, the threshold used throughout the
## paper. BH gives each phenotype x model its own cutoff, so each panel draws
## its own line rather than one shared line. In practice the nine cutoffs span
## -log10 6.05 to 6.97, under one log unit, so the panels stay comparable by
## eye. The Bonferroni and suggestive lines of the previous version are gone:
## neither threshold is used anywhere in the manuscript any more.
FDR_DIR <- file.path(ROOT, "final", "results", "goodman_buckler_amino_acid_gwas",
                     "by_model_fdr", "by_model")

fdr_cutoff <- function(trait, model) {
  f <- file.path(FDR_DIR, sprintf("%s_%s.RData", trait, model))
  if (!file.exists(f)) stop("no FDR file for ", trait, " / ", model, ": ", f)
  e <- new.env(parent = emptyenv()); nm <- load(f, envir = e)
  d <- as.data.table(get(nm[1], envir = e))
  if (!nrow(d[fdr_sig == TRUE])) stop(trait, "/", model, ": no FDR-significant SNP")
  max(d[fdr_sig == TRUE]$P.value)
}

GENE       <- "Zm00001eb393450"
GENE_LAB   <- "Lysine-rich arabinogalactan protein 17"
GENE_CHR   <- 9L
GENE_POS   <- 130359872       # gene midpoint, v5
PEAK_WINDOW <- 250e3          # anchor the arrow on the lead SNP within this

## columns of four, ordered by significance
## The locus is FDR-significant for NINE phenotypes. Eight are shown, in two
## columns of four, ordered by the strength of the association here. N.E is the
## ninth and weakest (P = 1.2e-07 against 4.5e-13 for EHPRQ) and is left out of
## the figure; the caption says so. To add it, append c("N.E", "FarmCPU") and
## put the thinned scan in final/results/manhattan_thin/.
COLUMNS <- list(
  list(panels = list(
    c("EHPRQ",   "FarmCPU"), c("LAV", "FarmCPU"),
    c("A",       "FarmCPU"), c("Total", "FarmCPU"))),
  list(panels = list(
    c("IVL",     "FarmCPU"), c("T", "FarmCPU"),
    c("E.Total", "BLINK"),   c("S", "FarmCPU")))
)

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
  f <- file.path(IN_DIR, sprintf("%s__%s.csv", trait, model))
  if (!file.exists(f)) stop("missing: ", f, "\n  run 32_thin_for_manhattan.R first")
  d <- fread(f, showProgress = FALSE)
  d[is.finite(P) & P > 0 & Chr >= 1 & Chr <= 10]
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
  d <- merge(read_scan(trait, model), chrlen[, .(Chr, offset)], by = "Chr")
  d[, gpos := Pos + offset][, logp := -log10(P)]
  d[, band := factor(Chr %% 2)]

  near <- d[abs(gpos - marker_g) <= PEAK_WINDOW]
  px <- if (nrow(near)) near$gpos[which.max(near$logp)] else marker_g
  py <- if (nrow(near)) max(near$logp) else 0
  ytop <- max(d$logp)
  a_lo <- py + 0.05 * ytop
  a_hi <- py + 0.20 * ytop
  head_room <- if (is.na(gene_lab)) 0.26 else 0.36

  ggplot(d, aes(gpos, logp, colour = band)) +
    geom_point(size = .8, alpha = .85) +
    geom_hline(yintercept = -log10(fdr_cutoff(trait, model)), linetype = "dashed",
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

cols_built <- lapply(seq_along(COLUMNS), function(ci) {
  cl <- COLUMNS[[ci]]
  ps <- lapply(seq_along(cl$panels), function(i) {
    p <- cl$panels[[i]]
    message("  ", p[1], " / ", p[2])
    panel(p[1], p[2],
          show_x   = (i == length(cl$panels)),
          ## label the gene once, on the very first panel -- with three columns
          ## the old `i == 1` rule printed it three times
          gene_lab = if (i == 1 && ci == 1) sprintf("%s  %s", GENE_LAB, GENE) else NA_character_)
  })
  wrap_plots(ps, ncol = 1)
})

fig <- wrap_plots(cols_built, nrow = 1)
dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)
ggsave(OUT_FILE, fig, width = 11, height = 11, dpi = 300, bg = "white", limitsize = FALSE)
message("Saved: ", OUT_FILE)
