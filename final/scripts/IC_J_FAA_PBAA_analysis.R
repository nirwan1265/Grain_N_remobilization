#!/usr/bin/env Rscript
################################################################################
## INDIAN CHIEF / JARVIS -- free and protein-bound amino acids, C0 vs C14
##
## Replaces scripts/10_Fig4_Fst_Tajima_Pi.R (FAA part), 13_IC_J_aa_plots.R,
## 14_IC_J_faa_derived_metrics.R and 18_SuppFig_AA_N_pools.R for the paper.
##
## Design: each biological sample was measured as 3 technical replicates; the
## replicates are averaged first, so every test is on BIOLOGICAL samples:
##   Indian Chief C0 n = 5, C14 n = 8;  Jarvis C0 n = 3, C14 n = 8.
## Test: two-sided Wilcoxon rank-sum (normal approximation, exact = FALSE, as in
## the original scripts), C0 vs C14 within each population.
## Smallest attainable P (complete separation): 5 vs 8 -> 0.0043; 3 vs 8 -> 0.019.
##
## Multiple testing (Benjamini-Hochberg, within population):
##   FAA metrics (Fig. 3, Table S14): 8 pre-specified metrics, all per mg of
##     sample -- total FAA per mg; Pro, Asn, Asp, Gln, Glu fractions of total
##     FAA; N-rich FAA (Asn+Gln+Arg+Lys+His) per mg; FAA-derived N per mg.
##     (The earlier 9-metric set counted the FAA-N pool twice, per sample and
##     per mg; with everything per mg they are one metric.)
##   Per-amino-acid relative composition (Table S15): BH over the 20 FAAs or
##     the 16 PBAAs.
##   PBAA metrics (Table S16): 6 metrics (total PBAA per mg; Pro, Asx, Glx
##     fractions; Asx+Glx+Arg+Lys+His per mg; estimated PBAA-derived N per mg).
##   Amino acid-derived N per mg (Supp. Fig., Table S17): 3 metrics (FAA, PBAA,
##     total).
## N atoms: Arg 4, His 3, Asn/Gln/Lys/Trp 2, others 1; PBAA Asx/Glx 1.5
## (unresolved Asp/Asn and Glu/Gln after acid hydrolysis).
##
## Inputs : data/FAA_IC_J.csv, data/PBAA_IC_J.csv
## Outputs: final/main_figs/Fig3_IC_J_FAA.{png,pdf}
##          final/supp_figs/SuppFig_IC_J_PBAA_composition.{png,pdf}
##          final/supp_figs/SuppFig_IC_J_AA_N_pools.{png,pdf}
##          final/supp_tables/SuppTable_S14..S17_IC_J_*.csv
################################################################################
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
gv   <- function(k, d) { v <- Sys.getenv(k, ""); if (nzchar(v)) v else d }
REPO <- gv("GN_REPO", "/Users/nirwantandukar/Documents/Github/Grain_N_remobilization")
F_FAA  <- gv("GN_FAA",  file.path(REPO, "data/FAA_IC_J.csv"))
F_PBAA <- gv("GN_PBAA", file.path(REPO, "data/PBAA_IC_J.csv"))
MAIN <- gv("GN_MAINDIR", file.path(REPO, "final/main_figs"))
SUPP <- gv("GN_SUPPDIR", file.path(REPO, "final/supp_figs"))
TAB  <- gv("GN_TABDIR",  file.path(REPO, "final/supp_tables"))
for (d in c(MAIN, SUPP, TAB)) dir.create(d, FALSE, TRUE)

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
GRP_COL <- c("Indian Chief C0" = "#92C5DE", "Indian Chief C14" = "#2166AC",
             "Jarvis C0" = "#FDB863", "Jarvis C14" = "#B35806")
legend_inside <- function(x, y, just) {
  if (utils::packageVersion("ggplot2") >= "3.5.0")
    theme(legend.position = "inside", legend.position.inside = c(x, y), legend.justification = just)
  else theme(legend.position = c(x, y), legend.justification = just)
}
stars <- function(p) fifelse(is.na(p), "", fifelse(p < 0.001, "***", fifelse(p < 0.01, "**",
                      fifelse(p < 0.05, "*", "ns"))))

## ---- read and average technical replicates -----------------------------------
bio_means <- function(f) {
  d <- fread(f, encoding = "UTF-8")
  setnames(d, 1:2, c("ID", "wt_mg"))
  stopifnot(all(grepl("^(IC|JV)(0|14)_[0-9]+_[0-9]+$", d$ID)))
  d[, `:=`(pop = fifelse(substr(ID, 1, 2) == "IC", "Indian Chief", "Jarvis"),
           cyc = fifelse(sub("^(IC|JV)(0|14)_.*", "\\2", ID) == "0", "C0", "C14"),
           bio = sub("^((IC|JV)(0|14)_[0-9]+)_[0-9]+$", "\\1", ID))]
  aa <- setdiff(names(d), c("ID", "wt_mg", "pop", "cyc", "bio"))
  b <- d[, lapply(.SD, mean, na.rm = TRUE), by = .(pop, cyc, bio), .SDcols = c("wt_mg", aa)]
  b[, group := paste(pop, cyc)]
  list(b = b[], aa = aa)
}
F <- bio_means(F_FAA); P <- bio_means(F_PBAA)
fa <- F$b; pa <- P$b; faa <- F$aa; paa <- P$aa
cat("biological samples:\n"); print(fa[, .N, by = group])

Natom <- c(Ala=1, Arg=4, Asn=2, Asp=1, Gln=2, Glu=1, Gly=1, His=3, Ile=1, Leu=1, Lys=2,
           Met=1, Phe=1, Pro=1, Ser=1, Trp=2, Thr=1, Tyr=1, Val=1, Cys=1, Asx=1.5, Glx=1.5)

fa[, Total := rowSums(.SD), .SDcols = faa]
fa[, `:=`(Total_FAA_per_mg = Total / wt_mg,
          Pro_fraction = Pro / Total, Asn_fraction = Asn / Total, Asp_fraction = Asp / Total,
          Gln_fraction = Gln / Total, Glu_fraction = Glu / Total,
          N_rich_FAA_per_mg = (Asn + Gln + Arg + Lys + His) / wt_mg,
          FAA_N_per_mg = as.numeric(as.matrix(.SD) %*% Natom[faa]) / wt_mg), .SDcols = faa]
pa[, Total := rowSums(.SD), .SDcols = paa]
pa[, `:=`(Total_PBAA_per_mg = Total / wt_mg,
          Pro_fraction = Pro / Total, Asx_fraction = Asx / Total, Glx_fraction = Glx / Total,
          Asx_Glx_basic_per_mg = (Asx + Glx + Arg + Lys + His) / wt_mg,
          PBAA_N_per_mg = as.numeric(as.matrix(.SD) %*% Natom[paa]) / wt_mg), .SDcols = paa]

## ---- tests -------------------------------------------------------------------
cliff <- function(x, y) mean(outer(x, y, ">")) - mean(outer(x, y, "<"))   # C14 vs C0
test_set <- function(d, metrics, desc) {
  out <- rbindlist(lapply(c("Indian Chief", "Jarvis"), function(pp) {
    r <- rbindlist(lapply(metrics, function(m) {
      x0 <- d[pop == pp & cyc == "C0"][[m]]; x1 <- d[pop == pp & cyc == "C14"][[m]]
      data.table(population = pp, metric = m, description = desc[m],
                 n_C0 = length(x0), n_C14 = length(x1),
                 mean_C0 = mean(x0), mean_C14 = mean(x1),
                 median_C0 = median(x0), median_C14 = median(x1),
                 P = wilcox.test(x0, x1, exact = FALSE)$p.value,
                 cliffs_delta = cliff(x1, x0))
    }))
    r[, P_BH := p.adjust(P, "BH")]
  }))
  out[, direction := fifelse(mean_C14 > mean_C0, "higher in C14", "lower in C14")]
  out[, significant_BH := P_BH < 0.05][]
}
faa_metrics <- c(Total_FAA_per_mg = "Total FAA per mg sample",
                 Pro_fraction = "Proline / total FAA", Asn_fraction = "Asparagine / total FAA",
                 Asp_fraction = "Aspartate / total FAA", Gln_fraction = "Glutamine / total FAA",
                 Glu_fraction = "Glutamate / total FAA",
                 N_rich_FAA_per_mg = "Asn + Gln + Arg + Lys + His per mg sample",
                 FAA_N_per_mg = "FAA-derived N per mg sample")
S13 <- test_set(fa, names(faa_metrics), faa_metrics)

relcomp <- function(d, aa, pool) {
  r <- copy(d); r[, (aa) := lapply(.SD, function(v) v / r$Total), .SDcols = aa]
  s <- test_set(r, aa, setNames(paste(aa, "/ total", pool), aa))
  s[, pool := pool][, metric := NULL]; setnames(s, "description", "amino_acid_fraction")
  s[, amino_acid := sub(" /.*", "", amino_acid_fraction)]
  list(stats = s, rel = r)
}
RF <- relcomp(fa, faa, "FAA"); RP <- relcomp(pa, paa, "PBAA")
S14 <- rbind(RF$stats, RP$stats)
setcolorder(S14, c("pool", "population", "amino_acid"))

pbaa_metrics <- c(Total_PBAA_per_mg = "Total PBAA per mg sample",
                  Pro_fraction = "Proline / total PBAA", Asx_fraction = "Asx / total PBAA",
                  Glx_fraction = "Glx / total PBAA",
                  Asx_Glx_basic_per_mg = "Asx + Glx + Arg + Lys + His per mg sample",
                  PBAA_N_per_mg = "Estimated PBAA-derived N per mg sample (Asx, Glx = 1.5 N)")
S15 <- test_set(pa, names(pbaa_metrics), pbaa_metrics)

nn <- merge(fa[, .(pop, cyc, bio, group, FAA_N_per_mg)], pa[, .(bio, PBAA_N_per_mg)], by = "bio")
nn[, Total_AA_N_per_mg := FAA_N_per_mg + PBAA_N_per_mg]
n_metrics <- c(FAA_N_per_mg = "FAA-derived N per mg sample",
               PBAA_N_per_mg = "Estimated PBAA-derived N per mg sample",
               Total_AA_N_per_mg = "Estimated total amino acid-derived N per mg sample")
S16 <- test_set(nn, names(n_metrics), n_metrics)

for (x in list(S13, S14, S15, S16)) {
  num <- names(x)[sapply(x, is.numeric) & !names(x) %in% c("n_C0", "n_C14")]
  x[, (num) := lapply(.SD, signif, 4), .SDcols = num]
}
fwrite(S13, file.path(TAB, "SuppTable_S14_IC_J_FAA_metrics.csv"))
fwrite(S14, file.path(TAB, "SuppTable_S15_IC_J_amino_acid_composition.csv"))
fwrite(S15, file.path(TAB, "SuppTable_S16_IC_J_PBAA_metrics.csv"))
fwrite(S16, file.path(TAB, "SuppTable_S17_IC_J_AA_N_per_mg.csv"))
cat("\nS13 FAA metrics:\n"); print(S13[, .(population, metric, n_C0, n_C14, mean_C0, mean_C14, P, P_BH)])
cat("\nS14 nominal (P < 0.05):\n"); print(S14[P < 0.05, .(pool, population, amino_acid, P, P_BH)])
cat("\nS15 PBAA metrics:\n"); print(S15[, .(population, metric, mean_C0, mean_C14, P, P_BH)])
cat("\nS16 N per mg:\n"); print(S16[, .(population, metric, mean_C0, mean_C14, P, P_BH)])

## ---- shared plotting helpers -------------------------------------------------
lvl <- names(GRP_COL)
box_panel <- function(d, ycol, ylab, title, stats, tag = NULL) {
  dd <- d[, .(group = factor(group, levels = lvl), v = get(ycol))]
  top <- max(dd$v) * 1.08
  br <- stats[metric == ycol][, .(x0 = match(paste(population, "C0"), lvl),
                                  lab = stars(P_BH))][, `:=`(x1 = x0 + 1, y = top)]
  ggplot(dd, aes(group, v, colour = group, fill = group)) +
    geom_boxplot(width = 0.6, alpha = 0.25, linewidth = 0.9, outlier.shape = NA) +
    geom_jitter(width = 0.12, size = 3.2, alpha = 0.85) +
    geom_segment(data = br, aes(x = x0, xend = x1, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.6) +
    geom_text(data = br, aes(x = (x0 + x1) / 2, y = y, label = lab), inherit.aes = FALSE,
              vjust = -0.3, size = 6, fontface = "bold") +
    scale_colour_manual(values = GRP_COL) + scale_fill_manual(values = GRP_COL) +
    scale_x_discrete(labels = function(z) sub(" C", "\nC", z)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.14))) +
    labs(x = NULL, y = ylab, title = title, tag = tag) + plot_theme +
    theme(legend.position = "none", axis.text.x = element_text(size = 14, face = "bold"),
          plot.title = element_text(size = 17, face = "bold", hjust = 0.5))
}

## ---- Figure 3 ----------------------------------------------------------------
pA <- box_panel(fa, "Total_FAA_per_mg", "Total FAA per mg sample", "Total free amino acids", S13)

frac <- c(Pro_fraction = "Pro", Asp_fraction = "Asp", Glu_fraction = "Glu",
          Asn_fraction = "Asn", Gln_fraction = "Gln")
fl <- melt(fa, id.vars = c("pop", "cyc", "group"), measure.vars = names(frac),
           variable.name = "metric", value.name = "v")
fl[, `:=`(aa = factor(frac[as.character(metric)], levels = frac),
          group = factor(group, levels = lvl))]
fs <- S13[metric %in% names(frac)][, aa := factor(frac[metric], levels = frac)]
ytop <- fl[, .(y = max(v)), by = .(aa, pop)]
fs <- merge(fs, ytop, by.x = c("aa", "population"), by.y = c("aa", "pop"))
fs[, lab := stars(P_BH)]
fs[, xpos := as.integer(aa) + fifelse(population == "Indian Chief", -0.2, 0.2)]
fs[, y := y + 0.025]
pB <- ggplot(fl, aes(aa, v, colour = group, fill = group)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.25,
               linewidth = 0.8, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.8),
             size = 2.4, alpha = 0.85) +
  geom_text(data = fs, aes(x = xpos, y = y, label = lab), inherit.aes = FALSE,
            size = 5.5, fontface = "bold") +
  scale_colour_manual(values = GRP_COL) + scale_fill_manual(values = GRP_COL) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.08))) +
  labs(x = NULL, y = "Fraction of total FAA", title = "Focal amino acid fractions") +
  plot_theme +
  legend_inside(0.98, 0.98, c(1, 1)) +
  theme(axis.text.x = element_text(size = 17, face = "bold"),
        plot.title = element_text(size = 17, face = "bold", hjust = 0.5)) +
  guides(colour = guide_legend(ncol = 2), fill = guide_legend(ncol = 2))

fig3 <- (pA | pB) + plot_layout(widths = c(1, 2.1)) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 28, face = "bold"))
for (ext in c("png", "pdf"))
  ggsave(file.path(MAIN, paste0("Fig3_IC_J_FAA.", ext)), fig3, width = 21, height = 8.5,
         dpi = 300, bg = "white")

## ---- Supp figure: amino acid-derived N per mg -----------------------------------
nA <- box_panel(nn, "FAA_N_per_mg", "Estimated N units per mg", "FAA-derived N", S16)
nB <- box_panel(nn, "PBAA_N_per_mg", "Estimated N units per mg", "Protein-bound amino acid-derived N", S16)
nC <- box_panel(nn, "Total_AA_N_per_mg", "Estimated N units per mg", "Total amino acid-derived N", S16)
figN <- (nA | nB | nC) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 28, face = "bold"))
for (ext in c("png", "pdf"))
  ggsave(file.path(SUPP, paste0("SuppFig_IC_J_AA_N_pools.", ext)), figN, width = 22, height = 8,
         dpi = 300, bg = "white")

## ---- Supp figure: PBAA composition, straight-edged radar ----------------------
## drawn in Cartesian coordinates so every edge (data, rings, spokes) is straight
radar <- function(pp, tag) {
  n <- length(paa); th <- pi / 2 - 2 * pi * (seq_len(n) - 1) / n
  g <- RP$rel[pop == pp, lapply(.SD, mean), by = .(group), .SDcols = paa]
  gl <- melt(g, id.vars = "group", variable.name = "aa", value.name = "v")
  gl[, i := match(as.character(aa), paa)][, `:=`(x = v * cos(th[i]), y = v * sin(th[i]))]
  setorder(gl, group, i)
  rmax <- max(gl$v) * 1.12
  rings <- pretty(c(0, rmax), 4); rings <- rings[rings > 0 & rings <= rmax]
  web <- rbindlist(lapply(rings, function(r) data.table(r = r, x = r * cos(c(th, th[1])),
                                                         y = r * sin(c(th, th[1])))))
  spokes <- data.table(x = 0, y = 0, xend = rmax * cos(th), yend = rmax * sin(th))
  sig <- S14[pool == "PBAA" & population == pp & P_BH < 0.05]$amino_acid
  lab <- data.table(x = 1.13 * rmax * cos(th), y = 1.13 * rmax * sin(th),
                    l = paste0(paa, fifelse(paa %in% sig, "*", "")))
  ggplot() +
    geom_path(data = web, aes(x, y, group = r), colour = "grey82", linewidth = 0.45) +
    geom_segment(data = spokes, aes(x, y, xend = xend, yend = yend), colour = "grey88", linewidth = 0.4) +
    annotate("text", x = 0.02 * rmax, y = rings, label = format(rings), size = 4,
             colour = "grey40", hjust = 0, vjust = -0.3) +
    geom_polygon(data = gl, aes(x, y, group = group, colour = group, fill = group),
                 alpha = 0.06, linewidth = 1.3) +
    geom_point(data = gl, aes(x, y, colour = group), size = 2.4) +
    geom_text(data = lab, aes(x, y, label = l), size = 5.8, fontface = "bold") +
    scale_colour_manual(values = GRP_COL) + scale_fill_manual(values = GRP_COL) +
    coord_equal(xlim = c(-1.25, 1.25) * rmax, ylim = c(-1.25, 1.25) * rmax) +
    labs(x = NULL, y = NULL, title = pp, tag = tag) + plot_theme +
    theme(axis.line = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), legend.position = "bottom", legend.text = element_text(size = 15),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
}
figP <- (radar("Indian Chief", "A") | radar("Jarvis", "B")) &
  theme(plot.tag = element_text(size = 28, face = "bold"))
for (ext in c("png", "pdf"))
  ggsave(file.path(SUPP, paste0("SuppFig_IC_J_PBAA_composition.", ext)), figP, width = 18,
         height = 9.5, dpi = 300, bg = "white")
cat("\nwritten: Fig3, two supp figs, S13-S16\n")
