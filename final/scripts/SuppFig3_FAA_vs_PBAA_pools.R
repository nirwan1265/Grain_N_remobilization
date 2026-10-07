#!/usr/bin/env Rscript
################################################################################
## FREE vs PROTEIN-BOUND AMINO ACID POOLS -- Goodman-Buckler panel
##
## Rebuilds the FAA/PBAA comparison, fixing a definition mismatch in
## scripts/24_Goodman_Buckler_FAA_PBAA_metrics.R.
##
## THE MISMATCH
##   The previous "nitrogen-rich fraction" used
##       FAA  : Asn + Gln + Arg + Lys + His
##       PBAA : Asx + Glx + Arg + Lys + His
##   but acid hydrolysis deamidates Asn to Asp and Gln to Glu, so Asx = Asp+Asn
##   and Glx = Glu+Gln. The PBAA numerator therefore silently carries aspartate
##   and glutamate, which have one N each and are not nitrogen-rich, while the
##   FAA numerator excludes them. The two fractions were not the same quantity.
##
##   This script reports both, clearly separated:
##     amide_basic_fraction  Asx + Glx + Arg + Lys + His, with FAA collapsed to
##                           (Asp+Asn) and (Glu+Gln) so both pools are computed
##                           identically. THIS is the comparable metric.
##     N_rich_strict         Asn + Gln + Arg + Lys + His. Well defined for FAA
##                           only; PBAA cannot resolve it, so it is reported as
##                           an FAA descriptive and never compared.
##
## THE N-ATOM ASSUMPTION
##   The N proxy weights each amino acid by its N atoms. Asx and Glx are
##   unresolved, and the previous script used 1.5 for both -- a 50/50 split of
##   Asp/Asn and Glu/Gln. That is a guess, and maize storage proteins are
##   glutamine-rich, so 1.5 likely understates PBAA nitrogen. The proxy is
##   therefore reported across the full range of the assumption (1.0, 1.5, 2.0)
##   rather than at a single unverifiable value.
##
## UNITS -- UNRESOLVED, read before using any absolute magnitude
##   FAA (own assay, 20 amino acids) and PBAA (Shrestha 2022, 15 amino acids,
##   acid hydrolysate) are different measurements and no unit metadata is
##   recorded in either file. Any FAA-vs-PBAA ratio of ABSOLUTE amounts
##   (the x-fold statements, the "% of grain amino acid N") is only meaningful
##   if both assays report in the same units, which is not established.
##   Fractions and CVs are within-pool and unaffected.
##
## Input : final/results/_stage/{FAA,PBAA}_Goodman_Buckler.csv
## Output: final/supp_figs/SuppFig3_FAA_vs_PBAA_pools.png / .pdf
##         final/supp_tables/SuppTable_S4_FAA_vs_PBAA_pools.csv
################################################################################

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

ROOT  <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "."
STAGE <- file.path(ROOT, "final", "results", "_stage")
FIG   <- file.path(ROOT, "final", "supp_figs")
TAB   <- file.path(ROOT, "final", "supp_tables")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB, showWarnings = FALSE, recursive = TRUE)

plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title   = element_text(size = 17, face = "bold", hjust = 0.5, margin = margin(b = 8)),
    plot.tag     = element_text(size = 22, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    axis.title.x = element_text(size = 17, face = "bold"),
    axis.title.y = element_text(size = 17, face = "bold"),
    axis.text.x  = element_text(size = 15, face = "bold", colour = "black"),
    axis.text.y  = element_text(size = 15, face = "bold", colour = "black"),
    axis.line    = element_line(colour = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
    legend.position = "none",
    plot.margin = margin(15, 15, 15, 15))

COL <- c(FAA = "#2166AC", PBAA = "#B35806")

FAA_AA  <- c("A","R","N","D","Q","E","G","H","I","L","K","M","F","P","S","W","T","Y","V","C")
PBAA_AA <- c("Ala","Arg","Asx","Glx","His","Ile","Leu","Lys","Met","Phe","Pro","Ser","Thr","Tyr","Val")

faa_n  <- c(A=1,R=4,N=2,D=1,Q=2,E=1,G=1,H=3,I=1,L=1,K=2,M=1,F=1,P=1,S=1,W=2,T=1,Y=1,V=1,C=1)

rd <- function(f) { d <- fread(f, showProgress = FALSE)
                    setnames(d, sub("^﻿", "", names(d))); d }

faa  <- rd(file.path(STAGE, "FAA_Goodman_Buckler.csv"))
pbaa <- rd(file.path(STAGE, "PBAA_Goodman_Buckler.csv"))

## ---- genotype means ----------------------------------------------------------
F <- faa[,  lapply(.SD, mean, na.rm = TRUE), by = taxa, .SDcols = FAA_AA]
PB <- pbaa[, lapply(.SD, mean, na.rm = TRUE), by = taxa, .SDcols = PBAA_AA]  # not P: F$P is proline

F[, Total := rowSums(.SD, na.rm = TRUE), .SDcols = FAA_AA]
PB[, Total := rowSums(.SD, na.rm = TRUE), .SDcols = PBAA_AA]

## collapse FAA to the hydrolysate's resolution so both are computed identically
F[, `:=`(Asx_equiv = D + N, Glx_equiv = E + Q)]

F[,  Proline_fraction := P / Total]      # P here is the proline column of F
PB[, Proline_fraction := Pro / Total]
F[, amide_basic_fraction := (Asx_equiv + Glx_equiv + R + K + H) / Total]
PB[, amide_basic_fraction := (Asx + Glx + Arg + Lys + His) / Total]
F[, N_rich_strict_fraction := (N + Q + R + K + H) / Total]   # FAA only

## N proxy under the three Asx/Glx assumptions
F[, N_proxy := as.matrix(.SD) %*% faa_n[FAA_AA], .SDcols = FAA_AA]
for (a in c(1.0, 1.5, 2.0)) {
  w <- c(Ala=1,Arg=4,Asx=a,Glx=a,His=3,Ile=1,Leu=1,Lys=2,Met=1,Phe=1,Pro=1,Ser=1,Thr=1,Tyr=1,Val=1)
  PB[[sprintf("N_proxy_%.1f", a)]] <- as.vector(as.matrix(PB[, ..PBAA_AA]) %*% w[PBAA_AA])
}

m <- merge(F[, .(taxa, FAA_total = Total, FAA_proline = Proline_fraction,
                 FAA_amide_basic = amide_basic_fraction,
                 FAA_N_rich_strict = N_rich_strict_fraction, FAA_N_proxy = N_proxy)],
           PB[, .(taxa, PBAA_total = Total, PBAA_proline = Proline_fraction,
                 PBAA_amide_basic = amide_basic_fraction,
                 PBAA_N_proxy_1.0 = N_proxy_1.0, PBAA_N_proxy_1.5 = N_proxy_1.5,
                 PBAA_N_proxy_2.0 = N_proxy_2.0)],
           by = "taxa")

cat(sprintf("matched taxa: %d\n\n", nrow(m)))

## ---- comparisons -------------------------------------------------------------
cmp <- function(fa, pb, label) {
  a <- m[[fa]]; b <- m[[pb]]
  ok <- is.finite(a) & is.finite(b)
  w <- suppressWarnings(wilcox.test(a[ok], b[ok], paired = TRUE))
  data.table(comparison = label, n_taxa = sum(ok),
             FAA_mean = mean(a[ok]),  FAA_median = median(a[ok]),
             FAA_sd = sd(a[ok]),      FAA_CV_pct = 100*sd(a[ok])/mean(a[ok]),
             PBAA_mean = mean(b[ok]), PBAA_median = median(b[ok]),
             PBAA_sd = sd(b[ok]),     PBAA_CV_pct = 100*sd(b[ok])/mean(b[ok]),
             fold_PBAA_over_FAA = mean(b[ok])/mean(a[ok]),
             paired_wilcoxon_p = w$p.value)
}

s4 <- rbindlist(list(
  cmp("FAA_total",       "PBAA_total",       "Total amino acid pool"),
  cmp("FAA_N_proxy",     "PBAA_N_proxy_1.5", "Amino acid-derived N proxy (Asx/Glx = 1.5 N)"),
  cmp("FAA_proline",     "PBAA_proline",     "Proline fraction"),
  cmp("FAA_amide_basic", "PBAA_amide_basic", "Asx + Glx + basic fraction (matched definition)")))
fwrite(s4, file.path(TAB, "SuppTable_S4_FAA_vs_PBAA_pools.csv"))

for (i in seq_len(nrow(s4))) { x <- s4[i]
  cat(sprintf("%-48s FAA %9.3f (CV %4.1f%%) | PBAA %9.3f (CV %4.1f%%) | %5.1fx | p=%.2g\n",
      x$comparison, x$FAA_mean, x$FAA_CV_pct, x$PBAA_mean, x$PBAA_CV_pct,
      x$fold_PBAA_over_FAA, x$paired_wilcoxon_p)) }

cat(sprintf("\nFAA-only strict N-rich fraction (Asn+Gln+Arg+Lys+His): %.1f%%\n",
            100*mean(m$FAA_N_rich_strict)))
cat("\nN proxy sensitivity to the Asx/Glx N-atom assumption:\n")
fa_mu <- mean(m$FAA_N_proxy, na.rm = TRUE)
for (a in c("1.0","1.5","2.0")) {
  pb_mu <- mean(m[[paste0("PBAA_N_proxy_", a)]], na.rm = TRUE)
  cat(sprintf("   Asx/Glx = %s N : PBAA proxy %8.1f | PBAA/FAA %5.1fx | FAA share %4.2f%%\n",
      a, pb_mu, pb_mu / fa_mu, 100 * fa_mu / (fa_mu + pb_mu)))
}
cat(sprintf("\n   (the FAA share of amino acid-derived N moves only %.2f-%.2f%% across the\n",
            100*fa_mu/(fa_mu+mean(m$PBAA_N_proxy_2.0, na.rm=TRUE)),
            100*fa_mu/(fa_mu+mean(m$PBAA_N_proxy_1.0, na.rm=TRUE))))
cat("    full range of the Asx/Glx assumption, so that number is robust to it --\n")
cat("    it is the UNIT question, not this one, that governs whether it is usable)\n")

## ---- figure ------------------------------------------------------------------
pan <- function(fa, pb, title, ylab, logy, tag) {
  d <- rbind(data.table(pool = "FAA",  v = m[[fa]]),
             data.table(pool = "PBAA", v = m[[pb]]))
  d[, pool := factor(pool, levels = c("FAA", "PBAA"))]
  p <- ggplot(d, aes(pool, v, colour = pool, fill = pool)) +
    geom_violin(alpha = .18, linewidth = .7) +
    geom_boxplot(width = .18, outlier.shape = NA, alpha = .5, linewidth = .6) +
    geom_jitter(width = .10, size = .9, alpha = .35) +
    scale_colour_manual(values = COL) + scale_fill_manual(values = COL) +
    labs(title = title, x = NULL, y = ylab, tag = tag) + plot_theme
  if (logy) p <- p + scale_y_log10() 
  p
}

fig <- (pan("FAA_total","PBAA_total","Total amino acid pool","Assay-scale abundance",TRUE,"A") |
        pan("FAA_N_proxy","PBAA_N_proxy_1.5","Amino acid-derived N proxy","Estimated N units",TRUE,"B")) /
       (pan("FAA_proline","PBAA_proline","Proline fraction","Fraction of pool",FALSE,"C") |
        pan("FAA_amide_basic","PBAA_amide_basic","Asx + Glx + basic fraction","Fraction of pool",FALSE,"D"))

ggsave(file.path(FIG, "SuppFig3_FAA_vs_PBAA_pools.png"), fig,
       width = 15, height = 13, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(FIG, "SuppFig3_FAA_vs_PBAA_pools.pdf"), fig,
       width = 15, height = 13, bg = "white", limitsize = FALSE)
message("Saved: ", file.path(FIG, "SuppFig3_FAA_vs_PBAA_pools.png"))
