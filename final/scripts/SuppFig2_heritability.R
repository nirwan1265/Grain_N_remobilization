#!/usr/bin/env Rscript
################################################################################
## SUPPLEMENTARY FIGURE -- heritability of grain amino acid phenotypes
##
## A  Broad-sense H2 (log scale, line-mean basis) for free and protein-bound
##    amino acids, ordered within each class. H2 comes from ~1,000 plot records
##    (2 years x 2 blocks) and is tight, so it is drawn as points without
##    intervals; the Cullis generalized H2 agrees to within 0.018 across every
##    trait and is not shown separately.
##
## B  Genomic h2 against broad-sense H2 with the y = x line. Points above the
##    line are NOT a contradiction: H2 is estimated from ~1,000 plots, h2 from
##    272 lines with a median 95% profile interval 0.52 wide, so the point
##    estimates cross freely. The panel therefore carries h2 error bars, and the
##    caption reports that H2 falls inside the h2 interval for 73% of traits.
##
## Traits whose PERMUTED phenotype returned a non-zero h2 (reliable = FALSE in
## H2_vs_h2_combined.csv; 15 of 116, all heavy-tailed PBAA ratio/contrast
## traits) are excluded from panel B and from the medians. Their H2 is
## unaffected and is still shown in panel A.
##
## Input : final/results/_stage/H2_FAA.csv, H2_PBAA.csv, H2_vs_h2_combined.csv
## Output: final/supp_figs/SuppFig2_heritability.png / .pdf
################################################################################

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

ROOT <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "."
IN   <- file.path(ROOT, "final", "results", "_stage")
OUT  <- file.path(ROOT, "final", "supp_figs")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- house theme (as in scripts/10_Fig4_Fst_Tajima_Pi.R) ---------------------
plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title   = element_text(size = 14, face = "bold", hjust = 0.5,
                                margin = margin(b = 10)),
    plot.tag     = element_text(size = 24, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x  = element_text(size = 16, colour = "black"),
    axis.text.y  = element_text(size = 16, face = "bold", colour = "black"),
    axis.line    = element_line(colour = "black"),
    panel.grid   = element_blank(),
    legend.background = element_rect(fill = "white", colour = "grey70", linewidth = 0.4),
    legend.title = element_blank(),
    legend.text  = element_text(size = 13),
    plot.margin  = margin(15, 15, 15, 15)
  )

COL_FAA  <- "#1B7837"
COL_PBAA <- "#762A83"

## ---- data --------------------------------------------------------------------
faa  <- fread(file.path(IN, "H2_FAA.csv"))[,  class := "Free (FAA)"]
pbaa <- fread(file.path(IN, "H2_PBAA.csv"))[, class := "Protein-bound (PBAA)"]
h2   <- rbind(faa, pbaa, fill = TRUE)
h2   <- h2[is.finite(H2_log)]

comb <- fread(file.path(IN, "H2_vs_h2_combined.csv"))
comb[, class := fifelse(source == "FAA", "Free (FAA)", "Protein-bound (PBAA)")]
comb_ok <- comb[reliable %in% c(TRUE, "True", "TRUE") & is.finite(h2_log) & is.finite(H2_log)]

cat(sprintf("panel A traits: %d FAA, %d PBAA\n",
            h2[class == "Free (FAA)", .N], h2[class == "Protein-bound (PBAA)", .N]))
cat(sprintf("panel B traits (reliable only): %d of %d\n", nrow(comb_ok), nrow(comb)))
cat(sprintf("median H2(log): FAA %.3f | PBAA %.3f\n",
            h2[class == "Free (FAA)", median(H2_log)],
            h2[class == "Protein-bound (PBAA)", median(H2_log)]))
inside <- comb_ok[H2_log >= h2_lo & H2_log <= h2_hi, .N]
cat(sprintf("H2 inside the h2 95%% interval: %d of %d (%.0f%%)\n",
            inside, nrow(comb_ok), 100 * inside / nrow(comb_ok)))

## ---- A: H2 per trait ---------------------------------------------------------
setorder(h2, class, H2_log)
h2[, trait_f := factor(trait, levels = trait)]

pA <- ggplot(h2, aes(H2_log, trait_f, colour = class)) +
  geom_segment(aes(x = 0, xend = H2_log, yend = trait_f), linewidth = .5, alpha = .45) +
  geom_point(size = 2.4) +
  facet_wrap(~ class, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = c("Free (FAA)" = COL_FAA,
                                 "Protein-bound (PBAA)" = COL_PBAA)) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, .04))) +
  labs(x = expression(bold(Broad*-sense~italic(H)^2~(log~scale))), y = NULL, tag = "A") +
  plot_theme +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 7, face = "plain"),
        strip.text  = element_text(size = 15, face = "bold"))

## ---- B: h2 vs H2 -------------------------------------------------------------
pB <- ggplot(comb_ok, aes(H2_log, h2_log, colour = class)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey35", linewidth = .7) +
  geom_errorbar(aes(ymin = h2_lo, ymax = h2_hi), width = 0, linewidth = .35, alpha = .28) +
  geom_point(size = 2.6, alpha = .9) +
  scale_colour_manual(values = c("Free (FAA)" = COL_FAA,
                                 "Protein-bound (PBAA)" = COL_PBAA)) +
  scale_x_continuous(limits = c(0, 1.02), expand = expansion(mult = .02)) +
  scale_y_continuous(limits = c(0, 1.02), expand = expansion(mult = .02)) +
  labs(x = expression(bold(Broad*-sense~italic(H)^2)),
       y = expression(bold(Genomic~italic(h)^2)), tag = "B") +
  plot_theme +
  ## bottom-right is the one empty quadrant: h2 is never high where H2 is
  ## high and the cloud sits left of centre, so the legend does not cover data
  theme(legend.position = c(0.98, 0.02),
        legend.justification = c(1, 0),
        legend.text = element_text(size = 14))

fig <- pA + pB + plot_layout(widths = c(1.35, 1))
ggsave(file.path(OUT, "SuppFig2_heritability.png"), fig,
       width = 18, height = 11, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(OUT, "SuppFig2_heritability.pdf"), fig,
       width = 18, height = 11, bg = "white", limitsize = FALSE)
message("Saved: ", file.path(OUT, "SuppFig2_heritability.png"))
