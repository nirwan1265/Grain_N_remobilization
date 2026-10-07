#!/usr/bin/env Rscript
################################################################################
### SUPPLEMENTARY FIGURE 1: FREE AMINO ACID COMPOSITION OF THE GOODMAN-BUCKLER
### PANEL
###
### Merges the Asn/Pro rank panels with the abundance/variability panel: all
### three describe the same thing -- which free amino acids dominate the kernel
### pool -- so they belong in one figure rather than three.
###
###   A  genotypes in which each FAA ranked first
###   B  genotypes in which Asn and Pro were the top two, either order
###   C  concentration (left) and coefficient of variation (right) per FAA
###
### Panel C carries one tag across its two sub-plots rather than splitting into
### C and D, following the convention used in the Fig3 Manhattan script: the two
### sub-plots are one statement about one set of traits, read left to right, and
### two letters would imply two findings.
###
### Everything is computed from the raw plot records. This replaces
### scripts/22_SuppFig_FAA_Asn_Pro_panel.R, in which the counts were typed in
### from the manuscript rather than derived from the data, and supersedes the
### separate SuppFig1_FAA_Asn_Pro_rank.R / SuppFig1_panelC_FAA_abundance.R pair.
###
###   Rscript SuppFig1_FAA_composition.R [path/to/FAA_Goodman_Buckler.csv] [outdir]
################################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(tidyr); library(ggplot2)
  library(patchwork); library(scales)
})

args    <- commandArgs(trailingOnly = TRUE)
in_csv  <- if (length(args) >= 1) args[1] else "data/FAA_Goodman_Buckler.csv"
out_dir <- if (length(args) >= 2) args[2] else "final/supp_figs"
tab_dir <- if (length(args) >= 3) args[3] else "final/supp_tables"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir, "SuppFig1_FAA_composition.png")

AA <- c("A","R","N","D","Q","E","G","H","I","L","K","M","F","P","S","W","T","Y","V","C")
FULL <- c(A="Ala", R="Arg", N="Asn", D="Asp", Q="Gln", E="Glu", G="Gly", H="His",
          I="Ile", L="Leu", K="Lys", M="Met", F="Phe", P="Pro", S="Ser", W="Trp",
          T="Thr", Y="Tyr", V="Val", C="Cys")

################################################################################
### COMPUTE EVERYTHING FROM THE DATA
################################################################################

raw <- read.csv(in_csv, check.names = FALSE)
names(raw) <- sub("^﻿", "", trimws(names(raw)))
stopifnot(all(c("taxa", AA) %in% names(raw)))

gmean <- raw %>%
  group_by(taxa) %>%
  summarise(across(all_of(AA), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

mat <- as.matrix(gmean[, AA]); rownames(mat) <- gmean$taxa
rk  <- t(apply(-mat, 1, rank, ties.method = "min", na.last = TRUE))
colnames(rk) <- colnames(mat)
top1 <- colnames(rk)[apply(rk, 1, which.min)]

n_genotypes <- nrow(mat)
n_asn  <- sum(top1 == "N", na.rm = TRUE)
n_pro  <- sum(top1 == "P", na.rm = TRUE)
n_oth  <- n_genotypes - n_asn - n_pro
n_both <- sum(rk[, "N"] <= 2 & rk[, "P"] <= 2, na.rm = TRUE)

cat("\n--- numbers for the Results text ---\n")
cat(sprintf("genotypes with FAA data            : %d\n", n_genotypes))
cat(sprintf("Asn ranked most abundant           : %d (%.0f%%)\n", n_asn, 100*n_asn/n_genotypes))
cat(sprintf("Pro ranked most abundant           : %d (%.0f%%)\n", n_pro, 100*n_pro/n_genotypes))
cat(sprintf("some other FAA ranked first        : %d (%.0f%%)\n", n_oth, 100*n_oth/n_genotypes))
cat(sprintf("Asn and Pro the top two (any order): %d (%.0f%%)\n", n_both, 100*n_both/n_genotypes))

## per-amino-acid abundance and variability across genotypes
long <- as_tibble(mat, rownames = "taxa") %>%
  pivot_longer(all_of(AA), names_to = "aa", values_to = "conc") %>%
  filter(is.finite(conc), conc > 0)

summ <- long %>%
  group_by(aa) %>%
  summarise(n_taxa = n(), median = median(conc),
            q25 = quantile(conc, .25), q75 = quantile(conc, .75),
            p05 = quantile(conc, .05), p95 = quantile(conc, .95),
            mean = mean(conc), sd = sd(conc), .groups = "drop") %>%
  mutate(CV = 100 * sd / mean, label = FULL[aa],
         hilite = ifelse(aa == "N", "Asn", ifelse(aa == "P", "Pro", "other"))) %>%
  arrange(median) %>%
  mutate(aa_f = factor(label, levels = label))

write.csv(summ %>% arrange(desc(median)) %>%
            select(amino_acid = label, code = aa, n_taxa, median, q25, q75,
                   p05, p95, mean, sd, CV),
          file.path(tab_dir, "SuppTable1b_FAA_abundance_summary.csv"), row.names = FALSE)

cat(sprintf("\nfold range, most to least abundant : %.0fx\n", max(summ$median)/min(summ$median)))
cat(sprintf("most variable FAA                  : %s (CV %.0f%%)\n",
            summ$label[which.max(summ$CV)], max(summ$CV)))

################################################################################
### FIGURE
################################################################################

plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title    = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 8)),
    plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 10)),
    plot.tag      = element_text(size = 22, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    axis.title.x  = element_text(size = 18, face = "bold"),
    axis.title.y  = element_text(size = 18, face = "bold"),
    axis.text.x   = element_text(size = 15, face = "bold", color = "black"),
    axis.text.y   = element_text(size = 15, face = "bold", color = "black"),
    axis.line     = element_line(color = "black"),
    axis.ticks    = element_line(color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.3),
    legend.position = "none",
    plot.margin = margin(15, 15, 15, 15)
  )

rank1_colors <- c("Asparagine" = "#2166AC", "Proline" = "#B35806",
                  "Other amino acids" = "#9E9E9E")
pair_colors  <- c("Asn + Pro top two" = "#6A3D9A", "Other top-two pairings" = "#BDBDBD")
## Panel C reuses the panel A colours exactly, so Asn is blue and Pro orange in
## every panel. Highlighting both in one colour would make Asn blue in A and
## orange in C, which reads as two different things being marked.
aa_colors    <- c("Asn" = "#2166AC", "Pro" = "#B35806", "other" = "#9E9E9E")

rank1_df <- tibble(
  category = factor(c("Asparagine", "Proline", "Other amino acids"),
                    levels = c("Asparagine", "Proline", "Other amino acids")),
  n_count  = c(n_asn, n_pro, n_oth)) %>% mutate(prop = n_count / n_genotypes)

top2_df <- tibble(
  category = factor(c("Asn + Pro top two", "Other top-two pairings"),
                    levels = c("Asn + Pro top two", "Other top-two pairings")),
  n_count  = c(n_both, n_genotypes - n_both)) %>% mutate(prop = n_count / n_genotypes)

bar <- function(df, cols, tag, title, subtitle, width) {
  ggplot(df, aes(x = category, y = n_count, fill = category)) +
    geom_col(width = width, alpha = 0.92) +
    geom_text(aes(label = paste0(n_count, " (", percent(prop, accuracy = 1), ")")),
              vjust = -0.6, size = 5) +
    scale_fill_manual(values = cols) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.10)), breaks = pretty_breaks(n = 5)) +
    labs(tag = tag, title = title, subtitle = subtitle, x = NULL, y = "Number of Genotypes") +
    plot_theme
}

p1 <- bar(rank1_df, rank1_colors, "A", "Most Abundant FAA by Genotype",
          "Counts of genotypes in which each amino acid ranked first", 0.68)
p2 <- bar(top2_df, pair_colors, "B", "Asparagine and Proline as the Top-Two FAA",
          "Counts of genotypes where Asn and Pro were the top two,\nregardless of order", 0.62)

## panel C, left: concentration. Log scale because medians span ~144-fold, so a
## linear axis would collapse the bottom two thirds of the amino acids onto zero.
c_left <- ggplot(summ, aes(y = aa_f, colour = hilite)) +
  geom_linerange(aes(xmin = p05, xmax = p95), linewidth = .7, alpha = .45) +
  geom_linerange(aes(xmin = q25, xmax = q75), linewidth = 2.2) +
  geom_point(aes(x = median), size = 3.1) +
  scale_colour_manual(values = aa_colors, guide = "none") +
  scale_x_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.4))) +
  labs(tag = "C", title = "Abundance and Variability of Each FAA",
       x = "Concentration (genotype mean, log scale)", y = NULL) +
  plot_theme +
  theme(axis.text.y = element_text(size = 13, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey88", linewidth = 0.3))

c_right <- ggplot(summ, aes(CV, aa_f, colour = hilite)) +
  geom_segment(aes(x = 0, xend = CV, yend = aa_f), linewidth = .5, alpha = .45) +
  geom_point(size = 3.1) +
  scale_colour_manual(values = aa_colors, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, .08))) +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.4))) +
  labs(title = " ", x = "CV across genotypes (%)", y = NULL) +
  plot_theme +
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey88", linewidth = 0.3))

panelC <- c_left + c_right + plot_layout(widths = c(1.3, 1))

fig <- (p1 + p2) / panelC + plot_layout(heights = c(1, 1.5))

ggsave(out_file, fig, width = 15, height = 15.5, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(sub("\\.png$", ".pdf", out_file), fig,
       width = 15, height = 15.5, bg = "white", limitsize = FALSE)
cat("Saved figure:\n  ", out_file, "\n", sep = "")
