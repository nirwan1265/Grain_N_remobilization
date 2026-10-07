################################################################################
### SUPPLEMENTARY FIGURE 1: ASPARAGINE / PROLINE FAA RANK COUNTS
###
### Replaces scripts/22_SuppFig_FAA_Asn_Pro_panel.R, in which the counts were
### typed in from the manuscript ("SOURCE COUNTS FROM MANUSCRIPT SUMMARY") rather
### than derived from the data. Everything below is computed from the raw plot
### records, and the counts that go in the text are printed to the console.
###
###   Rscript SuppFig1_FAA_Asn_Pro_rank.R [path/to/FAA_Goodman_Buckler.csv] [outdir]
################################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(ggplot2)
  library(patchwork); library(scales)
})

args    <- commandArgs(trailingOnly = TRUE)
in_csv  <- if (length(args) >= 1) args[1] else "data/FAA_Goodman_Buckler.csv"
out_dir <- if (length(args) >= 2) args[2] else "final/supp_figs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir, "SuppFig1_FAA_Asn_Pro_rank.png")

AA <- c("A","R","N","D","Q","E","G","H","I","L","K","M","F","P","S","W","T","Y","V","C")
FULL <- c(N = "Asparagine", P = "Proline")

################################################################################
### COMPUTE THE COUNTS FROM THE DATA
################################################################################

raw <- read.csv(in_csv, check.names = FALSE)
names(raw) <- sub("^﻿", "", trimws(names(raw)))
stopifnot(all(c("taxa", AA) %in% names(raw)))

# genotype means across replicates, then rank the 20 FAAs within each genotype
gmean <- raw %>%
  group_by(taxa) %>%
  summarise(across(all_of(AA), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

mat  <- as.matrix(gmean[, AA])
rownames(mat) <- gmean$taxa
# rank within genotype; a trait that is missing for every replicate of a
# genotype becomes NaN and is pushed to the back rather than breaking the rank
rk   <- t(apply(-mat, 1, rank, ties.method = "min", na.last = TRUE))
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
oth <- sort(table(top1[!top1 %in% c("N","P")]), decreasing = TRUE)
if (length(oth)) cat("  the other first-ranked FAAs:",
                     paste(sprintf("%s=%d", names(oth), as.integer(oth)), collapse = ", "), "\n")
reps <- table(table(raw$taxa))
cat("replicates per genotype            :",
    paste(sprintf("%s rep=%s", names(reps), as.integer(reps)), collapse = ", "), "\n\n")

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
    scale_y_continuous(expand = expansion(mult = c(0, 0.08)), breaks = pretty_breaks(n = 5)) +
    labs(tag = tag, title = title, subtitle = subtitle, x = NULL, y = "Number of Genotypes") +
    plot_theme
}

p1 <- bar(rank1_df, rank1_colors, "A", "Most Abundant FAA by Genotype",
          "Counts of genotypes in which each amino acid ranked first", 0.68)
p2 <- bar(top2_df, pair_colors, "B", "Asparagine and Proline as the Top-Two FAA",
          "Counts of genotypes where Asn and Pro were the top two,\nregardless of order", 0.62)

ggsave(out_file, p1 + p2 + plot_layout(widths = c(1, 1)),
       width = 13, height = 6.8, dpi = 300, bg = "white")
cat("Saved figure:\n  ", out_file, "\n", sep = "")
