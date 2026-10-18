#!/usr/bin/env Rscript
################################################################################
## SUPP FIGURE -- geography and soil N of the two focal soil-N SNPs
##
## Rows: S1_16128632 (omega-amidase, Zm00001eb005710) and S5_22338979 (aspartate kinase
## isozyme, Zm00001eb220040). 3,107 landraces with soil N.
## Left: collection sites; grey = homozygous major allele, white = heterozygous, black = homozygous minor.
## Right: soil N at 0-5 cm by genotype (boxplot; points shown for the minor genotypes), with a
## Kruskal-Wallis test across genotypes and Dunn post hoc tests (Holm correction).
## Input : final/results/soil_nitrogen_gwas/focal_SNP_genotypes_geography.csv (SNP, Taxa, alt_dosage,
##         N_0_5, lat, lon, region, country), built from all_maize2.vcf.gz, control_phenotypes_maize_landrace.csv
##         and crop_coordinates.csv
## Output: final/supp_figs/SuppFig_soilN_focal_allele_geography.{png,pdf}
################################################################################
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork); library(maps); library(rstatix) })
ROOT <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "."
IN  <- file.path(ROOT, "final", "results", "soil_nitrogen_gwas", "focal_SNP_genotypes_geography.csv")
OUT <- file.path(ROOT, "final", "supp_figs")

plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.tag     = element_text(size = 24, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x  = element_text(size = 15, colour = "black"),
    axis.text.y  = element_text(size = 15, colour = "black"),
    axis.line    = element_line(colour = "black"),
    panel.grid   = element_blank(),
    legend.position = "none",
    plot.margin  = margin(15, 15, 15, 15)
  )

d <- fread(IN)[!is.na(alt_dosage)]
info <- data.table(SNP = c("S1_16128632", "S5_22338979"), ref = c("G", "G"), alt = c("C", "T"),
                   label = c("S1_16128632, omega-amidase (Zm00001eb005710)",
                             "S5_22338979, aspartate kinase isozyme (Zm00001eb220040)"))
world <- as.data.table(map_data("world"))

panel <- function(snp) {
  x <- d[SNP == snp]; i <- info[SNP == snp]
  x[, carrier := alt_dosage > 0]
  ncar <- sum(x$carrier); nsite <- uniqueN(x[carrier == TRUE, .(round(lat, 3), round(lon, 3))])
  ncty <- uniqueN(x[carrier == TRUE, country])
  gt <- c(paste0(i$ref, i$ref), paste0(i$ref, i$alt), paste0(i$alt, i$alt))
  x[, geno := factor(gt[alt_dosage + 1], levels = gt)]
  fills <- setNames(c("grey70", "white", "black"), gt)
  setorder(x, alt_dosage)
  m <- ggplot() +
    geom_polygon(data = world, aes(long, lat, group = group), fill = "grey95", colour = "grey60", linewidth = 0.25) +
    geom_point(data = x[alt_dosage == 0], aes(lon, lat), colour = "grey70", size = 0.9) +
    geom_point(data = x[alt_dosage > 0], aes(lon, lat, fill = geno), colour = "black", shape = 21, size = 3.4, stroke = 0.8) +
    scale_fill_manual(values = fills[2:3], name = NULL) +
    annotate("label", x = -116, y = -49, hjust = 0, vjust = 0, size = 5, fill = "white", label.size = 0, lineheight = 0.95,
             label = sprintf("%s\n%d carriers, %d sites, %d countries", i$label, ncar, nsite, ncty)) +
    coord_quickmap(xlim = c(-118, -34), ylim = c(-52, 33), expand = FALSE) +
    labs(x = "Longitude", y = "Latitude") + plot_theme +
    theme(panel.border = element_rect(fill = NA, colour = "black"), axis.line = element_blank(),
          legend.position = c(0.86, 0.88), legend.background = element_rect(fill = "white", colour = "grey60"),
          legend.text = element_text(size = 15)) +
    guides(fill = guide_legend(override.aes = list(size = 5)))
  cnt <- x[, .N, by = geno]
  lev <- sprintf("%s\n(n = %s)", cnt[match(gt, geno), geno], trimws(format(cnt[match(gt, geno), N], big.mark = ",")))
  x[, geno_lab := factor(lev[as.integer(geno)], levels = lev)]
  ## Kruskal-Wallis across the three genotypes, Dunn post hoc tests with Holm correction
  df <- as.data.frame(x[, .(N_0_5, geno_lab)])
  kw <- kruskal_test(df, N_0_5 ~ geno_lab); es <- kruskal_effsize(df, N_0_5 ~ geno_lab)
  dn <- as.data.table(dunn_test(df, N_0_5 ~ geno_lab, p.adjust.method = "holm"))
  fmtp <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 1), formatC(p, format = "f", digits = 3))
  ymax <- max(x$N_0_5)
  dn[, `:=`(x1 = match(group1, lev), x2 = match(group2, lev))]
  dn <- dn[order(x2 - x1)]
  dn[, lab := paste0("P = ", fmtp(p.adj))]
  ds <- dn[p.adj < 0.05]                       # brackets for significant comparisons only
  ds[, y := ymax * (1.06 + 0.075 * (.I - 1))]
  kw_lab <- sprintf("Kruskal-Wallis chi^2 = %.1f, df = %d, P = %s, epsilon^2 = %.3f", kw$statistic, kw$df, fmtp(kw$p), es$effsize)
  stats_rows <<- rbind(stats_rows, data.table(SNP = snp, test = "Kruskal-Wallis", comparison = "all genotypes",
                       statistic = kw$statistic, P = kw$p, effect = es$effsize),
                       dn[, .(SNP = snp, test = "Dunn (Holm)", comparison = paste(sub("\n.*", "", group1), "vs", sub("\n.*", "", group2)),
                              statistic = statistic, P = p.adj, effect = NA_real_)])
  b <- ggplot(x, aes(geno_lab, N_0_5)) +
    geom_boxplot(fill = "white", colour = "black", outlier.shape = NA, width = 0.55, linewidth = 0.8) +
    geom_jitter(data = x[alt_dosage > 0], aes(fill = geno), width = 0.12, height = 0, size = 2.6, shape = 21, colour = "black", stroke = 0.6) +
    scale_fill_manual(values = fills) +
    geom_segment(data = ds, aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.6) +
    geom_segment(data = ds, aes(x = x1, xend = x1, y = y, yend = y - ymax * 0.015), inherit.aes = FALSE, linewidth = 0.6) +
    geom_segment(data = ds, aes(x = x2, xend = x2, y = y, yend = y - ymax * 0.015), inherit.aes = FALSE, linewidth = 0.6) +
    geom_text(data = ds, aes(x = (x1 + x2) / 2, y = y, label = lab), inherit.aes = FALSE, vjust = -0.4, size = 4.8) +
    annotate("text", x = 0.55, y = ymax * (1.06 + 0.075 * nrow(ds)) + ymax * 0.03, hjust = 0, size = 4.8,
             label = sprintf("Kruskal*'-'*Wallis~chi^2==%.1f*','~italic(P)=='%s'*','~epsilon^2=='%.3f'", kw$statistic, fmtp(kw$p), es$effsize), parse = TRUE) +
    scale_y_continuous(limits = c(NA, ymax * (1.06 + 0.075 * nrow(ds)) + ymax * 0.08)) +
    labs(x = "Genotype", y = "Soil N, 0-5 cm") + plot_theme
  list(m, b)
}
stats_rows <- data.table()
p1 <- panel("S1_16128632"); p2 <- panel("S5_22338979")
fwrite(stats_rows, file.path(dirname(IN), "focal_SNP_genotype_soilN_tests.csv"))
print(stats_rows)
fig <- (p1[[1]] | p1[[2]]) / (p2[[1]] | p2[[2]]) + plot_layout(widths = c(1.25, 1)) + plot_annotation(tag_levels = "A")
for (ext in c("png", "pdf"))
  ggsave(file.path(OUT, paste0("SuppFig_soilN_focal_allele_geography.", ext)), fig, width = 16, height = 17,
         dpi = 300, bg = "white", device = if (ext == "png") ragg::agg_png else "pdf")
print(d[, .(n = .N, carriers = sum(alt_dosage > 0), med_N_carrier = as.numeric(median(N_0_5[alt_dosage > 0])),
            med_N_non = as.numeric(median(N_0_5[alt_dosage == 0]))), by = SNP])
