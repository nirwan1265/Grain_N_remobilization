#!/usr/bin/env Rscript
################################################################################
## SUPPLEMENTARY FIGURE + SUPPLEMENTARY TABLES S11-S12 -- genotype-environment
## correlation test for the soil-N GWAS SNPs (FDR 5%).
##
## Input  geno_env_cor_9vars.csv.gz  (final/scripts/soilN_geno_env_correlation.py):
##        Pearson r between allele dosage and each of 9 collection-site variables
##        for all 495,755 SNPs (3,107 accessions), and each SNP's genome-wide
##        percentile of |r| (0.1 = top 0.1 % of the genome for that variable).
##        SuppTable_S11_soilN_SNPs_env_relations.csv : FDR soil-N SNPs, models,
##        lead genes, and which environmental GWAS each SNP is also significant in.
##
## A  the 15 SNPs significant in >= 2 models: tile = genome-wide percentile of
##    |r| with each variable; circle = SNP also FDR-significant in that
##    variable's GWAS (not drawn for soil N 0-5 cm, where every SNP is). Right: call.
## B  expectation test, simplified (2026-10-27): share of SNPs whose genotype is
##    in the top 1 % of the genome for latitude, longitude or elevation, for
##    (1) all SNPs, (2) the top 1 % of SNPs by genotype-soil N correlation alone
##    (no GWAS model), and (3) the soil-N GWAS SNPs (FDR). The 6-axis version and
##    the >= 2-model group are kept in Supplemental Table S13 only.
##
## Call (per SNP; comparators = 6 non-N variables)
##   linked to environment : top 1 % of the genome on >= 1 comparator
##   inconclusive          : the SNP's own genotype-soil N correlation is weak
##                           (outside the top 5 % of the genome), so there is no
##                           soil-N correlation to compare the other variables
##                           against; the test cannot separate it either way.
##                           Typical of SNPs whose GWAS signal appears only in
##                           BLINK/FarmCPU after adjusting for other loci.
##   not linked            : neither of the above
################################################################################
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
gv   <- function(k, d) { v <- Sys.getenv(k, ""); if (nzchar(v)) v else d }
REPO <- gv("GN_REPO", "/Users/nirwantandukar/Documents/Github/Grain_N_remobilization")
F_COR <- gv("GN_COR", file.path(REPO, "final/results/soil_nitrogen_gwas/geno_env_cor_9vars.csv.gz"))
F_S15 <- gv("GN_S15", file.path(REPO, "final/supp_tables/SuppTable_S11_soilN_SNPs_env_relations.csv"))
FIG   <- gv("GN_FIGDIR", file.path(REPO, "final/supp_figs"))
TAB   <- gv("GN_TABDIR", file.path(REPO, "final/supp_tables"))

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

VARS <- c("N_0_5", "N_5_15", "latitude", "longitude", "elevation", "pH", "CEC", "aridity")
LAB  <- c(N_0_5 = "Soil N\n0-5 cm", N_5_15 = "Soil N\n5-15 cm", latitude = "Latitude",
          longitude = "Longitude", elevation = "Elevation",
          pH = "pH", CEC = "CEC", aridity = "Aridity")
GEO <- c("latitude", "longitude", "elevation")
CMP <- c(GEO, "pH", "CEC", "aridity")

g <- fread(F_COR)
M <- nrow(g)
g[, geo_min := do.call(pmin, .SD), .SDcols = paste0("pct_", GEO)]
g[, cmp_min := do.call(pmin, .SD), .SDcols = paste0("pct_", CMP)]
g[, cmp_axis := CMP[apply(as.matrix(.SD), 1, which.min)], .SDcols = paste0("pct_", CMP)]

s <- fread(F_S15)
x <- merge(s[, .(SNP, Chr, Pos, MAF, n_models, min_P, env_hits, replicates_N_5_15,
                 lead_gene, lead_rel, lead_dist, lead_name)], g, by = "SNP")
x[, call := fifelse(cmp_min < 1, "linked to environment",
            fifelse(pct_N_0_5 >= 5, "inconclusive", "not linked"))]
cat(sprintf("%d FDR soil-N SNPs matched to correlations (of %d)\n", nrow(x), nrow(s)))
print(x[n_models >= 2, .N, by = call])

## ---- A: heatmap of the >= 2-model SNPs --------------------------------------
h <- x[n_models >= 2]
short <- function(z) { z <- sub(" domain-containing protein", "", z); z <- sub("Sec-independent protein translocase protein ", "", z)
  z <- sub(" chloroplastic", "", z, ignore.case = TRUE); ifelse(nchar(z) > 28, paste0(substr(z, 1, 26), "..."), z) }
h[, rowlab := sprintf("%s  %s", SNP, ifelse(lead_gene == "" | is.na(lead_gene), "no gene <=25 kb",
                                             paste0(lead_gene, " (", short(lead_name), ")")))]
h[lead_name == "" & lead_gene != "", rowlab := sprintf("%s  %s", SNP, lead_gene)]
## use the names given in Table 2 of the paper for the genes discussed there
NAME_FIX <- c(Zm00001eb220040 = "aspartate kinase isozyme", Zm00001eb130210 = "TTL1-like",
              Zm00001eb005710 = "omega-amidase")
h[lead_gene %in% names(NAME_FIX), rowlab := sprintf("%s  %s (%s)", SNP, lead_gene, NAME_FIX[lead_gene])]
setorder(h, cmp_min); h[, rowlab := factor(rowlab, levels = rev(rowlab))]

hm <- melt(h, id.vars = c("SNP", "rowlab", "env_hits", "replicates_N_5_15"),
           measure.vars = paste0("pct_", VARS), variable.name = "var", value.name = "pct")
hm[, var := sub("^pct_", "", var)]
hm[, bin := cut(pct, c(-Inf, 0.1, 1, 5, Inf), labels = c("top 0.1%", "top 1%", "top 5%", "> 5%"))]
hm[, hit := mapply(function(v, e, r) {
  if (v == "N_0_5") return(FALSE); if (v == "N_5_15") return(isTRUE(r))
  v %in% trimws(strsplit(e, ",")[[1]]) }, var, ifelse(is.na(env_hits), "", env_hits), replicates_N_5_15)]
hm[, var := factor(var, levels = VARS)]

callcol <- c(`linked to environment` = "#B2182B", `inconclusive` = "grey45",
             `not linked` = "#1B7837")
pA <- ggplot(hm, aes(var, rowlab)) +
  geom_tile(aes(fill = bin), colour = "white", linewidth = 1.2) +
  geom_point(data = hm[hit == TRUE], shape = 21, size = 4.2, fill = "white", colour = "black", stroke = 1) +
  geom_vline(xintercept = c(2.5, 5.5), linewidth = 0.8, colour = "grey20") +
  scale_fill_manual(values = c(`top 0.1%` = "#08306B", `top 1%` = "#2171B5", `top 5%` = "#9ECAE1", `> 5%` = "grey93"),
                    drop = FALSE, name = "Genome-wide rank of |r|") +
  geom_text(data = h, aes(x = length(VARS) + 0.75, y = rowlab, label = call, colour = call),
            hjust = 0, size = 5, fontface = "bold", inherit.aes = FALSE, show.legend = FALSE) +
  scale_colour_manual(values = callcol) +
  scale_x_discrete(labels = LAB, expand = expansion(add = c(0.6, 4.6))) +
  labs(x = NULL, y = NULL,
       caption = "Tile: genome-wide rank of |r(genotype, variable)|.  Circle: SNP also FDR-significant in that variable's GWAS.") +
  coord_cartesian(clip = "off") +
  plot_theme +
  theme(axis.line = element_blank(), axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 13, face = "bold", family = "mono"),
        legend.position = "top", legend.title = element_text(size = 14, face = "bold"),
        plot.caption = element_text(size = 13, colour = "grey30", hjust = 0))

## ---- B: expectation test ----------------------------------------------------
nN  <- g[pct_N_0_5 < 1]                     # SNPs correlated with soil N, no GWAS
grp <- list(`Whole genome` = g, `Top 1% |r| with soil N` = nN,
            `Soil-N GWAS SNPs (FDR)` = x, `Soil-N GWAS, >=2 models` = x[n_models >= 2])
B <- rbindlist(lapply(names(grp), function(k) { d <- grp[[k]]
  rbind(data.table(group = k, axis = "Geographic (3 axes)", n = nrow(d), k = sum(d$geo_min < 1)),
        data.table(group = k, axis = "Geographic + pH, CEC, aridity (6 axes)", n = nrow(d), k = sum(d$cmp_min < 1)))}))
B[, pct := 100 * k / n]
B[, group := factor(group, levels = names(grp))]
## binomial tests: GWAS SNPs vs genome (more than chance?) and vs soil-N-correlated SNPs (less?)
bt <- function(ax, against, alt) { b <- B[group == "Soil-N GWAS SNPs (FDR)" & axis == ax]
  p0 <- B[group == against & axis == ax]$k / B[group == against & axis == ax]$n
  binom.test(b$k, b$n, p0, alternative = alt)$p.value }
B[, P_vs_genome := NA_real_][, P_vs_Ncorrelated := NA_real_]
for (ax in unique(B$axis)) {
  B[group == "Soil-N GWAS SNPs (FDR)" & axis == ax, `:=`(
    P_vs_genome = bt(ax, "Whole genome", "greater"),
    P_vs_Ncorrelated = bt(ax, "Top 1% |r| with soil N", "less"))] }
print(B)

PB <- B[axis == "Geographic (3 axes)" & group != "Soil-N GWAS, >=2 models"]
PB[, group := droplevels(group)]
pB <- ggplot(PB, aes(group, pct, fill = group)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -1.9, size = 7, fontface = "bold") +
  geom_text(aes(label = sprintf("%s / %s", prettyNum(k, big.mark = ","), prettyNum(n, big.mark = ","))),
            vjust = -0.6, size = 4.6, colour = "grey30") +
  scale_fill_manual(values = c("grey70", "#6A51A3", "#B35806"), guide = "none") +
  scale_x_discrete(labels = c("All SNPs", "SNPs picked by\nsoil N correlation\nalone (no GWAS)",
                              "Soil N\nGWAS SNPs")) +
  scale_y_continuous(limits = c(0, 46), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "SNPs whose genotype follows latitude,\nlongitude or elevation (%)") +
  plot_theme +
  theme(axis.text.x = element_text(size = 14, face = "bold", lineheight = 0.95),
        axis.title.y = element_text(size = 16, face = "bold"))

fig <- (pA | pB) + plot_layout(widths = c(2.4, 1)) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 28, face = "bold"))
for (ext in c("png", "pdf"))
  ggsave(file.path(FIG, paste0("SuppFig_soilN_geno_env_correlation.", ext)), fig,
         width = 28, height = 12, dpi = 300, bg = "white",
         device = if (ext == "pdf") cairo_pdf else "png")

## ---- Supplementary Tables S12 and S12 -------------------------------------------------
setorder(x, -n_models, min_P)
s17 <- x[, c(list(SNP = SNP, Chr = Chr, Pos = Pos, MAF = round(MAF, 4), n_models = n_models,
                  min_P = signif(min_P, 3), lead_gene = lead_gene, lead_name = lead_name),
             setNames(lapply(VARS, function(v) round(get(paste0("r_", v)), 3)), paste0("r_", VARS)),
             setNames(lapply(VARS, function(v) signif(get(paste0("pct_", v)), 3)), paste0("pct_", VARS)),
             list(best_geographic_pct = signif(geo_min, 3), best_comparator_pct = signif(cmp_min, 3),
                  best_comparator = cmp_axis, env_GWAS_hits = env_hits,
                  also_N_5_15_GWAS = replicates_N_5_15,
                  call = fifelse(call == "inconclusive", "inconclusive (weak soil-N correlation)", call)))]
fwrite(s17, file.path(TAB, "SuppTable_S12_soilN_geno_env_correlation.csv"))
fwrite(B, file.path(TAB, "SuppTable_S13_soilN_geno_env_expectation.csv"))
cat("written figure + S11, S12\n")
