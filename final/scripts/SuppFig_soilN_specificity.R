#!/usr/bin/env Rscript
################################################################################
## SUPPLEMENTARY FIGURE + SUPPLEMENTARY TABLE S18 -- environmental specificity
## test for the soil-N GWAS signal, built from the control GWAS themselves.
##
## Each trait (soil N 0-5 / 5-15 cm, aridity, pH, CEC, elevation, latitude,
## longitude) was scanned with MLM, BLINK and FarmCPU (3 PCs,
## 495,755 SNPs); a SNP is significant for a trait if BH-FDR <= 0.05 in any
## model.
##
## A  share of each trait's significant SNPs that are also significant for a
##    geographic trait (union of latitude, longitude, elevation,
##    trait itself left out). Tick = share expected by chance.
## B  soil N (0-5 cm) against every other trait: observed shared SNPs vs the
##    hypergeometric expectation  E = n_soilN * n_trait / M.
##
## Test: one-sided hypergeometric, P(X >= observed), M = 495,755.
## It assumes the two SNP sets are drawn independently; for correlated traits
## (0-5 vs 5-15 cm, r = 0.84) enrichment is expected and reflects shared
## phenotype, not independent replication.
##
## Writes  final/supp_figs/SuppFig_soilN_specificity.{png,pdf}
##         final/supp_tables/SuppTable_S18_soilN_env_overlap_test.csv
################################################################################
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
gv   <- function(k, d) { v <- Sys.getenv(k, ""); if (nzchar(v)) v else d }
REPO <- gv("GN_REPO", "/Users/nirwantandukar/Documents/Github/Grain_N_remobilization")
D    <- gv("GN_ENV_DIR", file.path(REPO, "final", "results", "soil_nitrogen_gwas", "env_landrace"))
FIG  <- gv("GN_FIGDIR", file.path(REPO, "final", "supp_figs"))
TAB  <- gv("GN_TABDIR", file.path(REPO, "final", "supp_tables"))
dir.create(FIG, FALSE, TRUE); dir.create(TAB, FALSE, TRUE)

M <- 495755; MODELS <- c("MLM", "BLINK", "FarmCPU")
GEO <- c("latitude", "longitude", "elevation")
LAB <- c(N_0_5 = "Soil N (0-5 cm)", N_5_15 = "Soil N (5-15 cm)", aridity = "Aridity",
         pH = "Soil pH", CEC = "CEC", elevation = "Elevation", latitude = "Latitude",
         longitude = "Longitude")

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

.gg35 <- utils::packageVersion("ggplot2") >= "3.5.0"
legend_inside <- function(x, y, just) {
  if (.gg35) theme(legend.position = "inside", legend.position.inside = c(x, y),
                   legend.justification = just)
  else       theme(legend.position = c(x, y), legend.justification = just)
}

rdx <- function(tr, md) fread(file.path(D, sprintf("%s__%s.csv.gz", tr, md)),
                              showProgress = FALSE)[fdr_sig == TRUE, .(trait = tr, SNP)]
traits <- intersect(names(LAB), unique(sub("__.*", "", list.files(D, pattern = "\\.csv\\.gz$"))))
sig  <- unique(rbindlist(lapply(traits, function(t) rbindlist(lapply(MODELS, function(m) rdx(t, m))))))
sets <- split(sig$SNP, sig$trait)
hyp  <- function(o, a, b) phyper(o - 1, a, M - a, b, lower.tail = FALSE)

## ---- A: overlap with geographic union (leave-self-out) ----------------------
A <- rbindlist(lapply(traits, function(t) {
  s <- sets[[t]]; g <- unique(unlist(sets[setdiff(GEO, t)])); o <- length(intersect(s, g))
  e <- length(s) * length(g) / M
  data.table(trait = t, n = length(s), n_ref = length(g), ov = o, exp = e,
             pct = 100 * o / length(s), exp_pct = 100 * length(g) / M,
             P = hyp(o, length(g), length(s)))}))
A[, lab := LAB[trait]][, grp := fifelse(trait == "N_0_5", "focal",
                              fifelse(trait %in% GEO, "geography", "other soil"))]
setorder(A, pct); A[, lab := factor(lab, levels = lab)]

pA <- ggplot(A, aes(pct, lab, fill = grp)) +
  geom_col(width = .7) +
  geom_point(aes(x = exp_pct), shape = 124, size = 9, colour = "grey20", show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%  (%d/%d)", pct, ov, n)), hjust = -0.08,
            size = 5.2, fontface = "bold") +
  scale_fill_manual(values = c(focal = "#B35806", geography = "#2166AC", `other soil` = "#9E9E9E"),
                    breaks = c("focal", "geography", "other soil"),
                    labels = c("Soil N (focal)", "Geographic trait", "Other soil trait")) +
  scale_x_continuous(expand = expansion(mult = c(0, .35))) +
  labs(x = "SNPs also significant for a geographic trait (%)", y = NULL,
       caption = "| = overlap expected by chance") +
  plot_theme + legend_inside(.98, .04, c(1, 0)) +
  theme(plot.caption = element_text(size = 13, colour = "grey30"))

## ---- B: soil N vs each trait, observed vs expected --------------------------
N <- sets[["N_0_5"]]
B <- rbindlist(lapply(setdiff(traits, "N_0_5"), function(t) {
  o <- length(intersect(N, sets[[t]])); e <- length(N) * length(sets[[t]]) / M
  data.table(trait = t, n = length(sets[[t]]), ov = o, exp = e, P = hyp(o, length(sets[[t]]), length(N)))}))
B[, lab := LAB[trait]]; setorder(B, ov); B[, lab := factor(lab, levels = lab)]
B[, star := fifelse(P < 1e-20, "***", fifelse(P < 1e-5, "**", fifelse(P < 0.05, "*", "n.s.")))]

pB <- ggplot(B, aes(ov, lab)) +
  geom_segment(aes(x = exp, xend = ov, yend = lab), colour = "grey70", linewidth = 1.4) +
  geom_point(aes(x = exp, colour = "Expected"), size = 4.5) +
  geom_point(aes(colour = "Observed"), size = 6) +
  geom_text(aes(label = sprintf("%d  %s", ov, star)), hjust = -0.3, size = 5.2, fontface = "bold") +
  scale_colour_manual(values = c(Expected = "grey45", Observed = "#B35806")) +
  scale_x_continuous(expand = expansion(mult = c(.04, .28))) +
  labs(x = "Soil N (0-5 cm) SNPs shared with trait", y = NULL) +
  plot_theme + legend_inside(.98, .04, c(1, 0))

fig <- (pA | pB) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 28, face = "bold"))
for (ext in c("png", "pdf"))
  ggsave(file.path(FIG, paste0("SuppFig_soilN_specificity.", ext)), fig,
         width = 22, height = 8.5, dpi = 300, bg = "white")

## ---- Supplementary Table S18 ------------------------------------------------
fmtP <- function(p) ifelse(p < 1e-300, "<1e-300", formatC(p, format = "e", digits = 2))
s18a <- A[order(-as.integer(lab)), .(
  Test = "Trait vs geographic union", Trait = as.character(lab), Reference = "Geographic union (self excluded)",
  `n trait SNPs` = n, `n reference SNPs` = n_ref, `Observed shared` = ov,
  `Expected shared` = round(exp, 3), `Fold enrichment` = round(ov / exp, 1),
  `% of trait SNPs shared` = round(pct, 1), `Hypergeometric P` = fmtP(P))]
s18b <- B[order(-ov), .(
  Test = "Soil N (0-5 cm) vs trait", Trait = "Soil N (0-5 cm)", Reference = as.character(lab),
  `n trait SNPs` = length(N), `n reference SNPs` = n, `Observed shared` = ov,
  `Expected shared` = round(exp, 3), `Fold enrichment` = round(ov / exp, 1),
  `% of trait SNPs shared` = round(100 * ov / length(N), 1), `Hypergeometric P` = fmtP(P))]
s18 <- rbind(s18a, s18b)
fwrite(s18, file.path(TAB, "SuppTable_S18_soilN_env_overlap_test.csv"))
print(s18)
