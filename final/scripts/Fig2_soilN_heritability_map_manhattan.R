################################################################################
## Figure 2 — Soil nitrogen in maize landraces
##   A  soil N at 0-5 cm against 5-15 cm across the 3,107 accessions
##   B  geographic distribution of soil N across the sampled landraces
##   C  genomic h2 of both soil-N layers after removing 0 / 10 / 50 genotype
##      PCs, with the 95% profile CI at 0 PCs and the permutation null stated
##   D  cross-model Manhattan (MLM, BLINK, FarmCPU), full width, every SNP
##
## Three models only. MLMM is excluded: its scan is missing all 61,326
## chromosome-2 markers, so chromosome 2 cannot enter a four-model consensus.
## LFMM is not used here.
##
## THRESHOLD: Benjamini-Hochberg FDR 5%, per model, from the scan used in the
## Results text (all three models genome-wide, 3 PCs, matched 3,107 taxa).
## Panel D previously plotted the OLD run (genome-wide BLINK + per-chromosome
## FarmCPU) with Bonferroni/suggestive lines; that run is no longer reported.
## Panel C draws EVERY SNP from the complete scans (raw_GWAS_<MODEL>_3PC_N_0_5.txt,
## 495,755 markers x 3 models). No subsampling and no P cutoff. FDR is computed
## here by Benjamini-Hochberg over each complete scan, exactly as in
## collect_env_gwas.R, so the significant set matches the Results text.
##
## Run:  Rscript Fig2_soilN_heritability_map_manhattan.R
################################################################################

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(maps); library(viridis); library(scales)
})

## ---------------------------------------------------------------- paths ----
## Defaults are the paths on the analysis machine. Every one can be overridden
## with an environment variable, which is how the figure is regenerated
## elsewhere without editing this file.
gv <- function(k, d) { v <- Sys.getenv(k); if (nzchar(v)) v else d }

REPO <- gv("GN_REPO", "/Users/nirwantandukar/Documents/Github/Grain_N_remobilization")
SN   <- file.path(REPO, "final", "results", "soil_nitrogen_gwas")

## every input lives inside final/ -- nothing is read from outside the repo
DIR_PHENO <- gv("GN_DIR_PHENO", SN)
F_H2      <- gv("GN_F_H2",      file.path(SN, "h2_genomic_landrace_int.csv"))
DIR_FULL  <- gv("GN_DIR_FULL",  file.path(SN, "N_0_5_full"))   # complete scans, all 495,755 SNPs
F_GFF     <- gv("GN_GFF",       file.path(REPO, "final", "results", "annotation",
                                          "Zm-B73-REFERENCE-NAM-5.0_genes_only.gff3"))
OUTDIR    <- gv("GN_OUTDIR",    file.path(REPO, "final", "main_figs"))

## same phenotype file the GWAS ran on: 3,107 accessions with soil N
P_PHENO <- file.path(DIR_PHENO, "control_phenotypes_maize_landrace.csv")
P_H2    <- F_H2
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

MODELS   <- c("MLM", "BLINK", "FarmCPU")
MODEL_COL <- c(MLM = "#0072B2", BLINK = "#009E73", FarmCPU = "#D55E00")
LAYER_COL <- c(`0-5 cm` = "#1B7837", `5-15 cm` = "#762A83")


## ggplot2 >= 3.5 moved inside-panel legends to `legend.position.inside`.
## This keeps the script runnable under either version.
.gg35 <- utils::packageVersion("ggplot2") >= "3.5.0"
legend_inside <- function(x, y, just) {
  if (.gg35) theme(legend.position = "inside", legend.position.inside = c(x, y),
                   legend.justification = just)
  else       theme(legend.position = c(x, y), legend.justification = just)
}

## ---- house theme (as in final/scripts/SuppFig2_heritability.R and
## scripts/10_Fig4_Fst_Tajima_Pi.R) -- every panel uses it
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
base_thm <- plot_theme

################################################################################
## PANELS A and B — heritability
################################################################################
h2 <- fread(P_H2)
h2[, Layer := factor(fifelse(grepl("0-5", trait), "0-5 cm", "5-15 cm"),
                     levels = c("0-5 cm", "5-15 cm"))]

## ---- A: h2 with 95% CI vs the permutation null -----------------------------
perm_hi <- max(h2$perm_null_max)
pA <- ggplot(h2, aes(Layer, h2, colour = Layer)) +
  annotate("rect", xmin = 0.4, xmax = 2.6, ymin = 0, ymax = perm_hi,
           fill = "grey75", alpha = 0.55) +
  annotate("text", x = 2.58, y = perm_hi + 0.035, hjust = 1, size = 5,
           fontface = "italic", colour = "grey30",
           label = "permuted phenotypes (1,000x)") +
  geom_errorbar(aes(ymin = h2_lo, ymax = h2_hi), width = 0.10, linewidth = 0.6) +
  geom_point(size = 3) +
  scale_colour_manual(values = LAYER_COL, guide = "none") +
  scale_x_discrete(expand = c(0, 0.4)) +
  scale_y_continuous(limits = c(-0.02, 1.06), breaks = seq(0, 1, 0.25)) +
  labs(x = "Soil depth", y = expression(bold(paste("Genomic ", bolditalic(h)^2)))) +
  base_thm

## ---- B: h2 after removing genotype PCs from the phenotype ------------------
pB_dat <- melt(h2[, .(Layer, `0` = h2, `10` = h2_PC10adj, `50` = h2_PC50adj)],
               id.vars = "Layer", variable.name = "nPC", value.name = "h2")
pB_dat[, nPC := as.integer(as.character(nPC))]

perm_max <- max(h2$perm_null_max)
ci0 <- h2[, .(Layer, nPC = 0L, lo = h2_lo, hi = h2_hi)]
pB <- ggplot(pB_dat, aes(nPC, h2, colour = Layer, group = Layer)) +
  geom_hline(yintercept = 1, linewidth = 0.3, linetype = "dotted", colour = "grey60") +
  geom_errorbar(data = ci0, aes(x = nPC, ymin = lo, ymax = hi, colour = Layer),
                inherit.aes = FALSE, width = 1.6, linewidth = 1.0,
                position = position_dodge(width = 1.2)) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 4.5) +
  annotate("text", x = 50, y = 0.805, hjust = 1, vjust = 0, size = 5.2, colour = "grey30",
           parse = TRUE,
           label = sprintf('"permuted phenotypes (1,000x): "~italic(h)^2 <= %.2f', perm_max)) +
  scale_colour_manual(values = LAYER_COL, name = "Soil depth") +
  scale_x_continuous(breaks = c(0, 10, 50)) +
  scale_y_continuous(limits = c(0.80, 1.02), breaks = seq(0.80, 1.00, 0.05)) +
  labs(x = "Genotype PCs removed from phenotype",
       y = expression(bold(paste("Genomic ", bolditalic(h)^2)))) +
  base_thm +
  legend_inside(0.03, 0.30, c(0, 0)) +
  theme(legend.background = element_rect(fill = alpha("white", 0.85), colour = NA))

################################################################################
## PANEL A -- soil N at the two depths
################################################################################
dep <- fread(P_PHENO)[is.finite(N_0_5) & is.finite(N_5_15), .(N_0_5, N_5_15)]
rr  <- cor.test(dep$N_0_5, dep$N_5_15)
cat(sprintf("panel A: n = %d, Pearson r = %.3f (95%% CI %.3f-%.3f)\n",
            nrow(dep), rr$estimate, rr$conf.int[1], rr$conf.int[2]))
lim <- range(c(dep$N_0_5, dep$N_5_15))
pDep <- ggplot(dep, aes(N_0_5, N_5_15)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey55", linewidth = 0.4) +
  geom_point(size = 1.8, alpha = 0.35, shape = 16, colour = LAYER_COL[["0-5 cm"]]) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = 1.1) +
  annotate("text", x = lim[1], y = lim[2], hjust = 0, vjust = 1, size = 6.5,
           label = sprintf("italic(r) == %.2f", rr$estimate), parse = TRUE) +
  annotate("text", x = lim[1], y = lim[2] - 0.08 * diff(lim), hjust = 0, vjust = 1, size = 5.5,
           label = sprintf("n = %s", format(nrow(dep), big.mark = ","))) +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = "Soil N, 0-5 cm", y = "Soil N, 5-15 cm") +
  base_thm

################################################################################
## PANEL C — geographic distribution of soil N
################################################################################
geo <- fread(P_PHENO)[, .(Genotypes = Taxa, long = longitude, lat = latitude, soil_n = N_0_5)]
geo <- geo[is.finite(soil_n) & is.finite(long) & is.finite(lat)]
cat(sprintf("panel B: %d accessions, %d distinct soil-N values\n", nrow(geo), uniqueN(geo$soil_n)))

wmap <- as.data.table(map_data("world"))
wmap <- wmap[long >= -125 & long <= -28 & lat >= -56 & lat <= 42]

pC_map <- ggplot() +
  geom_polygon(data = wmap, aes(long, lat, group = group),
               fill = "grey95", colour = "grey70", linewidth = 0.18) +
  geom_point(data = geo, aes(long, lat, colour = soil_n), size = 1.3, alpha = 0.8) +
  coord_fixed(1.3, xlim = c(-118, -32), ylim = c(-48, 38)) +
  scale_colour_viridis_c(option = "D", name = "Soil N") +
  labs(x = "Longitude", y = "Latitude") +
  base_thm +
  legend_inside(0.97, 0.97, c(1, 1)) +
  theme(legend.title = element_text(size = 14, face = "bold"),
        legend.key.height = unit(1.6, "cm"), legend.key.width = unit(0.45, "cm"),
        legend.background = element_rect(fill = alpha("white", 0.8), colour = NA))

pC_hist <- ggplot(geo, aes(soil_n, fill = after_stat(x))) +
  geom_histogram(bins = 25, colour = "black", linewidth = 0.25) +
  scale_fill_viridis_c(option = "D", guide = "none") +
  labs(x = "Soil N", y = "Count") +
  theme_classic(base_size = 13) +
  theme(axis.text = element_text(colour = "black", size = 13),
        axis.title = element_text(size = 14, face = "bold"),
        plot.background = element_rect(fill = "white", colour = "grey60", linewidth = 0.3),
        plot.margin = margin(3, 4, 3, 3))

pC <- pC_map + inset_element(pC_hist, left = 0.00, bottom = 0.02,
                             right = 0.42, top = 0.40)

################################################################################
## PANEL D — cross-model Manhattan (3 models)
################################################################################
## raw files: header names 8 columns but FarmCPU rows leave MAF blank, so read by
## position after skipping the header (see collect_env_gwas.R for the details)
read_full <- function(m) {
  f <- file.path(DIR_FULL, sprintf("raw_GWAS_%s_3PC_N_0_5.txt", m))
  if (!file.exists(f)) stop("missing complete scan: ", f,
    "\n  copy raw_GWAS_{MLM,BLINK,FarmCPU}_3PC_N_0_5.txt from the HPC into ", DIR_FULL)
  d <- fread(f, skip = 1, header = FALSE, showProgress = FALSE)
  d <- d[, .(SNP = V1, Chr = as.integer(V2), Pos = as.numeric(V3), P = as.numeric(V4))]
  d <- d[is.finite(P) & P > 0 & !is.na(Chr) & !is.na(Pos)]
  mm <- nrow(d); setorder(d, P)
  d[, FDR := pmin(1, rev(cummin(rev(P * mm / seq_len(mm)))))]
  d[, `:=`(fdr_sig = FDR <= 0.05, Model = m)]
  cat(sprintf("  %-7s %d SNPs, %d FDR-significant\n", m, mm, d[, sum(fdr_sig)]))
  d[]
}
scan_all <- rbindlist(lapply(MODELS, read_full))

## per-model FDR cutoff = largest P still significant in that model
cut <- scan_all[fdr_sig == TRUE, .(cutoff = max(P)), by = Model]
cut[, logc := -log10(cutoff)]
print(cut)
## ONE line drawn: 5% FDR at the most lenient per-model cutoff, so every
## FDR-significant SNP (any model) sits on or above it
fdr_line <- min(cut$logc)
cat(sprintf("single FDR 5%% line at -log10(p) = %.2f\n", fdr_line))

## models in which each SNP is FDR-significant
nmod <- scan_all[fdr_sig == TRUE, .(n_models = uniqueN(Model)), by = SNP]

## one point per SNP; for significant SNPs, the best model among those calling it
sigrows <- scan_all[fdr_sig == TRUE]
best_sig <- sigrows[sigrows[, .I[which.min(P)], by = SNP]$V1]
best_all <- scan_all[scan_all[, .I[which.min(P)], by = SNP]$V1]
best <- rbind(best_sig, best_all[!SNP %in% best_sig$SNP])
best <- merge(best, nmod, by = "SNP", all.x = TRUE)
best[is.na(n_models), n_models := 0L]
best[, logp := -log10(P)]
cat(sprintf("panel D: %d SNPs drawn, %d FDR-significant, %d in >=2 models\n",
            nrow(best), best[n_models >= 1, .N], best[n_models >= 2, .N]))

## chromosome axis from the v5 gene annotation, not from the retained SNPs
gff <- fread(cmd = sprintf("awk -F'\\t' '$3==\"gene\"{print $1\"\\t\"$5}' %s", shQuote(F_GFF)),
             header = FALSE, col.names = c("chr", "end"))
gff[, Chr := suppressWarnings(as.integer(sub("^chr", "", chr)))]
chr_tbl <- gff[Chr %in% 1:10, .(len = max(end)), by = Chr][order(Chr)]
chr_tbl[, offset := cumsum(c(0, head(len, -1)))]
chr_tbl[, centre := offset + len / 2]
best <- merge(best, chr_tbl[, .(Chr, offset)], by = "Chr")
best[, x := Pos + offset]

bg  <- best[n_models == 0]          # every non-significant SNP, none dropped
sig <- best[n_models >= 1]
bg[, grp := fifelse(Chr %% 2 == 0, "even", "odd")]

## loci named in the Results
LABELS <- c(S5_22338979  = "Zm00001eb220040",
            S1_16128632  = "Zm00001eb005710",
            S8_138640842 = "Zm00001eb356170",
            S9_151742636 = "Zm00001eb399680")   # strongest; shared with pH, aridity
lab <- sig[SNP %in% names(LABELS)]
lab[, gene := LABELS[SNP]]
lab[, hj := c(S1_16128632 = 0, S5_22338979 = 0, S8_138640842 = 0.5, S9_151742636 = 1)[SNP]]

ymax <- max(best$logp) * 1.08
pD <- ggplot() +
  geom_point(data = bg, aes(x, logp, colour = grp),
             size = 0.9, alpha = 0.6, shape = 16, show.legend = FALSE) +
  geom_hline(yintercept = fdr_line, linewidth = 1.0, linetype = "dashed",
             colour = "red3") +
  annotate("text", x = min(best$x), y = fdr_line, label = "FDR 5%",
           hjust = 0, vjust = 1.6, size = 5.5, fontface = "bold", colour = "red3") +
  geom_point(data = sig, aes(x, logp, colour = Model), size = 3.6, alpha = 0.9) +
  geom_point(data = sig[n_models >= 2], aes(x, logp),
             shape = 21, size = 6.2, stroke = 1.0, colour = "black", fill = NA,
             show.legend = FALSE) +
  scale_colour_manual(values = c(odd = "grey68", even = "grey45", MODEL_COL),
                      breaks = MODELS, name = "Best model") +
  geom_text(data = lab, aes(x, logp + 0.7, label = gene, hjust = hj),
            size = 5.2, fontface = "italic") +
  scale_x_continuous(breaks = chr_tbl$centre, labels = chr_tbl$Chr,
                     expand = expansion(mult = 0.02)) +
  scale_y_continuous(limits = c(0, ymax), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Chromosome", y = expression(bold(-log[10](bolditalic(p))))) +
  base_thm +
  legend_inside(0.005, 0.995, c(0, 1)) +
  theme(legend.direction = "horizontal",
        legend.background = element_rect(fill = alpha("white", 0.85),
                                         colour = "grey70", linewidth = 0.3)) +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1)))

################################################################################
## ASSEMBLE
################################################################################
top <- (pDep | wrap_elements(full = pC) | pB) + plot_layout(widths = c(1, 1.05, 1.1))
fig <- (top / pD) +
  plot_layout(heights = c(1.15, 1.00)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 32, face = "bold"))

ggsave(file.path(OUTDIR, "Fig2_soilN.png"), fig,
       width = 26, height = 18, dpi = 300, limitsize = FALSE, bg = "white")
ggsave(file.path(OUTDIR, "Fig2_soilN.pdf"), fig,
       width = 26, height = 18, limitsize = FALSE, bg = "white")
cat("written:", file.path(OUTDIR, "Fig2_soilN.png"), "\n")
