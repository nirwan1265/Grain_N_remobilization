#!/usr/bin/env Rscript
################################################################################
## FIGURE 5 -- Cycle 0 vs Cycle 14 differentiation in Indian Chief and Jarvis
##
## One statistic only: F_ST between Cycle 0 and Cycle 14 within each population
## (ANGSD/realSFS site F_ST, summarised in GenWin spline windows; smoothness
## 500 kb). Tajima's D and pi are NOT used: both are computed from the same
## site-frequency spectrum, and both are biased by the deeper sequencing of
## Cycle 14.
##
## Top 5% = window mean F_ST >= the 95th percentile of that population's
## windows (109 of 2,162 Indian Chief windows, 108 of 2,151 Jarvis windows).
## These are empirical outliers, not significance calls: with N_e ~ 25-30,
## drift alone moves whole linked windows, so there is no neutral null here.
##
## The two populations have different window boundaries, so the genome is cut
## at every boundary of either population into segments; each segment carries
## the F_ST of the Indian Chief window and of the Jarvis window that contain it.
## One point = one segment.
##
## Overlap test: base pairs that are top 5% in both populations, compared with
## 10,000 circular shifts of the Jarvis outlier track along the concatenated
## genome (10-kb bins). The shift keeps the size and clustering of the outlier
## windows, so linkage is preserved.
##
## Genes: windows are intersected here with the B73 RefGen_v5 gene models (GFF3).
## The gene lists stored in the older window tables (tables/supplementary/
## SuppTable_7/8, GeneIDs column) are NOT used: they were attached by row
## position after a re-sort (17_spline_Fst_Tajima_Pi.R, annotate_windows_with_genes
## uses queryHits() as window_row_id), so 7% of Indian Chief and 62% of Jarvis
## gene-window pairs pointed at the wrong window.
##
## Inputs  (GN_WIN): SuppTable_7_IC_GenWin_windows_Fst_Tajima_pi_genes.csv,
##                   SuppTable_8_J_GenWin_windows_Fst_Tajima_pi_genes.csv (windows, F_ST)
##         (GN_GFF): Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3
##         (GN_GOU): GO_universe_genomewide.txt (PANTHER names, families, GO; all genes)
## Outputs: final/main_figs/Fig5_IC_J_Fst_scan.{png,pdf}
##          final/supp_tables/SuppTable_S17-S20 (windows and candidate genes)
##          final/results/ic_j_fst/IC_J_Fst_shared_regions.csv, _overlap_test.csv
################################################################################
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(ggrepel)})
gv   <- function(k, d) { v <- Sys.getenv(k, ""); if (nzchar(v)) v else d }
REPO <- gv("GN_REPO", "/Users/nirwantandukar/Documents/Github/Grain_N_remobilization")
WIN  <- gv("GN_WIN",  file.path(REPO, "tables/supplementary"))
GFF  <- gv("GN_GFF",  "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize/Maize.annotation/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3")
GOU  <- gv("GN_GOU",  file.path(REPO, "final/results/GO/GO_universe_genomewide.txt"))
OUT  <- gv("GN_MAINDIR", file.path(REPO, "final/main_figs"))
TAB  <- gv("GN_TABDIR",  file.path(REPO, "final/supp_tables"))
RES  <- gv("GN_RESDIR",  file.path(REPO, "final/results/ic_j_fst"))
for (d in c(OUT, TAB, RES)) dir.create(d, FALSE, TRUE)
set.seed(1)

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
## ggplot2 < 3.5 has no legend.position.inside
legend_inside <- function(x, y) {
  if (utils::packageVersion("ggplot2") >= "3.5.0")
    theme(legend.position = "inside", legend.position.inside = c(x, y), legend.justification = c(x, y))
  else theme(legend.position = c(x, y), legend.justification = c(x, y))
}

## ---- windows -------------------------------------------------------------------
rd <- function(f, pop) {
  d <- fread(file.path(WIN, f))
  d[, .(population = pop, chr, chr_num, window_id, start = window_start_bp, end = window_end_bp,
        n_sites = snp_count, mean_fst, top5 = fst_outlier == "Yes")]
}
W <- rbind(rd("SuppTable_7_IC_GenWin_windows_Fst_Tajima_pi_genes.csv", "Indian Chief"),
           rd("SuppTable_8_J_GenWin_windows_Fst_Tajima_pi_genes.csv", "Jarvis"))

## genes overlapping each window, from the GFF3 gene models
G <- fread(cmd = sprintf("awk -F'\t' '$3==\"gene\"' '%s'", GFF), header = FALSE, sep = "\t",
           select = c(1, 4, 5, 9), col.names = c("chr", "start", "end", "attr"))
G[, GeneID := sub("^ID=([^;]+).*", "\\1", attr)][, attr := NULL]
G <- G[grepl("^chr[0-9]+$", chr)]
setkey(G, chr, start, end)
ov <- foverlaps(W[, .(population, chr, window_id, start, end)], G,
                by.x = c("chr", "start", "end"), type = "any", nomatch = NULL)
ov <- ov[, .(gene_count = .N, GeneIDs = paste(GeneID, collapse = ";")), by = .(population, chr, window_id)]
W <- merge(W, ov, by = c("population", "chr", "window_id"), all.x = TRUE, sort = FALSE)
W[is.na(gene_count), `:=`(gene_count = 0L, GeneIDs = "")]
cat("genes assigned to windows:", W[, sum(gene_count)], "\n")
thr <- W[, .(thr = quantile(mean_fst, 0.95, na.rm = TRUE), n = .N, n_top = sum(top5)), by = population]
stopifnot(all(W[, all(top5 == (mean_fst >= thr[population == .BY$population, thr])), by = population]$V1))
print(thr)
TIC <- thr[population == "Indian Chief", thr]; TJ <- thr[population == "Jarvis", thr]

## ---- segments: cut the genome at every boundary of either population --------------
seg <- rbindlist(lapply(sort(unique(W$chr_num)), function(k) {
  wi <- W[population == "Indian Chief" & chr_num == k][order(start)]
  wj <- W[population == "Jarvis" & chr_num == k][order(start)]
  lo <- max(min(wi$start), min(wj$start)); hi <- min(max(wi$end), max(wj$end))
  br <- sort(unique(c(wi$start, wj$start, wi$end + 1, wj$end + 1)))
  br <- br[br >= lo & br <= hi + 1]
  s  <- data.table(chr_num = k, start = head(br, -1), end = tail(br, -1) - 1)
  s[, `:=`(i = findInterval(start, wi$start), j = findInterval(start, wj$start))]
  s[, `:=`(fst_IC = wi$mean_fst[i], fst_J = wj$mean_fst[j], top_IC = wi$top5[i], top_J = wj$top5[j],
           win_IC = wi$window_id[i], win_J = wj$window_id[j])]
  s[start <= wi$end[i] & start <= wj$end[j]]
}))
seg[, len := end - start + 1]
seg[, class := fifelse(top_IC & top_J, "Top 5% in both",
                fifelse(top_IC, "Top 5% Indian Chief only",
                fifelse(top_J, "Top 5% Jarvis only", "Other")))]
cat("segments:", nrow(seg), " genome covered (Mb):", round(sum(seg$len) / 1e6, 1), "\n")
rho <- cor(seg$fst_IC, seg$fst_J, method = "spearman")
rho_w <- cov.wt(cbind(rank(seg$fst_IC), rank(seg$fst_J)), wt = seg$len, cor = TRUE)$cor[1, 2]
cat(sprintf("Spearman r between IC and J window F_ST: %.3f (segment-length weighted %.3f)\n", rho, rho_w))

## ---- shared regions (adjacent shared segments merged) -----------------------------
sh <- seg[class == "Top 5% in both"][order(chr_num, start)]
sh[, grp := cumsum(c(TRUE, chr_num[-1] != chr_num[-.N] | start[-1] != end[-.N] + 1))]
genes_of <- function(pop, ids) {
  g <- W[population == pop & window_id %in% ids, GeneIDs]
  unique(unlist(strsplit(g[nzchar(g) & !is.na(g)], ";")))
}
REG <- sh[, {
  gi <- unique(unlist(lapply(unique(chr_num), function(k)
    W[population == "Indian Chief" & chr_num == k & window_id %in% win_IC, GeneIDs])))
  gj <- unique(unlist(lapply(unique(chr_num), function(k)
    W[population == "Jarvis" & chr_num == k & window_id %in% win_J, GeneIDs])))
  gi <- unique(unlist(strsplit(gi[!is.na(gi)], ";"))); gj <- unique(unlist(strsplit(gj[!is.na(gj)], ";")))
  .(chr = paste0("chr", chr_num[1]), start = min(start), end = max(end),
    fst_IC = max(fst_IC), fst_J = max(fst_J), shared_genes = length(intersect(gi, gj)),
    GeneIDs = paste(sort(intersect(gi, gj)), collapse = ";"))
}, by = grp]
REG[, region := sprintf("%s:%.1f-%.1f Mb", chr, start / 1e6, end / 1e6)]
REG[, length_kb := round((end - start + 1) / 1e3)]
print(REG[, .(region, length_kb, fst_IC = round(fst_IC, 3), fst_J = round(fst_J, 3), shared_genes)])
fwrite(REG[, .(region, chr, start, end, length_kb, fst_IC, fst_J, shared_genes, GeneIDs)],
       file.path(RES, "IC_J_Fst_shared_regions.csv"))

## ---- overlap test: circular shift of the Jarvis track (10-kb bins) ----------------
BIN <- 1e4
binvec <- function(top) {
  unlist(lapply(sort(unique(seg$chr_num)), function(k) {
    s <- seg[chr_num == k]; b <- seq(min(s$start), max(s$end), by = BIN)
    s[findInterval(b, s$start), get(top)]
  }))
}
bi <- binvec("top_IC"); bj <- binvec("top_J"); nb <- length(bi)
obs <- sum(bi & bj)
null <- vapply(sample.int(nb - 1, 10000, replace = TRUE), function(k)
  sum(bi & bj[c((k + 1):nb, 1:k)]), numeric(1))
TEST <- data.table(
  bins_10kb = nb, IC_top5_Mb = sum(bi) * BIN / 1e6, J_top5_Mb = sum(bj) * BIN / 1e6,
  shared_Mb = obs * BIN / 1e6, expected_independent_Mb = sum(bi) * sum(bj) / nb * BIN / 1e6,
  shift_mean_Mb = mean(null) * BIN / 1e6, shift_95pct_Mb = quantile(null, 0.95) * BIN / 1e6,
  fold_vs_shift = obs / mean(null), P_shift = (sum(null >= obs) + 1) / (length(null) + 1),
  n_shared_regions = nrow(REG), n_shared_genes = length(unique(unlist(strsplit(REG$GeneIDs, ";")))),
  spearman_r = rho)
print(TEST)
fwrite(TEST, file.path(RES, "IC_J_Fst_overlap_test.csv"))

## ---- figure ------------------------------------------------------------------------
COLS <- c("Other" = "grey78", "Top 5% Indian Chief only" = "#2166AC",
          "Top 5% Jarvis only" = "#B35806", "Top 5% in both" = "#B2182B")
seg[, class := factor(class, levels = names(COLS))]
lab <- REG[, .(fst_IC, fst_J, label = sprintf("%s %.1f Mb", chr, (start + end) / 2e6))]
p <- ggplot(seg[order(class)], aes(fst_IC, fst_J, colour = class)) +
  geom_vline(xintercept = TIC, linetype = "dashed", colour = "grey30", linewidth = 0.6) +
  geom_hline(yintercept = TJ, linetype = "dashed", colour = "grey30", linewidth = 0.6) +
  geom_point(aes(size = class == "Top 5% in both"), alpha = 0.85, stroke = 0) +
  scale_size_manual(values = c(`FALSE` = 1.6, `TRUE` = 3.2), guide = "none") +
  geom_text_repel(data = lab, aes(fst_IC, fst_J, label = label), inherit.aes = FALSE,
                  size = 4.2, colour = "#B2182B", min.segment.length = 0, box.padding = 0.5,
                  max.overlaps = Inf, seed = 1) +
  scale_colour_manual(values = COLS, breaks = names(COLS)[-1]) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  labs(x = expression(bold(paste(italic(F)[ST], "  Cycle 0 vs Cycle 14, Indian Chief"))),
       y = expression(bold(paste(italic(F)[ST], "  Cycle 0 vs Cycle 14, Jarvis")))) +
  coord_cartesian(clip = "off") +
  plot_theme + legend_inside(0.98, 0.98)
for (ext in c("png", "pdf"))
  ggsave(file.path(OUT, paste0("Fig5_IC_J_Fst_scan.", ext)), p, width = 10, height = 9, dpi = 300,
         bg = "white")

## ---- supplementary tables S17-S20 ----------------------------------------------------
shared_genes <- unique(unlist(strsplit(REG$GeneIDs, ";")))
for (pp in list(list("Indian Chief", "S17_IC_Fst_windows", "S19_IC_Fst_top5_genes"),
                list("Jarvis", "S18_J_Fst_windows", "S20_J_Fst_top5_genes"))) {
  w <- W[population == pp[[1]]][order(chr_num, start)]
  fwrite(w[, .(population, chr, window_id, start_bp = start, end_bp = end, n_sites,
               mean_fst = signif(mean_fst, 5), top5 = top5, gene_count, GeneIDs)],
         file.path(TAB, paste0("SuppTable_", pp[[2]], ".csv")))
  g <- w[top5 == TRUE, .(GeneID = unlist(strsplit(GeneIDs, ";"))),
         by = .(chr, window_id, start_bp = start, end_bp = end, mean_fst)]
  g <- g[, .SD[which.max(mean_fst)], by = GeneID]            # a gene spanning two windows: keep the higher
  a <- fread(GOU, quote = "\"")
  a[, GeneName := sub(";.*", "", GeneName)]
  g <- merge(g, a, by = "GeneID", all.x = TRUE, sort = FALSE)
  g[, shared_with_other_population := GeneID %in% shared_genes]
  setorder(g, -mean_fst, chr, start_bp, GeneID)
  g[, mean_fst := signif(mean_fst, 5)]
  setcolorder(g, c("GeneID", "GeneName", "chr", "window_id", "start_bp", "end_bp", "mean_fst",
                   "shared_with_other_population"))
  fwrite(g, file.path(TAB, paste0("SuppTable_", pp[[3]], ".csv")))
  cat(pp[[1]], ": windows", nrow(w), " top-5% genes", nrow(g), " annotated", sum(!is.na(g$GeneName)),
      " shared", sum(g$shared_with_other_population), "\n")
}
