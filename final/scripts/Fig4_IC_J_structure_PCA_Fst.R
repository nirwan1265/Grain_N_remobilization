#!/usr/bin/env Rscript
################################################################################
## FIGURE -- population structure, PCA and between-population differentiation
## of Indian Chief and Jarvis at Cycle 0 and Cycle 14 (skim sequencing, 379 plants)
##
## All three panels use GENOTYPE LIKELIHOODS, not hard-called genotypes. The
## samples are 0.6-0.9x and Cycle 14 was sequenced deeper than Cycle 0
## (IC 0.61x -> 0.88x, J 0.68x -> 0.76x); with hard calls, depth shows up as
## structure (it drove the old SNP-sharing panels and the SNPRelate PCA).
##
## A  NGSadmix (unsupervised), K = 2-4, seed with the highest likelihood.
##    K = 4 is the largest K at which independent seeds converge (K = 5 seeds
##    disagree and gain little likelihood).
## B  PCAngsd PCA (genotype-likelihood covariance). Low-depth samples carry less
##    information and are pulled toward the centre, so within a group the PC
##    score correlates with depth; the separation BETWEEN groups is the signal.
## C  Hudson F_ST between Indian Chief and Jarvis at Cycle 0 and at Cycle 14,
##    per chromosome, from the ANGSD/winsfs 2D site-frequency spectra
##    (ratio of sums over all sites; 96 diploid samples per group).
##
## Inputs
##   GN_STRUCT : dir with IC_J_all_sample_order.txt, IC_J_pcangsd.cov,
##               IC_J_K{K}_s{seed}.qopt / .log  (final/scripts/IC_J_structure_pipeline.sh)
##   GN_SFS    : dir with <chr>_IC01_J01.ml and <chr>_IC14_J14.ml
## Output: final/main_figs/Fig4_IC_J_structure.{png,pdf}
##         final/results/ic_j_structure/IC_J_between_pop_Fst.csv
################################################################################
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
gv   <- function(k, d) { v <- Sys.getenv(k, ""); if (nzchar(v)) v else d }
REPO <- gv("GN_REPO", "/Users/nirwantandukar/Documents/Github/Grain_N_remobilization")
STR  <- gv("GN_STRUCT", file.path(REPO, "final/results/ic_j_structure"))
SFS  <- gv("GN_SFS",    file.path(REPO, "data/Between_pop_Fst_IC_J"))
OUT  <- gv("GN_MAINDIR", file.path(REPO, "final/main_figs"))
RES  <- gv("GN_RESDIR",  STR)
dir.create(OUT, FALSE, TRUE); dir.create(RES, FALSE, TRUE)

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
GRP <- c(I01 = "Indian Chief C0", I14 = "Indian Chief C14", J01 = "Jarvis C0", J14 = "Jarvis C14")
COL <- c("Indian Chief C0" = "#92C5DE", "Indian Chief C14" = "#2166AC",
         "Jarvis C0" = "#FDB863", "Jarvis C14" = "#B35806")

ids <- readLines(file.path(STR, "IC_J_all_sample_order.txt"))
meta <- data.table(id = ids, group = factor(GRP[sub("_.*", "", ids)], levels = GRP))

## ---- A: NGSadmix ---------------------------------------------------------------
best_seed <- function(K) {
  lg <- Sys.glob(file.path(STR, sprintf("IC_J_K%d_s*.log", K)))
  ll <- sapply(lg, function(f) as.numeric(sub(".*best like=([-0-9.]+).*", "\\1",
                                               grep("best like", readLines(f), value = TRUE)[1])))
  sub("\\.log$", ".qopt", lg[which.max(ll)])
}
lik <- rbindlist(lapply(2:5, function(K) {
  lg <- Sys.glob(file.path(STR, sprintf("IC_J_K%d_s*.log", K)))
  ll <- sapply(lg, function(f) as.numeric(sub(".*best like=([-0-9.]+).*", "\\1",
                                               grep("best like", readLines(f), value = TRUE)[1])))
  data.table(K = K, n_seeds = length(ll), best_loglik = max(ll), spread = max(ll) - min(ll))
}))
lik[, gain := c(NA, diff(best_loglik))]
print(lik)
fwrite(lik, file.path(RES, "IC_J_ngsadmix_likelihoods.csv"))

## order samples: by group, then by their K = 3 membership of the Indian Chief C14 cluster
q4 <- as.matrix(fread(best_seed(4), header = FALSE))
q3 <- as.matrix(fread(best_seed(3), header = FALSE))
lead3 <- which.max(colMeans(q3[meta$group == "Indian Chief C14", , drop = FALSE]))
meta[, ord := order(order(group, -q3[, lead3]))]
adm <- rbindlist(lapply(2:4, function(K) {
  q <- as.matrix(fread(best_seed(K), header = FALSE))
  # colour each cluster by the group with the highest mean membership in it
  lab <- sapply(seq_len(K), function(k) names(which.max(tapply(q[, k], meta$group, mean))))
  d <- data.table(id = rep(ids, K), x = rep(meta$ord, K), group = rep(meta$group, K),
                  cluster = rep(seq_len(K), each = length(ids)), q = as.vector(q),
                  colour_group = rep(lab, each = length(ids)))
  d[, K := paste0("K = ", K)]
}))
bounds <- meta[, .(start = min(ord) - 0.5, end = max(ord) + 0.5, mid = mean(range(ord))), by = group]
pA <- ggplot(adm, aes(x, q, fill = colour_group, group = interaction(cluster, K))) +
  geom_col(width = 1, position = "stack", linewidth = 0) +
  geom_vline(xintercept = bounds$end[-nrow(bounds)], colour = "black", linewidth = 0.6) +
  facet_grid(K ~ ., switch = "y") +
  scale_fill_manual(values = COL, guide = "none") +
  scale_x_continuous(breaks = bounds$mid, labels = bounds$group, expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0.25, 0.5, 0.75), labels = c("", "0.5", ""), expand = c(0, 0)) +
  labs(x = NULL, y = "Ancestry proportion", title = "Unsupervised structure (NGSadmix)") +
  plot_theme +
  theme(strip.text.y.left = element_text(size = 15, face = "bold", angle = 0),
        strip.placement = "outside", panel.spacing.y = unit(10, "pt"),
        axis.text.x = element_text(size = 15, face = "bold"), axis.line.x = element_blank(),
        axis.text.y = element_text(size = 11), plot.title = element_text(size = 17, face = "bold"))

## ---- B: PCAngsd ----------------------------------------------------------------
C <- as.matrix(fread(file.path(STR, "IC_J_pcangsd.cov"), header = FALSE))
e <- eigen(C, symmetric = TRUE)
ve <- 100 * e$values / sum(e$values)
pc <- data.table(meta, PC1 = e$vectors[, 1], PC2 = e$vectors[, 2])
# orient: Indian Chief C14 positive on PC1, Cycle 14 negative on PC2 (sign is arbitrary)
if (pc[group == "Indian Chief C14", mean(PC1)] < 0) pc[, PC1 := -PC1]
if (pc[group %in% c("Indian Chief C14", "Jarvis C14"), mean(PC2)] > 0) pc[, PC2 := -PC2]
cat(sprintf("PCAngsd: PC1 %.1f%%, PC2 %.1f%%\n", ve[1], ve[2]))
pB <- ggplot(pc, aes(PC1, PC2, colour = group)) +
  geom_point(size = 2.6, alpha = 0.8) +
  stat_ellipse(level = 0.9, linewidth = 0.9, show.legend = FALSE) +
  scale_colour_manual(values = COL) +
  labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]),
       title = "Genotype-likelihood PCA (PCAngsd)") +
  plot_theme + theme(legend.position = "bottom", plot.title = element_text(size = 17, face = "bold")) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 5)))

## ---- C: between-population F_ST --------------------------------------------------
hudson <- function(f) {
  S <- matrix(scan(f, quiet = TRUE), nrow = sqrt(length(scan(f, quiet = TRUE))), byrow = TRUE)
  n1 <- nrow(S) - 1; n2 <- ncol(S) - 1
  p1 <- matrix((0:n1) / n1, nrow(S), ncol(S)); p2 <- matrix((0:n2) / n2, nrow(S), ncol(S), byrow = TRUE)
  num <- (p1 - p2)^2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
  den <- p1 * (1 - p2) + p2 * (1 - p1)
  c(num = sum(S * num), den = sum(S * den))
}
chrs <- sprintf("NC_%06d.1", 50096:50105)
fst <- rbindlist(lapply(seq_along(chrs), function(i) {
  a <- hudson(file.path(SFS, paste0(chrs[i], "_IC01_J01.ml")))
  b <- hudson(file.path(SFS, paste0(chrs[i], "_IC14_J14.ml")))
  data.table(chromosome = i, accession = chrs[i], Fst_C0 = a["num"] / a["den"], Fst_C14 = b["num"] / b["den"],
             num_C0 = a["num"], den_C0 = a["den"], num_C14 = b["num"], den_C14 = b["den"])
}))
gw <- fst[, .(C0 = sum(num_C0) / sum(den_C0), C14 = sum(num_C14) / sum(den_C14))]
cat(sprintf("genome-wide Hudson F_ST, Indian Chief vs Jarvis: C0 %.4f, C14 %.4f (%d/10 chromosomes higher at C14)\n",
            gw$C0, gw$C14, fst[, sum(Fst_C14 > Fst_C0)]))
fwrite(fst[, .(chromosome, accession, Fst_C0 = round(Fst_C0, 4), Fst_C14 = round(Fst_C14, 4))],
       file.path(RES, "IC_J_between_pop_Fst.csv"))
fl <- melt(fst[, .(chromosome, C0 = Fst_C0, C14 = Fst_C14)], id.vars = "chromosome",
           variable.name = "cycle", value.name = "Fst")
fl[, cycle := factor(fifelse(cycle == "C0", "Cycle 0", "Cycle 14"), levels = c("Cycle 0", "Cycle 14"))]
gwl <- data.table(cycle = factor(c("Cycle 0", "Cycle 14"), levels = c("Cycle 0", "Cycle 14")),
                  Fst = c(gw$C0, gw$C14))
pC <- ggplot(fl, aes(cycle, Fst)) +
  geom_line(aes(group = chromosome), colour = "grey65", linewidth = 0.8) +
  geom_point(aes(fill = cycle), shape = 21, size = 4.2, colour = "black", stroke = 0.4) +
  geom_point(data = gwl, shape = 23, size = 7, fill = "black") +
  geom_text(data = gwl, aes(label = sprintf("%.3f", Fst)), hjust = c(1.5, -0.5),
            size = 6, fontface = "bold") +
  scale_fill_manual(values = c("Cycle 0" = "grey80", "Cycle 14" = "grey35"), guide = "none") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = expression(bold(italic(F)[ST]*" (Indian Chief vs Jarvis)")),
       title = "Between-population differentiation") +
  plot_theme + theme(plot.title = element_text(size = 17, face = "bold"),
                     axis.text.x = element_text(size = 17, face = "bold"))

fig <- pA / (pB | pC) + plot_layout(heights = c(1, 1.3)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 28, face = "bold"))
for (ext in c("png", "pdf"))
  ggsave(file.path(OUT, paste0("Fig4_IC_J_structure.", ext)), fig, width = 20, height = 16,
         dpi = 300, bg = "white")
cat("written\n")
