#!/usr/bin/env Rscript
################################################################################
## Do Goodman-Buckler amino acid GWAS loci fall in Indian Chief / Jarvis top-5%
## F_ST windows more often than expected by chance?
##
## Locus-level, not gene-level: each GWAS locus (S7; significant SNPs merged
## within 25 kb) is placed at its lead SNP (B73 v5) and counted as inside or
## outside the top-5% windows of each population (S17, S18). Null: 10,000
## circular shifts of the top-5% window track along the concatenated genome,
## which keeps window sizes and clustering. Done for all loci and for loci with
## at least one aspartate-family phenotype (D, N, T, M, K, I, IMNTDK and their
## ratios).
##
## Output: final/results/ic_j_fst/GWAS_loci_in_Fst_windows.csv and the list of
##         Asp-family loci inside top-5% windows.
################################################################################
suppressPackageStartupMessages(library(data.table))
gv   <- function(k, d) { v <- Sys.getenv(k, ""); if (nzchar(v)) v else d }
REPO <- gv("GN_REPO", "/Users/nirwantandukar/Documents/Github/Grain_N_remobilization")
TAB  <- gv("GN_TABDIR", file.path(REPO, "final/supp_tables"))
RES  <- gv("GN_RESDIR", file.path(REPO, "final/results/ic_j_fst"))
set.seed(1); NSHIFT <- 10000

L <- fread(file.path(TAB, "SuppTable_S7_amino_gwas_loci.csv"))
L[, asp := vapply(strsplit(phenotypes, ", "), function(v)
  any(grepl("^(D|N|T|M|K|I)(\\.|$)", v) | v == "IMNTDK"), logical(1))]
cat("GWAS loci:", nrow(L), " with an Asp-family phenotype:", L[, sum(asp)], "\n")

run <- function(pop, f) {
  W <- fread(file.path(TAB, f))[, chr_num := as.integer(sub("chr", "", chr))]
  len <- W[, .(len = max(end_bp)), by = chr_num][order(chr_num)]
  len[, off := cumsum(c(0, head(len, -1)))]; G <- sum(len$len)
  W <- merge(W, len[, .(chr_num, off)], by = "chr_num")
  top <- W[top5 == TRUE, .(s = start_bp + off, e = end_bp + off)][order(s)]
  pos_all <- merge(L[, .(locus, chr_num = Chr, x = lead_Pos, asp)], len[, .(chr_num, off)], by = "chr_num")[, .(locus, g = x + off, asp)]
  inside <- function(g) { i <- findInterval(g, top$s); i > 0 & g <= top$e[pmax(i, 1)] }
  one <- function(sel) {
    g <- pos_all[sel, g]; obs <- sum(inside(g))
    null <- vapply(runif(NSHIFT, 0, G), function(k) sum(inside((g - k) %% G)), numeric(1))
    data.table(population = pop, loci = if (identical(sel, TRUE)) "all" else "Asp-family",
               n_loci = length(g), top5_genome_frac = round(top[, sum(e - s + 1)] / G, 4),
               observed_in_top5 = obs, expected_shift = round(mean(null), 2),
               P_shift = signif((sum(null >= obs) + 1) / (NSHIFT + 1), 3))
  }
  list(stats = rbind(one(TRUE), one(pos_all$asp)),
       hits = {ids <- pos_all[asp == TRUE][inside(g), locus]; h <- L[locus %in% ids]; h[, population := pop]; h})
}
r <- list(run("Indian Chief", "SuppTable_S17_IC_Fst_windows.csv"), run("Jarvis", "SuppTable_S18_J_Fst_windows.csv"))
ST <- rbindlist(lapply(r, `[[`, "stats")); print(ST)
H <- rbindlist(lapply(r, `[[`, "hits"))
print(H[, .(population, Chr, lead_Pos, phenotypes = substr(phenotypes, 1, 60), best_P = signif(best_P, 3),
            replicated, top_gene, family = substr(family, 1, 40))])
dir.create(RES, FALSE, TRUE)
fwrite(ST, file.path(RES, "GWAS_loci_in_Fst_windows.csv"))
fwrite(H[, .(population, locus, Chr, lead_SNP, lead_Pos, phenotypes, best_P, replicated, lead_models,
             genes, top_gene, family)], file.path(RES, "Asp_GWAS_loci_in_Fst_windows.csv"))
