#!/usr/bin/env Rscript
################################################################################
## Targeted look-up: do Asp-pathway genes found by the OTHER analyses (Indian
## Chief / Jarvis top-5% F_ST windows, soil-N GWAS) carry associations with
## Asp-family amino acid traits in the Goodman-Buckler amino acid GWAS?
##
## Genes (fixed before looking at the GWAS): enzymes of the aspartate hub
## (Asp -> Lys/Thr/Met/Ile, Asp <-> Glu/OAA, Asp -> Asn, PEPC anaplerosis,
## GS, IAA-Asp conjugation) present in S19, S20 or S8.
## Traits: all 23 aspartate-family phenotypes (D, N, T, M, K, I, the IMNTDK
## composite, their ratios to IMNTDK and to the total pool, N.D, N.E, N.Q,
## D.EHPRQ), in every model run for that trait (MLM, BLINK, FarmCPU).
##
## The by_model_fdr files keep only SNPs with P <= 1e-5, which is above the
## targeted threshold computed below, so "no SNP in the file" = not significant.
## Targeted Bonferroni: 0.05 / (SNPs in the gene windows x trait-model scans),
## with SNP density from the scan row counts (4,177,796 SNPs / 2.13 Gb).
##
## Input : final/results/goodman_buckler_amino_acid_gwas/by_model_fdr/by_model
## Output: final/results/ic_j_fst/Asp_pathway_targeted_lookup.csv
################################################################################
suppressPackageStartupMessages(library(data.table))
gv   <- function(k, d) { v <- Sys.getenv(k, ""); if (nzchar(v)) v else d }
REPO <- gv("GN_REPO", "/Users/nirwantandukar/Documents/Github/Grain_N_remobilization")
BYM  <- gv("GN_BYMODEL", file.path(REPO, "final/results/goodman_buckler_amino_acid_gwas/by_model_fdr/by_model"))
GFF  <- gv("GN_GFF", "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize/Maize.annotation/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3")
RES  <- gv("GN_RESDIR", file.path(REPO, "final/results/ic_j_fst"))
FLANK <- 25000; FLANK_WIDE <- 250000

GENES <- data.table(
  GeneID = c("Zm00001eb094670", "Zm00001eb152450", "Zm00001eb214360", "Zm00001eb294790",
             "Zm00001eb052750", "Zm00001eb346640", "Zm00001eb162900", "Zm00001eb151250",
             "Zm00001eb220040", "Zm00001eb005710"),
  enzyme = c("aspartate kinase", "aspartate aminotransferase (mitochondrial)",
             "aspartate-semialdehyde dehydrogenase", "threonine synthase",
             "methionine synthase (MetE)", "phosphoenolpyruvate carboxylase",
             "glutamine synthetase domain", "IAA-amido synthetase GH3.6",
             "aspartate kinase domain", "omega-amidase NIT2"),
  found_in = c("Jarvis F_ST", "Jarvis F_ST", "Indian Chief F_ST", "Indian Chief F_ST",
               "Indian Chief F_ST + soil-N GWAS", "Indian Chief F_ST", "Indian Chief F_ST",
               "Indian Chief F_ST", "soil-N GWAS", "soil-N GWAS"))
G <- fread(cmd = sprintf("awk -F'\\t' '$3==\"gene\"' '%s'", GFF), header = FALSE, sep = "\t",
           select = c(1, 4, 5, 9), col.names = c("chr", "start", "end", "attr"))
G[, GeneID := sub("^ID=([^;]+).*", "\\1", attr)]
GENES <- merge(GENES, G[, .(GeneID, Chr = as.integer(sub("chr", "", chr)), start, end)], by = "GeneID", sort = FALSE)

TRAITS <- c("D", "D.IMNTDK", "D.Total", "D.EHPRQ", "N", "N.IMNTDK", "N.Total", "N.D", "N.E", "N.Q",
            "T", "T.IMNTDK", "T.Total", "M", "M.IMNTDK", "M.Total", "K", "K.IMNTDK", "K.Total",
            "I", "I.IMNTDK", "I.Total", "IMNTDK")
files <- unlist(lapply(TRAITS, function(t) Sys.glob(file.path(BYM, paste0(t, "_*.RData")))))
files <- files[sub("_[A-Za-z]+\\.RData$", "", basename(files)) %in% TRAITS]
GW <- rbindlist(lapply(files, function(f) {
  e <- new.env(); load(f, e); x <- get(ls(e)[1], e)
  x[, .(trait = sub("_[A-Za-z]+\\.RData$", "", basename(f)), model, SNP, Chr, Pos, P.value, FDR)]
}))
n_scans <- length(files)

dens <- 4177796 / 2.13e9
n_snp_win <- GENES[, sum((end - start + 1 + 2 * FLANK) * dens)]
P_TARGET <- 0.05 / (n_snp_win * n_scans)
cat(sprintf("trait-model scans: %d | SNPs in gene +/- %d kb windows: ~%.0f | targeted Bonferroni P = %.1e\n",
            n_scans, FLANK / 1000, n_snp_win, P_TARGET))

hit <- function(g, fl) {
  x <- GW[Chr == g$Chr & Pos >= g$start - fl & Pos <= g$end + fl]
  if (!nrow(x)) return(data.table(n_snp_P_le_1e5 = 0L, min_P = NA_real_, trait = NA_character_,
                                    model = NA_character_, SNP = NA_character_, dist_kb = NA_real_))
  b <- x[which.min(P.value)]
  d <- if (b$Pos < g$start) g$start - b$Pos else if (b$Pos > g$end) b$Pos - g$end else 0
  data.table(n_snp_P_le_1e5 = nrow(x), min_P = b$P.value, trait = b$trait, model = b$model,
             SNP = b$SNP, dist_kb = round(d / 1000, 1))
}
OUT <- rbindlist(lapply(seq_len(nrow(GENES)), function(i) {
  g <- GENES[i]
  cbind(g[, .(GeneID, enzyme, found_in, chr = Chr, start, end)],
        setnames(hit(g, FLANK), function(n) paste0(n, "_25kb")),
        setnames(hit(g, FLANK_WIDE), function(n) paste0(n, "_250kb")))
}))
OUT[, significant_targeted := !is.na(min_P_25kb) & min_P_25kb <= P_TARGET]
OUT[, P_targeted_threshold := signif(P_TARGET, 3)]
print(OUT[, .(GeneID, enzyme, found_in, min_P_25kb, trait_25kb, model_25kb, dist_kb_25kb,
              min_P_250kb, trait_250kb, dist_kb_250kb, significant_targeted)])
dir.create(RES, FALSE, TRUE)
fwrite(OUT, file.path(RES, "Asp_pathway_targeted_lookup.csv"))
