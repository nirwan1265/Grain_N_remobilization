#!/usr/bin/env Rscript
################################################################################
## Goodman-Buckler amino acid GWAS: candidate genes (S6), loci (S7) and the
## main-text locus table, on B73 RefGen_v5 coordinates.
##
## WHY THIS SCRIPT EXISTS
##   The GWAS genotypes are HapMap 3.2.1 (hmp321_282_agpv4): the 'Pos' column of
##   every scan is AGPv4 and the SNP names carry the AGPv3 coordinate. Earlier
##   tables intersected that AGPv4 'Pos' directly with the v5 gene models, which
##   put most candidate genes on the wrong genes (v4 and v5 differ by -1 to +3.7
##   Mb along the chromosomes). Here every significant SNP is lifted AGPv4 -> v5
##   with the MaizeGDB/Ensembl chain file (liftover_v4_v5.R) before annotation.
##
## Significant SNPs: Benjamini-Hochberg FDR <= 0.05 per phenotype x model
##   (fdr_sig in the by_model_fdr files, MLM / BLINK / FarmCPU).
## Candidate genes: v5 gene models within 25 kb of the lifted SNP position.
## Loci: significant SNPs (any phenotype) within 25 kb of each other on v5.
## Annotation: PANTHER family, protein class and GO (GO_all_GWASs_LD.txt, with
##   GO_universe_genomewide.txt as fallback).
##
## Inputs  (env vars override the repo defaults)
##   GN_BYMODEL  final/results/goodman_buckler_amino_acid_gwas/by_model_fdr/by_model
##   GN_CHAIN    Sarah_amino_acid/chain/B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain.txt
##   GN_GFF      Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3
##   GN_GO       final/results/GO  (GO_all_GWASs_LD.txt, GO_universe_genomewide.txt)
## Outputs
##   final/supp_tables/SuppTable_S6_amino_gwas_candidate_genes.csv
##   final/supp_tables/SuppTable_S7_amino_gwas_loci.csv
##   final/main_tables/MainTable_amino_gwas_top_loci.tex
##   final/results/goodman_buckler_amino_acid_gwas/amino_gwas_v5_summary.txt
################################################################################
suppressPackageStartupMessages(library(data.table))
gv   <- function(k, d) { v <- Sys.getenv(k, ""); if (nzchar(v)) v else d }
REPO <- gv("GN_REPO", "/Users/nirwantandukar/Documents/Github/Grain_N_remobilization")
BYM  <- gv("GN_BYMODEL", file.path(REPO, "final/results/goodman_buckler_amino_acid_gwas/by_model_fdr/by_model"))
CHAIN<- gv("GN_CHAIN", "/Users/nirwantandukar/Documents/Research/data/Sarah_amino_acid/chain/B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain.txt")
GFF  <- gv("GN_GFF", "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize/Maize.annotation/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3")
GOD  <- gv("GN_GO", file.path(REPO, "final/results/GO"))
SUPP <- gv("GN_SUPP", file.path(REPO, "final/supp_tables"))
MAIN <- gv("GN_MAIN", file.path(REPO, "final/main_tables"))
RES  <- gv("GN_RES", file.path(REPO, "final/results/goodman_buckler_amino_acid_gwas"))
HERE <- gv("GN_SCRIPTS", file.path(REPO, "final/scripts"))
source(file.path(HERE, "liftover_v4_v5.R"))
WIN <- 25000; NTOP <- 3
log_ <- character(0); say <- function(...) { s <- sprintf(...); cat(s, "\n"); log_ <<- c(log_, s) }

## ---- significant SNPs, lifted to v5 --------------------------------------------
fs <- list.files(BYM, pattern = "\\.RData$", full.names = TRUE)
H <- rbindlist(lapply(fs, function(f) {
  e <- new.env(); load(f, e); d <- as.data.table(get(ls(e)[1], e))[fdr_sig == TRUE]
  if (!nrow(d)) return(NULL)
  d[, .(Phenotype = sub("_[A-Za-z]+\\.RData$", "", basename(f)), Model = model, SNP, Chr, Pos_v4 = Pos,
        P = P.value, FDR)]
}))
B <- read_chain_blocks(CHAIN)
u <- unique(H[, .(Chr, Pos_v4)]); u[, c("Pos_v5", "lift") := lift_v4_v5(Chr, Pos_v4, B)]
H <- merge(H, u, by = c("Chr", "Pos_v4"))
W <- H[, .(n_models_sig = uniqueN(Model), models_sig = paste(sort(unique(Model)), collapse = "+"),
           min_P = min(P), min_FDR = min(FDR)), by = .(Phenotype, SNP, Chr, Pos_v4, Pos_v5, lift)]
W[, replicated := n_models_sig >= 2]
## every FDR-significant phenotype x SNP (with or without a gene nearby)
fwrite(W[order(Phenotype, Chr, Pos_v5)], file.path(RES, "amino_gwas_v5_fdr_snps.csv"))
say("FDR-significant: %d SNP x phenotype rows, %d SNPs, %d phenotypes (lift: %s)",
    nrow(W), uniqueN(W$SNP), uniqueN(W$Phenotype), paste(names(table(u$lift)), table(u$lift), collapse = ", "))

## ---- genes within 25 kb (v5) ----------------------------------------------------
G <- fread(cmd = sprintf("awk -F'\\t' '$3==\"gene\"' '%s'", GFF), header = FALSE, sep = "\t",
           select = c(1, 4, 5, 9), col.names = c("chr", "GeneStart", "GeneEnd", "attr"))
G[, GeneID := sub("^ID=([^;]+).*", "\\1", attr)][, Chr := suppressWarnings(as.integer(sub("chr", "", chr)))]
G <- G[!is.na(Chr), .(Chr, GeneStart, GeneEnd, GeneID)]; setkey(G, Chr, GeneStart, GeneEnd)
q <- W[, .(Chr, s = pmax(1, Pos_v5 - WIN), e = Pos_v5 + WIN, idx = .I)]
ov <- foverlaps(q, G, by.x = c("Chr", "s", "e"), nomatch = NULL)
S6 <- cbind(W[ov$idx], ov[, .(GeneID, GeneStart, GeneEnd)])
S6[, Relation := fifelse(Pos_v5 >= GeneStart & Pos_v5 <= GeneEnd, "within", fifelse(Pos_v5 < GeneStart, "upstream", "downstream"))]
S6[, Distance_to_Gene_bp := fifelse(Relation == "within", 0, fifelse(Relation == "upstream", GeneStart - Pos_v5, Pos_v5 - GeneEnd))]
setorder(S6, min_P, Phenotype, Distance_to_Gene_bp)
setcolorder(S6, c("Phenotype", "SNP", "Chr", "Pos_v4", "Pos_v5", "lift", "n_models_sig", "models_sig", "min_P", "min_FDR", "replicated"))
fwrite(S6, file.path(SUPP, "SuppTable_S6_amino_gwas_candidate_genes.csv"))
pg <- S6[, .(np = uniqueN(Phenotype)), by = GeneID]
pt <- S6[, .(ng = uniqueN(GeneID)), by = Phenotype][order(-ng)]
say("S6: %d SNPs within 25 kb of a gene; %d candidate genes (%d phenotype-specific, %d shared) in %d phenotypes; most: %s",
    uniqueN(S6$SNP), nrow(pg), pg[np == 1, .N], pg[np > 1, .N], uniqueN(S6$Phenotype),
    paste(sprintf("%s=%d", head(pt, 4)$Phenotype, head(pt, 4)$ng), collapse = ", "))

## ---- annotation -----------------------------------------------------------------
ann <- unique(rbind(fread(file.path(GOD, "GO_all_GWASs_LD.txt"), quote = "\""),
                    fread(file.path(GOD, "GO_universe_genomewide.txt"), quote = "\""), fill = TRUE), by = "GeneID")
ann[, GeneName := sub(";.*", "", GeneName)]
ann <- ann[, .(GeneID, GeneName, Family_Subfamily, Protein_Class, GO_BP, GO_MF)]

## ---- families -------------------------------------------------------------------
AA_FAM <- c(D="Aspartate",N="Aspartate",T="Aspartate",M="Aspartate",K="Aspartate",I="Aspartate",
            E="Glutamate",Q="Glutamate",P="Glutamate",R="Glutamate",H="Glutamate",
            S="Serine",G="Serine",C="Serine", F="Aromatic",Y="Aromatic",W="Aromatic",
            A="Pyruvate",V="Pyruvate",L="Pyruvate")
COMPOSITE <- c(IMNTDK="Aspartate",EHPRQ="Glutamate",SGC="Serine",WFY="Aromatic",LAV="Pyruvate",IVL="Pyruvate")
FAM <- c("Aspartate","Glutamate","Serine","Aromatic","Pyruvate")
fam_of <- function(p) {
  if (grepl("^Total", p)) return("Whole pool")
  if (grepl("\\.", p)) { q <- strsplit(p, ".", fixed = TRUE)[[1]]
    if (q[2] %in% names(COMPOSITE)) return(unname(COMPOSITE[q[2]]))
    if (q[1] %in% names(COMPOSITE)) return(unname(COMPOSITE[q[1]]))
    if (q[1] %in% names(AA_FAM)) return(unname(AA_FAM[q[1]])); return("Whole pool") }
  if (p %in% names(COMPOSITE)) return(unname(COMPOSITE[p]))
  if (p %in% names(AA_FAM)) return(unname(AA_FAM[p])); "Whole pool"
}
W[, Family := vapply(Phenotype, fam_of, "")]; S6[, Family := vapply(Phenotype, fam_of, "")]
cnt <- merge(W[, .(phenotypes = uniqueN(Phenotype), SNPs_all = uniqueN(SNP)), by = Family],
             S6[, .(SNPs_near_gene = uniqueN(SNP), genes = uniqueN(GeneID)), by = Family], by = "Family", all.x = TRUE)
cnt[, genes_per_phenotype := round(genes / phenotypes, 1)]; setorder(cnt, -genes)
say("family counts:\n%s", paste(capture.output(print(cnt)), collapse = "\n"))

## ---- loci (all phenotypes) = S7 ---------------------------------------------------
mk_loci <- function(d) {
  s <- unique(d[, .(SNP, Chr, Pos_v5)]); setorder(s, Chr, Pos_v5)
  s[, locus := cumsum(c(TRUE, Chr[-1] != Chr[-.N] | diff(Pos_v5) > WIN))]
  d <- merge(d, s[, .(SNP, locus)], by = "SNP")
  L <- d[, .(Chr = Chr[1], start_v5 = min(Pos_v5), end_v5 = max(Pos_v5), n_SNPs = uniqueN(SNP),
             n_phenotypes = uniqueN(Phenotype), phenotypes = paste(sort(unique(Phenotype)), collapse = ", "),
             best_P = min(min_P), lead_SNP = SNP[which.min(min_P)], lead_Pos_v4 = Pos_v4[which.min(min_P)],
             lead_Pos_v5 = Pos_v5[which.min(min_P)], replicated = any(replicated),
             lead_models = models_sig[which.min(min_P)],
             models_any = paste(sort(unique(unlist(strsplit(models_sig, "\\+")))), collapse = "+")), by = locus]
  ## genes within 25 kb of any SNP of the locus, ordered by distance to the lead SNP
  q <- L[, .(Chr, s = pmax(1, start_v5 - WIN), e = end_v5 + WIN, locus)]
  o <- foverlaps(q, G, by.x = c("Chr", "s", "e"), nomatch = NULL)
  o <- merge(o, L[, .(locus, lead_Pos_v5)], by = "locus")
  o[, d := fifelse(lead_Pos_v5 >= GeneStart & lead_Pos_v5 <= GeneEnd, 0, pmin(abs(lead_Pos_v5 - GeneStart), abs(lead_Pos_v5 - GeneEnd)))]
  o <- merge(o, ann, by = "GeneID", all.x = TRUE); setorder(o, locus, d)
  o[, informative := !is.na(Family_Subfamily) & nzchar(Family_Subfamily)]
  gl <- o[, .(n_genes = .N, genes = paste(GeneID, collapse = ";"),
              gene_labels = paste(sprintf("%s (%s)", GeneID, fifelse(is.na(GeneName) | GeneName %in% c("", "Uncharacterized protein"), "uncharacterised", GeneName)), collapse = ", "),
              top_gene = GeneID[order(-informative, d)][1]), by = locus]
  gl <- merge(gl, ann[, .(top_gene = GeneID, top_gene_name = GeneName, family = Family_Subfamily, pclass = Protein_Class, go_bp = GO_BP, go_mf = GO_MF)], by = "top_gene", all.x = TRUE)
  L <- merge(L, gl, by = "locus", all.x = TRUE)
  ## nearest gene for loci with no gene within 25 kb
  nn <- L[is.na(n_genes), {
    g <- G[Chr == .BY$Chr]; dd <- pmin(abs(lead_Pos_v5 - g$GeneStart), abs(lead_Pos_v5 - g$GeneEnd))
    .(nearest_gene = g$GeneID[which.min(dd)], nearest_dist_bp = min(dd)) }, by = .(locus, Chr)]
  L <- merge(L, nn[, .(locus, nearest_gene, nearest_dist_bp)], by = "locus", all.x = TRUE)
  L[is.na(n_genes), n_genes := 0L]
  L
}
S7 <- mk_loci(W); setorder(S7, -n_phenotypes, best_P); S7[, locus := .I]
fwrite(S7, file.path(SUPP, "SuppTable_S7_amino_gwas_loci.csv"))
say("S7: %d loci (%d with a gene within 25 kb); top: %s", nrow(S7), S7[n_genes > 0, .N],
    paste(S7[1:5, sprintf("chr%d:%.2f Mb n=%d P=%.1e", Chr, lead_Pos_v5 / 1e6, n_phenotypes, best_P)], collapse = " | "))

## ---- main-text table: per family and across all phenotypes ------------------------
blocks <- rbindlist(lapply(c("All phenotypes", FAM), function(b) {
  d <- if (b == "All phenotypes") W else W[Family == b]
  L <- mk_loci(d)[n_phenotypes >= 2]; setorder(L, -n_phenotypes, best_P); L <- head(L, NTOP)
  if (nrow(L)) L[, Block := b]; L }), fill = TRUE)
fwrite(blocks, file.path(RES, "amino_gwas_v5_top_loci.csv"))
esc <- function(x) gsub("([&%#_$])", "\\\\\\1", x)
sci <- function(p) { e <- floor(log10(p)); sprintf("$%.1f\\times10^{%d}$", p / 10^e, e) }
phs <- function(s, k = 3) { v <- strsplit(s, ", ", fixed = TRUE)[[1]]
  if (length(v) <= k) paste(v, collapse = ", ") else paste0(paste(v[1:k], collapse = ", "), " (+", length(v) - k, ")") }
short_lab <- function(x, k = 4) { v <- strsplit(x, ", (?=Zm00001eb)", perl = TRUE)[[1]]
  if (length(v) > k) paste0(paste(v[1:k], collapse = ", "), sprintf(" (+%d more)", length(v) - k)) else paste(v, collapse = ", ") }
tex <- c("\\begin{table*}[t]", "\\centering",
  paste0("\\caption{Loci associated with the most amino acid phenotypes in the Goodman--Buckler FAA GWAS, ",
         "across all phenotypes and within each amino acid biosynthetic family. Significant SNPs (5\\% Benjamini--Hochberg FDR) ",
         "within 25~kb of each other form one locus; for each block the three loci significant for the most phenotypes are shown ",
         "(ties broken by the smallest $P$). Positions are B73 RefGen\\_v5, lifted from the AGPv4 coordinates of the genotype data, ",
         "and give the lead SNP. Genes within 25~kb of the locus are listed by distance from the lead SNP; ",
         "the full list is in Supplementary Table~S7.}"),
  "\\label{tab:amino_gwas_top_loci}", "\\footnotesize", "\\setlength{\\tabcolsep}{4pt}",
  "\\begin{tabularx}{\\textwidth}{@{}l l r r l >{\\raggedright\\arraybackslash}X@{}}", "\\hline",
  "Family & Locus (v5) & Best $P$ & $n$ phen. & Phenotype(s) & Candidate genes within 25 kb (putative function) \\\\", "\\hline")
bl <- unique(blocks$Block)
for (bi in seq_along(bl)) {
  d <- blocks[Block == bl[bi]]
  for (i in seq_len(nrow(d))) tex <- c(tex, sprintf("%s & chr%d:%.2f Mb & %s & %d & \\texttt{%s} & %s \\\\",
    if (i == 1) esc(bl[bi]) else "", d$Chr[i], d$lead_Pos_v5[i] / 1e6, sci(d$best_P[i]), d$n_phenotypes[i],
    esc(phs(d$phenotypes[i])),
    esc(if (d$n_genes[i] > 0) short_lab(d$gene_labels[i]) else sprintf("none; nearest %s (%.1f kb)", d$nearest_gene[i], d$nearest_dist_bp[i] / 1000))))
  if (bi < length(bl)) tex <- c(tex, "\\hline")
}
tex <- c(tex, "\\hline", "\\end{tabularx}", "\\end{table*}")
writeLines(tex, file.path(MAIN, "MainTable_amino_gwas_top_loci.tex"))
writeLines(log_, file.path(RES, "amino_gwas_v5_summary.txt"))
cat("written S6, S7, main table, summary\n")
