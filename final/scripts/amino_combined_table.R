################################################################################
## One main table: strongest loci within each amino acid biosynthetic family,
## plus an "All phenotypes" block ranked by pleiotropy (the old Table 1).
##
## Everything is at the LOCUS level -- significant SNPs within 25 kb are merged,
## the same window used for gene annotation -- so the family blocks and the
## all-phenotype block are directly comparable.
################################################################################
suppressPackageStartupMessages({library(data.table)})
gv <- function(k, d) { v <- Sys.getenv(k); if (nzchar(v)) v else d }
SUPP <- gv("GN_SUPP", "/Users/nirwantandukar/Documents/Github/Grain_N_remobilization/final/supp_tables")
GODIR<- gv("GN_GO",   "/Users/nirwantandukar/Documents/Github/Grain_N_remobilization/final/results/GO")
OUT  <- gv("GN_OUT",  "/Users/nirwantandukar/Documents/Github/Grain_N_remobilization/final/main_tables")
NTOP <- as.integer(gv("GN_NTOP", "5")); LOCUS_GAP <- 25000
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

AA_FAM <- c(D="Aspartate",N="Aspartate",T="Aspartate",M="Aspartate",K="Aspartate",I="Aspartate",
            E="Glutamate",Q="Glutamate",P="Glutamate",R="Glutamate",H="Glutamate",
            S="Serine",G="Serine",C="Serine",
            F="Aromatic",Y="Aromatic",W="Aromatic",
            A="Pyruvate",V="Pyruvate",L="Pyruvate")
COMPOSITE <- c(IMNTDK="Aspartate",EHPRQ="Glutamate",SGC="Serine",
               WFY="Aromatic",LAV="Pyruvate",IVL="Pyruvate")
FAM_LEVELS <- c("Aspartate","Glutamate","Serine","Aromatic","Pyruvate")

assign_family <- function(p) {
  if (grepl("^Total", p)) return("Whole pool")
  if (grepl("\\.", p)) {
    q <- strsplit(p, ".", fixed = TRUE)[[1]]
    if (q[2] %in% names(COMPOSITE)) return(unname(COMPOSITE[q[2]]))
    if (q[1] %in% names(COMPOSITE)) return(unname(COMPOSITE[q[1]]))
    if (nchar(q[1]) == 1 && q[1] %in% names(AA_FAM)) return(unname(AA_FAM[q[1]]))
    return("Whole pool")
  }
  if (p %in% names(COMPOSITE)) return(unname(COMPOSITE[p]))
  if (nchar(p) == 1 && p %in% names(AA_FAM)) return(unname(AA_FAM[p]))
  "Whole pool"
}

S5 <- fread(file.path(SUPP, "SuppTable_S5_amino_gwas_candidate_genes.csv"))
S5[, Family := vapply(Phenotype, assign_family, character(1))]

go <- rbindlist(lapply(c("GO_all_maize_1.txt","GO_all_maize_2.txt"),
        function(f) fread(file.path(GODIR,f), sep="\t", header=FALSE, fill=TRUE)), fill=TRUE)
setnames(go, c("long","GeneID","Desc","PantherFamily","ProtClass","GO_MF","GO_BP","GO_CC"))
go[, GeneName := tstrsplit(Desc, ";", fixed=TRUE)[[1]]]
go <- unique(go[, .(GeneID, GeneName)])[, .SD[1], by = GeneID]

## ---- merge SNPs into loci within a scope ("all" or one family) -------------
build_loci <- function(d, scope) {
  s <- unique(d[, .(SNP, Chr, Pos)])
  setorder(s, Chr, Pos)
  s[, newl := c(TRUE, (Pos[-1] - Pos[-.N] > LOCUS_GAP) | (Chr[-1] != Chr[-.N]))]
  s[, locus := cumsum(newl)]
  d <- merge(d, s[, .(SNP, locus)], by = "SNP")
  L <- d[, .(Scope = scope, Chr = Chr[1], start = min(Pos), end = max(Pos),
             n_snp = uniqueN(SNP), best_P = min(min_P),
             n_phen = uniqueN(Phenotype),
             phenos = paste(sort(unique(Phenotype)), collapse = ", "),
             models = models_sig[which.min(min_P)]), by = locus]
  G <- d[!is.na(GeneID) & GeneID != "",
         .(GeneID = unique(GeneID[order(Distance_to_Gene_bp)])), by = locus]
  G <- merge(G, go, by = "GeneID", all.x = TRUE, sort = FALSE)
  G[, lab := ifelse(is.na(GeneName) | GeneName == "" | GeneName == "Uncharacterized protein",
                    paste0(GeneID, " (uncharacterised)"),
                    paste0(GeneID, " (", GeneName, ")"))]
  G <- G[, .(n_genes = .N, genes = paste(lab, collapse = ", ")), by = locus]
  merge(L, G, by = "locus", all.x = TRUE)
}

fam <- rbindlist(lapply(FAM_LEVELS, function(f) {
  d <- build_loci(S5[Family == f], f); setorder(d, best_P); head(d, NTOP)
}))
all_l <- build_loci(S5, "All phenotypes")
setorder(all_l, -n_phen, best_P)
all_l <- head(all_l, NTOP)
tab <- rbind(fam, all_l)
fwrite(tab, file.path(OUT, "amino_family_and_all_top_loci.csv"))

cat("\n=== All-phenotype block (ranked by pleiotropy) ===\n")
print(all_l[, .(locus = sprintf("chr%d:%.1f Mb", Chr, start/1e6), n_phen,
                best_P = signif(best_P,3), n_snp, n_genes)])

## ---- LaTeX ------------------------------------------------------------------
esc <- function(x) gsub("([&%#_$])", "\\\\\\1", x)
sci <- function(p) { e <- floor(log10(p)); sprintf("$%.1f\\times10^{%d}$", p/10^e, e) }
phen_short <- function(s, k = 3) {
  v <- strsplit(s, ", ", fixed = TRUE)[[1]]
  if (length(v) <= k) paste(v, collapse = ", ")
  else paste0(paste(v[1:k], collapse = ", "), " (+", length(v) - k, ")")
}

lines <- c("\\begin{table*}[t]", "\\centering",
paste0("\\caption{Strongest loci in the Goodman--Buckler FAA GWAS, within each ",
"amino acid biosynthetic family and across all phenotypes. Phenotypes were ",
"assigned to a family by the amino acid they report on; for a ratio phenotype ",
"$X.Y$ the family is that of the denominator when the denominator is a family ",
"composite. Loci merge significant SNPs within 25 kb, the same window used for ",
"gene annotation; positions are B73 RefGen\\_v5 and give the first significant ",
"SNP in the locus. Family blocks list the five loci with the smallest $P$ ",
"($P \\le 1.20 \\times 10^{-8}$, Bonferroni); the final block lists the five ",
"loci significant for the most phenotypes, irrespective of family. All genes ",
"within 25 kb are given, ordered by distance from the lead SNP. Whole-pool ",
"phenotypes (Total, Total\\_N, Total\\_PBAA) belong to no single family and ",
"appear only in the final block.}"),
"\\label{tab:amino_gwas_top_loci}", "\\footnotesize",
"\\setlength{\\tabcolsep}{4pt}",
"\\begin{tabularx}{\\textwidth}{@{}l l r r l >{\\raggedright\\arraybackslash}X@{}}",
"\\hline",
"Family & Locus & Best $P$ & $n$ phen. & Phenotype(s) & Candidate genes within 25 kb (putative function) \\\\",
"\\hline")
blocks <- c(FAM_LEVELS, "All phenotypes")
for (bi in seq_along(blocks)) {
  b <- blocks[bi]; d <- tab[Scope == b]
  if (!nrow(d)) next
  if (b == "All phenotypes") setorder(d, -n_phen, best_P) else setorder(d, best_P)
  for (i in seq_len(nrow(d))) {
    lines <- c(lines, sprintf("%s & chr%d:%.1f Mb & %s & %d & %s & %s \\\\",
      if (i == 1) esc(b) else "", d$Chr[i], d$start[i]/1e6, sci(d$best_P[i]),
      d$n_phen[i], paste0("\\texttt{", esc(phen_short(d$phenos[i])), "}"),
      esc(ifelse(is.na(d$genes[i]), "no gene within 25 kb", d$genes[i]))))
  }
  if (bi < length(blocks)) lines <- c(lines, "\\hline")
}
lines <- c(lines, "\\hline", "\\end{tabularx}", "\\end{table*}")
writeLines(lines, file.path(OUT, "MainTable_amino_gwas_top_loci.tex"))
cat("\nwritten:", file.path(OUT, "MainTable_amino_gwas_top_loci.tex"), "\n")
