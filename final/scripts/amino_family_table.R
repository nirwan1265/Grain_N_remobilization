################################################################################
## Top loci per amino-acid biosynthetic family, Goodman-Buckler FAA GWAS.
##
## Every scanned phenotype is assigned to one of the six biosynthetic families
## by the amino acid it reports on.  For a ratio phenotype X.Y the family is Y
## when Y is a family composite (IMNTDK, EHPRQ, SGC, WFY, LAV, IVL) and the
## family of X otherwise; whole-pool phenotypes (Total, Total_N, Total_PBAA)
## are kept apart because they belong to no single family.
##
## Within each family the significant SNPs are collapsed (one row per SNP, best
## P across that family's phenotypes) and the strongest N are reported with the
## genes within 25 kb and their PANTHER descriptions.
################################################################################
suppressPackageStartupMessages({library(data.table)})

gv <- function(k, d) { v <- Sys.getenv(k); if (nzchar(v)) v else d }
SUPP <- gv("GN_SUPP", "/Users/nirwantandukar/Documents/Github/Grain_N_remobilization/final/supp_tables")
GODIR<- gv("GN_GO",   "/Users/nirwantandukar/Documents/Github/Grain_N_remobilization/final/results/GO")
OUT  <- gv("GN_OUT",  "/Users/nirwantandukar/Documents/Github/Grain_N_remobilization/final/main_tables")
N_PER_FAMILY <- as.integer(gv("GN_N_PER_FAMILY", "5"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

## ---- amino acid -> biosynthetic family -------------------------------------
AA_FAM <- c(D = "Aspartate", N = "Aspartate", T = "Aspartate",
            M = "Aspartate", K = "Aspartate", I = "Aspartate",
            E = "Glutamate", Q = "Glutamate", P = "Glutamate",
            R = "Glutamate", H = "Glutamate",
            S = "Serine", G = "Serine", C = "Serine",
            F = "Aromatic", Y = "Aromatic", W = "Aromatic",
            A = "Pyruvate", V = "Pyruvate", L = "Pyruvate")
COMPOSITE <- c(IMNTDK = "Aspartate", EHPRQ = "Glutamate", SGC = "Serine",
               WFY = "Aromatic", LAV = "Pyruvate", IVL = "Pyruvate")
FAM_LEVELS <- c("Aspartate", "Glutamate", "Serine", "Aromatic", "Pyruvate")

assign_family <- function(p) {
  if (grepl("^Total", p)) return("Whole pool")
  if (grepl("\\.", p)) {
    parts <- strsplit(p, ".", fixed = TRUE)[[1]]
    den <- parts[2]
    if (den %in% names(COMPOSITE)) return(unname(COMPOSITE[den]))
    num <- parts[1]
    if (num %in% names(COMPOSITE)) return(unname(COMPOSITE[num]))
    if (nchar(num) == 1 && num %in% names(AA_FAM)) return(unname(AA_FAM[num]))
    return("Whole pool")
  }
  if (p %in% names(COMPOSITE)) return(unname(COMPOSITE[p]))
  if (nchar(p) == 1 && p %in% names(AA_FAM)) return(unname(AA_FAM[p]))
  "Whole pool"
}

## ---- data -------------------------------------------------------------------
S5 <- fread(file.path(SUPP, "SuppTable_S5_amino_gwas_candidate_genes.csv"))
S5[, Family := vapply(Phenotype, assign_family, character(1))]

go <- rbindlist(lapply(c("GO_all_maize_1.txt", "GO_all_maize_2.txt"),
        function(f) fread(file.path(GODIR, f), sep = "\t", header = FALSE, fill = TRUE)),
        fill = TRUE)
setnames(go, c("long", "GeneID", "Desc", "PantherFamily", "ProtClass",
               "GO_MF", "GO_BP", "GO_CC"))
go[, GeneName := tstrsplit(Desc, ";", fixed = TRUE)[[1]]]
go <- unique(go[, .(GeneID, GeneName, PantherFamily)])[, .SD[1], by = GeneID]

## ---- counts per family ------------------------------------------------------
cnt <- S5[, .(phenotypes = uniqueN(Phenotype), SNPs = uniqueN(SNP),
              genes = uniqueN(GeneID)), by = Family]
cnt[, genes_per_phenotype := round(genes / phenotypes, 1)]
setorder(cnt, -genes)
cat("\n================ candidate genes by biosynthetic family ================\n")
print(cnt)
fwrite(cnt, file.path(OUT, "amino_family_counts.csv"))

## ---- top SNPs per family ----------------------------------------------------
snp <- S5[Family %in% FAM_LEVELS,
          .(min_P = min(min_P), n_pheno = uniqueN(Phenotype),
            phenos = paste(sort(unique(Phenotype)), collapse = ", "),
            models = models_sig[which.min(min_P)],
            n_models = n_models_sig[which.min(min_P)]),
          by = .(Family, SNP, Chr, Pos)]
setorder(snp, Family, min_P)
top <- snp[, head(.SD, N_PER_FAMILY), by = Family]

genes <- S5[, .(GeneID, Relation, Distance_to_Gene_bp), by = .(SNP, Family)]
genes <- unique(genes)[order(Distance_to_Gene_bp)]
genes <- merge(genes, go, by = "GeneID", all.x = TRUE)
genes[, lab := paste0(GeneID,
        ifelse(is.na(GeneName) | GeneName == "", "",
               paste0(" (", substr(GeneName, 1, 60), ")")))]
gsum <- genes[, .(n_genes = uniqueN(GeneID),
                  gene_list = paste(lab[order(Distance_to_Gene_bp)], collapse = "; ")),
              by = .(SNP, Family)]
top <- merge(top, gsum, by = c("SNP", "Family"), all.x = TRUE)
top[, Family := factor(Family, levels = FAM_LEVELS)]
setorder(top, Family, min_P)
fwrite(top, file.path(OUT, "amino_family_top_loci.csv"))

cat("\n================ top", N_PER_FAMILY, "SNPs per family ================\n")
for (f in FAM_LEVELS) {
  d <- top[Family == f]
  if (!nrow(d)) next
  cat("\n---- ", f, " family ----\n", sep = "")
  for (i in seq_len(nrow(d))) {
    cat(sprintf("%-14s chr%s:%s  P=%.2e  %d phenotype(s) [%s]\n",
                d$SNP[i], d$Chr[i], format(d$Pos[i], big.mark = ","),
                d$min_P[i], d$n_pheno[i], d$phenos[i]))
    cat(sprintf("   genes within 25 kb (%d): %s\n", d$n_genes[i],
                substr(d$gene_list[i], 1, 300)))
  }
}
cat("\nwritten:", file.path(OUT, "amino_family_top_loci.csv"), "\n")

################################################################################
## LaTeX main table, same style as MainTable1.
## SNP IDs are AGPv4 while Pos is v5, so only the v5 coordinate is printed --
## the SNP ID stays in the CSV.
################################################################################
esc <- function(x) { x <- gsub("([&%#_])", "\\\\\\1", x); gsub("~", "\\\\textasciitilde ", x) }
sci <- function(p) {
  e <- floor(log10(p)); m <- p / 10^e
  sprintf("$%.1f\\times10^{%d}$", m, e)
}
best_gene <- function(gl) {                     # first gene with a real description
  g <- strsplit(gl, "; ", fixed = TRUE)[[1]]
  named <- grep("\\(", g, value = TRUE)
  named <- named[!grepl("Uncharacterized protein", named)]
  pick <- if (length(named)) named[1] else g[1]
  id <- sub(" .*", "", pick)
  fn <- sub("^[^(]*\\((.*)\\)$", "\\1", pick); if (identical(fn, pick)) fn <- "---"
  c(id, fn, length(g))
}

lines <- c(
"\\begin{table*}[t]", "\\centering",
paste0("\\caption{Strongest associations within each amino acid biosynthetic ",
"family in the Goodman--Buckler FAA GWAS. Phenotypes were assigned to a family ",
"by the amino acid they report on; for a ratio phenotype $X.Y$ the family is ",
"that of the denominator when the denominator is a family composite. Within a ",
"family, SNPs were collapsed across phenotypes and the five with the smallest ",
"$P$ are shown ($P \\le 1.20 \\times 10^{-8}$, Bonferroni). Positions are B73 ",
"RefGen\\_v5. Where several genes fell within 25 kb, the first with an ",
"informative PANTHER description is shown and the total is given in ",
"parentheses. Whole-pool phenotypes (Total, Total\\_N, Total\\_PBAA) belong to ",
"no single family and are excluded.}"),
"\\label{tab:amino_family_top_loci}", "\\small",
"\\setlength{\\tabcolsep}{5pt}",
"\\begin{tabularx}{\\textwidth}{@{}l l r l l >{\\raggedright\\arraybackslash}X@{}}",
"\\hline",
"Family & Position & Best $P$ & Phenotype(s) & Candidate gene & Putative function \\\\",
"\\hline")
for (f in FAM_LEVELS) {
  d <- top[Family == f]; if (!nrow(d)) next
  for (i in seq_len(nrow(d))) {
    bg <- best_gene(d$gene_list[i])
    ph <- d$phenos[i]; if (nchar(ph) > 26) ph <- paste0(substr(ph, 1, 24), "\\ldots")
    lines <- c(lines, sprintf(
      "%s & chr%s:%.1f Mb & %s & \\texttt{%s} & \\texttt{%s%s} & %s \\\\",
      if (i == 1) f else "", d$Chr[i], d$Pos[i] / 1e6, sci(d$min_P[i]),
      esc(ph), bg[1], if (as.integer(bg[3]) > 1) sprintf(" (of %s)", bg[3]) else "",
      esc(bg[2])))
  }
  lines <- c(lines, "\\hline")
}
lines <- c(lines, "\\hline", "\\end{tabularx}", "\\end{table*}")
writeLines(lines, file.path(OUT, "MainTable2_amino_family_top_loci.tex"))
cat("written:", file.path(OUT, "MainTable2_amino_family_top_loci.tex"), "\n")
