#!/usr/bin/env Rscript
################################################################################
## MAIN TABLE + SUPPLEMENTARY TABLES -- Goodman-Buckler amino acid GWAS
##
## Rebuilds every number in Results paragraph 3 and the main candidate-locus
## table, starting from the per-model .RData files written by
## final/scripts/hpc/30_split_trait_by_model.R.
##
## WHY LOCI AND NOT GENES
##   A +/- 25 kb window around one SNP routinely catches several neighbouring
##   genes, so ranking genes by phenotype count produces near-duplicate rows
##   that are really one signal (e.g. four genes sharing a single SNP). Loci are
##   therefore defined by merging significant SNPs within LOCUS_GAP bp, and the
##   main table ranks loci. The per-gene table is still written for the
##   supplement.
##
## THRESHOLD
##   alpha/m with alpha = 0.05 and m = 4,177,796 markers.
##   NOTE: 4,177,796 is the value implied by the raw scan files
##   (8,355,592 rows / 2 models; 12,533,388 / 3 models). The manuscript Methods
##   previously read 4,117,796, a digit transposition.
##
## INPUTS
##   final/results/goodman_buckler_amino_acid_gwas/by_model/<trait>_<MODEL>.RData
##   <GFF v5>                                     B73 RefGen_v5 annotation
##   shiny/.../SuppTable1_GWAS_annotation_soilN_amino_with_GO.csv   PANTHER/GO
##
## OUTPUTS (final/supp_tables/ and final/main_tables/)
##   MainTable1_amino_gwas_top_loci.tex        main-text table, top 5 loci
##   SuppTable2_amino_gwas_candidate_genes.csv every Bonferroni gene
##   SuppTable2b_amino_gwas_genes_per_trait.csv
##   SuppTable3_amino_gwas_replicated_genes.csv
##   SuppTable4_amino_gwas_loci_ranked.csv     all loci, ranked, with GO
##
## Usage: Rscript amino_gwas_main_table.R [repo_root]
################################################################################

suppressPackageStartupMessages(library(data.table))

ROOT <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "."

BY_MODEL <- file.path(ROOT, "final", "results", "goodman_buckler_amino_acid_gwas", "by_model")
GFF      <- Sys.getenv("GFF_PATH", file.path(ROOT, "..", "..", "Maize.annotation",
                      "Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3"))
GO_FILE  <- file.path(ROOT, "shiny", "amino_gwas_app", "data", "supplementary",
                      "SuppTable1_GWAS_annotation_soilN_amino_with_GO.csv")
SUPP_DIR <- file.path(ROOT, "final", "supp_tables")
MAIN_DIR <- file.path(ROOT, "final", "main_tables")
dir.create(SUPP_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(MAIN_DIR, showWarnings = FALSE, recursive = TRUE)

M          <- 4177796
ALPHA      <- 0.05
P_BONF     <- ALPHA / M
WINDOW_BP  <- 25000
LOCUS_GAP  <- 100000     # merge significant SNPs closer than this into one locus
N_MAIN     <- 5          # rows in the main-text table

load1 <- function(p) { e <- new.env(parent = emptyenv()); n <- load(p, envir = e); get(n[1], envir = e) }

## ---- gather every per-model hit ---------------------------------------------

fs <- list.files(BY_MODEL, pattern = "\\.RData$", full.names = TRUE)
if (!length(fs)) stop("no .RData files under ", BY_MODEL)

meta <- rbindlist(lapply(fs, function(p) {
  stem <- sub("\\.RData$", "", basename(p)); mdl <- sub("^.*_", "", stem)
  data.table(Phenotype = sub(paste0("_", mdl, "$"), "", stem), Model = mdl, path = p)
}))

hits <- rbindlist(lapply(seq_len(nrow(meta)), function(i) {
  d <- load1(meta$path[i]); if (!nrow(d)) return(NULL)
  data.table(Phenotype = meta$Phenotype[i], Model = meta$Model[i],
             SNP = d$SNP, Chr = d$Chr, Pos = d$Pos, P = d$P.value)
}), fill = TRUE)

cat(sprintf("phenotypes with data : %d\n", uniqueN(meta$Phenotype)))
cat(sprintf("threshold            : P <= %.4e  (-log10 %.4f)\n", P_BONF, -log10(P_BONF)))

## one row per SNP x phenotype, carrying how many models called it
w <- hits[P <= P_BONF][, .(n_models_sig = uniqueN(Model),
                           models_sig   = paste(sort(unique(Model)), collapse = "+"),
                           min_P        = min(P)),
                       by = .(Phenotype, SNP, Chr, Pos)]
w[, replicated := n_models_sig >= 2]

## ---- genes within +/- 25 kb --------------------------------------------------

gff <- fread(cmd = sprintf("awk -F'\\t' '$3==\"gene\"' %s", shQuote(GFF)),
             header = FALSE, sep = "\t", showProgress = FALSE)
genes <- gff[, {
  id <- sub(".*ID=([^;]+).*", "\\1", V9)
  .(Chr = sub("^[Cc]hr", "", V1), Start = V4, End = V5, GeneID = id)
}]
genes <- genes[Chr %in% as.character(1:10)][, Chr := as.integer(Chr)]
setkey(genes, Chr, Start, End)

q <- w[, .(Chr = Chr, s = pmax(1, Pos - WINDOW_BP), e = Pos + WINDOW_BP, idx = .I)]
ov <- foverlaps(q, genes, by.x = c("Chr", "s", "e"), by.y = c("Chr", "Start", "End"), nomatch = 0L)
ann <- unique(cbind(w[ov$idx], ov[, .(GeneID, GeneStart = Start, GeneEnd = End)]))
ann[, Relation := fifelse(Pos >= GeneStart & Pos <= GeneEnd, "within",
                   fifelse(Pos <  GeneStart, "upstream", "downstream"))]
ann[, Distance_to_Gene_bp := fifelse(Relation == "within", 0,
      fifelse(Relation == "upstream", GeneStart - Pos, Pos - GeneEnd))]

## ---- headline numbers for Results paragraph 3 --------------------------------

say <- function(sn, a, lab) {
  pg <- a[, .(np = uniqueN(Phenotype)), by = GeneID]
  pt <- a[, .(ng = uniqueN(GeneID)), by = Phenotype][order(-ng)]
  cat(sprintf(paste0("\n%s\n  SNP x phenotype pairs %5d | distinct SNPs %5d | phenotypes %3d\n",
                     "  SNPs within 25 kb of a gene %4d | candidate genes %5d",
                     " (single %5d, multi %4d) in %d phenotypes\n  top: %s\n  Total_N %d | Total_PBAA %d\n"),
    lab, nrow(sn), uniqueN(sn$SNP), uniqueN(sn$Phenotype), uniqueN(a$SNP),
    nrow(pg), pg[np == 1, .N], pg[np > 1, .N], uniqueN(a$Phenotype),
    paste(sprintf("%s=%d", head(pt, 6)$Phenotype, head(pt, 6)$ng), collapse = ", "),
    {v <- pt[Phenotype == "Total_N", ng];    if (length(v)) v else 0L},
    {v <- pt[Phenotype == "Total_PBAA", ng]; if (length(v)) v else 0L}))
}
say(w, ann, "UNION (significant in >=1 of MLM / BLINK / FarmCPU)")
say(w[replicated == TRUE], ann[replicated == TRUE], "REPLICATED (>=2 of the 3 models)")

setorder(ann, min_P)
fwrite(ann, file.path(SUPP_DIR, "SuppTable2_amino_gwas_candidate_genes.csv"))
fwrite(ann[replicated == TRUE], file.path(SUPP_DIR, "SuppTable3_amino_gwas_replicated_genes.csv"))

per_trait <- merge(ann[, .(n_genes_union = uniqueN(GeneID)), by = Phenotype],
                   ann[replicated == TRUE, .(n_genes_replicated = uniqueN(GeneID)), by = Phenotype],
                   by = "Phenotype", all.x = TRUE)
per_trait[is.na(n_genes_replicated), n_genes_replicated := 0L]
setorder(per_trait, -n_genes_union)
fwrite(per_trait, file.path(SUPP_DIR, "SuppTable2b_amino_gwas_genes_per_trait.csv"))

## ---- loci --------------------------------------------------------------------

snps <- unique(w[, .(Chr, Pos, SNP)]); setorder(snps, Chr, Pos)
snps[, locus := cumsum(c(1L, as.integer(diff(Pos) > LOCUS_GAP | diff(Chr) != 0)))]
w <- merge(w, snps[, .(SNP, Chr, Pos, locus)], by = c("SNP", "Chr", "Pos"))

loc <- w[, .(Chr = Chr[1], start = min(Pos), end = max(Pos),
             n_SNPs = uniqueN(SNP), n_phenotypes = uniqueN(Phenotype),
             phenotypes = paste(sort(unique(Phenotype)), collapse = ", "),
             best_P = min(min_P), lead_SNP = SNP[which.min(min_P)],
             lead_Pos = Pos[which.min(min_P)], replicated = any(replicated),
             lead_models = models_sig[which.min(min_P)],
             models_any = paste(sort(unique(unlist(strsplit(models_sig, "\\+")))), collapse = "+")),
         by = locus]

go <- fread(GO_FILE, encoding = "UTF-8"); setnames(go, 1, "GeneID"); go <- unique(go, by = "GeneID")

lg <- unique(merge(w[, .(SNP, locus)], ann[, .(SNP, GeneID)], by = "SNP", allow.cartesian = TRUE)[, .(locus, GeneID)])
lg <- go[lg, on = "GeneID"]
lg[, informative := !is.na(Family_Subfamily) & Family_Subfamily != ""]
setorder(lg, locus, -informative)
loc_genes <- lg[, .(n_genes = uniqueN(GeneID), genes = paste(unique(GeneID), collapse = ";"),
                    top_gene = GeneID[1], family = Family_Subfamily[1],
                    pclass = Protein_Class[1], go_bp = GO_BO[1], go_mf = GO_MF[1]), by = locus]
loc <- merge(loc, loc_genes, by = "locus", all.x = TRUE)
setorder(loc, -n_phenotypes, best_P)
fwrite(loc, file.path(SUPP_DIR, "SuppTable4_amino_gwas_loci_ranked.csv"))
cat(sprintf("\nloci %d | with >=1 gene %d | replicated %d\n",
            nrow(loc), loc[!is.na(n_genes), .N], loc[replicated == TRUE, .N]))

## ---- main-text LaTeX table ---------------------------------------------------

tidy <- function(s) { s <- gsub("\\(PTHR[^)]*\\)", "", s); s <- gsub("\\(PC[0-9]+\\)", "", s)
  s <- gsub("_", "/", trimws(s)); s <- tolower(s); substr(s, 1, 1) <- toupper(substr(s, 1, 1)); s }

## The main table is drawn from ALL loci, not only "replicated" ones.
## Whether a locus can replicate is largely determined by which models ran for
## that phenotype (traits running both BLINK and FarmCPU replicate at 67%;
## traits missing one of them at 5%), so a replicated-only table selects on
## pipeline completeness. The Models column discloses the evidence instead.
top <- head(loc[!is.na(family) & family != ""], N_MAIN)

tex <- c("\\begin{table}[ht]", "\\centering", "\\small",
  paste0("\\caption{Top five loci from the Goodman--Buckler amino acid GWAS, ranked by the ",
         "number of phenotypes with a Bonferroni-significant association ",
         "($P \\le 1.20 \\times 10^{-8}$). Loci were defined by merging significant SNPs ",
         "within 100 kb. The Models column lists the models in which the lead SNP exceeded ",
         "the threshold; MLM returned results for 35 of 76 phenotypes and MLMM for none, so ",
         "the number of models available differs between phenotypes and cross-model agreement ",
         "is not used as a selection criterion. Candidate genes lie within 25 kb of a ",
         "significant SNP; where a window contained several genes, the gene with a PANTHER ",
         "family assignment is shown and the total is given in parentheses.}"),
  "\\label{tab:amino_gwas_top_loci}", "\\begin{tabular}{llrllll}", "\\hline",
  "Locus & Lead SNP & $n$ phen. & Best $P$ & Models & Candidate gene & Putative function \\\\", "\\hline")

for (i in seq_len(nrow(top))) {
  x  <- top[i]
  gl <- if (x$n_genes > 1) sprintf("%s (of %d)", x$top_gene, x$n_genes) else x$top_gene
  e  <- sprintf("%.1e", x$best_P)
  tex <- c(tex, sprintf("chr%d:%.1f Mb & %s & %d & $%s\\times10^{-%d}$ & %s & \\texttt{%s} & %s \\\\",
    x$Chr, x$lead_Pos / 1e6, x$lead_SNP, x$n_phenotypes,
    sub("e.*", "", e), abs(as.integer(sub(".*e", "", e))),
    gsub("\\+", ", ", x$lead_models), gl, tidy(x$family)))
}
tex <- c(tex, "\\hline", "\\end{tabular}", "\\end{table}")
writeLines(tex, file.path(MAIN_DIR, "MainTable1_amino_gwas_top_loci.tex"))

cat("\n", paste(tex, collapse = "\n"), "\n", sep = "")
cat("\nwrote tables to", SUPP_DIR, "and", MAIN_DIR, "\n")
