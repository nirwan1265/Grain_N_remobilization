#!/usr/bin/env Rscript
################################################################################
## Are Asp-related genes over-represented in the candidate gene sets?
##
## The Discussion argues that aspartate metabolism recurs across analyses. This
## asks how many Asp-related genes each candidate set would contain by chance.
##
## Asp-related genes (genome-wide, from the PANTHER annotation of all B73 v5
## genes): any of the GO terms
##   aspartate kinase (GO:0004072), L-aspartate:2-oxoglutarate transaminase
##   (GO:0004069), aspartate-tRNA ligase (GO:0004815), L-aspartate oxidase
##   (GO:0008734), aspartate-semialdehyde dehydrogenase (GO:0004073), aspartate
##   metabolic process (GO:0006531), de novo NAD+ from L-aspartate (GO:0034628),
##   L-threonine biosynthesis (GO:0009088), L-homoserine biosynthesis
##   (GO:0009090), aspartyl-tRNA aminoacylation (GO:0006422)
## or a gene / PANTHER family name containing aspartate, aspartokinase,
## homoserine or threonine synthase. Aspartic-type proteases and isoaspartate
## methyltransferases are excluded (they act on aspartyl residues in proteins,
## not on aspartate metabolism).
##
## Background: all annotated genes. Expected = n_set x K / N; P = hypergeometric
## P(X >= observed). Gene-level counts treat neighbouring genes as independent,
## which overstates support when two hits share one window or locus.
##
## Inputs : GO_universe_genomewide.txt, S6 (amino acid GWAS genes), S8 (soil-N
##          genes), S19 / S20 (Indian Chief / Jarvis top-5% F_ST genes)
## Output : final/results/IC_J_fst/Asp_related_genes_chance_test.csv
################################################################################
suppressPackageStartupMessages(library(data.table))
gv   <- function(k, d) { v <- Sys.getenv(k, ""); if (nzchar(v)) v else d }
REPO <- gv("GN_REPO", "/Users/nirwantandukar/Documents/Github/Grain_N_remobilization")
TAB  <- gv("GN_TABDIR", file.path(REPO, "final/supp_tables"))
GOU  <- gv("GN_GOU",  file.path(REPO, "final/results/GO/GO_universe_genomewide.txt"))
RES  <- gv("GN_RESDIR", file.path(REPO, "final/results/ic_j_fst"))
dir.create(RES, FALSE, TRUE)

U <- unique(fread(GOU, quote = "\""), by = "GeneID")
terms <- c("GO:0004072", "GO:0004069", "GO:0004815", "GO:0008734", "GO:0004073",
           "GO:0006531", "GO:0034628", "GO:0009088", "GO:0009090", "GO:0006422")
txt <- U[, paste(GeneName, Family_Subfamily, GO_MF, GO_BP)]
nm  <- U[, paste(GeneName, Family_Subfamily)]
by_go   <- Reduce(`|`, lapply(terms, function(t) grepl(t, txt, fixed = TRUE)))
by_name <- grepl("aspartate|aspartokinase|homoserine|threonine synthase", nm, ignore.case = TRUE) &
           !grepl("aspartic|isoaspartate|aspartyl protease", nm, ignore.case = TRUE)
ASP <- U[by_go | by_name, .(GeneID, GeneName = sub(";.*", "", GeneName))]
N <- nrow(U); K <- nrow(ASP)
cat("annotated genes:", N, " Asp-related:", K, "\n")

ids <- function(f, col = "GeneID") unique(fread(file.path(TAB, f))[[col]])
SETS <- list(
  "Amino acid GWAS (S6)"        = ids("SuppTable_S6_amino_gwas_candidate_genes.csv"),
  "Soil-N GWAS (S9)"            = ids("SuppTable_S9_soilN_candidate_genes.csv"),
  "Indian Chief top-5% F_ST (S20)" = ids("SuppTable_S20_IC_Fst_top5_genes.csv"),
  "Jarvis top-5% F_ST (S21)"    = ids("SuppTable_S21_J_Fst_top5_genes.csv"))
RESULT <- rbindlist(lapply(names(SETS), function(k) {
  s <- intersect(SETS[[k]], U$GeneID); hit <- ASP[GeneID %in% s]
  data.table(gene_set = k, n_genes = length(s), asp_related = nrow(hit),
             expected = round(length(s) * K / N, 2),
             P_hypergeometric = signif(phyper(nrow(hit) - 1, K, N - K, length(s), lower.tail = FALSE), 3),
             genes = paste(sprintf("%s (%s)", hit$GeneID, hit$GeneName), collapse = "; "))
}))
all_s <- intersect(unique(unlist(SETS)), U$GeneID); all_h <- ASP[GeneID %in% all_s]
RESULT <- rbind(RESULT, data.table(gene_set = "All four combined", n_genes = length(all_s),
  asp_related = nrow(all_h), expected = round(length(all_s) * K / N, 2),
  P_hypergeometric = signif(phyper(nrow(all_h) - 1, K, N - K, length(all_s), lower.tail = FALSE), 3),
  genes = ""))
print(RESULT[, -"genes"])
fwrite(RESULT, file.path(RES, "Asp_related_genes_chance_test.csv"))
fwrite(ASP, file.path(RES, "Asp_related_genes_genomewide.csv"))
