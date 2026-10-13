## Asp-pathway genes among LD-assigned GWAS candidates and enrichment vs the genome.
## narrow = Asp-family metabolism GO terms/names; broad = the Asp hub (+Asn, Lys branch, Met, PEPC, GH3, ASS, ATC, quinolinate).
## Inputs (copied from final/results): GO_all_GWASs_LD.txt, all_candidate_gene_sources.csv, amino LD outputs, soil-N LD outputs, amino_lead_traits.tsv
suppressPackageStartupMessages(library(data.table))
U <- unique(fread("/mnt/user-data/uploads/Grain_N_remobilization/final/results/GO/GO_universe_genomewide.txt", quote="\""), by="GeneID")
A <- fread("GO_all_GWASs_LD.txt", quote="\"")
txt <- function(d) d[, paste(GeneName, Family_Subfamily, GO_MF, GO_BP)]
excl <- "aspartic|isoaspartate|protein-aspartate|aspartyl protease|peptidase A1"
narrow_go <- c("GO:0004072","GO:0004069","GO:0004815","GO:0008734","GO:0004073","GO:0006531","GO:0034628","GO:0009088","GO:0009090","GO:0006422")
narrow_nm <- "aspartate|aspartokinase|homoserine|threonine synthase"
broad_nm  <- paste(narrow_nm, "asparagine synth|asparaginase|L-asparaginase|argininosuccinate synth|aspartate carbamoyltransferase|phosphoenolpyruvate carboxylase|quinolinate synth|prephenate aminotransferase|cystathionine|methionine synthase|homocysteine.*methyltransferase|dihydrodipicolinate|4-hydroxy-tetrahydrodipicolinate|diaminopimelate|indole-3-acetic acid-amido|GH3", sep="|")
flag <- function(d, nm) { t <- txt(d); (Reduce(`|`, lapply(narrow_go, function(g) grepl(g, t, fixed=TRUE))) | grepl(nm, t, ignore.case=TRUE)) & !grepl(excl, t, ignore.case=TRUE) }
U[, narrow := flag(U, narrow_nm)][, broad := flag(U, broad_nm)]
cat("genome: narrow", U[, sum(narrow)], " broad", U[, sum(broad)], " of", nrow(U), "\n")
S <- fread("all_candidate_gene_sources.csv")
sets <- c("aminoGWAS_LD","aminoGWAS_25kb","aminoGWAS_nearest","soilN_LD","soilN_25kb","IC_Fst","J_Fst")
res <- rbindlist(lapply(sets, function(s) {
  g <- intersect(S[get(s) == 1, GeneID], U$GeneID); n <- length(g)
  rbindlist(lapply(c("narrow","broad"), function(def) { K <- U[, sum(get(def))]; k <- U[GeneID %in% g, sum(get(def))]
    data.table(set = s, definition = def, n_genes = n, observed = k, expected = round(n * K / nrow(U), 2),
               P = signif(phyper(k - 1, K, nrow(U) - K, n, lower.tail = FALSE), 3)) }))
}))
print(res)
fwrite(res, "Asp_enrichment_LD_sets.csv")
## Asp-related genes assigned by LD, with their GWAS evidence
nm <- U[, .(GeneID, GeneName = sub(";.*", "", GeneName), Family = Family_Subfamily, narrow, broad)]
lg <- fread("amino_lead_genes.tsv")[in_ld == TRUE]; lv <- fread("amino_leads_v5.tsv"); tr <- fread("amino_lead_traits.tsv")
aa <- merge(merge(lg, nm[broad == TRUE], by = "GeneID"), lv[, .(lead, interval_kb, ld_found)], by = "lead")
aa <- merge(aa, tr, by.x = "lead", by.y = "SNP")
cat("\n== amino acid GWAS: Asp-pathway genes inside LD intervals ==\n")
print(aa[order(minP), .(GeneID, GeneName = substr(GeneName, 1, 40), narrow, lead, chr, pos_v5, interval_kb, n_traits, traits = substr(traits, 1, 60), models, minP = signif(minP, 3))])
sg <- fread("soilN_genes.tsv")[in_ld_block == TRUE]
ss <- merge(sg, nm[broad == TRUE], by = "GeneID")
cat("\n== soil N GWAS: Asp-pathway genes inside LD intervals ==\n"); print(ss)
fwrite(aa, "Asp_genes_amino_LD.csv"); fwrite(ss, "Asp_genes_soilN_LD.csv")
