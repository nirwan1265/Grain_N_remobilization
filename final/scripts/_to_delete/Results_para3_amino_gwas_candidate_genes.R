################################################################################
## Results paragraph 3 -- amino-acid GWAS candidate genes (Goodman-Buckler)
##
## Recomputes every number quoted in Results paragraph 3 from the annotated
## GWAS hit table, using an explicit Bonferroni threshold instead of the
## previous fixed -log10(p) >= 7 cutoff.
##
## Threshold:  alpha / m  with alpha = 0.05 and m = 4,117,796 SNPs
##             = 1.214 x 10^-8   (-log10 = 7.916)
##
## Input : shiny/amino_gwas_app/data/supplementary/
##           SuppTable_amino_gwas_gene_best_by_phenotype_25kb.csv
##         (one row per Gene x Phenotype, carrying the single best SNP across
##          the four GAPIT models; produced by scripts/12.2_make_amino_gene_summary.R)
## Output: final/supp_tables/SuppTable2_amino_gwas_candidate_genes_bonferroni.csv
##         final/supp_tables/SuppTable2b_amino_gwas_genes_per_trait.csv
################################################################################

library(data.table)

## Repo root: run from the repository root, or edit this path.
ROOT <- "."

IN_FILE <- file.path(ROOT, "shiny", "amino_gwas_app", "data", "supplementary",
                     "SuppTable_amino_gwas_gene_best_by_phenotype_25kb.csv")
OUT_DIR <- file.path(ROOT, "final", "supp_tables")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

N_SNP     <- 4117796
ALPHA     <- 0.05
P_BONF    <- ALPHA / N_SNP          # 1.2142e-08
LOG10_TH  <- -log10(P_BONF)         # 7.9157
P_ONE_FP  <- 1 / N_SNP              # 2.4285e-07 (one expected false positive)

cat(sprintf("Bonferroni threshold : P <= %.3e  (-log10 = %.3f)\n", P_BONF, LOG10_TH))
cat(sprintf("1/m  reference       : P <= %.3e  (-log10 = %.3f)\n", P_ONE_FP, -log10(P_ONE_FP)))

dt <- fread(IN_FILE, showProgress = FALSE)
dt[, log10_P := as.numeric(log10_P)]

sig <- dt[P.value <= P_BONF]

## ---- headline counts ---------------------------------------------------------
n_genes  <- uniqueN(sig$GeneID)
n_traits <- uniqueN(sig$Phenotype)

per_gene <- sig[, .(n_pheno = uniqueN(Phenotype)), by = GeneID]
n_single <- per_gene[n_pheno == 1, .N]
n_multi  <- per_gene[n_pheno  > 1, .N]

per_trait <- sig[, .(n_genes = uniqueN(GeneID)), by = Phenotype][order(-n_genes)]

cat(sprintf("\nUnique candidate genes          : %d\n", n_genes))
cat(sprintf("Traits with >=1 candidate gene   : %d (of %d phenotypes tested)\n",
            n_traits, uniqueN(dt$Phenotype)))
cat(sprintf("Genes in exactly one phenotype   : %d\n", n_single))
cat(sprintf("Genes in more than one phenotype : %d\n", n_multi))
cat("\nTop 10 traits by candidate-gene count:\n")
print(head(per_trait, 10))
cat("\nN proxy traits:\n")
print(per_trait[Phenotype %in% c("Total_N", "Total_PBAA")])

cat("\nModel recorded as best hit (NOT a model-overlap statistic -- the input\n",
    "table keeps only the single best SNP per Gene x Phenotype across models):\n", sep = "")
print(sig[, .N, by = Model][order(-N)])

## ---- outputs ----------------------------------------------------------------
setorder(sig, P.value)
fwrite(sig, file.path(OUT_DIR, "SuppTable2_amino_gwas_candidate_genes_bonferroni.csv"))
fwrite(per_trait, file.path(OUT_DIR, "SuppTable2b_amino_gwas_genes_per_trait.csv"))

cat("\nWrote:\n  ", file.path(OUT_DIR, "SuppTable2_amino_gwas_candidate_genes_bonferroni.csv"),
    "\n  ", file.path(OUT_DIR, "SuppTable2b_amino_gwas_genes_per_trait.csv"), "\n", sep = "")
