#!/usr/bin/env Rscript
################################################################################
## GO OVER-REPRESENTATION FOR THE AMINO ACID GWAS, ONE GENE SET PER AMINO ACID
##
## Design follows the lipid-GWAS approach: instead of pooling every candidate
## gene into one list, each focal amino acid gets its own set built from every
## phenotype in which that amino acid is the subject -- for asparagine that is
## N, N.D, N.E, N.Q, N.IMNTDK and N.Total. Family composites (LAV, IVL, IMNTDK,
## EHPRQ, SGC, WFY), Total, and the two N proxies are kept as separate sets.
##
## WHY THE PREVIOUS RUN RETURNED NOTHING
##   The 2026-10-05 offline run used the pooled candidate list (8,075 annotated
##   genes) as the BACKGROUND. With sets covering 14-25% of that background the
##   test reduces to "is this slice functionally different from the pool it
##   dominates", and the answer was no for all twelve sets: 0 of 2,182 tests
##   survived correction. The method was sound; the background was wrong.
##   This script uses g:Profiler's genome background for Zea mays, which is the
##   comparison the question actually asks.
##
## THIS SCRIPT MUST RUN ON A MACHINE WITH NETWORK ACCESS.
##   Neither the cloud sandbox nor the Cowork device bridge can reach
##   biit.cs.ut.ee (both return 403 from the egress proxy). Run it in your own
##   R session on the Mac, where g:Profiler is reachable.
##
##   install.packages("gprofiler2")
##   Rscript GO_amino_acid_genesets.R <repo_root>
##
## NON-INDEPENDENCE
##   Genes within one locus are not independent draws -- a single significant
##   SNP can contribute several neighbours within 25 kb. The script therefore
##   runs twice: once on all candidate genes, and once keeping a single
##   representative gene per locus. Terms that survive only the first run are
##   driven by gene clustering rather than by independent signals, so treat the
##   collapsed run as the conservative answer.
##
## Input : final/supp_tables/SuppTable_S5_amino_gwas_candidate_genes.csv
##         final/supp_tables/SuppTable_S6_amino_gwas_loci.csv   (for collapsing)
## Output: final/results/GO/genesets_amino_acid.csv       the sets themselves
##         final/results/GO/GO_enrichment_allgenes.csv
##         final/results/GO/GO_enrichment_locuscollapsed.csv
##         final/results/GO/GO_enrichment_summary.txt
################################################################################

suppressPackageStartupMessages(library(data.table))

ROOT <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "."
TAB  <- file.path(ROOT, "final", "supp_tables")
OUT  <- file.path(ROOT, "final", "results", "GO")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

ORGANISM <- "zmays"
SOURCES  <- c("GO:BP", "GO:MF", "GO:CC")   # add "KEGG","REAC" if you want pathways

AA <- c(A="Ala", R="Arg", N="Asn", D="Asp", Q="Gln", E="Glu", G="Gly", H="His",
        I="Ile", L="Leu", K="Lys", M="Met", F="Phe", P="Pro", S="Ser", W="Trp",
        T="Thr", Y="Tyr", V="Val", C="Cys")

## Biosynthetic families. These are not an outside classification imposed on the
## data -- they are the panel's own composite traits: IMNTDK is the aspartate
## family, EHPRQ the glutamate family, SGC the serine family, LAV the pyruvate
## family, WFY the aromatic family. A family set pools every phenotype whose
## subject is a member amino acid together with the composite trait itself, which
## is the amino acid equivalent of pooling every PC species into one lipid class.
FAMILY <- list(
  "Aspartate_family"  = c("I","M","N","T","D","K"),
  "Glutamate_family"  = c("E","H","P","R","Q"),
  "Serine_family"     = c("S","G","C"),
  "Pyruvate_family"   = c("L","A","V"),
  "Aromatic_family"   = c("W","F","Y"),
  "Branched_chain"    = c("I","V","L"))
FAMILY_TRAIT <- c(Aspartate_family = "IMNTDK", Glutamate_family = "EHPRQ",
                  Serine_family = "SGC", Pyruvate_family = "LAV",
                  Aromatic_family = "WFY", Branched_chain = "IVL")

## Amino acids carrying more than one nitrogen atom. The N-rich / N-poor split is
## the grouping this paper is actually about, and the two sets are nearly
## disjoint in the candidate data.
N_RICH <- c("N","Q","R","K","H","W")

genes <- fread(file.path(TAB, "SuppTable_S5_amino_gwas_candidate_genes.csv"))

## focal amino acid of a phenotype: the first token is its subject, so "N.E" is
## asparagine relative to glutamate and belongs to asparagine, not glutamate.
focal <- function(p) {
  if (p %in% c("Total", "Total_N", "Total_PBAA")) return(NA_character_)
  sub("\\..*$", "", p)
}
genes[, head_tok := vapply(Phenotype, focal, character(1))]

build <- function(keep, label, tier)
  data.table(set = label, tier = tier, GeneID = unique(genes[keep, GeneID]))

sets_all <- rbindlist(c(
  ## tier 1 -- one set per amino acid
  lapply(names(AA), function(a)
    build(genes$head_tok %in% a, unname(AA[a]), "amino_acid")),
  ## tier 2 -- biosynthetic family (the lipid-class analogue)
  lapply(names(FAMILY), function(f)
    build(genes$head_tok %in% c(FAMILY[[f]], FAMILY_TRAIT[[f]]), f, "biosynthetic_family")),
  ## tier 3 -- nitrogen content
  list(build(genes$head_tok %in% N_RICH, "N_rich_amino_acids", "nitrogen_content"),
       build(genes$head_tok %in% setdiff(names(AA), N_RICH), "N_poor_amino_acids", "nitrogen_content")),
  ## tier 4 -- whole-pool and proxy traits
  lapply(c("Total", "Total_N", "Total_PBAA"), function(t)
    build(genes$Phenotype == t, t, "pool_trait"))
))
sets_all <- sets_all[!is.na(GeneID)]
sets_all <- sets_all[, .(set, GeneID)][, head(.SD, .N), by = set]   # drop tier for downstream
tiers <- unique(rbindlist(c(
  lapply(names(AA), function(a) data.table(set = unname(AA[a]), tier = "amino_acid")),
  lapply(names(FAMILY), function(f) data.table(set = f, tier = "biosynthetic_family")),
  list(data.table(set = c("N_rich_amino_acids","N_poor_amino_acids"), tier = "nitrogen_content")),
  list(data.table(set = c("Total","Total_N","Total_PBAA"), tier = "pool_trait")))))



## ---- locus-collapsed version --------------------------------------------------
loci <- fread(file.path(TAB, "SuppTable_S6_amino_gwas_loci.csv"))
gene2locus <- loci[!is.na(genes) & genes != "",
                   .(GeneID = unlist(strsplit(genes, ";"))), by = locus]
gc <- merge(sets_all, genes[, .(GeneID, min_P)], by = "GeneID", allow.cartesian = TRUE)
gc <- merge(gc, gene2locus, by = "GeneID", allow.cartesian = TRUE)
setorder(gc, set, locus, min_P)
sets_collapsed <- gc[, .SD[1L], by = .(set, locus)][, .(set, GeneID)]

fwrite(rbind(sets_all[, .(set, GeneID, version = "all_genes")],
             sets_collapsed[, .(set, GeneID, version = "locus_collapsed")]),
       file.path(OUT, "genesets_amino_acid.csv"))

m <- merge(sets_all[, .(n_all = .N), by = set],
           sets_collapsed[, .(n_collapsed = .N), by = set], by = "set", all = TRUE)
m <- merge(tiers, m, by = "set"); setorder(m, tier, -n_all)
cat("gene set sizes\n"); print(m)
fwrite(m, file.path(OUT, "geneset_sizes.csv"))

## ---- enrichment ---------------------------------------------------------------
if (!requireNamespace("gprofiler2", quietly = TRUE))
  stop("install.packages(\"gprofiler2\") and run this where biit.cs.ut.ee is reachable")
library(gprofiler2)

run_sets <- function(dt, label) {
  qs <- split(dt$GeneID, dt$set)
  qs <- qs[vapply(qs, length, 1L) >= 10]        # <10 genes has no power
  cat(sprintf("\n[%s] querying %d sets\n", label, length(qs)))
  res <- gost(query = qs, organism = ORGANISM, sources = SOURCES,
              correction_method = "g_SCS",      # g:Profiler's own, not BH
              domain_scope = "annotated",       # genome background
              user_threshold = 0.05, multi_query = FALSE, evcodes = TRUE)
  if (is.null(res)) { cat("  no enrichment returned\n"); return(NULL) }
  r <- as.data.table(res$result)
  r[, gene_set := query]
  r[, .(gene_set, source, term_id, term_name, term_size, query_size,
        intersection_size, precision, recall, p_value,
        intersection = intersection)][order(gene_set, p_value)]
}

res_all <- run_sets(sets_all,       "all genes")
res_col <- run_sets(sets_collapsed, "locus collapsed")

if (!is.null(res_all)) fwrite(res_all, file.path(OUT, "GO_enrichment_allgenes.csv"))
if (!is.null(res_col)) fwrite(res_col, file.path(OUT, "GO_enrichment_locuscollapsed.csv"))

## ---- summary ------------------------------------------------------------------
lines <- c(sprintf("organism: %s | sources: %s | background: annotated genome",
                   ORGANISM, paste(SOURCES, collapse = ", ")),
           sprintf("correction: g:SCS, threshold 0.05"), "")
if (!is.null(res_all)) {
  lines <- c(lines, "significant terms per set, ALL GENES:",
             capture.output(print(res_all[, .N, by = gene_set][order(-N)])), "")
  if (!is.null(res_col)) {
    keep <- merge(res_all[, .(gene_set, term_id)], res_col[, .(gene_set, term_id)],
                  by = c("gene_set", "term_id"))
    lines <- c(lines,
      sprintf("terms surviving the locus-collapsed run: %d of %d (%.0f%%)",
              nrow(keep), nrow(res_all), 100 * nrow(keep) / max(1, nrow(res_all))),
      "terms that do NOT survive are driven by gene clustering within a locus.", "")
    lines <- c(lines, "significant terms per set, LOCUS COLLAPSED:",
               capture.output(print(res_col[, .N, by = gene_set][order(-N)])))
  }
} else lines <- c(lines, "no significant terms in any set.")

writeLines(lines, file.path(OUT, "GO_enrichment_summary.txt"))
cat("\n", paste(lines, collapse = "\n"), "\n", sep = "")
cat("\nwrote results to ", OUT, "\n", sep = "")
