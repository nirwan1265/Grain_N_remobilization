#!/usr/bin/env Rscript
################################################################################
## LD-AWARE GO OVER-REPRESENTATION FOR THE AMINO ACID GWAS
##
## Port of scripts/new_new_script/13_GO_enrichment_by_class.R from the SAP
## lipidomics pipeline. The statistical machinery -- interval collapsing,
## density-matched permutation null, two-stage permutation, BH by ontology --
## is kept as written. What changed is the biology layer: lipid classes become
## amino acid sets, and Condition x Layer x Class becomes tier x set.
##
## WHY A PERMUTATION NULL AND NOT A HYPERGEOMETRIC
##   A hypergeometric test treats every candidate gene as an independent
##   observation. After window- or LD-based mapping that is false: one
##   significant SNP recruits every gene in its window, and if that window holds
##   a tandem array of related genes they all carry the same annotation. The
##   test then reports a hugely enriched term when the evidence is one
##   association that landed on a cluster. In this dataset the effect is
##   measurable -- 192 asparagine candidate genes collapse to 101 loci -- so the
##   correction is not hypothetical.
##
##   The null therefore resamples LOCI, not genes. For the set being tested it
##   records which 250 kb intervals the real candidates came from and how many
##   genes each contributed, then every draw replaces each observed interval
##   with a random one FROM THE SAME GENE-DENSITY DECILE and pulls the same
##   number of genes. The fake set has the same number of loci, the same gene
##   count, and the same physical clustering; only the genomic location is
##   randomised. A term carried by one gene-rich window fails immediately,
##   because a single random window reproduces it just as easily.
##
##   Genes with no GO annotation stay in the draw and contribute nothing,
##   matching the observed set. Resampling them away would inflate enrichment.
##   For the same reason the draw is SIZED from all candidates in a locus, not
##   only the annotated ones -- see the note at the loop below.
##
## CALIBRATION
##   Feeding the test random loci matched to an observed set in number, gene
##   count and density decile should pass ~5% of terms at q<0.05. Measured over
##   12 null replicates: MIN_GENES=2 gives 2.8%, MIN_GENES=3 gives 5.5%,
##   MIN_GENES=4 gives 4.9%. MIN_GENES=3 is the default for that reason.
##
## BACKGROUND -- THE THING THAT BREAKS THIS ANALYSIS
##   GO_UNIVERSE must be a PANTHER export for the WHOLE GENOME. Two previous
##   attempts used a PANTHER export of the candidate list itself, which makes
##   the random intervals carry almost no annotation and every observed term
##   look enriched. The guard below refuses to call such a run a result; it
##   labels the output diagnostic-only. Do not remove it.
##
##   To produce the universe: paste final/results/GO/all_maize_v5_gene_ids.txt
##   (39,035 genes) into PANTHER and export the same seven columns as
##   GO_all_amino_acid.txt.
##
## WINDOW vs LD
##   Candidates here were defined by a +/-25 kb window, not by LD r^2 as in the
##   SAP pipeline. The null is agnostic to how candidates were recruited -- it
##   only needs the lead SNP position per candidate, which is the Pos column --
##   so this runs correctly on window-based candidates. Replacing the window
##   with PLINK --r2 partners would change the candidate set, not the test.
##
## Inputs  (env vars, all with defaults)
##   GO_CAND     final/supp_tables/SuppTable_S5_amino_gwas_candidate_genes.csv
##   GO_RANGE    final/results/GO/genes.range          chr start end GeneID
##   GO_ANNOT    final/results/GO/GO_all_amino_acid.txt   PANTHER, candidates
##   GO_UNIVERSE <genome-wide PANTHER export>          REQUIRED for a real result
##   GO_OUT      final/results/GO
##   N_PERM 2000 | N_PERM2 100000 | GO_ONTOLOGIES BP,MF
################################################################################

suppressPackageStartupMessages({ library(data.table) })

CAND     <- Sys.getenv("GO_CAND", "final/supp_tables/SuppTable_S5_amino_gwas_candidate_genes.csv")
RANGE    <- Sys.getenv("GO_RANGE", "final/results/GO/genes.range")
ANNOT    <- Sys.getenv("GO_ANNOT", "final/results/GO/GO_all_amino_acid.txt")
UNIVERSE <- Sys.getenv("GO_UNIVERSE", "")
OUT_DIR  <- Sys.getenv("GO_OUT", "final/results/GO")

N_PERM      <- as.integer(Sys.getenv("N_PERM",  "2000"))
N_PERM2     <- as.integer(Sys.getenv("N_PERM2", "100000"))
SCREEN_HITS <- 20L
INTERVAL    <- 250000
MIN_GENES   <- as.integer(Sys.getenv("MIN_GENES", "3"))
set.seed(20261008)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## ---- amino acid sets, replacing the lipid CLASSES vector ---------------------
AA <- c(A="Ala", R="Arg", N="Asn", D="Asp", Q="Gln", E="Glu", G="Gly", H="His",
        I="Ile", L="Leu", K="Lys", M="Met", F="Phe", P="Pro", S="Ser", W="Trp",
        T="Thr", Y="Tyr", V="Val", C="Cys")
FAMILY <- list(Aspartate_family = c("I","M","N","T","D","K"),
               Glutamate_family = c("E","H","P","R","Q"),
               Serine_family    = c("S","G","C"),
               Pyruvate_family  = c("L","A","V"),
               Aromatic_family  = c("W","F","Y"),
               ## Branched_chain overlaps Pyruvate (Leu, Val) and Aspartate (Ile)
               ## by construction. It is kept because IVL is a measured composite
               ## trait in its own right, but its hits are not independent of
               ## those two families and should not be counted as separate
               ## evidence.
               Branched_chain   = c("I","V","L"))
FAMILY_TRAIT <- c(Aspartate_family="IMNTDK", Glutamate_family="EHPRQ",
                  Serine_family="SGC", Pyruvate_family="LAV", Aromatic_family="WFY",
                  Branched_chain="IVL")
N_RICH <- c("N","Q","R","K","H","W")

## A phenotype's subject is its first token: N.E is asparagine relative to
## glutamate and belongs to asparagine. This replaces trait_class()/norm().
head_tok <- function(p) ifelse(p %in% c("Total","Total_N","Total_PBAA"), NA_character_,
                               sub("\\..*$", "", p))

## ---- 1. candidates ----------------------------------------------------------
a <- fread(CAND)
stopifnot(all(c("Phenotype","Chr","Pos","GeneID") %in% names(a)))
a[, tok := head_tok(Phenotype)]

assign_sets <- function(dt) {
  out <- list()
  for (x in names(AA))
    out[[length(out)+1L]] <- dt[tok == x, .(GeneID, Chr, Pos,
                                   Set = unname(AA[x]), Tier = "amino_acid")]
  for (f in names(FAMILY))
    out[[length(out)+1L]] <- dt[tok %in% c(FAMILY[[f]], FAMILY_TRAIT[[f]]),
                                .(GeneID, Chr, Pos, Set = f, Tier = "biosynthetic_family")]
  out[[length(out)+1L]] <- dt[tok %in% N_RICH,
                              .(GeneID, Chr, Pos, Set = "N_rich", Tier = "nitrogen_content")]
  out[[length(out)+1L]] <- dt[tok %in% setdiff(names(AA), N_RICH),
                              .(GeneID, Chr, Pos, Set = "N_poor", Tier = "nitrogen_content")]
  for (t in c("Total","Total_N","Total_PBAA"))
    out[[length(out)+1L]] <- dt[Phenotype == t, .(GeneID, Chr, Pos, Set = t, Tier = "pool_trait")]
  unique(rbindlist(out))
}
x <- assign_sets(a)
## Interval comes from the LEAD SNP, so every gene recruited by one association
## shares one interval id. That is the collapsing the null depends on.
x[, Interval := paste(Chr, Pos %/% INTERVAL, sep = ":")]
message("candidate rows: ", nrow(x), " across ", uniqueN(x$Set), " sets")

## ---- GO annotation ----------------------------------------------------------
ONTOLOGIES <- strsplit(Sys.getenv("GO_ONTOLOGIES", "BP,MF"), ",")[[1]]
go_col <- c(BP = "GO_BP", MF = "GO_MF", CC = "GO_CC")
stopifnot(all(ONTOLOGIES %in% names(go_col)))

parse_go <- function(dt, col, onto) {
  if (!col %in% names(dt)) return(NULL)
  d <- dt[nzchar(get(col)), .(GeneID, cell = get(col))]
  if (!nrow(d)) return(NULL)
  y <- d[, .(item = trimws(unlist(strsplit(cell, ";", fixed = TRUE)))), by = GeneID]
  y[, GO   := sub("^.*\\((GO:[0-9]+)\\)$", "\\1", item)]
  y[, Term := sub("\\s*\\(GO:[0-9]+\\)$", "", item)]
  unique(y[grepl("^GO:[0-9]+$", GO), .(GeneID, GO, Term, Ontology = onto)])
}

go   <- fread(ANNOT, quote = "", sep = "\t")
gmap <- rbindlist(lapply(ONTOLOGIES, function(o) parse_go(go, go_col[[o]], o)), use.names = TRUE)
message("candidate annotation: ", uniqueN(gmap$GeneID), " genes, ", uniqueN(gmap$GO), " terms")

if (nzchar(UNIVERSE)) {
  u <- fread(UNIVERSE, quote = "", sep = "\t")
  umap <- rbindlist(lapply(ONTOLOGIES, function(o) parse_go(u, go_col[[o]], o)), use.names = TRUE)
  message("genome universe: ", uniqueN(umap$GeneID), " genes, ", uniqueN(umap$GO), " terms")
  if (uniqueN(umap$GeneID) < 10000)
    warning("universe has < 10,000 genes - is it really genome-wide?")
  DIAGNOSTIC <- FALSE
} else {
  umap <- gmap; DIAGNOSTIC <- TRUE
  message("\n!! No GO_UNIVERSE supplied. Background falls back to the candidate pool.")
  message("!! Output is DIAGNOSTIC ONLY and must not be reported as a result.\n")
}

## ---- genome intervals and density deciles -----------------------------------
gr <- fread(RANGE, col.names = c("chr","start","end","GeneID"))
gr[, Interval := paste(chr, start %/% INTERVAL, sep = ":")]
gene_lvl <- sort(unique(gr$GeneID))
gr[, gidx := match(GeneID, gene_lvl)]
giv <- gr[, .(n_genes = .N, gidx = list(gidx)), by = Interval]
qs  <- unique(quantile(giv$n_genes, probs = seq(0, 1, 0.1)))
giv[, dens := if (length(qs) > 2) cut(n_genes, qs, include.lowest = TRUE, labels = FALSE) else 1L]
giv_by_dens <- split(seq_len(nrow(giv)), giv$dens)
message("genome intervals: ", nrow(giv), " of ", INTERVAL/1000, " kb, ",
        length(giv_by_dens), " density strata")

## ---- 3. enrichment ----------------------------------------------------------
bg_genes <- if (DIAGNOSTIC) intersect(unique(x$GeneID), gmap$GeneID) else unique(umap$GeneID)
bgm  <- umap[GeneID %in% bg_genes]
bg_n <- length(bg_genes)
gm   <- gmap[GeneID %in% intersect(unique(x$GeneID), gmap$GeneID)]
message("background = ", bg_n,
        if (DIAGNOSTIC) " annotated candidates (DIAGNOSTIC)" else " genes in the genome universe")

glob_key <- unique(paste(bgm$Ontology, bgm$GO)); n_glob <- length(glob_key)
bg2 <- bgm[GeneID %in% gene_lvl]
gt_glob <- vector("list", length(gene_lvl))
tmp <- split(match(paste(bg2$Ontology, bg2$GO), glob_key), match(bg2$GeneID, gene_lvl))
gt_glob[as.integer(names(tmp))] <- tmp

bgc <- bgm[, .(Background = uniqueN(GeneID)), by = .(Ontology, GO, Term)]
res <- list()

for (st in unique(x$Set)) {
  ## d_all keeps EVERY candidate, annotated or not; d keeps the annotated ones.
  ## The split matters: term counting uses annotated genes, but the null must
  ## draw the same TOTAL number of genes per locus, because only ~54% of the
  ## genome carries a GO term. Sizing the draw from annotated candidates only
  ## makes the null pull ~54% as many annotated genes as the observed set has,
  ## so every term looks enriched. That error inflated this analysis from 61 to
  ## 257 significant terms before it was caught by calibrating on null sets.
  d_all <- x[Set == st]
  d     <- d_all[GeneID %in% bg_genes]
  fg <- unique(d$GeneID)
  if (length(fg) < MIN_GENES) next
  tier <- d$Tier[1]

  obs <- merge(gm[GeneID %in% fg], unique(d[, .(GeneID, Interval)]),
               by = "GeneID", allow.cartesian = TRUE)[
         , .(Genes = uniqueN(GeneID), Intervals = uniqueN(Interval),
             Overlap_genes = paste(sort(unique(GeneID)), collapse = ";")),
         by = .(Ontology, GO, Term)][Genes >= MIN_GENES]
  if (!nrow(obs)) next
  obs <- merge(obs, bgc, by = c("Ontology","GO","Term"))
  obs[, `:=`(Set = st, Tier = tier,
             Fold = round((Genes / length(fg)) / (Background / bg_n), 2))]
  obs[, p_gene := phyper(Genes - 1, Background, bg_n - Background, length(fg), lower.tail = FALSE)]

  oc <- d_all[, .(n_cand = uniqueN(GeneID)), by = Interval]
  oc <- merge(oc, giv[, .(Interval, n_genes, dens)], by = "Interval", all.x = TRUE)
  oc[is.na(dens), dens := 1L]

  key <- paste(obs$Ontology, obs$GO)
  gi  <- match(key, glob_key); n_terms <- length(key); obs_n <- obs$Genes

  pools     <- lapply(oc$dens, function(dd) giv_by_dens[[as.character(dd)]])
  pool_flat <- unlist(pools, use.names = FALSE)
  pool_len  <- lengths(pools)
  pool_off  <- c(0L, cumsum(pool_len)[-length(pool_len)])
  ncand <- oc$n_cand; nr <- nrow(oc)

  draw <- function() {
    js <- pool_flat[pool_off + ceiling(runif(nr) * pool_len)]
    picked <- vector("list", nr)
    for (r in seq_len(nr)) {
      g <- giv$gidx[[ js[r] ]]
      k <- if (ncand[r] < length(g)) ncand[r] else length(g)
      picked[[r]] <- if (k == length(g)) g else g[sample.int(length(g), k)]
    }
    idx <- unlist(gt_glob[unique(unlist(picked, use.names = FALSE))], use.names = FALSE)
    if (is.null(idx)) return(integer(n_terms))
    tabulate(idx, nbins = n_glob)[gi]
  }

  hits <- integer(n_terms)
  for (i in seq_len(N_PERM)) hits <- hits + (draw() >= obs_n)
  deep <- which(hits < SCREEN_HITS)
  obs[, p_LD := (hits + 1) / (N_PERM + 1)][, N_perm := N_PERM]
  if (length(deep) && N_PERM2 > 0) {
    hits2 <- integer(n_terms)
    for (i in seq_len(N_PERM2)) hits2 <- hits2 + (draw() >= obs_n)
    obs$p_LD[deep]   <- (hits2[deep] + 1) / (N_PERM2 + 1)
    obs$N_perm[deep] <- N_PERM2
  }
  message(sprintf("  %-20s %4d genes, %3d loci, %4d terms tested",
                  st, length(fg), nrow(oc), nrow(obs)))
  res[[length(res) + 1L]] <- obs
}

out <- rbindlist(res, use.names = TRUE)
out[, q_gene := p.adjust(p_gene, "BH"), by = .(Tier, Ontology)]
out[, q_LD   := p.adjust(p_LD,   "BH"), by = .(Tier, Ontology)]
setcolorder(out, c("Tier","Set","Ontology","Term","GO","Genes","Background",
                   "Intervals","Fold","q_gene","q_LD","Overlap_genes"))
setorder(out, Tier, Set, Ontology, q_LD, -Fold)

fn <- if (DIAGNOSTIC) "GO_LDaware_DIAGNOSTIC_candidate_background.tsv" else "GO_LDaware_all.tsv"
fwrite(out[, .(Tier, Set, Ontology, Term, GO, Genes, Background, Intervals, Fold,
               q_gene = signif(q_gene, 3), q_LD = signif(q_LD, 3), N_perm, Overlap_genes)],
       file.path(OUT_DIR, fn), sep = "\t")

cat("\n-- tested terms --\n");  print(out[, .(terms = .N), by = .(Tier, Ontology)])
cat("\n-- passing q_LD < 0.05 --\n"); print(out[q_LD < 0.05, .(terms = .N), by = .(Tier, Ontology)])
cat("\n-- top by fold among those passing, >= 2 intervals --\n")
print(head(out[q_LD < 0.05 & Intervals >= 2][order(-Fold),
      .(Tier, Set, Term, Genes, Background, Intervals, Fold, q_LD)], 25), row.names = FALSE)
if (DIAGNOSTIC) cat("\n!! candidate-pool background, diagnostic only, NOT a result\n")
cat("\nSaved:", file.path(OUT_DIR, fn), "(", nrow(out), "rows )\n")
