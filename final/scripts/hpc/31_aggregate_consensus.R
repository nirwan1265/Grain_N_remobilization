#!/usr/bin/env Rscript
################################################################################
## 31 -- POOL THE PER-TRAIT CONSENSUS SETS AND ANNOTATE TO GENES (+/- 25 kb)
##
## Run after the 30_* array job has finished for every phenotype.
## Reports everything at BOTH thresholds: suggestive 1/m and Bonferroni 0.05/m.
##
## Input :  <BY_MODEL_DIR>/consensus/<trait>_all_model_common_SNP.RData
##          <BY_MODEL_DIR>/summary/<trait>_summary.csv
##          REFERENCE_GFF
##
## Output:  <OUT_DIR>/consensus_model_overlap_by_phenotype.csv
##          <OUT_DIR>/consensus_snps_all_models.csv
##          <OUT_DIR>/consensus_genes_all_models_25kb_suggestive.csv
##          <OUT_DIR>/consensus_genes_all_models_25kb_bonferroni.csv
##          <OUT_DIR>/consensus_summary.txt      <- the numbers for Results para 3
##
## Usage
##   Rscript 31_aggregate_consensus.R [by_model_dir] [out_dir]
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
  library(IRanges)
})

BY_MODEL_DIR  <- "/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/Sarah_N_grain/by_model"
REFERENCE_GFF <- "/rsstu/users/r/rrellan/sara/ref/gff3/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3"
WINDOW_BP     <- 25000

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L) BY_MODEL_DIR <- args[1]
OUT_DIR <- if (length(args) >= 2L) args[2] else file.path(BY_MODEL_DIR, "results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

load_one <- function(path) {
  env <- new.env(parent = emptyenv())
  nm  <- load(path, envir = env)
  get(nm[1], envir = env)
}

## ---- pooled per-trait summaries ---------------------------------------------

summ_files <- list.files(file.path(BY_MODEL_DIR, "summary"),
                         pattern = "_summary\\.csv$", full.names = TRUE)
if (length(summ_files) == 0L) stop("no per-trait summaries under ", BY_MODEL_DIR)

overlap <- rbindlist(lapply(summ_files, fread), fill = TRUE)
setorder(overlap, -N_bonf_all_models, -N_sugg_all_models, Phenotype)
fwrite(overlap, file.path(OUT_DIR, "consensus_model_overlap_by_phenotype.csv"))

cat("phenotypes summarised:", nrow(overlap), "\n")
if (uniqueN(overlap$Models_run) > 1L) {
  cat("NOTE: not every phenotype ran the same set of models --\n")
  print(overlap[, .N, by = Models_run])
  cat("  'common to all models' therefore means a different thing per phenotype.\n")
}

## ---- pooled consensus SNPs ---------------------------------------------------

cons_files <- list.files(file.path(BY_MODEL_DIR, "consensus"),
                         pattern = "_all_model_common_SNP\\.RData$", full.names = TRUE)

consensus <- rbindlist(Filter(Negate(is.null), lapply(cons_files, function(f) {
  d <- load_one(f)
  if (!is.data.table(d) || nrow(d) == 0L) return(NULL)
  d <- copy(d)
  d[, Phenotype := sub("_all_model_common_SNP\\.RData$", "", basename(f))]
  d
})), fill = TRUE)

if (nrow(consensus)) {
  setcolorder(consensus, "Phenotype")
  setorder(consensus, max_P)
}
fwrite(consensus, file.path(OUT_DIR, "consensus_snps_all_models.csv"))

## ---- annotate ----------------------------------------------------------------

annotate <- function(cons) {
  if (is.null(cons) || nrow(cons) == 0L) return(NULL)

  gff      <- import(REFERENCE_GFF)
  genes_gr <- gff[gff$type == "gene"]

  gene_label <- as.character(mcols(genes_gr)$ID)
  if ("Name" %in% colnames(mcols(genes_gr))) {
    nm <- as.character(mcols(genes_gr)$Name)
    nm[is.na(nm) | nm == ""] <- gene_label[is.na(nm) | nm == ""]
    gene_label <- nm
  }

  snps <- GRanges(seqnames = Rle(paste0("chr", cons$Chr)),
                  ranges   = IRanges(cons$Pos, cons$Pos))
  ext  <- snps
  start(ext) <- pmax(1L, start(snps) - WINDOW_BP)
  end(ext)   <- end(snps) + WINDOW_BP

  hits <- findOverlaps(genes_gr, ext, ignore.strand = TRUE)
  if (length(hits) == 0L) return(NULL)
  gi <- queryHits(hits); si <- subjectHits(hits)

  gs <- start(genes_gr)[gi]; ge <- end(genes_gr)[gi]; sp <- cons$Pos[si]
  rel <- ifelse(sp >= gs & sp <= ge, "within",
                ifelse(sp < gs, "upstream", "downstream"))

  out <- unique(data.table(
    GeneID     = as.character(mcols(genes_gr)$ID[gi]),
    GeneSymbol = gene_label[gi],
    GeneChr    = as.character(seqnames(genes_gr)[gi]),
    GeneStart  = as.integer(gs),
    GeneEnd    = as.integer(ge),
    Phenotype  = cons$Phenotype[si],
    SNP        = cons$SNP[si],
    Chr        = cons$Chr[si],
    SNP_Pos    = as.integer(sp),
    P_MLM      = cons$P_MLM[si],
    P_MLMM     = cons$P_MLMM[si],
    P_BLINK    = cons$P_BLINK[si],
    P_FarmCPU  = cons$P_FarmCPU[si],
    min_P      = cons$min_P[si],
    max_P      = cons$max_P[si],
    Relation   = rel,
    Distance_to_Gene_bp = as.integer(
      ifelse(rel == "within", 0L, ifelse(rel == "upstream", gs - sp, sp - ge)))))

  setorder(out, max_P)
  out
}

ann_sugg <- annotate(consensus)
ann_bonf <- annotate(if (nrow(consensus)) consensus[all_models_bonferroni == TRUE] else consensus)

if (!is.null(ann_sugg))
  fwrite(ann_sugg, file.path(OUT_DIR, "consensus_genes_all_models_25kb_suggestive.csv"))
if (!is.null(ann_bonf))
  fwrite(ann_bonf, file.path(OUT_DIR, "consensus_genes_all_models_25kb_bonferroni.csv"))

## ---- the numbers for Results paragraph 3 -------------------------------------

block <- function(ann, cons, label) {
  pg <- if (!is.null(ann)) ann[, .(n_pheno = uniqueN(Phenotype)), by = GeneID] else NULL
  pt <- if (!is.null(ann)) ann[, .(n_genes = uniqueN(GeneID)), by = Phenotype][order(-n_genes)] else NULL
  out <- c(
    sprintf("--- %s ---", label),
    sprintf("  consensus SNP x phenotype pairs : %d", nrow(cons)),
    sprintf("  distinct consensus SNPs         : %d", if (nrow(cons)) uniqueN(cons$SNP) else 0L),
    sprintf("  phenotypes with a consensus SNP : %d", if (nrow(cons)) uniqueN(cons$Phenotype) else 0L),
    sprintf("  unique candidate genes          : %d", if (!is.null(pg)) nrow(pg) else 0L),
    sprintf("    single-phenotype genes        : %d", if (!is.null(pg)) pg[n_pheno == 1, .N] else 0L),
    sprintf("    multi-phenotype  genes        : %d", if (!is.null(pg)) pg[n_pheno  > 1, .N] else 0L))
  if (!is.null(pt) && nrow(pt)) {
    out <- c(out, "  top traits by candidate-gene count:",
             sprintf("    %-12s %d", head(pt, 10)$Phenotype, head(pt, 10)$n_genes))
    for (tr in c("Total_N", "Total_PBAA")) {
      v <- pt[Phenotype == tr, n_genes]
      out <- c(out, sprintf("    %-12s %d", tr, if (length(v)) v else 0L))
    }
  }
  out
}

lines <- c(
  sprintf("m (markers assumed)              : %d", overlap$M_snps_assumed[1]),
  sprintf("suggestive  1/m                  : %.4e  (-log10 %.4f)",
          overlap$P_suggestive[1], -log10(overlap$P_suggestive[1])),
  sprintf("Bonferroni  0.05/m               : %.4e  (-log10 %.4f)",
          overlap$P_bonferroni[1], -log10(overlap$P_bonferroni[1])),
  sprintf("phenotypes scanned               : %d", nrow(overlap)),
  "",
  "SNPs significant in >=k models, summed over phenotypes:",
  sprintf("  suggestive  >=1 %8d  >=2 %8d  >=3 %8d  all %8d",
          sum(overlap$N_sugg_union), sum(overlap$N_sugg_2plus),
          sum(overlap$N_sugg_3plus), sum(overlap$N_sugg_all_models)),
  sprintf("  Bonferroni  >=1 %8d  >=2 %8d  >=3 %8d  all %8d",
          sum(overlap$N_bonf_union), sum(overlap$N_bonf_2plus),
          sum(overlap$N_bonf_3plus), sum(overlap$N_bonf_all_models)),
  "",
  "significant SNPs per model, summed over phenotypes:",
  sprintf("  %-8s suggestive %8d   Bonferroni %8d", c("MLM", "MLMM", "BLINK", "FarmCPU"),
          c(sum(overlap$N_sugg_MLM,     na.rm = TRUE), sum(overlap$N_sugg_MLMM,    na.rm = TRUE),
            sum(overlap$N_sugg_BLINK,   na.rm = TRUE), sum(overlap$N_sugg_FarmCPU, na.rm = TRUE)),
          c(sum(overlap$N_bonf_MLM,     na.rm = TRUE), sum(overlap$N_bonf_MLMM,    na.rm = TRUE),
            sum(overlap$N_bonf_BLINK,   na.rm = TRUE), sum(overlap$N_bonf_FarmCPU, na.rm = TRUE))),
  "",
  block(ann_sugg, consensus, "COMMON TO ALL MODELS at suggestive 1/m"),
  "",
  block(ann_bonf, if (nrow(consensus)) consensus[all_models_bonferroni == TRUE] else consensus,
        "COMMON TO ALL MODELS at Bonferroni 0.05/m"))

writeLines(lines, file.path(OUT_DIR, "consensus_summary.txt"))
cat("\n", paste(lines, collapse = "\n"), "\n", sep = "")
cat("\nwrote results to ", OUT_DIR, "\n", sep = "")
