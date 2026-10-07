################################################################################
## FOUR-MODEL CONSENSUS SNPs AND CANDIDATE GENES -- AMINO-ACID GWAS
##
## For each phenotype, keeps only SNPs that reach the Bonferroni threshold in
## EVERY GAPIT model run (MLM, MLMM, BLINK, FarmCPU), then annotates those SNPs
## to genes within +/- 25 kb.  Also reports the full model-overlap profile
## (how many SNPs are significant in 1, 2, 3 or 4 models) so the consensus can
## be reported honestly even if it turns out to be small or empty.
##
## MUST BE RUN WHERE THE FULL PER-MODEL SCANS LIVE.  The collapsed table in
## shiny/amino_gwas_app/data/supplementary/ cannot be used: it retains only the
## single best SNP per Gene x Phenotype across models, so model overlap is not
## recoverable from it.
##
##   INPUT_DIR/<trait>.csv  with columns: SNP, Chr, Pos, P.value, model
##
## Threshold: alpha / m, alpha = 0.05, m = 4,117,796 SNPs  ->  1.214e-08
##
## Outputs (written to OUTPUT_DIR):
##   consensus_snp_model_overlap_by_phenotype.csv  -- n SNPs per model-count class
##   consensus_snps_all_models.csv                 -- the consensus SNPs
##   consensus_genes_all_models_25kb.csv           -- annotated consensus genes
##   consensus_summary.txt                         -- headline numbers
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
  library(IRanges)
})

################################################################################
### CONFIGURATION
################################################################################

INPUT_DIR     <- "/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/Sarah_N_grain/split_by_trait"
REFERENCE_GFF <- "/rsstu/users/r/rrellan/sara/ref/gff3/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3"
OUTPUT_DIR    <- "consensus_four_model"

WINDOW_BP <- 25000
N_SNP     <- 4117796
ALPHA     <- 0.05
P_THRESH  <- ALPHA / N_SNP            # 1.2142e-08,  -log10 = 7.9157
MODELS    <- c("MLM", "MLMM", "BLINK", "FarmCPU")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

################################################################################
### HELPERS
################################################################################

normalize_model <- function(x) {
  x <- toupper(trimws(as.character(x)))
  c(MLM = "MLM", MLMM = "MLMM", BLINK = "BLINK", FARMCPU = "FarmCPU")[x]
}

get_gene_label <- function(genes_gr) {
  gene_id <- as.character(mcols(genes_gr)$ID)
  if ("Name" %in% colnames(mcols(genes_gr))) {
    nm <- as.character(mcols(genes_gr)$Name)
    nm[is.na(nm) | nm == ""] <- gene_id[is.na(nm) | nm == ""]
    return(nm)
  }
  gene_id
}

load_trait <- function(trait) {
  dt <- fread(file.path(INPUT_DIR, paste0(trait, ".csv")),
              select = c("SNP", "Chr", "Pos", "P.value", "model"),
              showProgress = FALSE)
  dt[, Chr := as.integer(Chr)]
  dt[, Pos := as.numeric(Pos)]
  dt[, P.value := as.numeric(P.value)]
  dt[, model := normalize_model(model)]
  dt[!is.na(Chr) & !is.na(Pos) & !is.na(model) &
       is.finite(P.value) & P.value > 0]
}

################################################################################
### PASS 1 -- model overlap per phenotype
################################################################################

traits <- sort(sub("\\.csv$", "", list.files(INPUT_DIR, pattern = "\\.csv$")))
cat("traits found:", length(traits), "\n")

overlap_rows  <- vector("list", length(traits))
consensus_rows <- vector("list", length(traits))

for (i in seq_along(traits)) {
  tr <- traits[i]
  dt <- load_trait(tr)

  models_run <- intersect(MODELS, sort(unique(dt$model)))
  sig <- dt[P.value <= P_THRESH]

  if (nrow(sig) == 0L) {
    overlap_rows[[i]] <- data.table(
      Phenotype = tr,
      Models_run = paste(models_run, collapse = "+"),
      N_models_run = length(models_run),
      N_sig_union = 0L, N_sig_2plus = 0L, N_sig_3plus = 0L, N_sig_all = 0L)
    next
  }

  per_snp <- sig[, .(n_models = uniqueN(model),
                     models = paste(sort(unique(model)), collapse = "+"),
                     max_P = max(P.value)), by = .(SNP, Chr, Pos)]

  overlap_rows[[i]] <- data.table(
    Phenotype    = tr,
    Models_run   = paste(models_run, collapse = "+"),
    N_models_run = length(models_run),
    N_sig_union  = nrow(per_snp),
    N_sig_2plus  = per_snp[n_models >= 2, .N],
    N_sig_3plus  = per_snp[n_models >= 3, .N],
    N_sig_all    = per_snp[n_models >= length(models_run), .N])

  cons <- per_snp[n_models >= length(models_run)]
  if (nrow(cons) > 0L) {
    cons[, Phenotype := tr]
    cons[, Models_run := paste(models_run, collapse = "+")]
    consensus_rows[[i]] <- cons
  }

  cat(sprintf("  %-18s models=%-28s sig=%6d  >=2:%5d  >=3:%5d  all:%5d\n",
              tr, paste(models_run, collapse = "+"),
              nrow(per_snp), per_snp[n_models >= 2, .N],
              per_snp[n_models >= 3, .N],
              per_snp[n_models >= length(models_run), .N]))
}

overlap <- rbindlist(overlap_rows, fill = TRUE)
fwrite(overlap, file.path(OUTPUT_DIR, "consensus_snp_model_overlap_by_phenotype.csv"))

consensus <- rbindlist(consensus_rows, fill = TRUE)
fwrite(consensus, file.path(OUTPUT_DIR, "consensus_snps_all_models.csv"))

################################################################################
### PASS 2 -- annotate consensus SNPs to genes within +/- 25 kb
################################################################################

ann <- NULL
if (nrow(consensus) > 0L) {
  gff <- import(REFERENCE_GFF)
  genes_gr <- gff[gff$type == "gene"]
  gene_label <- get_gene_label(genes_gr)

  snps <- GRanges(seqnames = Rle(paste0("chr", consensus$Chr)),
                  ranges   = IRanges(consensus$Pos, consensus$Pos))
  ext <- snps
  start(ext) <- pmax(1L, start(snps) - WINDOW_BP)
  end(ext)   <- end(snps) + WINDOW_BP

  hits <- findOverlaps(genes_gr, ext, ignore.strand = TRUE)
  gi <- queryHits(hits); si <- subjectHits(hits)

  gs <- start(genes_gr)[gi]; ge <- end(genes_gr)[gi]; sp <- consensus$Pos[si]
  rel <- ifelse(sp >= gs & sp <= ge, "within",
                ifelse(sp < gs, "upstream", "downstream"))

  ann <- unique(data.table(
    GeneID     = as.character(mcols(genes_gr)$ID[gi]),
    GeneSymbol = gene_label[gi],
    GeneChr    = as.character(seqnames(genes_gr)[gi]),
    GeneStart  = as.integer(gs),
    GeneEnd    = as.integer(ge),
    Phenotype  = consensus$Phenotype[si],
    SNP        = consensus$SNP[si],
    Chr        = consensus$Chr[si],
    SNP_Pos    = as.integer(sp),
    Max_P_across_models = consensus$max_P[si],
    Models     = consensus$models[si],
    Relation   = rel,
    Distance_to_Gene_bp = as.integer(ifelse(rel == "within", 0L,
                            ifelse(rel == "upstream", gs - sp, sp - ge)))))

  setorder(ann, Max_P_across_models)
  fwrite(ann, file.path(OUTPUT_DIR, "consensus_genes_all_models_25kb.csv"))
}

################################################################################
### SUMMARY
################################################################################

summary_lines <- c(
  sprintf("threshold                        : P <= %.4e (-log10 %.4f)", P_THRESH, -log10(P_THRESH)),
  sprintf("phenotypes scanned               : %d", length(traits)),
  sprintf("phenotypes with a consensus SNP  : %d", overlap[N_sig_all > 0, .N]),
  sprintf("consensus SNP x phenotype pairs  : %d", nrow(consensus)),
  sprintf("distinct consensus SNPs          : %d", if (nrow(consensus)) uniqueN(consensus$SNP) else 0L),
  sprintf("consensus candidate genes        : %d", if (!is.null(ann)) uniqueN(ann$GeneID) else 0L),
  sprintf("  single-phenotype genes         : %d",
          if (!is.null(ann)) ann[, uniqueN(Phenotype), by = GeneID][V1 == 1, .N] else 0L),
  sprintf("  multi-phenotype  genes         : %d",
          if (!is.null(ann)) ann[, uniqueN(Phenotype), by = GeneID][V1 > 1, .N] else 0L))

writeLines(summary_lines, file.path(OUTPUT_DIR, "consensus_summary.txt"))
cat("\n", paste(summary_lines, collapse = "\n"), "\n", sep = "")
