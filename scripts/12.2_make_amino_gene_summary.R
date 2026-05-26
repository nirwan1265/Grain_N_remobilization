################################################################################
### AMINO-ACID GWAS GENE SUMMARY (+/- 25 KB)
################################################################################

library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(IRanges)

################################################################################
### CONFIGURATION
################################################################################

INPUT_DIR <- "/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/Sarah_N_grain/split_by_trait"
REFERENCE_GFF <- "/rsstu/users/r/rrellan/sara/ref/gff3/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3"

WINDOW_BP <- 25000
P_CUTOFF <- 1e-7
MODEL_ORDER <- c("MLM", "MLMM", "BLINK", "FarmCPU")

OUTPUT_DIR <- "gene_count"
OUTPUT_BEST <- file.path(OUTPUT_DIR, "SuppTable_amino_gwas_gene_best_by_phenotype_25kb.csv")
OUTPUT_SUMMARY <- file.path(OUTPUT_DIR, "SuppTable_amino_gwas_gene_summary_25kb.csv")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

################################################################################
### HELPERS
################################################################################

normalize_model <- function(x) {
  x <- trimws(as.character(x))
  x[toupper(x) == "MLM"] <- "MLM"
  x[toupper(x) == "MLMM"] <- "MLMM"
  x[toupper(x) == "BLINK"] <- "BLINK"
  x[toupper(x) == "FARMCPU"] <- "FarmCPU"
  x
}

get_gene_label <- function(genes_gr) {
  gene_id <- as.character(mcols(genes_gr)$ID)

  if ("Name" %in% colnames(mcols(genes_gr))) {
    gene_name <- as.character(mcols(genes_gr)$Name)
    gene_name[is.na(gene_name) | gene_name == ""] <- gene_id[is.na(gene_name) | gene_name == ""]
    return(gene_name)
  }

  gene_id
}

load_trait_hits <- function(trait_name, p_cutoff = P_CUTOFF) {
  file_path <- file.path(INPUT_DIR, paste0(trait_name, ".csv"))

  dt <- fread(
    file_path,
    select = c("SNP", "Chr", "Pos", "P.value", "model"),
    showProgress = FALSE
  )

  setDT(dt)
  dt[, Chr := as.integer(Chr)]
  dt[, Pos := as.numeric(Pos)]
  dt[, P.value := as.numeric(P.value)]
  dt[, model := normalize_model(model)]

  dt <- dt[
    !is.na(Chr) &
      !is.na(Pos) &
      !is.na(P.value) &
      is.finite(P.value) &
      P.value > 0 &
      P.value <= p_cutoff
  ]

  if (nrow(dt) == 0) {
    return(NULL)
  }

  dt[, Phenotype := trait_name]
  dt[, log10_P := -log10(P.value)]
  dt[]
}

annotate_trait_hits <- function(hit_dt, genes_gr, window_bp = WINDOW_BP) {
  if (is.null(hit_dt) || nrow(hit_dt) == 0) {
    return(NULL)
  }

  snps <- GRanges(
    seqnames = Rle(paste0("chr", hit_dt$Chr)),
    ranges = IRanges(hit_dt$Pos, hit_dt$Pos)
  )
  mcols(snps)$row_idx <- seq_len(nrow(hit_dt))

  extended <- snps
  start(extended) <- pmax(1L, start(snps) - window_bp)
  end(extended) <- end(snps) + window_bp

  overlaps <- findOverlaps(genes_gr, extended, ignore.strand = TRUE)
  if (length(overlaps) == 0) {
    return(NULL)
  }

  gene_idx <- queryHits(overlaps)
  snp_idx <- mcols(extended)$row_idx[subjectHits(overlaps)]
  gene_label <- get_gene_label(genes_gr)
  snp_pos <- hit_dt$Pos[snp_idx]
  gene_start <- start(genes_gr)[gene_idx]
  gene_end <- end(genes_gr)[gene_idx]

  relation <- ifelse(
    snp_pos >= gene_start & snp_pos <= gene_end,
    "within",
    ifelse(snp_pos < gene_start, "upstream", "downstream")
  )

  distance_bp <- ifelse(
    relation == "within",
    0L,
    ifelse(
      relation == "upstream",
      as.integer(gene_start - snp_pos),
      as.integer(snp_pos - gene_end)
    )
  )

  out <- data.table(
    GeneID = as.character(mcols(genes_gr)$ID[gene_idx]),
    GeneSymbol = gene_label[gene_idx],
    GeneChr = as.character(seqnames(genes_gr)[gene_idx]),
    GeneStart = as.integer(gene_start),
    GeneEnd = as.integer(gene_end),
    Phenotype = hit_dt$Phenotype[snp_idx],
    Model = hit_dt$model[snp_idx],
    SNP = hit_dt$SNP[snp_idx],
    Chr = hit_dt$Chr[snp_idx],
    SNP_Pos = as.integer(snp_pos),
    P.value = hit_dt$P.value[snp_idx],
    log10_P = hit_dt$log10_P[snp_idx],
    Relation = relation,
    Distance_to_Gene_bp = distance_bp
  )

  unique(out)
}

collapse_best_by_gene_phenotype <- function(annotation_dt) {
  if (is.null(annotation_dt) || nrow(annotation_dt) == 0) {
    return(NULL)
  }

  dt <- copy(annotation_dt)
  setorder(dt, GeneID, Phenotype, P.value, Distance_to_Gene_bp, SNP)
  dt[, .SD[1L], by = .(GeneID, GeneSymbol, GeneChr, GeneStart, GeneEnd, Phenotype)]
}

make_gene_summary <- function(best_dt, phenotype_levels) {
  if (is.null(best_dt) || nrow(best_dt) == 0) {
    return(data.table(
      GeneID = character(),
      GeneSymbol = character(),
      GeneChr = character(),
      GeneStart = integer(),
      GeneEnd = integer(),
      Phenotypes = character(),
      Phenotypes_count = integer(),
      Best_Pvalues = character(),
      Best_log10P = character(),
      Best_Models = character(),
      Best_SNPs = character(),
      Top_Pvalue = numeric(),
      Top_log10P = numeric()
    ))
  }

  dt <- copy(best_dt)
  dt[, Phenotype := factor(Phenotype, levels = phenotype_levels)]
  setorder(dt, GeneSymbol, Phenotype)

  summary_dt <- dt[, .(
    Phenotypes = paste(as.character(Phenotype), collapse = ";"),
    Phenotypes_count = uniqueN(as.character(Phenotype)),
    Best_Pvalues = paste(formatC(P.value, digits = 3, format = "e"), collapse = ";"),
    Best_log10P = paste(formatC(log10_P, digits = 2, format = "f"), collapse = ";"),
    Best_Models = paste(Model, collapse = ";"),
    Best_SNPs = paste(SNP, collapse = ";"),
    Top_Pvalue = min(P.value),
    Top_log10P = max(log10_P)
  ), by = .(GeneID, GeneSymbol, GeneChr, GeneStart, GeneEnd)]

  setorder(summary_dt, -Phenotypes_count, Top_Pvalue, GeneSymbol)
  summary_dt[]
}

run_all_gene_summaries <- function() {
  trait_names <- sort(sub("\\.csv$", "", list.files(INPUT_DIR, pattern = "\\.csv$")))

  cat("Traits found:", length(trait_names), "\n")
  if (length(trait_names) == 0) {
    stop("No GWAS CSV files found in INPUT_DIR")
  }

  cat("\n=== Loading reference ===\n")
  reference_gr <- rtracklayer::import(REFERENCE_GFF)
  genes_gr <- reference_gr[mcols(reference_gr)$type == "gene"]
  cat("Genes in reference:", length(genes_gr), "\n")

  annotation_list <- vector("list", length(trait_names))
  names(annotation_list) <- trait_names

  for (trait_name in trait_names) {
    cat("\n=== Trait:", trait_name, "===\n")
    hit_dt <- load_trait_hits(trait_name)

    if (is.null(hit_dt) || nrow(hit_dt) == 0) {
      cat("No SNPs passed P.value <=", format(P_CUTOFF, scientific = TRUE), "\n")
      next
    }

    present_models <- MODEL_ORDER[MODEL_ORDER %in% unique(hit_dt$model)]
    cat("Significant SNP rows:", nrow(hit_dt), "\n")
    cat("Models:", paste(present_models, collapse = ", "), "\n")

    annotation_dt <- annotate_trait_hits(hit_dt, genes_gr, window_bp = WINDOW_BP)

    if (is.null(annotation_dt) || nrow(annotation_dt) == 0) {
      cat("Annotated gene hits: 0\n")
      next
    }

    annotation_list[[trait_name]] <- annotation_dt
    cat("Annotated gene hits:", nrow(annotation_dt), "\n")
    cat("Unique genes:", uniqueN(annotation_dt$GeneID), "\n")
  }

  all_annotations <- rbindlist(annotation_list, use.names = TRUE, fill = TRUE)
  best_by_gene_phenotype <- collapse_best_by_gene_phenotype(all_annotations)
  gene_summary <- make_gene_summary(best_by_gene_phenotype, phenotype_levels = trait_names)

  fwrite(best_by_gene_phenotype, OUTPUT_BEST)
  fwrite(gene_summary, OUTPUT_SUMMARY)

  cat("\nSaved files:\n")
  cat("  ", OUTPUT_BEST, "\n", sep = "")
  cat("  ", OUTPUT_SUMMARY, "\n", sep = "")
  cat("Rows (best by phenotype):", nrow(best_by_gene_phenotype), "\n")
  cat("Rows (gene summary):", nrow(gene_summary), "\n")

  invisible(list(
    best_by_gene_phenotype = best_by_gene_phenotype,
    gene_summary = gene_summary
  ))
}

if (sys.nframe() == 0) {
  run_all_gene_summaries()
}
