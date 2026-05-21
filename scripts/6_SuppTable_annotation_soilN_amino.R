################################################################################
### ANNOTATE GWAS: SOIL N + SINGLE AMINO ACIDS
################################################################################

library(tidyverse)
library(vroom)
library(rtracklayer)
library(GenomicRanges)

################################################################################
### CONFIGURATION
################################################################################

WINDOW_BP <- 25000
P_CUTOFF <- 1e-7

REFERENCE_GFF <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize/Maize.annotation/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3"

SOILN_MODELS <- c(
  MLM = "/Users/nirwantandukar/Documents/Research/results/GWAS/GAPIT/raw_GWAS_MLM_3PC_N.txt",
  MLMM = "/Users/nirwantandukar/Documents/Research/results/GWAS/GAPIT/raw_GWAS_MLMM_3PC_N.txt",
  BLINK = "/Users/nirwantandukar/Documents/Research/results/GWAS/GAPIT/raw_GWAS_BLINK_3PC_N.txt"
)
SOILN_FARMCPU_DIR <- "/Users/nirwantandukar/Documents/Research/results/GWAS/GAPIT/FarmCPU"

AMINO_DIR <- "/Users/nirwantandukar/Documents/Research/results/GWAS/Sarah_amino_acid/N_grain/Phenotypes_GWAS_Grain"
AMINO_TRAITS <- c("D", "E", "N", "P", "Q", "Total_N", "Total_PBAA")
PHENOTYPE_ORDER <- c("SoilN", AMINO_TRAITS)

OUTPUT_DIR <- "tables/supplementary"
OUTPUT_FILE <- file.path(OUTPUT_DIR, "SuppTable_GWAS_annotation_soilN_amino.csv")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

################################################################################
### HELPERS
################################################################################

get_gene_label <- function(genes_gr) {
  gene_id <- as.character(mcols(genes_gr)$ID)

  if ("Name" %in% colnames(mcols(genes_gr))) {
    gene_name <- as.character(mcols(genes_gr)$Name)
    ifelse(is.na(gene_name) | gene_name == "", gene_id, gene_name)
  } else {
    gene_id
  }
}

annotate_hits <- function(df, phenotype_name, model_name, genes_gr, window_bp = 25000) {
  if (is.null(df) || nrow(df) == 0) {
    return(NULL)
  }

  snps <- GRanges(
    seqnames = Rle(paste0("chr", df$Chr)),
    ranges = IRanges(df$Pos, df$Pos)
  )
  mcols(snps)$SNP <- df$SNP
  mcols(snps)$P.value <- df$P.value
  mcols(snps)$SNP_Pos <- df$Pos

  extended <- snps
  start(extended) <- pmax(1L, start(snps) - window_bp)
  end(extended) <- end(snps) + window_bp

  overlaps <- findOverlaps(genes_gr, extended, ignore.strand = TRUE)
  if (length(overlaps) == 0) {
    return(NULL)
  }

  gene_idx <- queryHits(overlaps)
  snp_idx <- subjectHits(overlaps)
  gene_label <- get_gene_label(genes_gr)
  snp_pos <- mcols(extended)$SNP_Pos[snp_idx]

  tibble(
    Phenotype = phenotype_name,
    Model = model_name,
    GeneID = as.character(mcols(genes_gr)$ID[gene_idx]),
    GeneName = gene_label[gene_idx],
    SNP = mcols(extended)$SNP[snp_idx],
    SNP_Pos = snp_pos,
    P.value = mcols(extended)$P.value[snp_idx],
    log10_P = -log10(mcols(extended)$P.value[snp_idx]),
    Relation = case_when(
      snp_pos >= start(genes_gr)[gene_idx] & snp_pos <= end(genes_gr)[gene_idx] ~ "within",
      snp_pos < start(genes_gr)[gene_idx] ~ "upstream",
      TRUE ~ "downstream"
    ),
    Distance_to_Gene = case_when(
      snp_pos >= start(genes_gr)[gene_idx] & snp_pos <= end(genes_gr)[gene_idx] ~ 0L,
      snp_pos < start(genes_gr)[gene_idx] ~ as.integer(start(genes_gr)[gene_idx] - snp_pos),
      TRUE ~ as.integer(snp_pos - end(genes_gr)[gene_idx])
    )
  )
}

collapse_best_per_phenotype <- function(annotation_df) {
  if (is.null(annotation_df) || nrow(annotation_df) == 0) {
    return(NULL)
  }

  annotation_df %>%
    dplyr::group_by(GeneID, GeneName, Phenotype) %>%
    dplyr::slice_min(P.value, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(log10_P = -log10(P.value)) %>%
    dplyr::select(GeneID, GeneName, Phenotype, log10_P)
}

make_summary_table <- function(best_hits_df) {
  best_hits_df %>%
    dplyr::mutate(Phenotype = factor(Phenotype, levels = PHENOTYPE_ORDER)) %>%
    dplyr::arrange(GeneName, Phenotype) %>%
    dplyr::group_by(GeneID, GeneName) %>%
    dplyr::summarise(
      Phenotypes = paste(as.character(Phenotype), collapse = ";"),
      Phenotypes_count = n(),
      pvalues = paste(formatC(log10_P, digits = 2, format = "f"), collapse = ";"),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(Phenotypes_count), GeneName) %>%
    dplyr::select(GeneName, Phenotypes, Phenotypes_count, pvalues)
}

################################################################################
### LOAD REFERENCE
################################################################################

cat("\n=== Loading reference ===\n")
reference_gr <- rtracklayer::import(REFERENCE_GFF)
genes_only <- reference_gr[mcols(reference_gr)$type == "gene"]
cat("Genes in reference:", length(genes_only), "\n")

################################################################################
### SOIL N
################################################################################

cat("\n=== Annotating SoilN ===\n")
soiln_annotations <- list()

for (model_name in names(SOILN_MODELS)) {
  model_path <- SOILN_MODELS[[model_name]]
  model_df <- vroom(model_path, delim = "\t", show_col_types = FALSE) %>%
    dplyr::select(SNP, Chr, Pos, P.value) %>%
    dplyr::filter(P.value <= P_CUTOFF)

  cat(model_name, "significant SNPs:", nrow(model_df), "\n")

  soiln_annotations[[model_name]] <- annotate_hits(
    df = model_df,
    phenotype_name = "SoilN",
    model_name = model_name,
    genes_gr = genes_only,
    window_bp = WINDOW_BP
  )
}

farmcpu_chunks <- list()
for (chr in 1:10) {
  farmcpu_path <- file.path(SOILN_FARMCPU_DIR, paste0("FarmCPU_TN_3PC_maize_chr", chr, ".rds"))
  farmcpu_obj <- readRDS(farmcpu_path)
  farmcpu_chunks[[chr]] <- farmcpu_obj$TN_maize$GWAS
}

farmcpu_df <- dplyr::bind_rows(farmcpu_chunks) %>%
  dplyr::select(SNP, Chromosome, Position, p.value) %>%
  dplyr::rename(Chr = Chromosome, Pos = Position, P.value = p.value) %>%
  dplyr::filter(P.value <= P_CUTOFF)

cat("FarmCPU significant SNPs:", nrow(farmcpu_df), "\n")

soiln_annotations[["FarmCPU"]] <- annotate_hits(
  df = farmcpu_df,
  phenotype_name = "SoilN",
  model_name = "FarmCPU",
  genes_gr = genes_only,
  window_bp = WINDOW_BP
)

soiln_best <- collapse_best_per_phenotype(dplyr::bind_rows(soiln_annotations))
cat("SoilN annotated genes:", ifelse(is.null(soiln_best), 0, nrow(soiln_best)), "\n")

################################################################################
### SINGLE AMINO ACIDS
################################################################################

cat("\n=== Annotating amino-acid traits ===\n")
amino_best_list <- list()

for (trait_name in AMINO_TRAITS) {
  trait_path <- file.path(AMINO_DIR, paste0(trait_name, ".csv"))
  trait_df <- vroom(trait_path, show_col_types = FALSE) %>%
    dplyr::filter(P.value <= P_CUTOFF) %>%
    dplyr::select(SNP, Chr, Pos, P.value, model)

  cat("\nTrait:", trait_name, "\n")
  cat("Significant SNP rows:", nrow(trait_df), "\n")

  trait_annotations <- list()
  for (model_name in unique(trait_df$model)) {
    model_df <- trait_df %>%
      dplyr::filter(model == model_name) %>%
      dplyr::select(SNP, Chr, Pos, P.value)

    trait_annotations[[model_name]] <- annotate_hits(
      df = model_df,
      phenotype_name = trait_name,
      model_name = model_name,
      genes_gr = genes_only,
      window_bp = WINDOW_BP
    )
  }

  amino_best_list[[trait_name]] <- collapse_best_per_phenotype(dplyr::bind_rows(trait_annotations))
  cat("Annotated genes:", ifelse(is.null(amino_best_list[[trait_name]]), 0, nrow(amino_best_list[[trait_name]])), "\n")
}

################################################################################
### FINAL SUPPLEMENTARY TABLE
################################################################################

cat("\n=== Building final table ===\n")
combined_best <- dplyr::bind_rows(soiln_best, dplyr::bind_rows(amino_best_list))

supp_table <- make_summary_table(combined_best)
write_csv(supp_table, OUTPUT_FILE)

cat("Saved CSV:\n")
cat("  ", OUTPUT_FILE, "\n", sep = "")
cat("Rows:", nrow(supp_table), "\n")
cat("\nPreview:\n")
print(head(supp_table, 10))
