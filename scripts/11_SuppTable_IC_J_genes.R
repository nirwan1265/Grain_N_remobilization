################################################################################
### SUPPLEMENTARY TABLE: INDIAN CHIEF AND JARVIS FIGURE 4 WINDOWS WITH GENES
################################################################################

library(tidyverse)
library(GenomicRanges)
library(rtracklayer)

################################################################################
### CONFIGURATION
################################################################################

fst_dir_ic <- "/Users/nirwantandukar/Documents/Research/results/Indian_Jarvis/ANGSD_Fst/Indian_chief/sliding_window"
fst_dir_j <- "/Users/nirwantandukar/Documents/Research/results/Indian_Jarvis/ANGSD_Fst/Jarvis/sliding_window"
tajima_dir_ic <- "/Users/nirwantandukar/Documents/Research/results/Indian_Jarvis/ANGSD_TajimaD/Indian_Chief/sliding_window"
tajima_dir_j <- "/Users/nirwantandukar/Documents/Research/results/Indian_Jarvis/ANGSD_TajimaD/Jarvis/sliding_window"

gff_file <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize/Maize.annotation/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3"
output_dir <- "tables/supplementary"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

window_size_bp <- 250000
fst_percentile <- 0.95
tail_percentile <- 0.05

chr_map <- tibble(
  nc_id = c(
    "NC_050096.1", "NC_050097.1", "NC_050098.1", "NC_050099.1", "NC_050100.1",
    "NC_050101.1", "NC_050102.1", "NC_050103.1", "NC_050104.1", "NC_050105.1"
  ),
  chr = paste0("chr", 1:10),
  chr_num = 1:10
)

################################################################################
### HELPERS
################################################################################

read_fst_files <- function(fst_dir, population_name) {
  files <- list.files(fst_dir, pattern = "\\.txt$", full.names = TRUE)

  if (length(files) == 0) {
    stop("No FST files found in: ", fst_dir)
  }

  purrr::map_dfr(files, function(path) {
    readr::read_tsv(
      path,
      show_col_types = FALSE,
      col_names = c("region", "Chr", "WinCenter", "Nsites", "Fst"),
      skip = 1
    ) %>%
      left_join(chr_map, by = c("Chr" = "nc_id")) %>%
      filter(!is.na(chr)) %>%
      transmute(
        Population = population_name,
        chr,
        chr_num,
        WinCenter = as.numeric(WinCenter),
        Nsites_fst = as.numeric(Nsites),
        Fst = as.numeric(Fst)
      )
  }) %>%
    arrange(chr_num, WinCenter)
}

read_tajima_files <- function(tajima_dir, generation_pattern) {
  files <- list.files(tajima_dir, pattern = generation_pattern, full.names = TRUE)

  if (length(files) == 0) {
    stop("No Tajima's D files found in: ", tajima_dir, " for pattern ", generation_pattern)
  }

  purrr::map_dfr(files, function(path) {
    readr::read_tsv(
      path,
      comment = "#",
      show_col_types = FALSE,
      col_names = c(
        "index", "Chr", "WinCenter", "tW", "tP", "tF", "tH", "tL",
        "Tajima", "fuf", "fud", "fayh", "zeng", "nSites"
      )
    ) %>%
      left_join(chr_map, by = c("Chr" = "nc_id")) %>%
      filter(!is.na(chr)) %>%
      transmute(
        chr,
        chr_num,
        WinCenter = as.numeric(WinCenter),
        Tajima = as.numeric(Tajima),
        tP = as.numeric(tP),
        nSites = as.numeric(nSites)
      )
  }) %>%
    arrange(chr_num, WinCenter)
}

prepare_population_windows <- function(fst_dir, tajima_dir, gen0_pattern, gen14_pattern, population_name) {
  fst_df <- read_fst_files(fst_dir, population_name)

  tajima_gen0 <- read_tajima_files(tajima_dir, gen0_pattern) %>%
    transmute(
      chr,
      chr_num,
      WinCenter,
      Tajima_Gen0 = Tajima,
      pi_Gen0 = tP / nSites,
      nSites_Gen0 = nSites
    )

  tajima_gen14 <- read_tajima_files(tajima_dir, gen14_pattern) %>%
    transmute(
      chr,
      chr_num,
      WinCenter,
      Tajima_Gen14 = Tajima,
      pi_Gen14 = tP / nSites,
      nSites_Gen14 = nSites
    )

  merged_df <- fst_df %>%
    inner_join(tajima_gen0, by = c("chr", "chr_num", "WinCenter")) %>%
    inner_join(tajima_gen14, by = c("chr", "chr_num", "WinCenter")) %>%
    mutate(
      window_start = pmax(1, as.integer(WinCenter - window_size_bp / 2)),
      window_end = as.integer(WinCenter + window_size_bp / 2),
      deltaD = Tajima_Gen14 - Tajima_Gen0,
      pi_ratio = pi_Gen14 / pi_Gen0,
      log_pi_ratio = log2(pi_ratio)
    ) %>%
    filter(
      is.finite(Fst),
      is.finite(deltaD),
      is.finite(pi_ratio),
      is.finite(log_pi_ratio)
    ) %>%
    arrange(chr_num, WinCenter)

  fst_threshold <- quantile(merged_df$Fst, fst_percentile, na.rm = TRUE)
  deltaD_threshold <- quantile(merged_df$deltaD, tail_percentile, na.rm = TRUE)
  log_pi_threshold <- quantile(merged_df$log_pi_ratio, tail_percentile, na.rm = TRUE)

  merged_df %>%
    mutate(
      fst_outlier = Fst >= fst_threshold,
      deltaD_outlier = deltaD <= deltaD_threshold,
      pi_outlier = log_pi_ratio <= log_pi_threshold,
      all_three_signals = fst_outlier & deltaD_outlier & pi_outlier
    )
}

annotate_windows_with_genes <- function(df, genes_gr) {
  win_gr <- GRanges(
    seqnames = df$chr,
    ranges = IRanges(start = df$window_start, end = df$window_end),
    window_id = seq_len(nrow(df))
  )

  hits <- findOverlaps(win_gr, genes_gr, ignore.strand = TRUE)

  gene_summary <- if (length(hits) == 0) {
    tibble(
      window_id = integer(),
      gene_count = integer(),
      GeneIDs = character(),
      GeneNames = character(),
      GeneBiotypes = character()
    )
  } else {
    gene_meta <- mcols(genes_gr)

    tibble(
      window_id = queryHits(hits),
      GeneID = as.character(gene_meta$ID[subjectHits(hits)]),
      GeneName = dplyr::coalesce(
        as.character(gene_meta$Name[subjectHits(hits)]),
        as.character(gene_meta$ID[subjectHits(hits)])
      ),
      GeneBiotype = as.character(gene_meta$biotype[subjectHits(hits)])
    ) %>%
      group_by(window_id) %>%
      summarise(
        gene_count = n_distinct(GeneID),
        GeneIDs = paste(unique(GeneID), collapse = ";"),
        GeneNames = paste(unique(GeneName), collapse = ";"),
        GeneBiotypes = paste(unique(stats::na.omit(GeneBiotype)), collapse = ";"),
        .groups = "drop"
      )
  }

  df %>%
    mutate(window_id = row_number()) %>%
    left_join(gene_summary, by = "window_id") %>%
    mutate(
      gene_count = replace_na(gene_count, 0L),
      GeneIDs = replace_na(GeneIDs, ""),
      GeneNames = replace_na(GeneNames, ""),
      GeneBiotypes = replace_na(GeneBiotypes, ""),
      threshold_label = "5pct",
      all_three_signals = if_else(all_three_signals, "Yes", "No"),
      fst_outlier = if_else(fst_outlier, "Yes", "No"),
      deltaD_outlier = if_else(deltaD_outlier, "Yes", "No"),
      pi_outlier = if_else(pi_outlier, "Yes", "No")
    ) %>%
    select(
      Population,
      chr,
      chr_num,
      WinCenter,
      window_start,
      window_end,
      threshold_label,
      Fst,
      Tajima_Gen0,
      Tajima_Gen14,
      deltaD,
      pi_Gen0,
      pi_Gen14,
      pi_ratio,
      log_pi_ratio,
      fst_outlier,
      deltaD_outlier,
      pi_outlier,
      all_three_signals,
      gene_count,
      GeneIDs,
      GeneNames,
      GeneBiotypes,
      Nsites_fst,
      nSites_Gen0,
      nSites_Gen14
    ) %>%
    arrange(chr_num, WinCenter)
}

################################################################################
### LOAD GENE ANNOTATION
################################################################################

cat("\n=== Loading maize gene annotation ===\n")
reference_gr <- rtracklayer::import(gff_file)
genes_only <- reference_gr[mcols(reference_gr)$type == "gene"]
cat("Genes in reference:", length(genes_only), "\n")

################################################################################
### BUILD TABLES
################################################################################

cat("\n=== Processing Indian Chief windows ===\n")
ic_windows <- prepare_population_windows(
  fst_dir = fst_dir_ic,
  tajima_dir = tajima_dir_ic,
  gen0_pattern = "IC01",
  gen14_pattern = "IC14",
  population_name = "Indian Chief"
)
ic_annotated <- annotate_windows_with_genes(ic_windows, genes_only)

cat("\n=== Processing Jarvis windows ===\n")
j_windows <- prepare_population_windows(
  fst_dir = fst_dir_j,
  tajima_dir = tajima_dir_j,
  gen0_pattern = "J01",
  gen14_pattern = "J14",
  population_name = "Jarvis"
)
j_annotated <- annotate_windows_with_genes(j_windows, genes_only)

################################################################################
### SAVE OUTPUTS
################################################################################

ic_file <- file.path(output_dir, "SuppTable_Fig4_windows_IndianChief_genes_5pct.csv")
j_file <- file.path(output_dir, "SuppTable_Fig4_windows_Jarvis_genes_5pct.csv")

readr::write_csv(ic_annotated, ic_file)
readr::write_csv(j_annotated, j_file)

cat("\nSaved supplementary tables to:\n")
cat("  ", ic_file, "\n", sep = "")
cat("  ", j_file, "\n", sep = "")

cat("\nSummary:\n")
cat(
  "  Indian Chief rows:", nrow(ic_annotated),
  "| rows with genes:", sum(ic_annotated$gene_count > 0),
  "| all three signals:", sum(ic_annotated$all_three_signals == "Yes"), "\n"
)
cat(
  "  Jarvis rows:", nrow(j_annotated),
  "| rows with genes:", sum(j_annotated$gene_count > 0),
  "| all three signals:", sum(j_annotated$all_three_signals == "Yes"), "\n"
)
