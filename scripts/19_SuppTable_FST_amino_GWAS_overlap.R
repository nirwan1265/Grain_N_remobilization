################################################################################
### SUPPLEMENTARY TABLE: FST x AMINO-ACID GWAS OVERLAP
################################################################################

library(tidyverse)

################################################################################
### CONFIGURATION
################################################################################

gwas_file <- "tables/supplementary/SuppTable_amino_gwas_gene_best_by_phenotype_25kb.csv"
ic_file <- "tables/supplementary/SuppTable_IC_GenWin_windows_Fst_Tajima_pi_genes.csv"
j_file <- "tables/supplementary/SuppTable_J_GenWin_windows_Fst_Tajima_pi_genes.csv"
output_dir <- "tables/supplementary"

csv_file <- file.path(
  output_dir,
  "SuppTable_IC_J_FST_amino_GWAS_overlap_all_traits.csv"
)
tex_file <- file.path(
  output_dir,
  "SuppTable_IC_J_FST_amino_GWAS_overlap_all_traits.tex"
)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
### HELPERS
################################################################################

normalize_flag <- function(x) {
  tolower(trimws(as.character(x))) %in% c("yes", "true", "1")
}

split_gene_ids <- function(x) {
  x %>%
    as.character() %>%
    strsplit(";", fixed = TRUE) %>%
    purrr::flatten_chr() %>%
    trimws() %>%
    discard(~ .x == "" || is.na(.x)) %>%
    unique() %>%
    sort()
}

collapse_or_none <- function(x) {
  if (length(x) == 0) {
    return("None")
  }
  paste(x, collapse = "; ")
}

latex_escape <- function(x) {
  x %>%
    gsub("\\\\", "\\\\textbackslash{}", ., perl = TRUE) %>%
    gsub("_", "\\\\_", ., fixed = TRUE) %>%
    gsub("%", "\\\\%", ., fixed = TRUE) %>%
    gsub("&", "\\\\&", ., fixed = TRUE) %>%
    gsub("#", "\\\\#", ., fixed = TRUE)
}

extract_fst_gene_set <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE)

  yes_rows <- df %>%
    filter(normalize_flag(fst_outlier))

  split_gene_ids(yes_rows$GeneIDs)
}

################################################################################
### LOAD
################################################################################

gwas_df <- readr::read_csv(gwas_file, show_col_types = FALSE)
ic_genes <- extract_fst_gene_set(ic_file)
j_genes <- extract_fst_gene_set(j_file)

################################################################################
### BUILD OVERLAP TABLE
################################################################################

overlap_df <- gwas_df %>%
  group_by(Phenotype) %>%
  summarise(
    gwas_genes = list(sort(unique(GeneID))),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    IC_overlap_genes = list(intersect(gwas_genes, ic_genes)),
    J_overlap_genes = list(intersect(gwas_genes, j_genes)),
    Shared_IC_J_genes = list(intersect(IC_overlap_genes, J_overlap_genes)),
    IC_overlap_n = length(IC_overlap_genes),
    J_overlap_n = length(J_overlap_genes),
    Shared_IC_J_n = length(Shared_IC_J_genes),
    IC_overlap_genes_txt = collapse_or_none(IC_overlap_genes),
    J_overlap_genes_txt = collapse_or_none(J_overlap_genes),
    Shared_IC_J_genes_txt = collapse_or_none(Shared_IC_J_genes)
  ) %>%
  ungroup() %>%
  filter(IC_overlap_n > 0 | J_overlap_n > 0) %>%
  arrange(desc(IC_overlap_n + J_overlap_n), Phenotype) %>%
  select(
    Trait = Phenotype,
    IC_overlap_n,
    IC_overlap_genes = IC_overlap_genes_txt,
    J_overlap_n,
    J_overlap_genes = J_overlap_genes_txt,
    Shared_IC_J_n,
    Shared_IC_J_genes = Shared_IC_J_genes_txt
  )

################################################################################
### WRITE CSV
################################################################################

readr::write_csv(overlap_df, csv_file)

################################################################################
### WRITE LATEX LONGTABLE
################################################################################

caption_text <- paste(
  "Overlap between amino-acid GWAS candidate genes and spline-based",
  "$F_{ST}$ outlier genes in Indian Chief (IC) and Jarvis (J).",
  "$F_{ST}$ genes were defined from spline windows with",
  "\\texttt{fst\\_outlier = Yes}.",
  "Traits are ordered by the total number of overlapping IC and J genes."
)

tex_lines <- c(
  "\\begin{longtable}{p{1.9cm}p{0.9cm}p{4.8cm}p{0.9cm}p{4.8cm}p{0.9cm}p{3.2cm}}",
  paste0("\\caption{", caption_text, "}\\\\"),
  "\\label{tab:fst_amino_gwas_overlap_all}\\\\",
  "\\hline",
  "Trait & IC $n$ & IC overlap genes & J $n$ & J overlap genes & Shared $n$ & Shared IC/J genes \\\\",
  "\\hline",
  "\\endfirsthead",
  "\\hline",
  "Trait & IC $n$ & IC overlap genes & J $n$ & J overlap genes & Shared $n$ & Shared IC/J genes \\\\",
  "\\hline",
  "\\endhead",
  "\\hline",
  "\\endfoot"
)

row_lines <- overlap_df %>%
  mutate(
    Trait = latex_escape(Trait),
    IC_overlap_genes = latex_escape(IC_overlap_genes),
    J_overlap_genes = latex_escape(J_overlap_genes),
    Shared_IC_J_genes = latex_escape(Shared_IC_J_genes)
  ) %>%
  transmute(
    line = paste0(
      Trait, " & ",
      IC_overlap_n, " & ",
      IC_overlap_genes, " & ",
      J_overlap_n, " & ",
      J_overlap_genes, " & ",
      Shared_IC_J_n, " & ",
      Shared_IC_J_genes, " \\\\"
    )
  ) %>%
  pull(line)

tex_lines <- c(
  tex_lines,
  row_lines,
  "\\hline",
  "\\end{longtable}"
)

writeLines(tex_lines, tex_file)

cat("\nSaved overlap CSV to:\n")
cat("  ", csv_file, "\n", sep = "")
cat("\nSaved overlap LaTeX table to:\n")
cat("  ", tex_file, "\n", sep = "")
cat("\nRows with any overlap:", nrow(overlap_df), "\n")
cat("Traits with shared IC/J genes:", sum(overlap_df$Shared_IC_J_n > 0), "\n")
