################################################################################
### GOODMAN-BUCKLER FAA / PBAA PANEL METRICS
################################################################################

library(tidyverse)

################################################################################
### CONFIGURATION
################################################################################

faa_file <- "data/FAA_Goodman_Buckler.csv"
pbaa_file <- "data/PBAA_Goodman_Buckler.csv"
tables_dir <- "tables/supplementary"

taxa_means_file <- file.path(
  tables_dir,
  "SuppTable_12_Goodman_Buckler_FAA_PBAA_taxa_means.csv"
)
summary_file <- file.path(
  tables_dir,
  "SuppTable_13_Goodman_Buckler_FAA_PBAA_metric_summary.csv"
)
paired_tests_file <- file.path(
  tables_dir,
  "SuppTable_14_Goodman_Buckler_FAA_PBAA_paired_tests.csv"
)

dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
### HELPERS
################################################################################

significance_label <- function(p_value) {
  dplyr::case_when(
    is.na(p_value) ~ "",
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

summarise_metric <- function(df, metric_col, metric_label, source_label) {
  vals <- df[[metric_col]]
  tibble(
    source = source_label,
    metric = metric_col,
    metric_label = metric_label,
    n_taxa = sum(!is.na(vals)),
    mean = mean(vals, na.rm = TRUE),
    median = median(vals, na.rm = TRUE),
    sd = sd(vals, na.rm = TRUE),
    min = min(vals, na.rm = TRUE),
    q1 = quantile(vals, probs = 0.25, na.rm = TRUE),
    q3 = quantile(vals, probs = 0.75, na.rm = TRUE),
    max = max(vals, na.rm = TRUE)
  )
}

################################################################################
### LOAD
################################################################################

faa_raw <- readr::read_csv(faa_file, show_col_types = FALSE)
pbaa_raw <- readr::read_csv(pbaa_file, show_col_types = FALSE) %>%
  select(-any_of("Unnamed: 79"))

faa_cols <- c(
  "A", "R", "N", "D", "Q", "E", "G", "H", "I", "L",
  "K", "M", "F", "P", "S", "W", "T", "Y", "V", "C"
)

pbaa_cols <- c(
  "Ala", "Arg", "Asx", "Glx", "His", "Ile", "Leu", "Lys",
  "Met", "Phe", "Pro", "Ser", "Thr", "Tyr", "Val"
)

faa_n_atoms <- c(
  A = 1, R = 4, N = 2, D = 1, Q = 2, E = 1, G = 1, H = 3, I = 1, L = 1,
  K = 2, M = 1, F = 1, P = 1, S = 1, W = 2, T = 1, Y = 1, V = 1, C = 1
)

# Asx and Glx are unresolved Asp/Asn and Glu/Gln hydrolysis pools in PBAA.
pbaa_n_atoms_est <- c(
  Ala = 1, Arg = 4, Asx = 1.5, Glx = 1.5, His = 3, Ile = 1, Leu = 1, Lys = 2,
  Met = 1, Phe = 1, Pro = 1, Ser = 1, Thr = 1, Tyr = 1, Val = 1
)

################################################################################
### BUILD GENOTYPE MEANS
################################################################################

faa_taxa <- faa_raw %>%
  group_by(taxa) %>%
  summarise(
    n_faa_reps = n(),
    across(all_of(faa_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    Total_FAA = sum(c_across(all_of(faa_cols)), na.rm = TRUE),
    Proline_fraction_FAA = if_else(Total_FAA > 0, P / Total_FAA, NA_real_),
    Asn_fraction_FAA = if_else(Total_FAA > 0, N / Total_FAA, NA_real_),
    Asp_fraction_FAA = if_else(Total_FAA > 0, D / Total_FAA, NA_real_),
    Gln_fraction_FAA = if_else(Total_FAA > 0, Q / Total_FAA, NA_real_),
    Glu_fraction_FAA = if_else(Total_FAA > 0, E / Total_FAA, NA_real_),
    N_rich_FAA_pool = N + Q + R + K + H,
    N_rich_FAA_fraction = if_else(Total_FAA > 0, N_rich_FAA_pool / Total_FAA, NA_real_),
    FAA_N_proxy = sum(c_across(all_of(faa_cols)) * faa_n_atoms[faa_cols], na.rm = TRUE)
  ) %>%
  ungroup()

pbaa_taxa <- pbaa_raw %>%
  group_by(taxa) %>%
  summarise(
    n_pbaa_reps = n(),
    across(all_of(pbaa_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    Total_PBAA = sum(c_across(all_of(pbaa_cols)), na.rm = TRUE),
    Proline_fraction_PBAA = if_else(Total_PBAA > 0, Pro / Total_PBAA, NA_real_),
    Asx_fraction_PBAA = if_else(Total_PBAA > 0, Asx / Total_PBAA, NA_real_),
    Glx_fraction_PBAA = if_else(Total_PBAA > 0, Glx / Total_PBAA, NA_real_),
    N_rich_PBAA_pool = Asx + Glx + Arg + Lys + His,
    N_rich_PBAA_fraction = if_else(Total_PBAA > 0, N_rich_PBAA_pool / Total_PBAA, NA_real_),
    PBAA_N_proxy = sum(c_across(all_of(pbaa_cols)) * pbaa_n_atoms_est[pbaa_cols], na.rm = TRUE)
  ) %>%
  ungroup()

combined_taxa <- faa_taxa %>%
  full_join(pbaa_taxa, by = "taxa") %>%
  mutate(
    matched_FAA_PBAA = !is.na(n_faa_reps) & !is.na(n_pbaa_reps),
    Total_AA_N_proxy = if_else(matched_FAA_PBAA, FAA_N_proxy + PBAA_N_proxy, NA_real_)
  ) %>%
  arrange(taxa)

################################################################################
### DESCRIPTIVE SUMMARIES
################################################################################

summary_df <- bind_rows(
  summarise_metric(combined_taxa, "Total_FAA", "Broad free amino acid pool", "FAA"),
  summarise_metric(combined_taxa, "FAA_N_proxy", "Estimated FAA-derived nitrogen proxy", "FAA"),
  summarise_metric(combined_taxa, "Proline_fraction_FAA", "Proline / Total FAA", "FAA"),
  summarise_metric(combined_taxa, "N_rich_FAA_fraction", "Nitrogen-rich FAA pool / Total FAA", "FAA"),
  summarise_metric(combined_taxa, "Total_PBAA", "Broad protein-bound amino acid pool", "PBAA"),
  summarise_metric(combined_taxa, "PBAA_N_proxy", "Estimated PBAA-derived nitrogen proxy", "PBAA"),
  summarise_metric(combined_taxa, "Proline_fraction_PBAA", "Proline / Total PBAA", "PBAA"),
  summarise_metric(combined_taxa, "N_rich_PBAA_fraction", "Nitrogen-rich PBAA pool / Total PBAA", "PBAA"),
  summarise_metric(combined_taxa, "Total_AA_N_proxy", "Estimated total amino-acid-derived nitrogen proxy", "Combined")
)

paired_defs <- tribble(
  ~faa_metric, ~pbaa_metric, ~comparison_label,
  "Total_FAA", "Total_PBAA", "Broad amino acid pool: FAA vs PBAA",
  "FAA_N_proxy", "PBAA_N_proxy", "Estimated amino-acid-derived nitrogen proxy: FAA vs PBAA",
  "Proline_fraction_FAA", "Proline_fraction_PBAA", "Proline fraction: FAA vs PBAA",
  "N_rich_FAA_fraction", "N_rich_PBAA_fraction", "Nitrogen-rich fraction: FAA vs PBAA"
)

paired_tests_df <- paired_defs %>%
  rowwise() %>%
  mutate(
    pair_data = list(
      combined_taxa %>%
        select(taxa, faa_value = all_of(faa_metric), pbaa_value = all_of(pbaa_metric)) %>%
        drop_na()
    ),
    n_taxa = nrow(pair_data),
    mean_faa = mean(pair_data$faa_value, na.rm = TRUE),
    mean_pbaa = mean(pair_data$pbaa_value, na.rm = TRUE),
    median_faa = median(pair_data$faa_value, na.rm = TRUE),
    median_pbaa = median(pair_data$pbaa_value, na.rm = TRUE),
    p_value = tryCatch(
      wilcox.test(pair_data$faa_value, pair_data$pbaa_value, paired = TRUE, exact = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    sig = significance_label(p_value)
  ) %>%
  ungroup() %>%
  select(
    comparison_label,
    faa_metric,
    pbaa_metric,
    n_taxa,
    mean_faa,
    mean_pbaa,
    median_faa,
    median_pbaa,
    p_value,
    sig
  )

################################################################################
### WRITE
################################################################################

readr::write_csv(combined_taxa, taxa_means_file)
readr::write_csv(summary_df, summary_file)
readr::write_csv(paired_tests_df, paired_tests_file)

cat("\nSaved Goodman-Buckler FAA/PBAA taxa means to:\n")
cat("  ", taxa_means_file, "\n", sep = "")

cat("\nSaved Goodman-Buckler FAA/PBAA metric summary to:\n")
cat("  ", summary_file, "\n", sep = "")

cat("\nSaved Goodman-Buckler FAA/PBAA paired tests to:\n")
cat("  ", paired_tests_file, "\n", sep = "")

cat("\nDone.\n")
