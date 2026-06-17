################################################################################
### SELECTED FAA / PBAA METRICS AND AA-DERIVED N FOR IC / J
################################################################################

library(tidyverse)
library(effsize)

################################################################################
### CONFIGURATION
################################################################################

faa_file <- "data/FAA_IC_J.csv"
pbaa_file <- "data/PBAA_IC_J.csv"
tables_dir <- "tables/supplementary"

faa_samples_file <- file.path(
  tables_dir,
  "SuppTable_IC_J_FAA_biological_means_selected_metrics.csv"
)
faa_summary_file <- file.path(
  tables_dir,
  "SuppTable_IC_J_FAA_selected_metrics_summary.csv"
)
pbaa_samples_file <- file.path(
  tables_dir,
  "SuppTable_IC_J_PBAA_biological_means_selected_metrics.csv"
)
pbaa_summary_file <- file.path(
  tables_dir,
  "SuppTable_IC_J_PBAA_selected_metrics_summary.csv"
)
aa_n_samples_file <- file.path(
  tables_dir,
  "SuppTable_IC_J_AA_N_biological_means.csv"
)
aa_n_summary_file <- file.path(
  tables_dir,
  "SuppTable_IC_J_AA_N_summary.csv"
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
    TRUE ~ ""
  )
}

parse_sample_ids <- function(ids) {
  tibble(ID = ids) %>%
    tidyr::extract(
      ID,
      into = c("population_code", "generation_code", "bio_id", "tech_id"),
      regex = "^(IC|JV)(0|14)_([0-9]+)_([0-9]+)$",
      remove = FALSE
    ) %>%
    mutate(
      population = dplyr::recode(
        population_code,
        "IC" = "Indian Chief",
        "JV" = "Jarvis"
      ),
      generation = dplyr::recode(
        generation_code,
        "0" = "G0",
        "14" = "G14"
      ),
      group_label = paste(population, generation),
      biological_sample = paste0(population_code, generation_code, "_", bio_id)
    )
}

summarize_result <- function(mean_group0, mean_group14, p_value) {
  if (is.na(p_value) || is.na(mean_group0) || is.na(mean_group14)) {
    return(NA_character_)
  }

  if (p_value < 0.05) {
    if (mean_group14 > mean_group0) {
      return("Significantly increased in G14")
    }
    return("Significantly decreased in G14")
  }

  if (mean_group14 > mean_group0) {
    return("Numerically increased in G14")
  }
  if (mean_group14 < mean_group0) {
    return("Numerically decreased in G14")
  }
  "No mean change"
}

build_biological_means <- function(path, amino_cols) {
  raw_df <- readr::read_csv(path, show_col_types = FALSE)

  raw_df %>%
    left_join(parse_sample_ids(raw_df$ID), by = "ID") %>%
    group_by(population, generation, group_label, biological_sample, bio_id) %>%
    summarise(
      weight_mg = mean(`Wt (mg)`, na.rm = TRUE),
      across(all_of(amino_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    arrange(population, generation, bio_id)
}

compute_metric_stats <- function(samples_df, metric_names) {
  purrr::map_dfr(unique(samples_df$population), function(pop_name) {
    purrr::map_dfr(metric_names, function(metric_name) {
      sub_df <- samples_df %>%
        filter(population == pop_name)

      g0_vals <- sub_df %>%
        filter(generation == "G0") %>%
        pull(all_of(metric_name))

      g14_vals <- sub_df %>%
        filter(generation == "G14") %>%
        pull(all_of(metric_name))

      p_val <- tryCatch(
        wilcox.test(g0_vals, g14_vals, exact = FALSE)$p.value,
        error = function(e) NA_real_
      )

      cliff_obj <- tryCatch(
        effsize::cliff.delta(g14_vals, g0_vals),
        error = function(e) NULL
      )

      tibble(
        population = pop_name,
        comparison = if_else(pop_name == "Indian Chief", "IC0_vs_IC14", "JV0_vs_JV14"),
        metric = metric_name,
        n_group0 = length(g0_vals),
        n_group14 = length(g14_vals),
        mean_group0 = mean(g0_vals, na.rm = TRUE),
        mean_group14 = mean(g14_vals, na.rm = TRUE),
        median_group0 = median(g0_vals, na.rm = TRUE),
        median_group14 = median(g14_vals, na.rm = TRUE),
        p_value = p_val,
        cliffs_delta = if (is.null(cliff_obj)) NA_real_ else unname(cliff_obj$estimate),
        cliffs_magnitude = if (is.null(cliff_obj)) NA_character_ else as.character(cliff_obj$magnitude)
      )
    }) %>%
      mutate(
        p_adj = p.adjust(p_value, method = "BH"),
        raw_sig = vapply(p_value, significance_label, character(1)),
        bh_sig = vapply(p_adj, significance_label, character(1))
      )
  })
}

build_summary_table <- function(stats_df, metric_descriptions) {
  metric_names <- names(metric_descriptions)

  stats_df %>%
    filter(metric %in% metric_names) %>%
    mutate(
      metric_description = dplyr::recode(metric, !!!metric_descriptions),
      result_summary = purrr::pmap_chr(
        list(mean_group0, mean_group14, p_value),
        summarize_result
      )
    ) %>%
    select(
      population,
      comparison,
      metric,
      metric_description,
      n_group0,
      n_group14,
      mean_group0,
      mean_group14,
      median_group0,
      median_group14,
      p_value,
      p_adj,
      cliffs_delta,
      cliffs_magnitude,
      raw_sig,
      bh_sig,
      result_summary
    )
}

################################################################################
### FAA
################################################################################

cat("\n=== Building FAA, PBAA, and AA-derived N metrics ===\n")

faa_cols <- c(
  "Ala", "Arg", "Asn", "Asp", "Gln", "Glu", "Gly", "His", "Ile", "Leu",
  "Lys", "Met", "Phe", "Pro", "Ser", "Trp", "Thr", "Tyr", "Val", "Cys"
)

faa_n_atoms <- c(
  Ala = 1, Arg = 4, Asn = 2, Asp = 1, Gln = 2, Glu = 1, Gly = 1, His = 3,
  Ile = 1, Leu = 1, Lys = 2, Met = 1, Phe = 1, Pro = 1, Ser = 1, Trp = 2,
  Thr = 1, Tyr = 1, Val = 1, Cys = 1
)

faa_samples_df <- build_biological_means(faa_file, faa_cols) %>%
  rowwise() %>%
  mutate(
    Total_FAA = sum(c_across(all_of(faa_cols)), na.rm = TRUE),
    Proline_fraction = if_else(Total_FAA > 0, Pro / Total_FAA, NA_real_),
    Asn_fraction = if_else(Total_FAA > 0, Asn / Total_FAA, NA_real_),
    Asp_fraction = if_else(Total_FAA > 0, Asp / Total_FAA, NA_real_),
    Gln_fraction = if_else(Total_FAA > 0, Gln / Total_FAA, NA_real_),
    Glu_fraction = if_else(Total_FAA > 0, Glu / Total_FAA, NA_real_),
    N_rich_FAA_pool = Asn + Gln + Arg + Lys + His,
    FAA_N_pool = sum(c_across(all_of(faa_cols)) * faa_n_atoms[faa_cols], na.rm = TRUE),
    FAA_N_per_mg = if_else(weight_mg > 0, FAA_N_pool / weight_mg, NA_real_)
  ) %>%
  ungroup() %>%
  select(
    population,
    generation,
    biological_sample,
    bio_id,
    weight_mg,
    all_of(faa_cols),
    Total_FAA,
    Proline_fraction,
    Asn_fraction,
    Asp_fraction,
    Gln_fraction,
    Glu_fraction,
    N_rich_FAA_pool,
    FAA_N_pool,
    FAA_N_per_mg
  )

faa_metric_descriptions <- c(
  "Total_FAA" = "Broad amino acid pool",
  "Proline_fraction" = "Proline / Total FAA",
  "Asn_fraction" = "Asn / Total FAA",
  "Asp_fraction" = "Asp / Total FAA",
  "Gln_fraction" = "Gln / Total FAA",
  "Glu_fraction" = "Glu / Total FAA",
  "N_rich_FAA_pool" = "Nitrogen-rich amino acid pool",
  "FAA_N_pool" = "Free amino-acid-derived nitrogen pool",
  "FAA_N_per_mg" = "Free amino-acid-derived nitrogen per mg sample"
)

faa_stats_df <- compute_metric_stats(faa_samples_df, names(faa_metric_descriptions))
faa_summary_df <- build_summary_table(faa_stats_df, faa_metric_descriptions)

################################################################################
### PBAA
################################################################################

pbaa_cols <- c(
  "Ala", "Arg", "Asx", "Glx", "Gly", "His", "Ile", "Leu",
  "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Tyr", "Val"
)

# Asx and Glx are unresolved Asp/Asn and Glu/Gln hydrolysis pools, so we use
# 1.5 N atoms as an explicit midpoint estimate for PBAA-derived N.
pbaa_n_atoms_est <- c(
  Ala = 1, Arg = 4, Asx = 1.5, Glx = 1.5, Gly = 1, His = 3, Ile = 1, Leu = 1,
  Lys = 2, Met = 1, Phe = 1, Pro = 1, Ser = 1, Thr = 1, Tyr = 1, Val = 1
)

pbaa_samples_df <- build_biological_means(pbaa_file, pbaa_cols) %>%
  rowwise() %>%
  mutate(
    Total_PBAA = sum(c_across(all_of(pbaa_cols)), na.rm = TRUE),
    Proline_fraction = if_else(Total_PBAA > 0, Pro / Total_PBAA, NA_real_),
    Asx_fraction = if_else(Total_PBAA > 0, Asx / Total_PBAA, NA_real_),
    Glx_fraction = if_else(Total_PBAA > 0, Glx / Total_PBAA, NA_real_),
    N_rich_PBAA_pool = Asx + Glx + Arg + Lys + His,
    PBAA_N_pool_est = sum(c_across(all_of(pbaa_cols)) * pbaa_n_atoms_est[pbaa_cols], na.rm = TRUE),
    PBAA_N_per_mg_est = if_else(weight_mg > 0, PBAA_N_pool_est / weight_mg, NA_real_)
  ) %>%
  ungroup() %>%
  select(
    population,
    generation,
    biological_sample,
    bio_id,
    weight_mg,
    all_of(pbaa_cols),
    Total_PBAA,
    Proline_fraction,
    Asx_fraction,
    Glx_fraction,
    N_rich_PBAA_pool,
    PBAA_N_pool_est,
    PBAA_N_per_mg_est
  )

pbaa_metric_descriptions <- c(
  "Total_PBAA" = "Broad protein-bound amino acid pool",
  "Proline_fraction" = "Proline / Total PBAA",
  "Asx_fraction" = "Asx / Total PBAA",
  "Glx_fraction" = "Glx / Total PBAA",
  "N_rich_PBAA_pool" = "Nitrogen-rich protein-bound amino acid pool",
  "PBAA_N_pool_est" = "Estimated protein-bound amino-acid-derived nitrogen pool",
  "PBAA_N_per_mg_est" = "Estimated protein-bound amino-acid-derived nitrogen per mg sample"
)

pbaa_stats_df <- compute_metric_stats(pbaa_samples_df, names(pbaa_metric_descriptions))
pbaa_summary_df <- build_summary_table(pbaa_stats_df, pbaa_metric_descriptions)

################################################################################
### COMBINED AA-DERIVED N
################################################################################

aa_n_samples_df <- faa_samples_df %>%
  select(
    population,
    generation,
    biological_sample,
    bio_id,
    FAA_weight_mg = weight_mg,
    Total_FAA,
    FAA_N_pool,
    FAA_N_per_mg
  ) %>%
  left_join(
    pbaa_samples_df %>%
      select(
        population,
        generation,
        biological_sample,
        bio_id,
        PBAA_weight_mg = weight_mg,
        Total_PBAA,
        PBAA_N_pool_est,
        PBAA_N_per_mg_est
      ),
    by = c("population", "generation", "biological_sample", "bio_id")
  ) %>%
  mutate(
    Total_AA_N_pool_est = FAA_N_pool + PBAA_N_pool_est,
    Total_AA_N_per_mg_est = FAA_N_per_mg + PBAA_N_per_mg_est
  ) %>%
  arrange(population, generation, bio_id)

aa_n_metric_descriptions <- c(
  "FAA_N_pool" = "Free amino-acid-derived nitrogen pool",
  "FAA_N_per_mg" = "Free amino-acid-derived nitrogen per mg sample",
  "PBAA_N_pool_est" = "Estimated protein-bound amino-acid-derived nitrogen pool",
  "PBAA_N_per_mg_est" = "Estimated protein-bound amino-acid-derived nitrogen per mg sample",
  "Total_AA_N_pool_est" = "Estimated total amino-acid-derived nitrogen pool",
  "Total_AA_N_per_mg_est" = "Estimated total amino-acid-derived nitrogen per mg sample"
)

aa_n_stats_df <- compute_metric_stats(aa_n_samples_df, names(aa_n_metric_descriptions))
aa_n_summary_df <- build_summary_table(aa_n_stats_df, aa_n_metric_descriptions)

################################################################################
### WRITE
################################################################################

readr::write_csv(faa_samples_df, faa_samples_file)
readr::write_csv(faa_summary_df, faa_summary_file)
readr::write_csv(pbaa_samples_df, pbaa_samples_file)
readr::write_csv(pbaa_summary_df, pbaa_summary_file)
readr::write_csv(aa_n_samples_df, aa_n_samples_file)
readr::write_csv(aa_n_summary_df, aa_n_summary_file)

cat("\nSaved FAA biological means to:\n")
cat("  ", faa_samples_file, "\n", sep = "")

cat("\nSaved FAA summary to:\n")
cat("  ", faa_summary_file, "\n", sep = "")

cat("\nSaved PBAA biological means to:\n")
cat("  ", pbaa_samples_file, "\n", sep = "")

cat("\nSaved PBAA summary to:\n")
cat("  ", pbaa_summary_file, "\n", sep = "")

cat("\nSaved combined AA-derived N biological means to:\n")
cat("  ", aa_n_samples_file, "\n", sep = "")

cat("\nSaved combined AA-derived N summary to:\n")
cat("  ", aa_n_summary_file, "\n", sep = "")

cat("\nDone.\n")
