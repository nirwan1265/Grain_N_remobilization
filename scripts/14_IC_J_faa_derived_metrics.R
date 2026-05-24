################################################################################
### SELECTED FAA AND PBAA METRICS FOR IC / J
################################################################################

library(tidyverse)

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
        p_value = p_val
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
      raw_sig,
      bh_sig,
      result_summary
    )
}

################################################################################
### FAA
################################################################################

cat("\n=== Building selected FAA and PBAA metrics ===\n")

faa_cols <- c(
  "Ala", "Arg", "Asn", "Asp", "Gln", "Glu", "Gly", "His", "Ile", "Leu",
  "Lys", "Met", "Phe", "Pro", "Ser", "Trp", "Thr", "Tyr", "Val", "Cys"
)

faa_samples_df <- build_biological_means(faa_file, faa_cols) %>%
  mutate(
    Total_FAA = Ala + Arg + Asn + Asp + Gln + Glu + Gly + His + Ile + Leu + Lys +
      Met + Phe + Pro + Ser + Trp + Thr + Tyr + Val + Cys,
    Proline_fraction = if_else(Total_FAA > 0, Pro / Total_FAA, NA_real_),
    Asn_fraction = if_else(Total_FAA > 0, Asn / Total_FAA, NA_real_),
    Asp_fraction = if_else(Total_FAA > 0, Asp / Total_FAA, NA_real_),
    Gln_fraction = if_else(Total_FAA > 0, Gln / Total_FAA, NA_real_),
    Glu_fraction = if_else(Total_FAA > 0, Glu / Total_FAA, NA_real_),
    N_rich_FAA_pool = Asn + Gln + Arg + Lys + His
  ) %>%
  select(
    population,
    generation,
    biological_sample,
    bio_id,
    all_of(faa_cols),
    Total_FAA,
    Proline_fraction,
    Asn_fraction,
    Asp_fraction,
    Gln_fraction,
    Glu_fraction,
    N_rich_FAA_pool
  )

faa_metric_descriptions <- c(
  "Total_FAA" = "Broad amino acid pool",
  "Proline_fraction" = "Proline / Total FAA",
  "Asn_fraction" = "Asn / Total FAA",
  "Asp_fraction" = "Asp / Total FAA",
  "Gln_fraction" = "Gln / Total FAA",
  "Glu_fraction" = "Glu / Total FAA",
  "N_rich_FAA_pool" = "Nitrogen-rich amino acid pool"
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

pbaa_samples_df <- build_biological_means(pbaa_file, pbaa_cols) %>%
  mutate(
    Total_PBAA = Ala + Arg + Asx + Glx + Gly + His + Ile + Leu +
      Lys + Met + Phe + Pro + Ser + Thr + Tyr + Val,
    Proline_fraction = if_else(Total_PBAA > 0, Pro / Total_PBAA, NA_real_),
    Asx_fraction = if_else(Total_PBAA > 0, Asx / Total_PBAA, NA_real_),
    Glx_fraction = if_else(Total_PBAA > 0, Glx / Total_PBAA, NA_real_),
    N_rich_PBAA_pool = Asx + Glx + Arg + Lys + His
  ) %>%
  select(
    population,
    generation,
    biological_sample,
    bio_id,
    all_of(pbaa_cols),
    Total_PBAA,
    Proline_fraction,
    Asx_fraction,
    Glx_fraction,
    N_rich_PBAA_pool
  )

pbaa_metric_descriptions <- c(
  "Total_PBAA" = "Broad protein-bound amino acid pool",
  "Proline_fraction" = "Proline / Total PBAA",
  "Asx_fraction" = "Asx / Total PBAA",
  "Glx_fraction" = "Glx / Total PBAA",
  "N_rich_PBAA_pool" = "Nitrogen-rich protein-bound amino acid pool"
)

pbaa_stats_df <- compute_metric_stats(pbaa_samples_df, names(pbaa_metric_descriptions))
pbaa_summary_df <- build_summary_table(pbaa_stats_df, pbaa_metric_descriptions)

################################################################################
### WRITE
################################################################################

readr::write_csv(faa_samples_df, faa_samples_file)
readr::write_csv(faa_summary_df, faa_summary_file)
readr::write_csv(pbaa_samples_df, pbaa_samples_file)
readr::write_csv(pbaa_summary_df, pbaa_summary_file)

cat("\nSaved FAA biological means to:\n")
cat("  ", faa_samples_file, "\n", sep = "")

cat("\nSaved FAA summary to:\n")
cat("  ", faa_summary_file, "\n", sep = "")

cat("\nSaved PBAA biological means to:\n")
cat("  ", pbaa_samples_file, "\n", sep = "")

cat("\nSaved PBAA summary to:\n")
cat("  ", pbaa_summary_file, "\n", sep = "")
