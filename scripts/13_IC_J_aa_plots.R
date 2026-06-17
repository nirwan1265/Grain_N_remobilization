################################################################################
### IC / J AMINO ACID PROFILE PLOTS
### FAA radial panels for the main figure and PBAA radial panels for supplement
################################################################################

library(tidyverse)
library(ggplot2)
library(patchwork)
library(scales)
library(effsize)

################################################################################
### CONFIGURATION
################################################################################

faa_file <- "data/FAA_IC_J.csv"
pbaa_file <- "data/PBAA_IC_J.csv"

supp_dir <- "Figs/Supplementary"
tables_dir <- "tables/supplementary"

supp_plot_file <- file.path(supp_dir, "SuppFig_IC_J_amino_profiles.png")
stats_file <- file.path(tables_dir, "SuppTable_IC_J_amino_profile_stats.csv")

dir.create(supp_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

group_colors <- c(
  "Indian Chief G0" = "#92C5DE",
  "Indian Chief G14" = "#2166AC",
  "Jarvis G0" = "#FDB863",
  "Jarvis G14" = "#B35806"
)

amino_plot_theme <- theme_minimal(base_size = 22) +
  theme(
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 8)
    ),
    plot.tag = element_text(size = 24, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    legend.background = element_rect(fill = "white", color = "grey70", linewidth = 0.4),
    legend.text = element_text(size = 18, color = "black"),
    legend.key.width = unit(26, "pt"),
    legend.key.height = unit(14, "pt"),
    legend.spacing.x = unit(10, "pt"),
    legend.margin = margin(4, 8, 4, 8),
    plot.margin = margin(6, 8, 6, 8)
  )

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
      group_label = paste(population, generation)
    )
}

prepare_pool_data <- function(path, pool_name) {
  raw_df <- readr::read_csv(path, show_col_types = FALSE)
  amino_order <- setdiff(names(raw_df), c("ID", "Wt (mg)"))

  sample_info <- parse_sample_ids(raw_df$ID)

  long_df <- raw_df %>%
    pivot_longer(
      cols = all_of(amino_order),
      names_to = "amino_acid",
      values_to = "value"
    ) %>%
    left_join(sample_info, by = "ID") %>%
    mutate(
      pool = pool_name,
      amino_acid = factor(amino_acid, levels = amino_order)
    )

  bio_means <- long_df %>%
    group_by(pool, population, generation, group_label, bio_id, amino_acid) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

  group_means <- bio_means %>%
    group_by(pool, population, generation, group_label, amino_acid) %>%
    summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")

  list(
    long = long_df,
    bio_means = bio_means,
    group_means = group_means,
    amino_order = amino_order
  )
}

convert_to_relative_pool <- function(pool_obj) {
  rel_bio_means <- pool_obj$bio_means %>%
    group_by(pool, population, generation, group_label, bio_id) %>%
    mutate(total_value = sum(value, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(value = dplyr::if_else(total_value > 0, value / total_value, NA_real_)) %>%
    select(-total_value)

  rel_group_means <- rel_bio_means %>%
    group_by(pool, population, generation, group_label, amino_acid) %>%
    summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")

  list(
    long = pool_obj$long,
    bio_means = rel_bio_means,
    group_means = rel_group_means,
    amino_order = pool_obj$amino_order
  )
}

compute_generation_stats <- function(bio_means, pool_name) {
  populations <- unique(bio_means$population)

  purrr::map_dfr(populations, function(pop_name) {
    amino_levels <- levels(droplevels(filter(bio_means, population == pop_name)$amino_acid))

    purrr::map_dfr(amino_levels, function(amino_name) {
      sub_df <- bio_means %>%
        filter(population == pop_name, amino_acid == amino_name)

      g0_vals <- sub_df %>%
        filter(generation == "G0") %>%
        pull(value)

      g14_vals <- sub_df %>%
        filter(generation == "G14") %>%
        pull(value)

      if (length(g0_vals) == 0 || length(g14_vals) == 0) {
        return(NULL)
      }

      p_val <- tryCatch(
        wilcox.test(g0_vals, g14_vals, exact = FALSE)$p.value,
        error = function(e) NA_real_
      )

      cliff_obj <- tryCatch(
        effsize::cliff.delta(g14_vals, g0_vals),
        error = function(e) NULL
      )

      comparison_label <- if (pop_name == "Indian Chief") {
        "IC0_vs_IC14"
      } else {
        "JV0_vs_JV14"
      }

      tibble(
        pool = pool_name,
        population = pop_name,
        comparison = comparison_label,
        amino_acid = amino_name,
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

build_radar_line_df <- function(group_means, amino_order, population_name) {
  group_means %>%
    filter(population == population_name) %>%
    mutate(amino_acid = factor(amino_acid, levels = amino_order)) %>%
    arrange(group_label, amino_acid) %>%
    group_by(group_label) %>%
    reframe({
      tibble(
        x = c(seq_along(amino_order), 1),
        y = c(mean_value, mean_value[1])
      )
    }) %>%
    ungroup()
}

build_label_df <- function(amino_order, stats_df, population_name, label_radius) {
  n_amino <- length(amino_order)

  base_df <- tibble(
    amino_acid = factor(amino_order, levels = amino_order),
    x = seq_len(n_amino)
  ) %>%
    mutate(
      angle = 90 - 360 * (x - 1) / n_amino,
      hjust = if_else(angle < -90, 1, 0),
      angle = if_else(angle < -90, angle + 180, angle),
      y = label_radius,
      label = as.character(amino_acid)
    )

  sig_df <- base_df %>%
    left_join(
      stats_df %>%
        filter(population == population_name) %>%
        select(amino_acid, raw_sig),
      by = "amino_acid"
    ) %>%
    filter(!is.na(raw_sig), raw_sig != "") %>%
    mutate(
      y = label_radius * 1.05,
      label = raw_sig
    )

  list(base = base_df, sig = sig_df)
}

make_radar_plot <- function(
  group_means,
  stats_df,
  amino_order,
  population_name,
  plot_title,
  base_theme = NULL,
  amino_label_size = 3.3,
  sig_label_size = 3.3,
  grid_label_size = 2.8,
  trace_linewidth = 1.2,
  amino_label_face = "plain",
  sig_label_face = "bold",
  plot_title_size = 14
) {
  line_df <- build_radar_line_df(group_means, amino_order, population_name)

  max_value <- max(line_df$y, na.rm = TRUE)
  radius_limit <- max_value * 1.22
  label_radius <- radius_limit * 0.92

  grid_breaks <- pretty(c(0, radius_limit), n = 4)
  grid_breaks <- grid_breaks[grid_breaks > 0 & grid_breaks < radius_limit]

  label_dfs <- build_label_df(
    amino_order = amino_order,
    stats_df = stats_df,
    population_name = population_name,
    label_radius = label_radius
  )

  grid_df <- tibble(
    x = 1,
    y = grid_breaks,
    label = format(grid_breaks, trim = TRUE, scientific = FALSE)
  )

  plot_colors <- if (population_name == "Indian Chief") {
    group_colors[c("Indian Chief G0", "Indian Chief G14")]
  } else {
    group_colors[c("Jarvis G0", "Jarvis G14")]
  }

  spoke_color <- "grey78"
  ring_color <- "grey76"

  applied_theme <- if (is.null(base_theme)) amino_plot_theme else base_theme

  ggplot() +
    geom_hline(yintercept = grid_breaks, color = ring_color, linewidth = 0.6) +
    geom_vline(
      xintercept = seq_along(amino_order),
      color = spoke_color,
      linewidth = 0.5
    ) +
    geom_polygon(
      data = line_df,
      aes(x = x, y = y, group = group_label, fill = group_label),
      alpha = 0.05,
      color = NA,
      show.legend = FALSE
    ) +
    geom_path(
      data = line_df,
      aes(x = x, y = y, group = group_label, color = group_label),
      linewidth = trace_linewidth,
      lineend = "round"
    ) +
    geom_text(
      data = label_dfs$base,
      aes(x = x, y = y, label = label, angle = angle, hjust = hjust),
      color = "black",
      size = amino_label_size,
      fontface = amino_label_face,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = label_dfs$sig,
      aes(x = x, y = y, label = label, angle = angle, hjust = hjust),
      color = "firebrick",
      size = sig_label_size,
      fontface = sig_label_face,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = grid_df,
      aes(x = x, y = y, label = label),
      color = "grey35",
      size = grid_label_size,
      hjust = -0.1,
      inherit.aes = FALSE
    ) +
    coord_polar(start = -pi / 2, clip = "off") +
    scale_x_continuous(limits = c(0.5, length(amino_order) + 0.5), breaks = NULL) +
    scale_y_continuous(limits = c(0, radius_limit), breaks = NULL) +
    scale_color_manual(values = plot_colors) +
    scale_fill_manual(values = plot_colors) +
    labs(title = plot_title, x = NULL, y = NULL) +
    applied_theme +
    theme(
      plot.title = element_text(size = plot_title_size, face = "bold", hjust = 0.5, margin = margin(b = 8)),
      axis.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      panel.grid = element_blank(),
      legend.position = "top",
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", color = "grey70", linewidth = 0.4),
      legend.text = element_text(size = 24, color = "black"),
      legend.key.width = unit(38, "pt"),
      legend.key.height = unit(20, "pt"),
      legend.spacing.x = unit(16, "pt"),
      legend.margin = margin(7, 12, 7, 12),
      plot.margin = margin(6, 8, 6, 8)
    )
}

build_amino_profile_outputs <- function(
  base_theme = NULL,
  faa_amino_label_size = 4.0,
  faa_sig_label_size = 4.5,
  faa_grid_label_size = 3.2,
  faa_trace_linewidth = 1.55,
  faa_amino_label_face = "plain",
  faa_sig_label_face = "bold",
  faa_plot_title_size = 14
) {
  faa_obj <- prepare_pool_data(faa_file, "FAA")
  faa_relative_obj <- convert_to_relative_pool(faa_obj)
  pbaa_obj <- prepare_pool_data(pbaa_file, "PBAA")
  pbaa_relative_obj <- convert_to_relative_pool(pbaa_obj)

  stats_df <- bind_rows(
    compute_generation_stats(faa_relative_obj$bio_means, "FAA"),
    compute_generation_stats(pbaa_relative_obj$bio_means, "PBAA")
  )

  faa_ic_plot <- make_radar_plot(
    group_means = faa_relative_obj$group_means,
    stats_df = stats_df %>% filter(pool == "FAA"),
    amino_order = faa_relative_obj$amino_order,
    population_name = "Indian Chief",
    plot_title = "Indian Chief relative FAA composition, C0 vs C14",
    base_theme = base_theme,
    amino_label_size = faa_amino_label_size,
    sig_label_size = faa_sig_label_size,
    grid_label_size = faa_grid_label_size,
    trace_linewidth = faa_trace_linewidth,
    amino_label_face = faa_amino_label_face,
    sig_label_face = faa_sig_label_face,
    plot_title_size = faa_plot_title_size
  )

  faa_j_plot <- make_radar_plot(
    group_means = faa_relative_obj$group_means,
    stats_df = stats_df %>% filter(pool == "FAA"),
    amino_order = faa_relative_obj$amino_order,
    population_name = "Jarvis",
    plot_title = "Jarvis relative FAA composition, C0 vs C14",
    base_theme = base_theme,
    amino_label_size = faa_amino_label_size,
    sig_label_size = faa_sig_label_size,
    grid_label_size = faa_grid_label_size,
    trace_linewidth = faa_trace_linewidth,
    amino_label_face = faa_amino_label_face,
    sig_label_face = faa_sig_label_face,
    plot_title_size = faa_plot_title_size
  )

  pbaa_ic_plot <- make_radar_plot(
    group_means = pbaa_relative_obj$group_means,
    stats_df = stats_df %>% filter(pool == "PBAA"),
    amino_order = pbaa_relative_obj$amino_order,
    population_name = "Indian Chief",
    plot_title = "Indian Chief relative protein-bound amino acid composition, C0 vs C14"
  )

  pbaa_j_plot <- make_radar_plot(
    group_means = pbaa_relative_obj$group_means,
    stats_df = stats_df %>% filter(pool == "PBAA"),
    amino_order = pbaa_relative_obj$amino_order,
    population_name = "Jarvis",
    plot_title = "Jarvis relative protein-bound amino acid composition, C0 vs C14"
  )

  pbaa_figure <- (pbaa_ic_plot + pbaa_j_plot) +
    plot_layout(ncol = 2) +
    plot_annotation(tag_levels = "A")

  list(
    faa_ic_plot = faa_ic_plot,
    faa_j_plot = faa_j_plot,
    pbaa_ic_plot = pbaa_ic_plot,
    pbaa_j_plot = pbaa_j_plot,
    pbaa_figure = pbaa_figure,
    stats = stats_df
  )
}

save_amino_profile_outputs <- function(outputs) {
  readr::write_csv(outputs$stats, stats_file)
  ggsave(
    supp_plot_file,
    outputs$pbaa_figure,
    width = 14,
    height = 7.8,
    dpi = 300,
    bg = "white"
  )
}

################################################################################
### RUN
################################################################################

if (sys.nframe() == 0) {
  cat("\n=== Building IC / J amino-acid figures ===\n")
  amino_outputs <- build_amino_profile_outputs()
  save_amino_profile_outputs(amino_outputs)

  cat("\nSaved supplementary figure to:\n")
  cat("  ", supp_plot_file, "\n", sep = "")

  cat("\nSaved stats table to:\n")
  cat("  ", stats_file, "\n", sep = "")
}
