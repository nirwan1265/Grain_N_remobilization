#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)

input_file <- if (length(args) >= 1) {
  args[[1]]
} else {
  "/Users/nirwantandukar/Documents/Research/data/Cell_designer/Grain_N/simulations/base_values.txt"
}

output_dir <- if (length(args) >= 2) {
  args[[2]]
} else {
  "data/cell_designer_simulations/base_values"
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_section <- function(lines, marker) {
  marker_idx <- match(marker, lines)
  if (is.na(marker_idx)) {
    stop(sprintf("Could not find section marker: %s", marker))
  }

  start_idx <- marker_idx + 1
  end_candidates <- which(seq_along(lines) > marker_idx & grepl("^#", lines))
  end_idx <- if (length(end_candidates) > 0) {
    min(end_candidates) - 1
  } else {
    length(lines)
  }

  section_lines <- lines[start_idx:end_idx]
  section_lines <- section_lines[nzchar(section_lines)]

  read.table(
    text = paste(section_lines, collapse = "\n"),
    header = TRUE,
    check.names = FALSE
  )
}

plot_facets <- function(df, title, subtitle, outfile, ncol = 4) {
  ggplot(df, aes(x = time, y = value)) +
    geom_line(linewidth = 0.5, color = "#1f78b4") +
    facet_wrap(~ variable, scales = "free_y", ncol = ncol) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Time",
      y = "Value"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "#e8f1fb"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )

  ggsave(outfile, width = 14, height = 9, dpi = 300)
}

lines <- readLines(input_file, warn = FALSE)
species_df <- read_section(lines, "#species")
params_df <- read_section(lines, "#parameters")

species_long <- species_df %>%
  pivot_longer(-time, names_to = "variable", values_to = "value") %>%
  mutate(module = sub("_.*", "", variable))

params_long <- params_df %>%
  pivot_longer(-time, names_to = "variable", values_to = "value") %>%
  mutate(module = sub("_.*", "", variable))

species_summary <- species_long %>%
  group_by(variable, module) %>%
  summarise(
    start_value = first(value),
    end_value = last(value),
    min_value = min(value, na.rm = TRUE),
    max_value = max(value, na.rm = TRUE),
    abs_change = end_value - start_value,
    value_range = max_value - min_value,
    .groups = "drop"
  ) %>%
  arrange(desc(value_range), desc(abs(abs_change)))

params_summary <- params_long %>%
  group_by(variable, module) %>%
  summarise(
    start_value = first(value),
    end_value = last(value),
    min_value = min(value, na.rm = TRUE),
    max_value = max(value, na.rm = TRUE),
    abs_change = end_value - start_value,
    value_range = max_value - min_value,
    .groups = "drop"
  ) %>%
  arrange(desc(value_range), desc(abs(abs_change)))

write_tsv(species_summary, file.path(output_dir, "species_summary.tsv"))
write_tsv(params_summary, file.path(output_dir, "parameter_summary.tsv"))

core_species <- c(
  "pyr_pyruvate", "pyr_AcCoA",
  "tca_oaa", "tca_akg", "tca_suc",
  "nh4_NH4", "nh4_GLU", "nh4_GLN",
  "asp_Asp", "asp_Thr", "asp_Lys", "asp_Hser", "asp_Ile"
)

core_species_present <- intersect(core_species, unique(species_long$variable))

core_plot_df <- species_long %>%
  filter(variable %in% core_species_present)

top_species <- species_summary %>%
  filter(value_range > 0) %>%
  slice_head(n = 16) %>%
  pull(variable)

top_species_df <- species_long %>%
  filter(variable %in% top_species)

parameter_plot_df <- params_long

plot_facets(
  core_plot_df,
  title = "Base Simulation: Core C-N Bridge Species",
  subtitle = basename(input_file),
  outfile = file.path(output_dir, "base_values_core_species.png"),
  ncol = 4
)

plot_facets(
  top_species_df,
  title = "Base Simulation: Most Dynamic Species",
  subtitle = "Top species ranked by value range across the time course",
  outfile = file.path(output_dir, "base_values_top_dynamic_species.png"),
  ncol = 4
)

plot_facets(
  parameter_plot_df,
  title = "Base Simulation: Dynamic Parameters",
  subtitle = basename(input_file),
  outfile = file.path(output_dir, "base_values_parameters.png"),
  ncol = 4
)

message("Wrote plots and summaries to: ", normalizePath(output_dir))
