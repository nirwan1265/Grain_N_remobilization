#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(purrr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: plot_n_sink_retune_check.R <input_dir_with_txt_exports> <output_dir>")
}

input_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_section <- function(lines, marker) {
  marker_idx <- match(marker, lines)
  if (is.na(marker_idx)) stop(sprintf("Could not find section marker: %s", marker))
  start_idx <- marker_idx + 1
  end_candidates <- which(seq_along(lines) > marker_idx & grepl("^#", lines))
  end_idx <- if (length(end_candidates) > 0) min(end_candidates) - 1 else length(lines)
  section_lines <- lines[start_idx:end_idx]
  section_lines <- section_lines[nzchar(section_lines)]
  read.table(text = paste(section_lines, collapse = "\n"), header = TRUE, check.names = FALSE)
}

parse_run_label <- function(path) {
  base <- tools::file_path_sans_ext(basename(path))
  recode(
    base,
    "base_new" = "Baseline",
    "baseline_new" = "Baseline",
    "baseline" = "Baseline",
    "base_values" = "Baseline",
    "n_sink_glu_new" = "GLU sink x0.75",
    "n_sink_glu_0p75" = "GLU sink x0.75",
    "n_sink_gln_new" = "GLN sink x0.75",
    "n_sink_gln_0p75" = "GLN sink x0.75",
    "n_sink_both_0p75" = "Both sinks x0.75",
    "n_sink_both_0p5" = "Both sinks x0.5",
    .default = base
  )
}

load_run <- function(path) {
  label <- parse_run_label(path)
  lines <- readLines(path, warn = FALSE)
  species_df <- read_section(lines, "#species") %>%
    pivot_longer(-time, names_to = "variable", values_to = "value") %>%
    mutate(run = label, source_file = basename(path))
  species_df
}

txt_files <- list.files(input_dir, pattern = "\\.txt$", full.names = TRUE)
all_species <- map_dfr(txt_files, load_run)

if (!"Baseline" %in% all_species$run) {
  stop("No baseline file found; expected base_new.txt, baseline.txt, or base_values.txt")
}

plot_theme <- theme_minimal(base_size = 20) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "white", color = "grey80"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 11, face = "bold", color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.35),
    legend.position = "bottom"
  )

focus_vars <- c("nh4_NH4", "nh4_GLU", "nh4_GLN", "tca_akg", "asp_Asp", "asp_Thr")
display_labels <- c(
  "nh4_NH4" = "NH4",
  "nh4_GLU" = "Glu",
  "nh4_GLN" = "Gln",
  "tca_akg" = "2-OG",
  "asp_Asp" = "Asp",
  "asp_Thr" = "Thr"
)

palette <- c(
  "Baseline" = "#1B1B1B",
  "GLU sink x0.75" = "#D95F02",
  "GLN sink x0.75" = "#7570B3",
  "Both sinks x0.75" = "#1B9E77",
  "Both sinks x0.5" = "#1F78B4"
)

timecourse_df <- all_species %>%
  filter(variable %in% focus_vars) %>%
  mutate(
    display = factor(display_labels[variable], levels = unname(display_labels[focus_vars])),
    run = factor(run, levels = names(palette)[names(palette) %in% unique(run)])
  )

timecourse_plot <- ggplot(timecourse_df, aes(x = time, y = value, color = run)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ display, scales = "free_y", ncol = 3) +
  scale_color_manual(values = palette, drop = FALSE) +
  labs(
    title = "Conservative N-sink retune check",
    subtitle = "Baseline versus mild glutamate/glutamine demand reductions",
    x = "Time",
    y = "Concentration / model value",
    color = NULL
  ) +
  plot_theme

ggsave(
  file.path(output_dir, "n_sink_retune_timecourses.png"),
  timecourse_plot,
  width = 15,
  height = 10,
  dpi = 300
)

endpoint_summary <- all_species %>%
  filter(variable %in% focus_vars) %>%
  group_by(run, variable) %>%
  arrange(time, .by_group = TRUE) %>%
  summarise(
    start_value = first(value),
    end_value = last(value),
    min_value = min(value, na.rm = TRUE),
    max_value = max(value, na.rm = TRUE),
    .groups = "drop"
  )

baseline_end <- endpoint_summary %>%
  filter(run == "Baseline") %>%
  select(variable, baseline_end = end_value)

endpoint_summary <- endpoint_summary %>%
  left_join(baseline_end, by = "variable") %>%
  mutate(
    delta_end = end_value - baseline_end,
    log2fc_end = log2((end_value + 1e-12) / (baseline_end + 1e-12)),
    display = display_labels[variable]
  )

write_tsv(endpoint_summary, file.path(output_dir, "n_sink_retune_endpoint_summary.tsv"))

heatmap_df <- endpoint_summary %>%
  filter(run != "Baseline") %>%
  mutate(
    run = factor(run, levels = c("GLU sink x0.75", "GLN sink x0.75", "Both sinks x0.75", "Both sinks x0.5")),
    display = factor(display, levels = unname(display_labels[focus_vars]))
  )

heatmap_plot <- ggplot(heatmap_df, aes(x = display, y = run, fill = log2fc_end)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(low = "#2166AC", mid = "#FFF2CC", high = "#B2182B", midpoint = 0) +
  labs(
    title = "Endpoint effect of conservative N-sink retunes",
    subtitle = "Relative to the new baseline",
    x = NULL,
    y = NULL,
    fill = "log2FC"
  ) +
  plot_theme +
  theme(axis.text.x = element_text(angle = 0))

ggsave(
  file.path(output_dir, "n_sink_retune_endpoint_heatmap.png"),
  heatmap_plot,
  width = 10,
  height = 5.5,
  dpi = 300
)

message("Wrote figures and summaries to: ", normalizePath(output_dir))
