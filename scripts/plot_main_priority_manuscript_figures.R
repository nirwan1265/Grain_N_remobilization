#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(tidyr)
  library(purrr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
full_args <- commandArgs(trailingOnly = FALSE)

if (length(args) != 2) {
  stop("Usage: plot_main_priority_manuscript_figures.R <input_dir_with_txt_exports> <output_dir>")
}

input_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

script_flag <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(script_flag) > 0) sub("^--file=", "", script_flag[[1]]) else "scripts/plot_main_priority_manuscript_figures.R"
root_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
map_path <- file.path(root_dir, "data", "cell_designer_simulations", "grain_n_candidate_model_map.tsv")

plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 10)),
    plot.tag = element_text(size = 24, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.7),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.35),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "grey94", color = NA),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    plot.margin = margin(15, 15, 15, 15)
  )

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
  case_when(
    base %in% c("base_values", "baseline") ~ "baseline",
    TRUE ~ base
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
if (length(txt_files) < 2) {
  stop("Expected at least baseline plus one perturbation .txt file in the input directory")
}

all_species <- map_dfr(txt_files, load_run)
scenario_map <- read_tsv(map_path, show_col_types = FALSE)

run_display_map <- tribble(
  ~run, ~scenario_label,
  "baseline", "Baseline",
  "asp_kinase_x2", "Asp kinase x2",
  "asadh_x2", "ASADH x2",
  "asadh_half", "ASADH x0.5",
  "aspat_capacity_x2", "AspAT x2",
  "aspat_capacity_half", "AspAT x0.5",
  "gs_gogat_half", "GS/GOGAT x0.5",
  "akg_pool_x2", "2-OG x2",
  "akg_pool_half", "2-OG x0.5",
  "kdh_capacity_half", "KDH x0.5",
  "kbiosyn_x2", "2-OG sink x2"
)

variable_meta <- tribble(
  ~variable,      ~display,        ~module,
  "asp_Asp",      "Asp",           "Asp-family",
  "asp_ASA",      "ASA",           "Asp-family",
  "asp_Hser",     "Hser",          "Asp-family",
  "asp_PHser",    "PHser",         "Asp-family",
  "asp_Thr",      "Thr",           "Asp-family",
  "asp_Lys",      "Lys",           "Asp-family",
  "asp_Ile",      "Ile",           "Asp-family",
  "nh4_GLU",      "Glu",           "N assimilation",
  "nh4_GLN",      "Gln",           "N assimilation",
  "tca_akg",      "2-OG",          "TCA / 2-OG",
  "tca_oaa",      "OAA",           "TCA / 2-OG",
  "pyr_AcCoA",    "AcCoA",         "Pyruvate / entry",
  "pyr_pyruvate", "Pyruvate",      "Pyruvate / entry"
)

n_weights <- tribble(
  ~variable,     ~n_atoms, ~metric,
  "nh4_NH4",     1,        "Inorganic N",
  "nh4_GLU",     1,        "Glu/Gln N",
  "nh4_GLN",     2,        "Glu/Gln N",
  "asp_Asp",     1,        "Asp-family N",
  "asp_ASA",     1,        "Asp-family N",
  "asp_Hser",    1,        "Asp-family N",
  "asp_PHser",   1,        "Asp-family N",
  "asp_Thr",     1,        "Asp-family N",
  "asp_Ile",     1,        "Asp-family N",
  "asp_Lys",     2,        "Asp-family N"
)

target_runs <- run_display_map$run
all_species <- all_species %>%
  filter(run %in% target_runs)

if (!"baseline" %in% all_species$run) {
  stop("No baseline.txt file found in input directory")
}

baseline_df <- all_species %>%
  filter(run == "baseline")

baseline_vars <- c(
  "pyr_pyruvate", "pyr_AcCoA",
  "tca_akg", "tca_oaa",
  "nh4_GLU", "nh4_GLN",
  "asp_Asp", "asp_ASA",
  "asp_PHser", "asp_Thr",
  "asp_Lys", "asp_Ile"
)

baseline_plot_df <- baseline_df %>%
  filter(variable %in% baseline_vars) %>%
  left_join(variable_meta, by = "variable") %>%
  mutate(display = factor(display, levels = variable_meta$display[match(baseline_vars, variable_meta$variable)]))

baseline_panel <- ggplot(baseline_plot_df, aes(x = time, y = value)) +
  geom_line(linewidth = 0.8, color = "#2166AC") +
  facet_wrap(~ display, scales = "free_y", ncol = 4) +
  labs(
    title = "Baseline simulation supports a functioning C-N bridge and Asp-family branch",
    subtitle = "Core metabolites and downstream Asp-family products from the baseline model",
    x = "Time",
    y = "Concentration / model value"
  ) +
  plot_theme

proxy_timeseries <- bind_rows(
  baseline_df %>%
    inner_join(n_weights, by = "variable") %>%
    mutate(weighted_n = value * n_atoms) %>%
    group_by(time, metric) %>%
    summarise(value = sum(weighted_n, na.rm = TRUE), .groups = "drop"),
  baseline_df %>%
    filter(variable %in% c("asp_Thr", "asp_Ile", "asp_Lys")) %>%
    left_join(select(n_weights, variable, n_atoms), by = "variable") %>%
    mutate(weighted_n = value * n_atoms, metric = "Downstream product N") %>%
    group_by(time, metric) %>%
    summarise(value = sum(weighted_n, na.rm = TRUE), .groups = "drop")
) %>%
  distinct()

proxy_panel <- ggplot(proxy_timeseries, aes(x = time, y = value, color = metric)) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = c(
      "Asp-family N" = "#1B9E77",
      "Glu/Gln N" = "#D95F02",
      "Inorganic N" = "#7570B3",
      "Downstream product N" = "#E7298A"
    )
  ) +
  labs(
    title = "Baseline nitrogen-proxy pools remain interpretable across the time course",
    x = "Time",
    y = "Weighted N proxy",
    color = NULL
  ) +
  plot_theme +
  theme(legend.position = "bottom")

baseline_figure <- baseline_panel / proxy_panel + plot_layout(heights = c(3, 1.3))

ggsave(
  file.path(output_dir, "supp_baseline_model_behavior.png"),
  baseline_figure,
  width = 16,
  height = 13,
  dpi = 300
)

focus_vars <- c("asp_Asp", "asp_Thr", "asp_Lys", "asp_Ile", "nh4_GLU", "nh4_GLN", "tca_akg", "tca_oaa", "pyr_AcCoA")
eps <- 1e-12

endpoint_summary <- all_species %>%
  filter(variable %in% focus_vars) %>%
  group_by(run, source_file, variable) %>%
  arrange(time, .by_group = TRUE) %>%
  summarise(
    start_value = first(value),
    end_value = last(value),
    .groups = "drop"
  ) %>%
  left_join(variable_meta, by = "variable")

baseline_end <- endpoint_summary %>%
  filter(run == "baseline") %>%
  select(variable, baseline_end = end_value)

endpoint_summary <- endpoint_summary %>%
  left_join(baseline_end, by = "variable") %>%
  mutate(log2fc_end = log2((end_value + eps) / (baseline_end + eps))) %>%
  left_join(run_display_map, by = "run") %>%
  left_join(select(scenario_map, scenario, pathway_module), by = c("run" = "scenario")) %>%
  mutate(
    pathway_module = coalesce(pathway_module, "Baseline"),
    display = factor(display, levels = variable_meta$display[match(focus_vars, variable_meta$variable)]),
    run_display = factor(scenario_label, levels = rev(run_display_map$scenario_label))
  )

endpoint_heatmap <- ggplot(endpoint_summary, aes(x = display, y = run_display, fill = log2fc_end)) +
  geom_tile(color = "white", linewidth = 0.4) +
  facet_grid(pathway_module ~ module, scales = "free", space = "free") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(
    title = "Perturbations shift Asp-family, nitrogen-assimilation, and 2-OG/TCA states in distinct ways",
    subtitle = "Endpoint response relative to the baseline simulation",
    x = NULL,
    y = NULL,
    fill = "log2FC"
  ) +
  plot_theme +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    legend.position = "right"
  )

nitrogen_summary <- bind_rows(
  all_species %>%
    inner_join(n_weights, by = "variable") %>%
    mutate(weighted_n = value * n_atoms) %>%
    group_by(run, source_file, time, metric) %>%
    summarise(value = sum(weighted_n, na.rm = TRUE), .groups = "drop"),
  all_species %>%
    filter(variable %in% c("asp_Thr", "asp_Ile", "asp_Lys")) %>%
    left_join(select(n_weights, variable, n_atoms), by = "variable") %>%
    mutate(weighted_n = value * n_atoms, metric = "Downstream product N") %>%
    group_by(run, source_file, time, metric) %>%
    summarise(value = sum(weighted_n, na.rm = TRUE), .groups = "drop")
) %>%
  group_by(run, source_file, metric) %>%
  arrange(time, .by_group = TRUE) %>%
  summarise(
    start_value = first(value),
    end_value = last(value),
    .groups = "drop"
  )

baseline_proxy <- nitrogen_summary %>%
  filter(run == "baseline") %>%
  select(metric, baseline_end = end_value)

nitrogen_summary <- nitrogen_summary %>%
  left_join(baseline_proxy, by = "metric") %>%
  mutate(log2fc_end = log2((end_value + eps) / (baseline_end + eps))) %>%
  left_join(run_display_map, by = "run") %>%
  left_join(select(scenario_map, scenario, pathway_module), by = c("run" = "scenario")) %>%
  mutate(
    pathway_module = coalesce(pathway_module, "Baseline"),
    run_display = factor(scenario_label, levels = rev(run_display_map$scenario_label))
  )

proxy_heatmap <- ggplot(nitrogen_summary, aes(x = metric, y = run_display, fill = log2fc_end)) +
  geom_tile(color = "white", linewidth = 0.4) +
  facet_grid(pathway_module ~ ., scales = "free", space = "free") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(
    title = "Scenario effects on model-side nitrogen proxies",
    subtitle = "Weighted concentration proxies, not direct elemental nitrogen fluxes",
    x = NULL,
    y = NULL,
    fill = "log2FC"
  ) +
  plot_theme +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
    legend.position = "right"
  )

main_figure <- endpoint_heatmap / proxy_heatmap + plot_layout(heights = c(2.6, 1.5))

ggsave(
  file.path(output_dir, "main_priority_perturbation_response.png"),
  main_figure,
  width = 16,
  height = 13,
  dpi = 300
)

write_tsv(
  endpoint_summary %>%
    transmute(run, source_file, scenario_label, pathway_module, variable, display,
              module, start_value, end_value, baseline_end, log2fc_end),
  file.path(output_dir, "main_priority_endpoint_summary.tsv")
)

write_tsv(
  nitrogen_summary %>%
    transmute(run, source_file, scenario_label, pathway_module, metric,
              start_value, end_value, baseline_end, log2fc_end),
  file.path(output_dir, "main_priority_nitrogen_proxy_summary.tsv")
)

message("Wrote manuscript-style baseline and perturbation figures to: ", normalizePath(output_dir))
