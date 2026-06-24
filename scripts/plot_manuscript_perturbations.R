#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(forcats)
})

args <- commandArgs(trailingOnly = TRUE)
full_args <- commandArgs(trailingOnly = FALSE)

if (length(args) != 2) {
  stop("Usage: plot_manuscript_perturbations.R <input_dir_with_txt_exports> <output_dir>")
}

input_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

script_flag <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(script_flag) > 0) sub("^--file=", "", script_flag[[1]]) else "scripts/plot_manuscript_perturbations.R"
root_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
map_path <- file.path(root_dir, "data", "cell_designer_simulations", "grain_n_candidate_model_map.tsv")

`%||%` <- function(x, y) if (is.null(x)) y else x

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
run_manifest <- tibble(
  run = map_chr(txt_files, parse_run_label),
  source_file = basename(txt_files)
) %>%
  distinct()

if (!"baseline" %in% run_manifest$run) {
  stop("No baseline export found. Name the baseline file base_values.txt or baseline.txt")
}

scenario_map <- read_tsv(map_path, show_col_types = FALSE)

run_order <- c("baseline", scenario_map$scenario)
run_manifest <- run_manifest %>%
  mutate(run = factor(run, levels = unique(run_order))) %>%
  arrange(run)

variable_meta <- tribble(
  ~variable,      ~display,        ~module,
  "asp_Asp",      "Asp",           "Asp-family",
  "asp_ASA",      "ASA",           "Asp-family",
  "asp_Hser",     "Hser",          "Asp-family",
  "asp_PHser",    "PHser",         "Asp-family",
  "asp_Thr",      "Thr",           "Asp-family",
  "asp_Lys",      "Lys",           "Asp-family",
  "asp_Ile",      "Ile",           "Asp-family",
  "nh4_NH4",      "NH4",           "N assimilation",
  "nh4_GLU",      "Glu",           "N assimilation",
  "nh4_GLN",      "Gln",           "N assimilation",
  "tca_akg",      "2-OG",          "TCA / 2-OG",
  "tca_cit",      "Citrate",       "TCA / 2-OG",
  "tca_mal",      "Malate",        "TCA / 2-OG",
  "tca_oaa",      "OAA",           "TCA / 2-OG",
  "tca_suc",      "Succinate",     "TCA / 2-OG",
  "pyr_AcCoA",    "AcCoA",         "Pyruvate / entry",
  "pyr_pyruvate", "Pyruvate",      "Pyruvate / entry"
)

focus_vars <- variable_meta$variable
eps <- 1e-12

endpoint_summary <- all_species %>%
  filter(variable %in% focus_vars) %>%
  arrange(run, variable, time) %>%
  group_by(run, source_file, variable) %>%
  summarise(
    start_value = first(value),
    end_value = last(value),
    min_value = min(value, na.rm = TRUE),
    max_value = max(value, na.rm = TRUE),
    tail_mean = mean(tail(value, pmax(3, floor(length(value) * 0.1))), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(variable_meta, by = "variable")

baseline_end <- endpoint_summary %>%
  filter(run == "baseline") %>%
  select(variable, baseline_end = end_value, baseline_tail = tail_mean)

endpoint_summary <- endpoint_summary %>%
  left_join(baseline_end, by = "variable") %>%
  mutate(
    delta_end = end_value - baseline_end,
    log2fc_end = log2((end_value + eps) / (baseline_end + eps))
  ) %>%
  left_join(scenario_map, by = c("run" = "scenario")) %>%
  mutate(
    scenario_label = if_else(as.character(run) == "baseline", "baseline", coalesce(candidate_label, as.character(run))),
    pathway_module = coalesce(pathway_module, "Baseline")
  )

heatmap_df <- endpoint_summary %>%
  mutate(
    run_chr = as.character(run),
    run_chr = if_else(run_chr == "baseline", "baseline", run_chr),
    run_chr = factor(run_chr, levels = rev(unique(c("baseline", scenario_map$scenario))))
  )

heatmap_plot <- ggplot(heatmap_df, aes(x = display, y = run_chr, fill = log2fc_end)) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_grid(pathway_module ~ module, scales = "free", space = "free") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
  labs(
    title = "Scenario endpoint response relative to baseline",
    subtitle = "Tile color shows log2 fold-change in final value relative to baseline",
    x = NULL,
    y = NULL,
    fill = "log2FC"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "#e8f1fb"),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(output_dir, "manuscript_endpoint_heatmap.png"),
  heatmap_plot,
  width = 15,
  height = 9,
  dpi = 300
)

n_weights <- tribble(
  ~variable,     ~n_atoms, ~proxy_group,
  "nh4_NH4",     1,        "Inorganic N",
  "nh4_GLU",     1,        "Glu/Gln pool",
  "nh4_GLN",     2,        "Glu/Gln pool",
  "asp_Asp",     1,        "Asp-family pool",
  "asp_ASA",     1,        "Asp-family pool",
  "asp_Hser",    1,        "Asp-family pool",
  "asp_PHser",   1,        "Asp-family pool",
  "asp_Thr",     1,        "Asp-family pool",
  "asp_Ile",     1,        "Asp-family pool",
  "asp_Lys",     2,        "Asp-family pool"
)

nitrogen_long <- all_species %>%
  inner_join(n_weights, by = "variable") %>%
  mutate(weighted_n = value * n_atoms)

proxy_timeseries <- bind_rows(
  nitrogen_long %>%
    group_by(run, source_file, time) %>%
    summarise(metric = "modeled_total_N_proxy", value = sum(weighted_n, na.rm = TRUE), .groups = "drop"),
  nitrogen_long %>%
    filter(proxy_group == "Inorganic N") %>%
    group_by(run, source_file, time) %>%
    summarise(metric = "inorganic_N_proxy", value = sum(weighted_n, na.rm = TRUE), .groups = "drop"),
  nitrogen_long %>%
    filter(proxy_group == "Glu/Gln pool") %>%
    group_by(run, source_file, time) %>%
    summarise(metric = "glugln_N_proxy", value = sum(weighted_n, na.rm = TRUE), .groups = "drop"),
  nitrogen_long %>%
    filter(proxy_group == "Asp-family pool") %>%
    group_by(run, source_file, time) %>%
    summarise(metric = "asp_family_N_proxy", value = sum(weighted_n, na.rm = TRUE), .groups = "drop"),
  all_species %>%
    filter(variable %in% c("asp_Thr", "asp_Ile", "asp_Lys")) %>%
    left_join(select(n_weights, variable, n_atoms), by = "variable") %>%
    mutate(weighted_n = value * n_atoms) %>%
    group_by(run, source_file, time) %>%
    summarise(metric = "downstream_product_N_proxy", value = sum(weighted_n, na.rm = TRUE), .groups = "drop")
) %>%
  distinct()

nitrogen_summary <- proxy_timeseries %>%
  arrange(run, metric, time) %>%
  group_by(run, metric) %>%
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
  mutate(
    delta_end = end_value - baseline_end,
    log2fc_end = log2((end_value + eps) / (baseline_end + eps))
  ) %>%
  left_join(scenario_map, by = c("run" = "scenario")) %>%
  mutate(pathway_module = coalesce(pathway_module, "Baseline"))

proxy_plot <- nitrogen_summary %>%
  mutate(run = factor(as.character(run), levels = rev(unique(c("baseline", scenario_map$scenario))))) %>%
  ggplot(aes(x = metric, y = run, fill = log2fc_end)) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_grid(pathway_module ~ ., scales = "free", space = "free") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
  labs(
    title = "Nitrogen proxy response relative to baseline",
    subtitle = "These are model-side concentration proxies, not direct elemental N flux measurements",
    x = NULL,
    y = NULL,
    fill = "log2FC"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "#eef5e2"),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(output_dir, "manuscript_nitrogen_proxy_heatmap.png"),
  proxy_plot,
  width = 12,
  height = 8,
  dpi = 300
)

timecourse_vars <- c("asp_Asp", "asp_Thr", "asp_Lys", "nh4_GLU", "nh4_GLN", "tca_akg", "tca_oaa", "pyr_AcCoA")
max_runs_for_lines <- 8
selected_runs <- c("baseline", intersect(scenario_map$scenario, unique(as.character(run_manifest$run))))
selected_runs <- selected_runs[seq_len(min(length(selected_runs), max_runs_for_lines))]

timecourse_df <- all_species %>%
  filter(variable %in% timecourse_vars, run %in% selected_runs) %>%
  left_join(variable_meta, by = "variable") %>%
  mutate(run = factor(run, levels = selected_runs))

timecourse_plot <- ggplot(timecourse_df, aes(x = time, y = value, color = run)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ display, scales = "free_y", ncol = 4) +
  labs(
    title = "Core time courses for baseline and leading perturbations",
    x = "Time",
    y = "Value",
    color = "Run"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "#f3edf8"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(output_dir, "manuscript_core_timecourses.png"),
  timecourse_plot,
  width = 14,
  height = 10,
  dpi = 300
)

write_tsv(
  endpoint_summary %>%
    select(run, source_file = source_file, pathway_module, candidate_gene_id, candidate_label,
           variable, display, module, start_value, end_value, tail_mean, delta_end, log2fc_end),
  file.path(output_dir, "manuscript_endpoint_summary.tsv")
)

write_tsv(
  nitrogen_summary %>%
    select(run, pathway_module, candidate_gene_id, candidate_label,
           metric, start_value, end_value, delta_end, log2fc_end),
  file.path(output_dir, "manuscript_nitrogen_proxy_summary.tsv")
)

write_tsv(run_manifest %>% mutate(run = as.character(run)), file.path(output_dir, "scenario_file_manifest.tsv"))

message("Wrote manuscript perturbation outputs to: ", normalizePath(output_dir))
