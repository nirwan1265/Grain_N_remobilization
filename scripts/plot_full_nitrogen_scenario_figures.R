#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggbreak)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(forcats)
})

args <- commandArgs(trailingOnly = TRUE)
full_args <- commandArgs(trailingOnly = FALSE)

if (length(args) != 2) {
  stop("Usage: plot_full_nitrogen_scenario_figures.R <input_dir_with_txt_exports> <output_dir>")
}

input_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

script_flag <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(script_flag) > 0) sub("^--file=", "", script_flag[[1]]) else "scripts/plot_full_nitrogen_scenario_figures.R"
root_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
map_path <- file.path(root_dir, "data", "cell_designer_simulations", "grain_n_candidate_model_map.tsv")

plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "white", color = "grey80"),
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
    base == "ask_kinase_x2" ~ "asp_kinase_x2",
    TRUE ~ base
  )
}

simple_run_label <- function(run) {
  dplyr::case_when(
    run == "baseline" ~ "Baseline",
    run == "nh4_half" ~ "External N x0.5",
    run == "nh4_x2" ~ "External N x2",
    run == "nh4_import_half" ~ "NH4 import x0.5",
    run == "nh4_import_x2" ~ "NH4 import x2",
    run == "gs_half" ~ "GS x0.5",
    run == "gs_x2" ~ "GS x2",
    run == "gogat_half" ~ "GOGAT x0.5",
    run == "gogat_x2" ~ "GOGAT x2",
    run == "gs_gogat_half" ~ "GS + GOGAT x0.5",
    run == "gs_gogat_x2" ~ "GS + GOGAT x2",
    run == "gdh_half" ~ "GDH x0.5",
    run == "gdh_x2" ~ "GDH x2",
    run == "asp_kinase_half" ~ "Asp kinase x0.5",
    run == "asp_kinase_x2" ~ "Asp kinase x2",
    run == "asadh_half" ~ "ASADH x0.5",
    run == "asadh_x2" ~ "ASADH x2",
    run == "aspat_capacity_half" ~ "AspAT x0.5",
    run == "aspat_capacity_x2" ~ "AspAT x2",
    run == "akg_pool_half" ~ "2-OG pool x0.5",
    run == "akg_pool_x2" ~ "2-OG pool x2",
    run == "kdh_capacity_half" ~ "KDH drain x0.5",
    run == "kdh_capacity_x2" ~ "KDH drain x2",
    run == "pc_bridge_half" ~ "PC bridge x0.5",
    run == "pc_bridge_x2" ~ "PC bridge x2",
    run == "pepck_bridge_half" ~ "PEPCK bridge x0.5",
    run == "pepck_bridge_x2" ~ "PEPCK bridge x2",
    run == "kbiosyn_half" ~ "2-OG sink x0.5",
    run == "kbiosyn_x2" ~ "2-OG sink x2",
    run == "glu_half" ~ "Glu pool x0.5",
    run == "glu_x2" ~ "Glu pool x2",
    run == "gln_half" ~ "Gln pool x0.5",
    run == "gln_x2" ~ "Gln pool x2",
    run == "asns_half" ~ "ASNS x0.5",
    run == "asns_x2" ~ "ASNS x2",
    run == "omega_half" ~ "OMEGA/NIT2 x0.5",
    run == "omega_x2" ~ "OMEGA/NIT2 x2",
    run == "p5cs_half" ~ "P5CS x0.5",
    run == "p5cs_x2" ~ "P5CS x2",
    run == "p5cr_half" ~ "P5CR x0.5",
    run == "p5cr_x2" ~ "P5CR x2",
    run == "n_sink_both_0p75" ~ "Reduced N demand x0.75",
    run == "n_sink_both_0p5" ~ "Reduced N demand x0.5",
    run == "n_sink_gln_0p75" ~ "Reduced GLN demand x0.75",
    run == "n_sink_glu_0p75" ~ "Reduced GLU demand x0.75",
    run == "n_demand_x1p25" ~ "Increased N demand x1.25",
    run == "asp_pool_half" ~ "Asp pool x0.5",
    run == "asp_pool_x2" ~ "Asp pool x2",
    TRUE ~ run
  )
}

scenario_group <- c(
  "KDH drain x0.5" = "TCA / 2-OG hub",
  "KDH drain x2" = "TCA / 2-OG hub",
  "2-OG pool x0.5" = "TCA / 2-OG hub",
  "2-OG pool x2" = "TCA / 2-OG hub",
  "2-OG sink x0.5" = "TCA / 2-OG hub",
  "2-OG sink x2" = "TCA / 2-OG hub",
  "PC bridge x0.5" = "TCA / 2-OG hub",
  "PC bridge x2" = "TCA / 2-OG hub",
  "PEPCK bridge x0.5" = "TCA / 2-OG hub",
  "PEPCK bridge x2" = "TCA / 2-OG hub",
  "AspAT x0.5" = "C-N interface (AspAT/GOT)",
  "AspAT x2" = "C-N interface (AspAT/GOT)",
  "ASNS x0.5" = "Asp / Asn branch",
  "ASNS x2" = "Asp / Asn branch",
  "GS x0.5" = "Glu/Gln assimilation",
  "GS x2" = "Glu/Gln assimilation",
  "GOGAT x0.5" = "Glu/Gln assimilation",
  "GOGAT x2" = "Glu/Gln assimilation",
  "GS + GOGAT x0.5" = "Glu/Gln assimilation",
  "GS + GOGAT x2" = "Glu/Gln assimilation",
  "GDH x0.5" = "Glu/Gln assimilation",
  "GDH x2" = "Glu/Gln assimilation",
  "Glu pool x0.5" = "Glu/Gln assimilation",
  "Glu pool x2" = "Glu/Gln assimilation",
  "Gln pool x0.5" = "Glu/Gln assimilation",
  "Gln pool x2" = "Glu/Gln assimilation",
  "P5CS x0.5" = "Proline axis",
  "P5CS x2" = "Proline axis",
  "P5CR x0.5" = "Proline axis",
  "P5CR x2" = "Proline axis",
  "Reduced GLU demand x0.75" = "N demand (biosynthetic withdrawal)",
  "Reduced GLN demand x0.75" = "N demand (biosynthetic withdrawal)",
  "Reduced N demand x0.75" = "N demand (biosynthetic withdrawal)",
  "Reduced N demand x0.5" = "N demand (biosynthetic withdrawal)",
  "Increased N demand x1.25" = "N demand (biosynthetic withdrawal)",
  "Asp kinase x0.5" = "Aspartate-family",
  "Asp kinase x2" = "Aspartate-family",
  "ASADH x0.5" = "Aspartate-family",
  "ASADH x2" = "Aspartate-family",
  "Asp pool x0.5" = "Aspartate-family",
  "Asp pool x2" = "Aspartate-family",
  "OMEGA/NIT2 x0.5" = "TCA / 2-OG hub",
  "OMEGA/NIT2 x2" = "TCA / 2-OG hub",
  "External N x0.5" = "N supply / transport",
  "External N x2" = "N supply / transport",
  "NH4 import x0.5" = "N supply / transport",
  "NH4 import x2" = "N supply / transport"
)

scenario_group_levels <- c(
  "TCA / 2-OG hub",
  "C-N interface (AspAT/GOT)",
  "Asp / Asn branch",
  "Glu/Gln assimilation",
  "Proline axis",
  "N demand (biosynthetic withdrawal)",
  "Aspartate-family",
  "N supply / transport"
)

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
all_species <- all_species %>%
  filter(time <= 100)
if (!"baseline" %in% all_species$run) {
  stop("No baseline.txt or base_values.txt file found in input directory")
}

scenario_map <- read_tsv(map_path, show_col_types = FALSE)

preferred_run_order <- c(
  "baseline",
  "nh4_half", "nh4_x2", "nh4_import_half", "nh4_import_x2",
  "asns_half", "asns_x2",
  "gs_half", "gs_x2", "gogat_half", "gogat_x2", "gs_gogat_half", "gs_gogat_x2", "gdh_half", "gdh_x2",
  "omega_half", "omega_x2", "p5cs_half", "p5cs_x2", "p5cr_half", "p5cr_x2",
  "n_sink_glu_0p75", "n_sink_gln_0p75", "n_sink_both_0p75", "n_sink_both_0p5", "n_demand_x1p25",
  "glu_half", "glu_x2", "gln_half", "gln_x2",
  "asp_kinase_half", "asp_kinase_x2",
  "asadh_half", "asadh_x2",
  "asp_pool_half", "asp_pool_x2",
  "aspat_capacity_half", "aspat_capacity_x2",
  "kdh_capacity_half", "kdh_capacity_x2", "pc_bridge_half", "pc_bridge_x2",
  "akg_pool_half", "akg_pool_x2",
  "kbiosyn_half", "kbiosyn_x2"
)

run_display_map <- tibble(run = unique(all_species$run)) %>%
  left_join(select(scenario_map, scenario, candidate_label, pathway_module), by = c("run" = "scenario")) %>%
  mutate(
    scenario_label = simple_run_label(run),
    pathway_module = case_when(
      run == "baseline" ~ "Baseline",
      scenario_label %in% names(scenario_group) ~ unname(scenario_group[scenario_label]),
      TRUE ~ pathway_module
    )
  ) %>%
  mutate(
    pathway_module = factor(
      pathway_module,
      levels = c("Baseline", scenario_group_levels)
    ),
    run_order = match(run, preferred_run_order),
    run_order = if_else(is.na(run_order), length(preferred_run_order) + row_number(), run_order),
    scenario_label_unique = make.unique(scenario_label, sep = " | ")
  ) %>%
  arrange(run_order)

variable_meta <- tribble(
  ~variable,      ~display,        ~module,
  "asp_Asp",      "Asp",           "Asp processes",
  "asn_c",        "Asn",           "Storage / remobilization",
  "asp_ASA",      "ASA",           "Asp processes",
  "asp_Hser",     "Hser",          "Asp processes",
  "asp_PHser",    "PHser",         "Asp processes",
  "asp_Thr",      "Thr",           "Asp processes",
  "asp_Lys",      "Lys",           "Asp processes",
  "asp_Ile",      "Ile",           "Asp processes",
  "nh4_NH4",      "NH4",           "N assimilation",
  "nh4_GLU",      "Glu",           "N assimilation",
  "nh4_GLN",      "Gln",           "N assimilation",
  "pro_c",        "Pro",           "Storage / remobilization",
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
  "asn_c",       2,        "Storage / stress N",
  "pro_c",       1,        "Storage / stress N",
  "asp_Asp",     1,        "Asp-family N",
  "asp_ASA",     1,        "Asp-family N",
  "asp_Hser",    1,        "Asp-family N",
  "asp_PHser",   1,        "Asp-family N",
  "asp_Thr",     1,        "Asp-family N",
  "asp_Ile",     1,        "Asp-family N",
  "asp_Lys",     2,        "Asp-family N"
)

baseline_df <- all_species %>% filter(run == "baseline")
baseline_vars <- c(
  "pyr_AcCoA",
  "tca_akg",
  "nh4_NH4", "nh4_GLU", "nh4_GLN",
  "asp_Asp", "asn_c", "pro_c", "asp_PHser", "asp_Thr", "asp_Lys", "asp_Ile"
)

baseline_plot_df <- baseline_df %>%
  filter(variable %in% baseline_vars) %>%
  left_join(variable_meta, by = "variable") %>%
  mutate(display = factor(display, levels = variable_meta$display[match(baseline_vars, variable_meta$variable)]))

baseline_panel <- ggplot(baseline_plot_df, aes(x = time, y = value)) +
  geom_line(linewidth = 0.8, color = "#2166AC") +
  facet_wrap(~ display, scales = "free_y", ncol = 4) +
  labs(
    title = "Baseline model remains stable with explicit nitrogen support",
    subtitle = "Core C-N bridge metabolites and Asp-family outputs",
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
      "Storage / stress N" = "#66A61E",
      "Downstream product N" = "#E7298A"
    )
  ) +
  labs(
    title = "Baseline nitrogen-proxy pools remain interpretable",
    x = "Time",
    y = "Weighted N proxy",
    color = NULL
  ) +
  plot_theme +
  theme(legend.position = "bottom")

baseline_figure <- baseline_panel / proxy_panel + plot_layout(heights = c(3, 1.3))

ggsave(
  file.path(output_dir, "supp_baseline_with_nitrogen.png"),
  baseline_figure,
  width = 16,
  height = 13,
  dpi = 300
)

focus_vars <- c("asp_Asp", "asn_c", "pro_c", "asp_Thr", "asp_Lys", "asp_Ile", "nh4_NH4", "nh4_GLU", "nh4_GLN", "tca_akg")
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
  mutate(
    display = factor(display, levels = variable_meta$display[match(focus_vars, variable_meta$variable)]),
    run_display = factor(scenario_label_unique, levels = rev(run_display_map$scenario_label_unique))
  )

endpoint_heatmap_df <- endpoint_summary %>%
  filter(run != "baseline")

endpoint_heatmap <- ggplot(endpoint_heatmap_df, aes(x = display, y = run_display, fill = log2fc_end)) +
  geom_tile(color = "white", linewidth = 0.4) +
  facet_grid(pathway_module ~ module, scales = "free", space = "free", switch = "y") +
  scale_fill_gradientn(
    colours = c("#2166AC", "#92C5DE", "#D9D9D9", "#FDAE61", "#B2182B"),
    values = scales::rescale(c(-1.5, -0.4, 0, 0.4, 1.5)),
    limits = c(-1.5, 1.5),
    oob = scales::squish,
    breaks = c(-1, 0, 1),
    labels = c("Low", "Baseline", "High"),
    guide = guide_colorbar(
      direction = "horizontal",
      barheight = grid::unit(7, "mm"),
      barwidth = grid::unit(125, "mm"),
      ticks = FALSE
    )
  ) +
  labs(
    title = "Updated scenario panel shows distinct carbon-nitrogen responses",
    subtitle = "Endpoint response relative to baseline; near-baseline scenarios appear pale yellow",
    x = NULL,
    y = NULL,
    fill = "Direction"
  ) +
  plot_theme +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 10),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, size = 11, face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 11, lineheight = 0.9),
    legend.position = "bottom"
  )

ggsave(
  file.path(output_dir, "all_scenarios_endpoint_heatmap.png"),
  endpoint_heatmap,
  width = 18,
  height = 10,
  dpi = 300
)

focused_runs <- c(
  "nh4_half", "nh4_x2",
  "nh4_import_half", "nh4_import_x2",
  "gs_half", "gs_x2", "gogat_half", "gogat_x2", "gs_gogat_half", "gs_gogat_x2",
  "gdh_half", "gdh_x2",
  "asns_half", "asns_x2",
  "omega_half", "omega_x2",
  "p5cs_half", "p5cs_x2", "p5cr_half", "p5cr_x2",
  "n_sink_glu_0p75", "n_sink_gln_0p75",
  "n_sink_both_0p75", "n_sink_both_0p5",
  "n_demand_x1p25",
  "glu_half", "glu_x2", "gln_half", "gln_x2",
  "asp_kinase_half", "asp_kinase_x2",
  "asp_pool_half", "asp_pool_x2",
  "aspat_capacity_half", "aspat_capacity_x2",
  "asadh_half", "asadh_x2",
  "kdh_capacity_half", "kdh_capacity_x2",
  "pc_bridge_half", "pc_bridge_x2"
)

focused_endpoint_df <- endpoint_summary %>%
  filter(run %in% focused_runs) %>%
  mutate(
    run_display = factor(
      scenario_label_unique,
      levels = rev(run_display_map$scenario_label_unique[run_display_map$run %in% focused_runs])
    )
  )

focused_heatmap <- ggplot(focused_endpoint_df, aes(x = display, y = run_display, fill = log2fc_end)) +
  geom_tile(color = "white", linewidth = 0.4) +
  facet_grid(pathway_module ~ module, scales = "free", space = "free", switch = "y") +
  scale_fill_gradientn(
    colours = c("#2166AC", "#92C5DE", "#D9D9D9", "#FDAE61", "#B2182B"),
    values = scales::rescale(c(-1.5, -0.4, 0, 0.4, 1.5)),
    limits = c(-1.5, 1.5),
    oob = scales::squish,
    breaks = c(-1, 0, 1),
    labels = c("Low", "Baseline", "High"),
    guide = guide_colorbar(
      direction = "horizontal",
      barheight = grid::unit(7, "mm"),
      barwidth = grid::unit(125, "mm"),
      ticks = FALSE
    )
  ) +
  labs(
    title = "Focused perturbations summarize the full carbon-nitrogen scaffold",
    subtitle = "Shown scenarios cover N supply, N assimilation, Asp/Asn branching, Pro remodeling, C-N interface, and contrasting KDH-driven 2-OG/TCA competition",
    x = NULL,
    y = NULL,
    fill = "Direction"
  ) +
  plot_theme +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 11),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, size = 11, face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 11, lineheight = 0.9),
    legend.position = "bottom"
  )

ggsave(
  file.path(output_dir, "focused_relevant_endpoint_heatmap.png"),
  focused_heatmap,
  width = 16,
  height = 8,
  dpi = 300
)

focused_timecourse_runs <- c(
  "baseline",
  "nh4_half", "nh4_x2",
  "gs_half", "gs_x2", "gogat_half", "gogat_x2",
  "gs_gogat_half", "gs_gogat_x2",
  "n_sink_both_0p75", "n_sink_both_0p5",
  "asp_kinase_x2",
  "aspat_capacity_half", "aspat_capacity_x2",
  "asadh_half", "asadh_x2",
  "kdh_capacity_half", "kdh_capacity_x2"
)
focused_timecourse_labels <- run_display_map %>%
  filter(run %in% focused_timecourse_runs) %>%
  select(run, scenario_label)

focused_timecourse_vars <- c("nh4_NH4", "nh4_GLU", "nh4_GLN", "tca_akg", "asp_Asp", "asp_Thr")
focused_timecourse_df <- all_species %>%
  filter(run %in% focused_timecourse_runs, variable %in% focused_timecourse_vars) %>%
  left_join(variable_meta, by = "variable") %>%
  left_join(focused_timecourse_labels, by = "run") %>%
  mutate(
    display = factor(display, levels = variable_meta$display[match(focused_timecourse_vars, variable_meta$variable)]),
    scenario_label = factor(scenario_label, levels = focused_timecourse_labels$scenario_label)
  )

focused_palette <- c(
  "Baseline" = "#1B1B1B",
  "External N x0.5" = "#D95F02",
  "External N x2" = "#1B9E77",
  "GS x0.5" = "#E6AB02",
  "GS x2" = "#A6761D",
  "GOGAT x0.5" = "#666666",
  "GOGAT x2" = "#4DAF4A",
  "GS + GOGAT x0.5" = "#66A61E",
  "GS + GOGAT x2" = "#A6D854",
  "Reduced N demand x0.75" = "#A6761D",
  "Reduced N demand x0.5" = "#E6AB02",
  "Asp kinase x2" = "#7570B3",
  "AspAT x0.5" = "#FB9A99",
  "AspAT x2" = "#6A3D9A",
  "ASADH x0.5" = "#E7298A",
  "ASADH x2" = "#B15928",
  "KDH drain x0.5" = "#1F78B4",
  "KDH drain x2" = "#A6CEE3"
)

make_focused_panel <- function(df, var_id, panel_title, break_range = NULL) {
  p <- df %>%
    filter(variable == var_id) %>%
    ggplot(aes(x = time, y = value, color = scenario_label)) +
    geom_line(linewidth = 0.85) +
    scale_color_manual(values = focused_palette, drop = FALSE) +
    labs(title = panel_title, x = "Time", y = NULL, color = NULL) +
    plot_theme +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "none"
    )

  if (!is.null(break_range)) {
    p <- p + ggbreak::scale_y_break(break_range, scales = 0.8)
  }

  p
}

nh4_focus_panel <- make_focused_panel(focused_timecourse_df, "nh4_NH4", "NH4")
glu_focus_panel <- make_focused_panel(focused_timecourse_df, "nh4_GLU", "Glu", c(0.02, 0.95))
gln_focus_panel <- make_focused_panel(focused_timecourse_df, "nh4_GLN", "Gln", c(0.002, 0.95))
akg_focus_panel <- make_focused_panel(focused_timecourse_df, "tca_akg", "2-OG")
asp_focus_panel <- make_focused_panel(focused_timecourse_df, "asp_Asp", "Asp")
thr_focus_panel <- make_focused_panel(focused_timecourse_df, "asp_Thr", "Thr")

legend_df <- tibble(
  time = 0,
  value = 0,
  scenario_label = factor(names(focused_palette), levels = names(focused_palette))
)

legend_plot <- ggplot(legend_df, aes(x = time, y = value, color = scenario_label)) +
  geom_line(linewidth = 1, alpha = 1) +
  scale_color_manual(values = focused_palette, drop = FALSE) +
  labs(color = NULL) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 11),
    plot.margin = margin(0, 0, 0, 0)
  )

focused_timecourse_plot <- (
  nh4_focus_panel | glu_focus_panel | gln_focus_panel
) / (
  akg_focus_panel | asp_focus_panel | thr_focus_panel
) / (
  legend_plot
) +
  plot_annotation(
    title = "Focused transient time courses at t = 100",
    subtitle = "Selected nitrogen, Asp-branch, C-N interface, and 2-OG/TCA scenarios"
  ) +
  plot_layout(heights = c(1, 1, 0.22))

ggsave(
  file.path(output_dir, "focused_relevant_timecourses.png"),
  focused_timecourse_plot,
  width = 16,
  height = 10,
  dpi = 300
)

write_tsv(endpoint_summary, file.path(output_dir, "all_scenarios_endpoint_summary.tsv"))

message("Wrote figures and summaries to: ", normalizePath(output_dir))
