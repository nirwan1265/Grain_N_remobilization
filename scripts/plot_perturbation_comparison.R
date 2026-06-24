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

if (length(args) < 3) {
  stop("Usage: plot_perturbation_comparison.R <output_dir> <label=txt> <label=txt> [...]")
}

output_dir <- args[[1]]
inputs <- args[-1]
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

load_run <- function(spec) {
  parts <- str_split_fixed(spec, "=", 2)
  label <- parts[, 1]
  file <- parts[, 2]
  lines <- readLines(file, warn = FALSE)
  species_df <- read_section(lines, "#species") %>%
    pivot_longer(-time, names_to = "variable", values_to = "value") %>%
    mutate(run = label)
  species_df
}

all_species <- map_dfr(inputs, load_run)

focus <- c(
  "asp_Asp", "asp_Thr", "asp_Lys", "asp_Ile", "asp_PHser",
  "nh4_GLU", "nh4_GLN", "tca_akg", "tca_cit", "tca_mal", "tca_oaa",
  "pyr_AcCoA", "pyr_pyruvate"
)

plot_df <- all_species %>%
  filter(variable %in% focus)

p <- ggplot(plot_df, aes(x = time, y = value, color = run)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ variable, scales = "free_y", ncol = 4) +
  labs(
    title = "Perturbation Comparison",
    x = "Time",
    y = "Value",
    color = "Run"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "#e8f1fb"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(output_dir, "perturbation_comparison.png"), p, width = 15, height = 10, dpi = 300)

summary_df <- plot_df %>%
  group_by(run, variable) %>%
  summarise(
    start_value = first(value),
    end_value = last(value),
    min_value = min(value, na.rm = TRUE),
    max_value = max(value, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(summary_df, file.path(output_dir, "perturbation_comparison_summary.tsv"))
message("Wrote comparison outputs to: ", normalizePath(output_dir))
