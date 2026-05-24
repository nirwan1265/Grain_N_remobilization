################################################################################
### FIGURE 5: SPLINE FST, DELTA TAJIMA'S D, AND NUCLEOTIDE DIVERSITY
################################################################################

library(tidyverse)
library(ggplot2)
library(patchwork)
library(scales)

################################################################################
### CONFIGURATION
################################################################################

ic_windows_file <- "tables/supplementary/SuppTable_IC_GenWin_windows_Fst_Tajima_pi_genes.csv"
j_windows_file <- "tables/supplementary/SuppTable_J_GenWin_windows_Fst_Tajima_pi_genes.csv"

main_dir <- "Figs/main"
main_file <- file.path(main_dir, "Fig5.png")

dir.create(main_dir, showWarnings = FALSE, recursive = TRUE)

min_display_window_bp <- 12000000
candidate_gap_bp <- 500000

plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 10)
    ),
    plot.tag = element_text(size = 24, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, face = "bold", color = "black"),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    legend.background = element_rect(fill = "white", color = "grey70", linewidth = 0.4),
    legend.title = element_blank(),
    legend.text = element_text(size = 13),
    plot.margin = margin(15, 15, 15, 15)
  )

################################################################################
### HELPERS
################################################################################

check_required_files <- function(paths) {
  missing <- paths[!file.exists(paths)]

  if (length(missing) > 0) {
    stop(
      "Missing required spline workflow files:\n  ",
      paste(missing, collapse = "\n  "),
      "\nRun scripts/16_spling_Fst.R and scripts/17_spline_Fst_Tajima_Pi.R first."
    )
  }
}

normalize_flag <- function(x) {
  if (is.logical(x)) {
    return(replace_na(x, FALSE))
  }

  tolower(as.character(x)) %in% c("true", "yes", "1")
}

running_median <- function(x, k = 5) {
  n <- length(x)
  out <- rep(NA_real_, n)
  half_k <- floor(k / 2)

  for (i in seq_len(n)) {
    left <- max(1, i - half_k)
    right <- min(n, i + half_k)
    out[i] <- median(x[left:right], na.rm = TRUE)
  }

  out
}

merge_candidate_regions <- function(df, candidate_gap_bp) {
  candidate_df <- df %>%
    filter(all_three_signals)

  if (nrow(candidate_df) == 0) {
    return(tibble(
      chr = character(),
      chr_num = integer(),
      region_id = integer(),
      n_windows = integer(),
      region_start = numeric(),
      region_end = numeric(),
      region_center = numeric(),
      peak_Fst = numeric(),
      min_deltaD = numeric(),
      min_log_pi_ratio = numeric(),
      region_length_bp = numeric()
    ))
  }

  candidate_df %>%
    arrange(chr_num, window_start_bp, window_end_bp) %>%
    mutate(
      new_region = if_else(
        row_number() == 1L,
        TRUE,
        chr_num != lag(chr_num) | (window_start_bp - lag(window_end_bp)) > candidate_gap_bp
      ),
      region_id = cumsum(new_region)
    ) %>%
    group_by(chr, chr_num, region_id) %>%
    summarise(
      n_windows = n(),
      region_start = min(window_start_bp, na.rm = TRUE),
      region_end = max(window_end_bp, na.rm = TRUE),
      region_center = mean(c(region_start, region_end)),
      peak_Fst = max(mean_fst, na.rm = TRUE),
      min_deltaD = min(deltaD, na.rm = TRUE),
      min_log_pi_ratio = min(log_pi_ratio, na.rm = TRUE),
      region_length_bp = region_end - region_start,
      .groups = "drop"
    ) %>%
    arrange(chr_num, region_start)
}

add_cumulative_positions <- function(df, chr_info) {
  df %>%
    left_join(chr_info %>% select(chr_num, offset, center), by = "chr_num") %>%
    mutate(pos_cum = window_mid_bp + offset)
}

add_region_positions <- function(regions_df, chr_info) {
  if (nrow(regions_df) == 0) {
    return(regions_df)
  }

  regions_df %>%
    left_join(chr_info %>% select(chr_num, offset), by = "chr_num") %>%
    mutate(
      region_start_cum = region_start + offset,
      region_end_cum = region_end + offset,
      region_center_cum = region_center + offset,
      display_width_bp = pmax(region_end - region_start, min_display_window_bp),
      display_start_cum = region_center_cum - display_width_bp / 2,
      display_end_cum = region_center_cum + display_width_bp / 2
    )
}

plot_population_panel <- function(df, regions_df, chr_info, pop_label) {
  df <- add_cumulative_positions(df, chr_info)
  regions_df <- add_region_positions(regions_df, chr_info)

  thresholds <- list(
    fst = quantile(df$mean_fst, 0.95, na.rm = TRUE),
    deltaD = quantile(df$deltaD, 0.05, na.rm = TRUE),
    log_pi = quantile(df$log_pi_ratio, 0.05, na.rm = TRUE)
  )

  shade_fill <- "#F4B942"
  shade_edge <- "#9A3412"
  outlier_color <- "#B22222"
  fst_color <- "#163A5F"
  delta_color <- "#8E3B8A"
  pi_color <- "#127475"

  base_scale <- scale_x_continuous(
    breaks = chr_info$center,
    labels = chr_info$chr_num,
    expand = expansion(mult = c(0.005, 0.005))
  )

  p_fst <- ggplot(df, aes(x = pos_cum)) +
    {if (nrow(regions_df) > 0)
      geom_rect(
        data = regions_df,
        aes(
          xmin = display_start_cum,
          xmax = display_end_cum,
          ymin = -Inf,
          ymax = Inf
        ),
        fill = shade_fill,
        color = shade_edge,
        linewidth = 0.8,
        alpha = 0.42,
        inherit.aes = FALSE
      )
    } +
    geom_point(aes(y = mean_fst), color = "grey82", size = 0.8, alpha = 0.5) +
    geom_line(aes(y = Fst_smooth, group = chr_num), color = fst_color, linewidth = 1.05) +
    geom_point(
      data = df %>% filter(fst_outlier),
      aes(y = mean_fst),
      color = outlier_color,
      size = 1.3,
      alpha = 0.95
    ) +
    geom_hline(
      yintercept = thresholds$fst,
      linetype = "dashed",
      color = outlier_color,
      linewidth = 0.55
    ) +
    labs(title = pop_label, x = NULL, y = expression(F[ST])) +
    base_scale +
    plot_theme +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(10, 15, 0, 15)
    )

  p_delta <- ggplot(df, aes(x = pos_cum)) +
    {if (nrow(regions_df) > 0)
      geom_rect(
        data = regions_df,
        aes(
          xmin = display_start_cum,
          xmax = display_end_cum,
          ymin = -Inf,
          ymax = Inf
        ),
        fill = shade_fill,
        color = shade_edge,
        linewidth = 0.8,
        alpha = 0.42,
        inherit.aes = FALSE
      )
    } +
    geom_point(aes(y = deltaD), color = "grey82", size = 0.8, alpha = 0.5) +
    geom_line(aes(y = deltaD_smooth, group = chr_num), color = delta_color, linewidth = 1.05) +
    geom_point(
      data = df %>% filter(deltaD_outlier),
      aes(y = deltaD),
      color = outlier_color,
      size = 1.3,
      alpha = 0.95
    ) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    geom_hline(
      yintercept = thresholds$deltaD,
      linetype = "dashed",
      color = outlier_color,
      linewidth = 0.55
    ) +
    labs(x = NULL, y = expression(Delta * "D")) +
    base_scale +
    plot_theme +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(0, 15, 0, 15)
    )

  p_pi <- ggplot(df, aes(x = pos_cum)) +
    {if (nrow(regions_df) > 0)
      geom_rect(
        data = regions_df,
        aes(
          xmin = display_start_cum,
          xmax = display_end_cum,
          ymin = -Inf,
          ymax = Inf
        ),
        fill = shade_fill,
        color = shade_edge,
        linewidth = 0.8,
        alpha = 0.42,
        inherit.aes = FALSE
      )
    } +
    geom_point(aes(y = log_pi_ratio), color = "grey82", size = 0.8, alpha = 0.5) +
    geom_line(aes(y = log_pi_smooth, group = chr_num), color = pi_color, linewidth = 1.05) +
    geom_point(
      data = df %>% filter(pi_outlier),
      aes(y = log_pi_ratio),
      color = outlier_color,
      size = 1.3,
      alpha = 0.95
    ) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    geom_hline(
      yintercept = thresholds$log_pi,
      linetype = "dashed",
      color = outlier_color,
      linewidth = 0.55
    ) +
    labs(
      x = "Chromosome",
      y = expression(log[2](pi[14] / pi[0]))
    ) +
    base_scale +
    plot_theme +
    theme(plot.margin = margin(0, 15, 10, 15))

  p_fst / p_delta / p_pi +
    plot_layout(heights = c(1, 1, 1))
}

################################################################################
### LOAD DATA
################################################################################

cat("\n=== Building Fig5 spline selection-scan panels ===\n")

check_required_files(c(
  ic_windows_file,
  j_windows_file
))

ic_df <- readr::read_csv(ic_windows_file, show_col_types = FALSE) %>%
  mutate(
    fst_outlier = normalize_flag(fst_outlier),
    deltaD_outlier = normalize_flag(deltaD_outlier),
    pi_outlier = normalize_flag(pi_outlier),
    all_three_signals = normalize_flag(all_three_signals)
  ) %>%
  arrange(chr_num, window_id) %>%
  group_by(chr_num) %>%
  mutate(
    Fst_smooth = running_median(mean_fst, k = 5),
    deltaD_smooth = running_median(deltaD, k = 5),
    log_pi_smooth = running_median(log_pi_ratio, k = 5)
  ) %>%
  ungroup()

j_df <- readr::read_csv(j_windows_file, show_col_types = FALSE) %>%
  mutate(
    fst_outlier = normalize_flag(fst_outlier),
    deltaD_outlier = normalize_flag(deltaD_outlier),
    pi_outlier = normalize_flag(pi_outlier),
    all_three_signals = normalize_flag(all_three_signals)
  ) %>%
  arrange(chr_num, window_id) %>%
  group_by(chr_num) %>%
  mutate(
    Fst_smooth = running_median(mean_fst, k = 5),
    deltaD_smooth = running_median(deltaD, k = 5),
    log_pi_smooth = running_median(log_pi_ratio, k = 5)
  ) %>%
  ungroup()

ic_regions <- merge_candidate_regions(ic_df, candidate_gap_bp = candidate_gap_bp)
j_regions <- merge_candidate_regions(j_df, candidate_gap_bp = candidate_gap_bp)

combined_df <- bind_rows(ic_df, j_df)

chr_info <- combined_df %>%
  group_by(chr_num) %>%
  summarise(chr_len = max(window_end_bp, na.rm = TRUE), .groups = "drop") %>%
  arrange(chr_num) %>%
  mutate(
    offset = lag(cumsum(chr_len), default = 0),
    center = offset + chr_len / 2
  )

################################################################################
### BUILD FIGURE
################################################################################

panel_ic <- plot_population_panel(
  ic_df,
  ic_regions,
  chr_info,
  "Indian Chief selection scan: FST, Delta D, log2(pi14/pi0)"
)

panel_j <- plot_population_panel(
  j_df,
  j_regions,
  chr_info,
  "Jarvis selection scan: FST, Delta D, log2(pi14/pi0)"
)

fig5 <- (wrap_elements(panel = panel_ic) + wrap_elements(panel = panel_j)) +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A")

ggsave(main_file, fig5, width = 18, height = 10.5, dpi = 300, bg = "white")

cat("\nSaved main figure to:\n")
cat("  ", main_file, "\n", sep = "")
