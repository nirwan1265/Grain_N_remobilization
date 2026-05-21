################################################################################
### FIGURE 4: FST, DELTA TAJIMA'S D, AND NUCLEOTIDE DIVERSITY
################################################################################

library(tidyverse)
library(ggplot2)
library(patchwork)
library(scales)

################################################################################
### CONFIGURATION
################################################################################

fst_dir_ic <- "/Users/nirwantandukar/Documents/Research/results/Indian_Jarvis/ANGSD_Fst/Indian_chief/sliding_window"
fst_dir_j <- "/Users/nirwantandukar/Documents/Research/results/Indian_Jarvis/ANGSD_Fst/Jarvis/sliding_window"
tajima_dir_ic <- "/Users/nirwantandukar/Documents/Research/results/Indian_Jarvis/ANGSD_TajimaD/Indian_Chief/sliding_window"
tajima_dir_j <- "/Users/nirwantandukar/Documents/Research/results/Indian_Jarvis/ANGSD_TajimaD/Jarvis/sliding_window"

main_dir <- "Figs/main"

main_file <- file.path(main_dir, "Fig4.png")

dir.create(main_dir, showWarnings = FALSE, recursive = TRUE)

window_size_bp <- 250000
candidate_gap_bp <- 500000
min_display_window_bp <- 12000000
fst_percentile <- 0.95
tail_percentile <- 0.05

chr_map <- tibble(
  nc_id = c(
    "NC_050096.1", "NC_050097.1", "NC_050098.1", "NC_050099.1", "NC_050100.1",
    "NC_050101.1", "NC_050102.1", "NC_050103.1", "NC_050104.1", "NC_050105.1"
  ),
  chr = paste0("chr", 1:10),
  chr_num = 1:10
)

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

read_fst_files <- function(fst_dir, population_name) {
  files <- list.files(fst_dir, pattern = "\\.txt$", full.names = TRUE)

  if (length(files) == 0) {
    stop("No FST files found in: ", fst_dir)
  }

  purrr::map_dfr(files, function(path) {
    readr::read_tsv(
      path,
      show_col_types = FALSE,
      col_names = c("region", "Chr", "WinCenter", "Nsites", "Fst"),
      skip = 1
    ) %>%
      left_join(chr_map, by = c("Chr" = "nc_id")) %>%
      filter(!is.na(chr)) %>%
      transmute(
        Population = population_name,
        chr,
        chr_num,
        WinCenter = as.numeric(WinCenter),
        Nsites_fst = as.numeric(Nsites),
        Fst = as.numeric(Fst)
      )
  }) %>%
    arrange(chr_num, WinCenter)
}

read_tajima_files <- function(tajima_dir, generation_pattern, generation_name) {
  files <- list.files(tajima_dir, pattern = generation_pattern, full.names = TRUE)

  if (length(files) == 0) {
    stop("No Tajima's D files found in: ", tajima_dir, " for pattern ", generation_pattern)
  }

  purrr::map_dfr(files, function(path) {
    readr::read_tsv(
      path,
      comment = "#",
      show_col_types = FALSE,
      col_names = c(
        "index", "Chr", "WinCenter", "tW", "tP", "tF", "tH", "tL",
        "Tajima", "fuf", "fud", "fayh", "zeng", "nSites"
      )
    ) %>%
      left_join(chr_map, by = c("Chr" = "nc_id")) %>%
      filter(!is.na(chr)) %>%
      transmute(
        Generation = generation_name,
        chr,
        chr_num,
        WinCenter = as.numeric(WinCenter),
        Tajima = as.numeric(Tajima),
        tP = as.numeric(tP),
        nSites = as.numeric(nSites)
      )
  }) %>%
    arrange(chr_num, WinCenter)
}

merge_candidate_regions <- function(df, window_size_bp, candidate_gap_bp) {
  if (nrow(df) == 0) {
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
      min_log_pi_ratio = numeric()
    ))
  }

  half_window <- window_size_bp / 2

  df %>%
    arrange(chr_num, WinCenter) %>%
    mutate(
      new_region = dplyr::if_else(
        row_number() == 1,
        TRUE,
        chr_num != lag(chr_num) | (WinCenter - lag(WinCenter)) > candidate_gap_bp
      ),
      region_id = cumsum(new_region)
    ) %>%
    group_by(chr, chr_num, region_id) %>%
    summarise(
      n_windows = n(),
      region_start = min(WinCenter) - half_window,
      region_end = max(WinCenter) + half_window,
      region_center = mean(c(region_start, region_end)),
      peak_Fst = max(Fst, na.rm = TRUE),
      min_deltaD = min(deltaD, na.rm = TRUE),
      min_log_pi_ratio = min(log_pi_ratio, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      region_start = pmax(1, region_start),
      region_length_bp = region_end - region_start
    ) %>%
    arrange(chr_num, region_start)
}

prepare_population_data <- function(fst_dir, tajima_dir, gen0_pattern, gen14_pattern, population_name) {
  fst_df <- read_fst_files(fst_dir, population_name)

  tajima_gen0 <- read_tajima_files(tajima_dir, gen0_pattern, "Gen0") %>%
    transmute(
      chr,
      chr_num,
      WinCenter,
      Tajima_Gen0 = Tajima,
      pi_Gen0 = tP / nSites,
      nSites_Gen0 = nSites
    )

  tajima_gen14 <- read_tajima_files(tajima_dir, gen14_pattern, "Gen14") %>%
    transmute(
      chr,
      chr_num,
      WinCenter,
      Tajima_Gen14 = Tajima,
      pi_Gen14 = tP / nSites,
      nSites_Gen14 = nSites
    )

  merged_df <- fst_df %>%
    inner_join(tajima_gen0, by = c("chr", "chr_num", "WinCenter")) %>%
    inner_join(tajima_gen14, by = c("chr", "chr_num", "WinCenter")) %>%
    mutate(
      deltaD = Tajima_Gen14 - Tajima_Gen0,
      pi_ratio = pi_Gen14 / pi_Gen0,
      log_pi_ratio = log2(pi_ratio)
    ) %>%
    filter(
      is.finite(Fst),
      is.finite(deltaD),
      is.finite(pi_ratio),
      is.finite(log_pi_ratio)
    ) %>%
    arrange(chr_num, WinCenter)

  fst_threshold <- quantile(merged_df$Fst, fst_percentile, na.rm = TRUE)
  deltaD_threshold <- quantile(merged_df$deltaD, tail_percentile, na.rm = TRUE)
  log_pi_threshold <- quantile(merged_df$log_pi_ratio, tail_percentile, na.rm = TRUE)

  flagged_df <- merged_df %>%
    mutate(
      fst_outlier = Fst >= fst_threshold,
      deltaD_outlier = deltaD <= deltaD_threshold,
      pi_outlier = log_pi_ratio <= log_pi_threshold,
      is_candidate = fst_outlier & deltaD_outlier & pi_outlier
    ) %>%
    group_by(chr_num) %>%
    mutate(
      Fst_smooth = running_median(Fst, k = 5),
      deltaD_smooth = running_median(deltaD, k = 5),
      log_pi_smooth = running_median(log_pi_ratio, k = 5)
    ) %>%
    ungroup()

  candidate_regions <- merge_candidate_regions(
    flagged_df %>% filter(is_candidate),
    window_size_bp = window_size_bp,
    candidate_gap_bp = candidate_gap_bp
  )

  list(
    data = flagged_df,
    regions = candidate_regions,
    thresholds = list(
      fst = fst_threshold,
      deltaD = deltaD_threshold,
      log_pi = log_pi_threshold
    )
  )
}

add_cumulative_positions <- function(df, chr_info) {
  df %>%
    left_join(chr_info %>% select(chr_num, offset, center), by = "chr_num") %>%
    mutate(pos_cum = WinCenter + offset)
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

plot_population_panel <- function(pop_obj, pop_label, chr_info) {
  df <- add_cumulative_positions(pop_obj$data, chr_info)
  regions_df <- add_region_positions(pop_obj$regions, chr_info)

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
    geom_point(aes(y = Fst), color = "grey82", size = 0.22, alpha = 0.3) +
    geom_line(aes(y = Fst_smooth, group = chr_num), color = fst_color, linewidth = 1.05) +
    geom_point(
      data = df %>% filter(fst_outlier),
      aes(y = Fst),
      color = outlier_color,
      size = 1.15,
      alpha = 0.95
    ) +
    geom_hline(
      yintercept = pop_obj$thresholds$fst,
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
    geom_point(aes(y = deltaD), color = "grey82", size = 0.22, alpha = 0.3) +
    geom_line(aes(y = deltaD_smooth, group = chr_num), color = delta_color, linewidth = 1.05) +
    geom_point(
      data = df %>% filter(deltaD_outlier),
      aes(y = deltaD),
      color = outlier_color,
      size = 1.15,
      alpha = 0.95
    ) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    geom_hline(
      yintercept = pop_obj$thresholds$deltaD,
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
    geom_point(aes(y = log_pi_ratio), color = "grey82", size = 0.22, alpha = 0.3) +
    geom_line(aes(y = log_pi_smooth, group = chr_num), color = pi_color, linewidth = 1.05) +
    geom_point(
      data = df %>% filter(pi_outlier),
      aes(y = log_pi_ratio),
      color = outlier_color,
      size = 1.15,
      alpha = 0.95
    ) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    geom_hline(
      yintercept = pop_obj$thresholds$log_pi,
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

cat("\n=== Loading FST, Tajima's D, and pi data ===\n")

ic_obj <- prepare_population_data(
  fst_dir = fst_dir_ic,
  tajima_dir = tajima_dir_ic,
  gen0_pattern = "IC01",
  gen14_pattern = "IC14",
  population_name = "Indian Chief"
)

j_obj <- prepare_population_data(
  fst_dir = fst_dir_j,
  tajima_dir = tajima_dir_j,
  gen0_pattern = "J01",
  gen14_pattern = "J14",
  population_name = "Jarvis"
)

combined_df <- bind_rows(ic_obj$data, j_obj$data)

chr_info <- combined_df %>%
  group_by(chr_num) %>%
  summarise(chr_len = max(WinCenter, na.rm = TRUE), .groups = "drop") %>%
  arrange(chr_num) %>%
  mutate(
    offset = lag(cumsum(chr_len), default = 0),
    center = offset + chr_len / 2
  )

################################################################################
### BUILD FIGURES
################################################################################

cat("\n=== Building figure panels ===\n")

panel_ic <- plot_population_panel(ic_obj, "Indian Chief", chr_info)
panel_j <- plot_population_panel(j_obj, "Jarvis", chr_info)

fig4 <- (wrap_elements(panel = panel_ic) / wrap_elements(panel = panel_j)) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A")

################################################################################
### SAVE OUTPUTS
################################################################################

ggsave(main_file, fig4, width = 15, height = 16, dpi = 300, bg = "white")

################################################################################
### REPORT
################################################################################

cat("\nSaved main figure to:\n")
cat("  ", main_file, "\n", sep = "")

cat("\nSummary:\n")
cat(
  "  Indian Chief windows:", nrow(ic_obj$data),
  "| candidate regions:", nrow(ic_obj$regions), "\n"
)
cat(
  "  Jarvis windows:", nrow(j_obj$data),
  "| candidate regions:", nrow(j_obj$regions), "\n"
)
