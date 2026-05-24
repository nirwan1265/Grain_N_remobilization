################################################################################
### GENWIN SPLINE FST WINDOWS FROM ANGSD realSFS OUTPUT
################################################################################

required_packages <- c(
  "GenWin",
  "dplyr",
  "ggplot2",
  "pspline",
  "purrr",
  "readr",
  "stringr",
  "tibble"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them with install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))"
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
})

################################################################################
### CONFIGURATION
################################################################################

population_code <- toupper(Sys.getenv("POPULATION_CODE", unset = "IC"))

if (population_code == "IC") {
  gen0_label <- "IC01"
  gen14_label <- "IC14"
  default_fst_idx_dir <- "/Users/nirwantandukar/Documents/Research/data/Indian_Jarvis/Fst/Individual/IC"
} else if (population_code == "J") {
  gen0_label <- "J01"
  gen14_label <- "J14"
  default_fst_idx_dir <- "/Users/nirwantandukar/Documents/Research/data/Indian_Jarvis/Fst/Individual/J"
} else {
  stop("Unsupported POPULATION_CODE: ", population_code, ". Use IC or J.")
}

comparison_label <- Sys.getenv(
  "COMPARISON_LABEL",
  unset = paste0(gen0_label, "_vs_", gen14_label)
)

fst_idx_dir <- Sys.getenv("FST_IDX_DIR", unset = default_fst_idx_dir)

# If realSFS is not on PATH, set REALSFS_BIN or replace this with the full path.
realSFS_bin <- Sys.getenv("REALSFS_BIN", unset = "realSFS")

# This matches the 500 kb smoothing scale from the earlier maize GenWin workflow.
smoothness_bp <- as.numeric(Sys.getenv("SMOOTHNESS_BP", unset = "500000"))
output_tag <- Sys.getenv("OUTPUT_TAG", unset = "")

if (!is.finite(smoothness_bp) || smoothness_bp <= 0) {
  stop("SMOOTHNESS_BP must be a positive number. Got: ", smoothness_bp)
}

genwin_method <- 4
min_sites_per_chr <- 25
plot_point_cap <- 25000

table_dir <- "tables/supplementary"
fig_dir <- "Figs/Supplementary"

smoothness_label <- paste0(as.integer(round(smoothness_bp / 1000)), "kb")
file_suffix <- if (nzchar(output_tag)) {
  paste0("_", gsub("[^A-Za-z0-9._-]+", "_", output_tag))
} else {
  ""
}

windows_file <- file.path(
  table_dir,
  paste0("SuppTable_", population_code, "_GenWin_spline_Fst_windows", file_suffix, ".csv")
)
breaks_file <- file.path(
  table_dir,
  paste0("SuppTable_", population_code, "_GenWin_spline_Fst_breaks", file_suffix, ".csv")
)
summary_file <- file.path(
  table_dir,
  paste0("SuppTable_", population_code, "_GenWin_spline_Fst_summary", file_suffix, ".csv")
)
plot_file <- file.path(
  fig_dir,
  paste0("SuppFig_", population_code, "_GenWin_spline_Fst", file_suffix, ".pdf")
)

dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

chr_map <- tibble::tribble(
  ~nc_id,         ~chr,   ~chr_num,
  "NC_050096.1",  "chr1", 1L,
  "NC_050097.1",  "chr2", 2L,
  "NC_050098.1",  "chr3", 3L,
  "NC_050099.1",  "chr4", 4L,
  "NC_050100.1",  "chr5", 5L,
  "NC_050101.1",  "chr6", 6L,
  "NC_050102.1",  "chr7", 7L,
  "NC_050103.1",  "chr8", 8L,
  "NC_050104.1",  "chr9", 9L,
  "NC_050105.1",  "chr10", 10L
)

################################################################################
### HELPERS
################################################################################

resolve_realSFS_bin <- function(realSFS_bin) {
  if (file.exists(realSFS_bin)) {
    return(normalizePath(realSFS_bin, mustWork = TRUE))
  }

  located <- Sys.which(realSFS_bin)

  if (nzchar(located)) {
    return(located)
  }

  stop(
    "Could not find realSFS.\n",
    "Put it on PATH or set REALSFS_BIN / `realSFS_bin` to the executable path."
  )
}

normalize_names <- function(x) {
  tolower(gsub("[^a-z0-9]+", "", x))
}

match_first_column <- function(norm_names, candidates) {
  idx <- which(norm_names %in% candidates)

  if (length(idx) == 0) {
    return(NA_integer_)
  }

  idx[[1]]
}

first_line_is_data <- function(line) {
  tokens <- strsplit(stringr::str_squish(line), "\\s+")[[1]]

  if (length(tokens) < 4) {
    return(FALSE)
  }

  probe <- suppressWarnings(as.numeric(tokens[2:4]))
  all(!is.na(probe))
}

root_positions <- function(second_derivative_values) {
  if (length(second_derivative_values) < 2) {
    return(numeric())
  }

  left_vals <- c(NA_real_, second_derivative_values[1:(length(second_derivative_values) - 1)])
  right_vals <- second_derivative_values

  pos_to_neg <- which(left_vals > 0 & right_vals < 0) - 0.5
  neg_to_pos <- which(left_vals < 0 & right_vals > 0) - 0.5
  exact_zero <- which(second_derivative_values == 0)

  sort(c(pos_to_neg, neg_to_pos, exact_zero))
}

parse_realSFS_print <- function(path, source_idx_path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!stringr::str_detect(stringr::str_trim(lines), "^$|^->")]

  if (length(lines) == 0) {
    stop("No usable output from realSFS fst print for: ", source_idx_path)
  }

  has_header <- !first_line_is_data(lines[[1]])

  raw_df <- readr::read_table(
    I(paste(lines, collapse = "\n")),
    col_names = has_header,
    progress = FALSE
  )

  if (!has_header) {
    names(raw_df)[seq_len(min(4, ncol(raw_df)))] <- c("chr", "position", "A", "B")[seq_len(min(4, ncol(raw_df)))]
  }

  fallback_chr <- stringr::str_extract(basename(source_idx_path), "NC_[0-9]+\\.[0-9]+")
  norm_names <- normalize_names(names(raw_df))

  chr_col <- match_first_column(
    norm_names,
    c("chr", "chrom", "chromosome", "region", "contig", "scaffold")
  )
  pos_col <- match_first_column(
    norm_names,
    c("pos", "position", "bp", "site", "coordinate")
  )
  a_col <- match_first_column(
    norm_names,
    c("a", "alpha", "reynoldsalpha", "numerator")
  )
  b_col <- match_first_column(
    norm_names,
    c("b", "beta", "alphabeta", "denominator")
  )
  fst_col <- match_first_column(
    norm_names,
    c("fst", "weightedfst", "fstweight", "fstvaoverb")
  )

  if (is.na(chr_col) && ncol(raw_df) >= 1) {
    chr_col <- 1L
  }

  if (is.na(pos_col) && ncol(raw_df) >= 2) {
    pos_col <- 2L
  }

  if (is.na(a_col) && ncol(raw_df) >= 3) {
    a_col <- 3L
  }

  if (is.na(b_col) && ncol(raw_df) >= 4) {
    b_col <- 4L
  }

  chr_values <- rep(fallback_chr, nrow(raw_df))
  if (!is.na(chr_col)) {
    chr_values <- as.character(raw_df[[chr_col]])
  }

  pos_values <- rep(NA_real_, nrow(raw_df))
  if (!is.na(pos_col)) {
    pos_values <- suppressWarnings(as.numeric(raw_df[[pos_col]]))
  }

  if (all(is.na(pos_values)) && !is.na(chr_col)) {
    split_region <- stringr::str_match(chr_values, "^(.+?):([0-9]+)$")

    if (ncol(split_region) == 3 && any(!is.na(split_region[, 1]))) {
      matched_rows <- !is.na(split_region[, 1])
      chr_values[matched_rows] <- split_region[matched_rows, 2]
      pos_values[matched_rows] <- as.numeric(split_region[matched_rows, 3])
    }
  }

  a_values <- rep(NA_real_, nrow(raw_df))
  b_values <- rep(NA_real_, nrow(raw_df))

  if (!is.na(a_col)) {
    a_values <- suppressWarnings(as.numeric(raw_df[[a_col]]))
  }

  if (!is.na(b_col)) {
    b_values <- suppressWarnings(as.numeric(raw_df[[b_col]]))
  }

  fst_values <- rep(NA_real_, nrow(raw_df))

  if (!is.na(fst_col)) {
    fst_values <- suppressWarnings(as.numeric(raw_df[[fst_col]]))
  } else if (!is.na(a_col) && !is.na(b_col)) {
    fst_values <- ifelse(b_values == 0, NA_real_, a_values / b_values)
  } else {
    stop(
      "Could not identify FST columns in realSFS output for: ",
      source_idx_path
    )
  }

  tibble::tibble(
    SourceFile = basename(source_idx_path),
    Chr = chr_values,
    Position = pos_values,
    A = a_values,
    B = b_values,
    Fst = fst_values
  ) %>%
    dplyr::filter(
      !is.na(Chr),
      is.finite(Position),
      is.finite(Fst)
    )
}

read_realSFS_print <- function(idx_path, realSFS_bin) {
  out_path <- tempfile(fileext = ".txt")
  err_path <- tempfile(fileext = ".log")
  on.exit(unlink(c(out_path, err_path)), add = TRUE)

  status <- system2(
    command = realSFS_bin,
    args = c("fst", "print", idx_path),
    stdout = out_path,
    stderr = err_path
  )

  err_lines <- readLines(err_path, warn = FALSE)

  if (!identical(status, 0L)) {
    stop(
      "realSFS fst print failed for: ", idx_path, "\n",
      paste(err_lines, collapse = "\n")
    )
  }

  parse_realSFS_print(out_path, idx_path)
}

load_all_fst_sites <- function(fst_idx_dir, realSFS_bin) {
  idx_files <- list.files(
    fst_idx_dir,
    pattern = "\\.fst\\.idx$",
    full.names = TRUE
  )

  if (length(idx_files) == 0) {
    stop("No .fst.idx files found in: ", fst_idx_dir)
  }

  idx_files <- sort(idx_files)

  purrr::map_dfr(idx_files, function(idx_path) {
    cat("  extracting:", basename(idx_path), "\n")
    read_realSFS_print(idx_path, realSFS_bin)
  }) %>%
    dplyr::left_join(chr_map, by = c("Chr" = "nc_id")) %>%
    dplyr::mutate(
      chr = dplyr::coalesce(chr, Chr),
      chr_num = dplyr::coalesce(
        chr_num,
        as.integer(factor(chr, levels = unique(chr)))
      )
    ) %>%
    dplyr::arrange(chr_num, Position)
}

build_single_window_table <- function(chr_df, global_mean, global_var) {
  mean_fst <- mean(chr_df$Fst, na.rm = TRUE)
  snp_count <- nrow(chr_df)

  if (is.na(global_var) || global_var <= 0 || snp_count == 0) {
    w_stat <- NA_real_
  } else {
    w_stat <- (mean_fst - global_mean) / sqrt(global_var / snp_count)
  }

  tibble::tibble(
    WindowStart = min(chr_df$Position, na.rm = TRUE),
    WindowStop = max(chr_df$Position, na.rm = TRUE),
    SNPcount = snp_count,
    MeanY = mean_fst,
    Wstat = w_stat
  )
}

sample_sites_for_plot <- function(site_df, plot_point_cap) {
  site_df %>%
    dplyr::group_by(chr, chr_num) %>%
    dplyr::group_modify(function(.x, .y) {
      if (nrow(.x) <= plot_point_cap) {
        .x
      } else {
        dplyr::slice_sample(.x, n = plot_point_cap)
      }
    }) %>%
    dplyr::ungroup()
}

run_genwin_for_chr <- function(
  chr_df,
  population_code,
  comparison_label,
  smoothness_bp,
  genwin_method,
  global_mean,
  global_var
) {
  chr_df <- chr_df %>%
    dplyr::group_by(Chr, chr, chr_num, Position) %>%
    dplyr::summarise(
      A = sum(A, na.rm = TRUE),
      B = sum(B, na.rm = TRUE),
      Fst = mean(Fst, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(Fst)) %>%
    dplyr::arrange(Position)

  if (nrow(chr_df) < min_sites_per_chr) {
    stop(
      "Too few sites for ", unique(chr_df$chr),
      ": found ", nrow(chr_df), ", need at least ", min_sites_per_chr
    )
  }

  fit <- GenWin::splineAnalyze(
    Y = chr_df$Fst,
    map = chr_df$Position,
    smoothness = smoothness_bp,
    s2 = global_var,
    mean = global_mean,
    plotRaw = FALSE,
    plotWindows = FALSE,
    method = genwin_method
  )

  windows_df <- tibble::as_tibble(fit$windowData) %>%
    dplyr::mutate(
      Population = population_code,
      Comparison = comparison_label,
      nc_id = unique(chr_df$Chr),
      chr = unique(chr_df$chr),
      chr_num = unique(chr_df$chr_num),
      window_id = dplyr::row_number(),
      # These disjoint bp bounds are the ones to reuse for pi and Tajima's D.
      window_start_bp = dplyr::if_else(
        window_id == 1L,
        pmax(1, floor(WindowStart)),
        floor(WindowStart) + 1
      ),
      window_end_bp = dplyr::if_else(
        window_id == dplyr::n(),
        ceiling(WindowStop),
        floor(WindowStop)
      ),
      window_end_bp = pmax(window_end_bp, window_start_bp),
      window_mid_bp = round((window_start_bp + window_end_bp) / 2),
      smoothness_bp = smoothness_bp,
      genwin_method = genwin_method,
      n_input_sites = nrow(chr_df)
    ) %>%
    dplyr::rename(
      window_start_raw = WindowStart,
      window_end_raw = WindowStop,
      snp_count = SNPcount,
      mean_fst = MeanY
    ) %>%
    dplyr::select(
      Population,
      Comparison,
      nc_id,
      chr,
      chr_num,
      window_id,
      window_start_bp,
      window_end_bp,
      window_mid_bp,
      window_start_raw,
      window_end_raw,
      snp_count,
      mean_fst,
      Wstat,
      n_input_sites,
      smoothness_bp,
      genwin_method
    )

  if (length(fit$breaks) == 0) {
    breaks_df <- tibble::tibble(
      Population = character(),
      Comparison = character(),
      nc_id = character(),
      chr = character(),
      chr_num = integer(),
      break_id = integer(),
      break_position = numeric()
    )
  } else {
    breaks_df <- tibble::tibble(
      Population = population_code,
      Comparison = comparison_label,
      nc_id = unique(chr_df$Chr),
      chr = unique(chr_df$chr),
      chr_num = unique(chr_df$chr_num),
      break_id = seq_along(fit$breaks),
      break_position = fit$breaks
    )
  }

  grid_pos <- seq(
    min(chr_df$Position, na.rm = TRUE),
    max(chr_df$Position, na.rm = TRUE),
    by = smoothness_bp
  )

  spline_df <- tibble::tibble(
    nc_id = unique(chr_df$Chr),
    chr = unique(chr_df$chr),
    chr_num = unique(chr_df$chr_num),
    Position = grid_pos,
    Fst_spline = as.numeric(stats::predict(fit$rawSpline, grid_pos))
  )

  summary_df <- tibble::tibble(
    Population = population_code,
    Comparison = comparison_label,
    nc_id = unique(chr_df$Chr),
    chr = unique(chr_df$chr),
    chr_num = unique(chr_df$chr_num),
    n_sites = nrow(chr_df),
    n_windows = nrow(windows_df),
    n_breaks = length(fit$breaks),
    chr_start_bp = min(chr_df$Position, na.rm = TRUE),
    chr_end_bp = max(chr_df$Position, na.rm = TRUE),
    mean_fst_chr = mean(chr_df$Fst, na.rm = TRUE)
  )

  list(
    sites = chr_df,
    windows = windows_df,
    breaks = breaks_df,
    spline = spline_df,
    summary = summary_df
  )
}

build_spline_plot <- function(site_df, spline_df, breaks_df, population_code, smoothness_label) {
  sampled_sites <- sample_sites_for_plot(site_df, plot_point_cap = plot_point_cap)

  ggplot2::ggplot() +
    ggplot2::geom_point(
      data = sampled_sites,
      ggplot2::aes(x = Position, y = Fst),
      color = "grey65",
      alpha = 0.25,
      size = 0.16
    ) +
    ggplot2::geom_line(
      data = spline_df,
      ggplot2::aes(x = Position, y = Fst_spline),
      color = "#163A5F",
      linewidth = 0.45
    ) +
    ggplot2::geom_vline(
      data = breaks_df,
      ggplot2::aes(xintercept = break_position),
      color = "#B22222",
      linewidth = 0.28,
      alpha = 0.7
    ) +
    ggplot2::facet_wrap(~ chr, scales = "free_x", ncol = 2) +
    ggplot2::labs(
      title = paste0(population_code, " GenWin spline FST windows (", smoothness_label, ")"),
      subtitle = paste0(
        "Points = per-site FST from `realSFS fst print`; ",
        "red lines = spline breakpoints"
      ),
      x = "Position (bp)",
      y = expression(F[ST])
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black"),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10)
    )
}

################################################################################
### RUN
################################################################################

cat("\n=== Running GenWin spline FST workflow ===\n")

realSFS_bin <- resolve_realSFS_bin(realSFS_bin)

cat("\nInput directory:\n  ", fst_idx_dir, "\n", sep = "")
cat("realSFS executable:\n  ", realSFS_bin, "\n", sep = "")

cat("\n=== Extracting per-site FST from ANGSD indices ===\n")
all_sites <- load_all_fst_sites(
  fst_idx_dir = fst_idx_dir,
  realSFS_bin = realSFS_bin
)

global_mean_fst <- mean(all_sites$Fst, na.rm = TRUE)
global_var_fst <- stats::var(all_sites$Fst, na.rm = TRUE)

cat("\nGlobal FST mean:", signif(global_mean_fst, 5), "\n")
cat("Global FST variance:", signif(global_var_fst, 5), "\n")

cat("\n=== Fitting spline windows chromosome by chromosome ===\n")

chr_results <- all_sites %>%
  split(.$chr) %>%
  purrr::map(function(chr_df) {
    cat("  fitting:", unique(chr_df$chr), "\n")
    run_genwin_for_chr(
      chr_df = chr_df,
      population_code = population_code,
      comparison_label = comparison_label,
      smoothness_bp = smoothness_bp,
      genwin_method = genwin_method,
      global_mean = global_mean_fst,
      global_var = global_var_fst
    )
  })

window_table <- purrr::map_dfr(chr_results, "windows") %>%
  dplyr::arrange(chr_num, window_id)

break_table <- purrr::map_dfr(chr_results, "breaks") %>%
  dplyr::arrange(chr_num, break_id)

summary_table <- purrr::map_dfr(chr_results, "summary") %>%
  dplyr::arrange(chr_num) %>%
  dplyr::mutate(
    global_mean_fst = global_mean_fst,
    global_var_fst = global_var_fst,
    smoothness_bp = smoothness_bp,
    genwin_method = genwin_method
  )

spline_table <- purrr::map_dfr(chr_results, "spline")
site_table <- purrr::map_dfr(chr_results, "sites")

cat("\n=== Writing outputs ===\n")

readr::write_csv(window_table, windows_file)
readr::write_csv(break_table, breaks_file)
readr::write_csv(summary_table, summary_file)

plot_obj <- build_spline_plot(
  site_df = site_table,
  spline_df = spline_table,
  breaks_df = break_table,
  population_code = population_code,
  smoothness_label = smoothness_label
)

ggplot2::ggsave(
  filename = plot_file,
  plot = plot_obj,
  width = 14,
  height = 11,
  dpi = 300,
  bg = "white"
)

cat("\nSaved files:\n")
cat("  ", windows_file, "\n", sep = "")
cat("  ", breaks_file, "\n", sep = "")
cat("  ", summary_file, "\n", sep = "")
cat("  ", plot_file, "\n", sep = "")

cat("\nWindow summary by chromosome:\n")
print(summary_table)
