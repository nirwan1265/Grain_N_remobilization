################################################################################
### GENWIN SPLINE FST WINDOWS WITH TAJIMA'S D, PI, FIGURE, AND CANDIDATES
################################################################################

required_packages <- c(
  "data.table",
  "dplyr",
  "GenomicRanges",
  "ggplot2",
  "IRanges",
  "patchwork",
  "purrr",
  "readr",
  "rtracklayer",
  "scales",
  "stringr",
  "tibble",
  "tidyr"
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
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(ggplot2)
  library(IRanges)
  library(patchwork)
  library(purrr)
  library(readr)
  library(rtracklayer)
  library(scales)
  library(stringr)
  library(tibble)
  library(tidyr)
})

################################################################################
### CONFIGURATION
################################################################################

population_code <- toupper(Sys.getenv("POPULATION_CODE", unset = "IC"))

if (population_code == "IC") {
  population_name <- "Indian Chief"
  gen0_label <- "IC01"
  gen14_label <- "IC14"
  default_theta_dir <- "/Users/nirwantandukar/Documents/Research/data/Indian_Jarvis/TajimaD_NucleotideDiversity/IC"
  default_pestpg_dir <- "/Users/nirwantandukar/Documents/Research/results/Indian_Jarvis/ANGSD_TajimaD/Indian_Chief/sliding_window"
} else if (population_code == "J") {
  population_name <- "Jarvis"
  gen0_label <- "J01"
  gen14_label <- "J14"
  default_theta_dir <- "/Users/nirwantandukar/Documents/Research/data/Indian_Jarvis/TajimaD_NucleotideDiversity/J"
  default_pestpg_dir <- "/Users/nirwantandukar/Documents/Research/results/Indian_Jarvis/ANGSD_TajimaD/Jarvis/sliding_window"
} else {
  stop("Unsupported POPULATION_CODE: ", population_code, ". Use IC or J.")
}

comparison_label <- Sys.getenv(
  "COMPARISON_LABEL",
  unset = paste0(gen0_label, "_vs_", gen14_label)
)

theta_dir <- Sys.getenv("THETA_DIR", unset = default_theta_dir)
pestpg_dir <- Sys.getenv("PESTPG_DIR", unset = default_pestpg_dir)
theta_mode <- tolower(Sys.getenv("THETA_MODE", unset = "auto"))
theta_stat_bin <- Sys.getenv("THETASTAT_BIN", unset = "thetaStat")
output_tag <- Sys.getenv("OUTPUT_TAG", unset = "")

gff_file <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize/Maize.annotation/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3"

file_suffix <- if (nzchar(output_tag)) {
  paste0("_", gsub("[^A-Za-z0-9._-]+", "_", output_tag))
} else {
  ""
}

windows_file <- file.path(
  "tables/supplementary",
  paste0("SuppTable_", population_code, "_GenWin_spline_Fst_windows", file_suffix, ".csv")
)

window_table_file <- file.path(
  "tables/supplementary",
  paste0("SuppTable_", population_code, "_GenWin_windows_Fst_Tajima_pi", file_suffix, ".csv")
)
window_gene_file <- file.path(
  "tables/supplementary",
  paste0("SuppTable_", population_code, "_GenWin_windows_Fst_Tajima_pi_genes", file_suffix, ".csv")
)
candidate_regions_file <- file.path(
  "tables/supplementary",
  paste0("SuppTable_", population_code, "_GenWin_candidate_regions", file_suffix, ".csv")
)
candidate_windows_file <- file.path(
  "tables/supplementary",
  paste0("SuppTable_", population_code, "_GenWin_candidate_windows", file_suffix, ".csv")
)
figure_file <- file.path(
  "Figs/Supplementary",
  paste0("SuppFig_", population_code, "_GenWin_Fst_Tajima_pi", file_suffix, ".pdf")
)

dir.create(dirname(window_table_file), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(figure_file), showWarnings = FALSE, recursive = TRUE)

fst_percentile <- 0.95
tail_percentile <- 0.05
candidate_gap_bp <- 500000
min_display_window_bp <- 12000000

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

resolve_binary <- function(bin_path, label) {
  if (file.exists(bin_path)) {
    return(normalizePath(bin_path, mustWork = TRUE))
  }

  located <- Sys.which(bin_path)

  if (nzchar(located)) {
    return(located)
  }

  stop("Could not find ", label, " at: ", bin_path)
}

resolve_pestpg_file <- function(pestpg_dir, nc_id, generation_label) {
  matches <- Sys.glob(file.path(pestpg_dir, paste0(nc_id, "_", generation_label, "*.pestPG")))
  matches <- sort(unique(matches))

  if (length(matches) == 0) {
    stop(
      "Could not find a pestPG summary file for ",
      nc_id, " ", generation_label, " in: ", pestpg_dir
    )
  }

  matches[[1]]
}

detect_theta_mode <- function(theta_mode, theta_dir, pestpg_dir, chr_ids, gen0_label, gen14_label) {
  if (!theta_mode %in% c("auto", "raw", "pestpg")) {
    stop("Unsupported THETA_MODE: ", theta_mode, ". Use auto, raw, or pestpg.")
  }

  raw_paths <- c(
    file.path(theta_dir, paste0(chr_ids, "_", gen0_label, ".thetas.idx")),
    file.path(theta_dir, paste0(chr_ids, "_", gen14_label, ".thetas.idx"))
  )
  raw_ready <- all(file.exists(raw_paths))

  pestpg_ready <- all(vapply(chr_ids, function(chr_id) {
    gen0_hits <- Sys.glob(file.path(pestpg_dir, paste0(chr_id, "_", gen0_label, "*.pestPG")))
    gen14_hits <- Sys.glob(file.path(pestpg_dir, paste0(chr_id, "_", gen14_label, "*.pestPG")))
    length(gen0_hits) > 0 && length(gen14_hits) > 0
  }, logical(1)))

  if (theta_mode == "raw") {
    if (!raw_ready) {
      stop(
        "THETA_MODE=raw was requested, but not all theta index files exist in:\n  ",
        theta_dir
      )
    }
    return("raw")
  }

  if (theta_mode == "pestpg") {
    if (!pestpg_ready) {
      stop(
        "THETA_MODE=pestpg was requested, but not all pestPG files exist in:\n  ",
        pestpg_dir
      )
    }
    return("pestpg")
  }

  if (raw_ready) {
    return("raw")
  }

  if (pestpg_ready) {
    return("pestpg")
  }

  stop(
    "Could not detect usable theta inputs for ", paste(chr_ids, collapse = ", "), ".\n",
    "Checked raw theta indexes in:\n  ", theta_dir, "\n",
    "and pestPG summaries in:\n  ", pestpg_dir
  )
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

a1f <- function(n) {
  sum(1 / seq_len(n - 1))
}

a2f <- function(n) {
  sum(1 / (seq_len(n - 1) ^ 2))
}

b1f <- function(n) {
  (n + 1) / (3 * (n - 1))
}

b2f <- function(n) {
  (2 * (n * n + n + 3)) / (9 * n * (n - 1))
}

c1f <- function(a1, b1) {
  b1 - (1 / a1)
}

c2f <- function(n, a1, a2, b2) {
  b2 - ((n + 2) / (a1 * n)) + (a2 / (a1 * a1))
}

e1f <- function(a1, c1) {
  c1 / a1
}

e2f <- function(a1, a2, c2) {
  c2 / ((a1 * a1) + a2)
}

tajd <- function(n, thetaW, sumk) {
  a1 <- a1f(n)
  segsites <- thetaW * a1

  if (!is.finite(segsites) || segsites == 0) {
    return(0)
  }

  a2 <- a2f(n)
  b1 <- b1f(n)
  b2 <- b2f(n)
  c1 <- c1f(a1, b1)
  c2 <- c2f(n, a1, a2, b2)
  e1 <- e1f(a1, c1)
  e2 <- e2f(a1, a2, c2)

  denom <- sqrt((e1 * segsites) + ((e2 * segsites) * (segsites - 1)))

  if (!is.finite(denom) || denom == 0) {
    return(0)
  }

  (sumk - thetaW) / denom
}

cn_fun <- function(n) {
  (2 * n * a1f(n) - 4 * (n - 1)) / ((n - 1) * (n - 2))
}

vd_fun <- function(n) {
  led1 <- (a1f(n) * a1f(n)) / (a2f(n) + a1f(n) * a1f(n))
  led2 <- cn_fun(n) - (n + 1) / (n - 1)
  1 + led1 * led2
}

ud_fun <- function(n) {
  a1f(n) - 1 - vd_fun(n)
}

vf_fun <- function(n) {
  top <- cn_fun(n) + 2 * (n * n + n + 3) / (9 * n * (n - 1)) - 2 / (n - 1)
  bot <- a1f(n) * a1f(n) + a2f(n)
  top / bot
}

uf_fun <- function(n) {
  top <- 1 + (n + 1) / (3 * (n - 1)) -
    4 * ((n + 1) / ((n - 1) * (n - 1))) * (a1f(n + 1) - (2 * n) / (n + 1))
  (top / a1f(n)) - vf_fun(n)
}

fulid <- function(n, thetaW, thetaFL) {
  S <- thetaW * a1f(n)
  L <- thetaFL * a1f(n)
  top <- S - L
  bot <- ud_fun(n) * S + vd_fun(n) * S * S

  if (!is.finite(bot) || bot <= 0) {
    return(0)
  }

  top / sqrt(bot)
}

fulif <- function(n, thetaW, thetaFL, thetaPi) {
  S <- thetaW * a1f(n)
  top <- thetaPi - thetaFL
  bot <- uf_fun(n) * S + vf_fun(n) * S * S

  if (!is.finite(bot) || bot <= 0) {
    return(0)
  }

  top / sqrt(bot)
}

H_var <- function(n) {
  top <- 18 * n * n * (3 * n + 2) * a2f(n + 1) - (88 * n * n * n + 9 * n * n - 13 * n + 6)
  bot <- 9 * n * (n - 1) * (n - 1)
  top / bot
}

fayh_fun <- function(n, thetaW, thetaH, thetaPi) {
  S <- thetaW * a1f(n)
  led1 <- (n - 2) / (6 * (n - 1)) * S
  led2 <- H_var(n) * S * S
  denom <- sqrt(led1 + led2)

  if (!is.finite(denom) || denom == 0) {
    return(0)
  }

  (thetaPi - thetaH) / denom
}

E_var1 <- function(n) {
  n / (2 * (n - 1)) - 1 / a1f(n)
}

E_var2 <- function(n) {
  led1 <- a2f(n) / (a1f(n) * a1f(n))
  tmp <- n / (n - 1)
  led2 <- 2 * tmp * tmp * a2f(n)
  led3 <- (2 * (n * a2f(n) - n + 1)) / ((n - 1) * a1f(n))
  led4 <- (3 * n + 1) / (n - 1)
  led1 + led2 - led3 - led4
}

zenge_fun <- function(n, thetaW, thetaL) {
  S <- thetaW * a1f(n)
  bot <- E_var1(n) * S + E_var2(n) * S * S

  if (!is.finite(bot) || bot <= 0) {
    return(0)
  }

  (thetaL - thetaW) / sqrt(bot)
}

parse_nchr_from_stderr <- function(stderr_lines, idx_path) {
  match <- stringr::str_match(stderr_lines, "pc\\.nChr=([0-9]+)")
  vals <- suppressWarnings(as.integer(match[, 2]))
  vals <- vals[!is.na(vals)]

  if (length(vals) == 0) {
    stop("Could not parse nChr from thetaStat stderr for: ", idx_path)
  }

  vals[[1]]
}

run_theta_print <- function(idx_path, theta_stat_bin) {
  out_path <- tempfile(fileext = ".txt")
  err_path <- tempfile(fileext = ".log")
  on.exit(unlink(c(out_path, err_path)), add = TRUE)

  status <- system2(
    command = theta_stat_bin,
    args = c("print", idx_path),
    stdout = out_path,
    stderr = err_path
  )

  stderr_lines <- readLines(err_path, warn = FALSE)

  if (!identical(status, 0L)) {
    stop(
      "thetaStat print failed for: ", idx_path, "\n",
      paste(stderr_lines, collapse = "\n")
    )
  }

  dt <- data.table::fread(out_path, showProgress = FALSE)
  n_chr <- parse_nchr_from_stderr(stderr_lines, idx_path)

  names(dt) <- names(dt) |>
    gsub("^#", "", x = _) |>
    trimws()

  col_map <- c(
    "Chromo" = "Chr",
    "Chr" = "Chr",
    "Pos" = "Position",
    "Position" = "Position",
    "Watterson" = "log_tW",
    "Pairwise" = "log_tP",
    "thetaSingleton" = "log_tF",
    "thetaH" = "log_tH",
    "thetaL" = "log_tL"
  )

  matched_old <- intersect(names(col_map), names(dt))

  if (length(matched_old) > 0) {
    setnames(dt, old = matched_old, new = unname(col_map[matched_old]))
  }

  required_cols <- c("Chr", "Position", "log_tW", "log_tP", "log_tF", "log_tH", "log_tL")
  missing_cols <- setdiff(required_cols, names(dt))

  if (length(missing_cols) > 0) {
    stop(
      "thetaStat output is missing expected columns for: ", idx_path, "\n",
      "Found columns: ", paste(names(dt), collapse = ", "), "\n",
      "Missing: ", paste(missing_cols, collapse = ", ")
    )
  }

  dt[, `:=`(
    tW = exp(log_tW),
    tP = exp(log_tP),
    tF = exp(log_tF),
    tH = exp(log_tH),
    tL = exp(log_tL)
  )]

  list(data = dt, nChr = n_chr)
}

read_pestpg_summary <- function(pestpg_path) {
  dt <- data.table::fread(pestpg_path, showProgress = FALSE)

  names(dt) <- names(dt) |>
    gsub("^#", "", x = _) |>
    trimws()

  meta_col <- names(dt)[1]

  if (!all(c("PestStart", "PestEnd") %in% names(dt))) {
    win_match <- stringr::str_match(as.character(dt[[meta_col]]), ".*\\(([0-9]+),([0-9]+)\\)$")
    dt[, `:=`(
      PestStart = suppressWarnings(as.numeric(win_match[, 2])),
      PestEnd = suppressWarnings(as.numeric(win_match[, 3]))
    )]
  }

  required_cols <- c("Chr", "WinCenter", "PestStart", "PestEnd", "tW", "tP", "Tajima", "nSites")
  missing_cols <- setdiff(required_cols, names(dt))

  if (length(missing_cols) > 0) {
    stop(
      "pestPG file is missing expected columns for: ", pestpg_path, "\n",
      "Found columns: ", paste(names(dt), collapse = ", "), "\n",
      "Missing: ", paste(missing_cols, collapse = ", ")
    )
  }

  dt[, `:=`(
    WinCenter = as.numeric(WinCenter),
    PestStart = as.numeric(PestStart),
    PestEnd = as.numeric(PestEnd),
    tW = as.numeric(tW),
    tP = as.numeric(tP),
    Tajima = as.numeric(Tajima),
    nSites = as.numeric(nSites)
  )]

  dt
}

summarise_theta_for_windows <- function(theta_dt, windows_dt, n_chr, generation_label) {
  theta_dt <- data.table::as.data.table(copy(theta_dt))
  windows_dt <- data.table::as.data.table(copy(windows_dt))

  setkey(theta_dt, Position, Position)

  if (!"window_row_id" %in% names(windows_dt)) {
    windows_dt[, window_row_id := .I]
  }

  setkey(windows_dt, window_start_bp, window_end_bp)

  joined <- foverlaps(
    x = theta_dt[, .(
      Position,
      PositionEnd = Position,
      tW,
      tP,
      tF,
      tH,
      tL
    )],
    y = windows_dt[, .(
      window_row_id,
      window_start_bp,
      window_end_bp
    )],
    by.x = c("Position", "PositionEnd"),
    by.y = c("window_start_bp", "window_end_bp"),
    type = "within",
    nomatch = 0L
  )

  if (nrow(joined) == 0) {
    stop("No theta sites overlapped the GenWin windows for ", generation_label)
  }

  summary_dt <- joined[
    ,
    .(
      n_theta_sites = .N,
      tW = sum(tW, na.rm = TRUE),
      tP = sum(tP, na.rm = TRUE),
      tF = sum(tF, na.rm = TRUE),
      tH = sum(tH, na.rm = TRUE),
      tL = sum(tL, na.rm = TRUE)
    ),
    by = window_row_id
  ]

  summary_dt[, `:=`(
    Tajima = vapply(seq_len(.N), function(i) tajd(n_chr, tW[i], tP[i]), numeric(1)),
    fuf = vapply(seq_len(.N), function(i) fulif(n_chr, tW[i], tF[i], tP[i]), numeric(1)),
    fud = vapply(seq_len(.N), function(i) fulid(n_chr, tW[i], tF[i]), numeric(1)),
    fayh = vapply(seq_len(.N), function(i) fayh_fun(n_chr, tW[i], tH[i], tP[i]), numeric(1)),
    zeng = vapply(seq_len(.N), function(i) zenge_fun(n_chr, tW[i], tL[i]), numeric(1)),
    pi = ifelse(n_theta_sites > 0, tP / n_theta_sites, NA_real_),
    thetaW_per_site = ifelse(n_theta_sites > 0, tW / n_theta_sites, NA_real_)
  )]

  summary_dt[, `:=`(
    support_windows = NA_integer_,
    mean_window_nSites = NA_real_
  )]

  summary_tbl <- as_tibble(summary_dt) %>%
    dplyr::rename_with(~ paste0(., "_", generation_label), c(
      "n_theta_sites",
      "support_windows",
      "mean_window_nSites",
      "tW",
      "tP",
      "tF",
      "tH",
      "tL",
      "Tajima",
      "fuf",
      "fud",
      "fayh",
      "zeng",
      "pi",
      "thetaW_per_site"
    ))

  n_chr_col <- paste0("nChr_", generation_label)
  summary_tbl[[n_chr_col]] <- n_chr
  summary_tbl
}

summarise_pestpg_for_windows <- function(pestpg_dt, windows_dt, generation_label) {
  pestpg_dt <- data.table::as.data.table(copy(pestpg_dt))
  windows_dt <- data.table::as.data.table(copy(windows_dt))

  if (!"window_row_id" %in% names(windows_dt)) {
    windows_dt[, window_row_id := .I]
  }

  setkey(pestpg_dt, PestStart, PestEnd)
  setkey(windows_dt, window_start_bp, window_end_bp)

  joined <- foverlaps(
    x = pestpg_dt[, .(
      PestStart,
      PestEnd,
      tW,
      tP,
      Tajima,
      nSites
    )],
    y = windows_dt[, .(
      window_row_id,
      window_start_bp,
      window_end_bp
    )],
    by.x = c("PestStart", "PestEnd"),
    by.y = c("window_start_bp", "window_end_bp"),
    type = "any",
    nomatch = 0L
  )

  if (nrow(joined) == 0) {
    stop("No pestPG summary windows overlapped the GenWin windows for ", generation_label)
  }

  joined[, `:=`(
    overlap_bp = pmin(PestEnd, window_end_bp) - pmax(PestStart, window_start_bp) + 1,
    pest_window_bp = PestEnd - PestStart + 1
  )]

  joined <- joined[is.finite(overlap_bp) & overlap_bp > 0 & is.finite(pest_window_bp) & pest_window_bp > 0]
  joined[, overlap_fraction := overlap_bp / pest_window_bp]
  joined[, `:=`(
    effective_sites = nSites * overlap_fraction,
    effective_tW = tW * overlap_fraction,
    effective_tP = tP * overlap_fraction
  )]

  summary_dt <- joined[
    ,
    .(
      n_theta_sites = NA_real_,
      support_windows = .N,
      mean_window_nSites = mean(nSites, na.rm = TRUE),
      tW = sum(effective_tW, na.rm = TRUE),
      tP = sum(effective_tP, na.rm = TRUE),
      tF = NA_real_,
      tH = NA_real_,
      tL = NA_real_,
      Tajima = stats::weighted.mean(Tajima, w = pmax(effective_sites, 1e-9), na.rm = TRUE),
      fuf = NA_real_,
      fud = NA_real_,
      fayh = NA_real_,
      zeng = NA_real_,
      pi = ifelse(sum(effective_sites, na.rm = TRUE) > 0, sum(effective_tP, na.rm = TRUE) / sum(effective_sites, na.rm = TRUE), NA_real_),
      thetaW_per_site = ifelse(sum(effective_sites, na.rm = TRUE) > 0, sum(effective_tW, na.rm = TRUE) / sum(effective_sites, na.rm = TRUE), NA_real_)
    ),
    by = window_row_id
  ]

  summary_tbl <- as_tibble(summary_dt) %>%
    dplyr::rename_with(~ paste0(., "_", generation_label), c(
      "n_theta_sites",
      "support_windows",
      "mean_window_nSites",
      "tW",
      "tP",
      "tF",
      "tH",
      "tL",
      "Tajima",
      "fuf",
      "fud",
      "fayh",
      "zeng",
      "pi",
      "thetaW_per_site"
    ))

  n_chr_col <- paste0("nChr_", generation_label)
  summary_tbl[[n_chr_col]] <- NA_integer_
  summary_tbl
}

merge_candidate_regions <- function(df, candidate_gap_bp) {
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

  df %>%
    arrange(chr_num, window_start_bp) %>%
    mutate(
      new_region = if_else(
        row_number() == 1,
        TRUE,
        chr_num != lag(chr_num) | (window_start_bp - lag(window_end_bp)) > candidate_gap_bp
      ),
      region_id = cumsum(new_region)
    ) %>%
    group_by(chr, chr_num, region_id) %>%
    summarise(
      n_windows = n(),
      region_start = min(window_start_bp),
      region_end = max(window_end_bp),
      region_center = mean(c(region_start, region_end)),
      peak_Fst = max(mean_fst, na.rm = TRUE),
      min_deltaD = min(deltaD, na.rm = TRUE),
      min_log_pi_ratio = min(log_pi_ratio, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(region_length_bp = region_end - region_start) %>%
    arrange(chr_num, region_start)
}

annotate_windows_with_genes <- function(df, genes_gr) {
  win_gr <- GRanges(
    seqnames = df$chr,
    ranges = IRanges(start = df$window_start_bp, end = df$window_end_bp),
    window_row_id = df$window_row_id
  )

  hits <- findOverlaps(win_gr, genes_gr, ignore.strand = TRUE)

  gene_summary <- if (length(hits) == 0) {
    tibble(
      window_row_id = integer(),
      gene_count = integer(),
      GeneIDs = character(),
      GeneNames = character(),
      GeneBiotypes = character()
    )
  } else {
    gene_meta <- mcols(genes_gr)

    tibble(
      window_row_id = queryHits(hits),
      GeneID = as.character(gene_meta$ID[subjectHits(hits)]),
      GeneName = coalesce(
        as.character(gene_meta$Name[subjectHits(hits)]),
        as.character(gene_meta$ID[subjectHits(hits)])
      ),
      GeneBiotype = as.character(gene_meta$biotype[subjectHits(hits)])
    ) %>%
      group_by(window_row_id) %>%
      summarise(
        gene_count = n_distinct(GeneID),
        GeneIDs = paste(unique(GeneID), collapse = ";"),
        GeneNames = paste(unique(GeneName), collapse = ";"),
        GeneBiotypes = paste(unique(stats::na.omit(GeneBiotype)), collapse = ";"),
        .groups = "drop"
      )
  }

  df %>%
    left_join(gene_summary, by = "window_row_id") %>%
    mutate(
      gene_count = replace_na(gene_count, 0L),
      GeneIDs = replace_na(GeneIDs, ""),
      GeneNames = replace_na(GeneNames, ""),
      GeneBiotypes = replace_na(GeneBiotypes, "")
    )
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

plot_population_panel <- function(df, regions_df, thresholds, chr_info, pop_label) {
  df <- add_cumulative_positions(df, chr_info)
  regions_df <- add_region_positions(regions_df, chr_info)

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
### RUN
################################################################################

cat("\n=== Running spline-window Tajima and pi workflow ===\n")

if (!file.exists(windows_file)) {
  stop(
    "Spline window file not found: ", windows_file, "\n",
    "Run scripts/16_spling_Fst.R first."
  )
}

windows_tbl <- readr::read_csv(windows_file, show_col_types = FALSE) %>%
  mutate(window_row_id = row_number()) %>%
  arrange(chr_num, window_id)

theta_input_mode <- detect_theta_mode(
  theta_mode = theta_mode,
  theta_dir = theta_dir,
  pestpg_dir = pestpg_dir,
  chr_ids = unique(windows_tbl$nc_id),
  gen0_label = gen0_label,
  gen14_label = gen14_label
)

if (theta_input_mode == "raw") {
  theta_stat_bin <- resolve_binary(theta_stat_bin, "thetaStat")
  cat("Theta mode:\n  raw theta indexes\n", sep = "")
  cat("thetaStat executable:\n  ", theta_stat_bin, "\n", sep = "")
  cat("Theta directory:\n  ", theta_dir, "\n", sep = "")
} else {
  cat("Theta mode:\n  pestPG window summaries\n", sep = "")
  cat("pestPG directory:\n  ", pestpg_dir, "\n", sep = "")
}

cat("Spline windows:\n  ", windows_file, "\n", sep = "")

cat("\n=== Summarising theta files into spline windows ===\n")

theta_summaries <- windows_tbl %>%
  group_split(chr, .keep = TRUE) %>%
  purrr::map(function(chr_windows) {
    nc_id <- unique(chr_windows$nc_id)

    cat("  chromosome:", unique(chr_windows$chr), "\n")

    if (theta_input_mode == "raw") {
      idx_gen0 <- file.path(theta_dir, paste0(nc_id, "_", gen0_label, ".thetas.idx"))
      idx_gen14 <- file.path(theta_dir, paste0(nc_id, "_", gen14_label, ".thetas.idx"))

      gen0_obj <- run_theta_print(idx_gen0, theta_stat_bin)
      gen14_obj <- run_theta_print(idx_gen14, theta_stat_bin)

      gen0_summary <- summarise_theta_for_windows(
        theta_dt = gen0_obj$data,
        windows_dt = chr_windows,
        n_chr = gen0_obj$nChr,
        generation_label = "Gen0"
      )

      gen14_summary <- summarise_theta_for_windows(
        theta_dt = gen14_obj$data,
        windows_dt = chr_windows,
        n_chr = gen14_obj$nChr,
        generation_label = "Gen14"
      )
    } else {
      pestpg_gen0 <- resolve_pestpg_file(pestpg_dir, nc_id, gen0_label)
      pestpg_gen14 <- resolve_pestpg_file(pestpg_dir, nc_id, gen14_label)

      gen0_summary <- summarise_pestpg_for_windows(
        pestpg_dt = read_pestpg_summary(pestpg_gen0),
        windows_dt = chr_windows,
        generation_label = "Gen0"
      )

      gen14_summary <- summarise_pestpg_for_windows(
        pestpg_dt = read_pestpg_summary(pestpg_gen14),
        windows_dt = chr_windows,
        generation_label = "Gen14"
      )
    }

    chr_windows %>%
      left_join(gen0_summary, by = "window_row_id") %>%
      left_join(gen14_summary, by = "window_row_id")
  })

merged_tbl <- bind_rows(theta_summaries) %>%
  mutate(
    theta_summary_mode = theta_input_mode,
    deltaD = Tajima_Gen14 - Tajima_Gen0,
    pi_ratio = pi_Gen14 / pi_Gen0,
    log_pi_ratio = log2(pi_ratio)
  ) %>%
  filter(
    is.finite(mean_fst),
    is.finite(Tajima_Gen0),
    is.finite(Tajima_Gen14),
    is.finite(deltaD),
    is.finite(pi_Gen0),
    is.finite(pi_Gen14),
    is.finite(pi_ratio),
    is.finite(log_pi_ratio)
  ) %>%
  arrange(chr_num, window_id)

fst_threshold <- quantile(merged_tbl$mean_fst, fst_percentile, na.rm = TRUE)
deltaD_threshold <- quantile(merged_tbl$deltaD, tail_percentile, na.rm = TRUE)
log_pi_threshold <- quantile(merged_tbl$log_pi_ratio, tail_percentile, na.rm = TRUE)

merged_tbl <- merged_tbl %>%
  mutate(
    fst_outlier = mean_fst >= fst_threshold,
    deltaD_outlier = deltaD <= deltaD_threshold,
    pi_outlier = log_pi_ratio <= log_pi_threshold,
    all_three_signals = fst_outlier & deltaD_outlier & pi_outlier
  ) %>%
  group_by(chr_num) %>%
  mutate(
    Fst_smooth = running_median(mean_fst, k = 5),
    deltaD_smooth = running_median(deltaD, k = 5),
    log_pi_smooth = running_median(log_pi_ratio, k = 5)
  ) %>%
  ungroup()

candidate_regions <- merge_candidate_regions(
  merged_tbl %>% filter(all_three_signals),
  candidate_gap_bp = candidate_gap_bp
)

cat("\n=== Loading gene annotation ===\n")
reference_gr <- rtracklayer::import(gff_file)
genes_only <- reference_gr[mcols(reference_gr)$type == "gene"]
cat("Genes in reference:", length(genes_only), "\n")

annotated_tbl <- annotate_windows_with_genes(merged_tbl, genes_only) %>%
  mutate(
    Population = population_name,
    Comparison = comparison_label,
    theta_summary_mode = theta_summary_mode,
    fst_outlier = if_else(fst_outlier, "Yes", "No"),
    deltaD_outlier = if_else(deltaD_outlier, "Yes", "No"),
    pi_outlier = if_else(pi_outlier, "Yes", "No"),
    all_three_signals = if_else(all_three_signals, "Yes", "No"),
    threshold_label = "5pct"
  ) %>%
  select(
    Population,
    Comparison,
    theta_summary_mode,
    chr,
    chr_num,
    nc_id,
    window_id,
    window_start_bp,
    window_end_bp,
    window_mid_bp,
    mean_fst,
    Tajima_Gen0,
    Tajima_Gen14,
    deltaD,
    pi_Gen0,
    pi_Gen14,
    pi_ratio,
    log_pi_ratio,
    fst_outlier,
    deltaD_outlier,
    pi_outlier,
    all_three_signals,
    gene_count,
    GeneIDs,
    GeneNames,
    GeneBiotypes,
    snp_count,
    n_theta_sites_Gen0,
    n_theta_sites_Gen14,
    support_windows_Gen0,
    support_windows_Gen14,
    mean_window_nSites_Gen0,
    mean_window_nSites_Gen14,
    nChr_Gen0,
    nChr_Gen14
  ) %>%
  arrange(chr_num, window_id)

candidate_windows_tbl <- annotated_tbl %>%
  filter(all_three_signals == "Yes") %>%
  arrange(chr_num, window_id)

cat("\n=== Writing tables ===\n")
readr::write_csv(merged_tbl, window_table_file)
readr::write_csv(annotated_tbl, window_gene_file)
readr::write_csv(candidate_regions, candidate_regions_file)
readr::write_csv(candidate_windows_tbl, candidate_windows_file)

chr_info <- merged_tbl %>%
  group_by(chr_num) %>%
  summarise(chr_len = max(window_end_bp, na.rm = TRUE), .groups = "drop") %>%
  arrange(chr_num) %>%
  mutate(
    offset = lag(cumsum(chr_len), default = 0),
    center = offset + chr_len / 2
  )

thresholds <- list(
  fst = fst_threshold,
  deltaD = deltaD_threshold,
  log_pi = log_pi_threshold
)

panel_ic <- plot_population_panel(
  df = merged_tbl,
  regions_df = candidate_regions,
  thresholds = thresholds,
  chr_info = chr_info,
  pop_label = paste0(population_name, " spline selection scan: FST, Delta D, log2(pi14/pi0)")
)

ggsave(figure_file, panel_ic, width = 9.8, height = 10.5, dpi = 300, bg = "white")

cat("\nSaved files:\n")
cat("  ", window_table_file, "\n", sep = "")
cat("  ", window_gene_file, "\n", sep = "")
cat("  ", candidate_regions_file, "\n", sep = "")
cat("  ", candidate_windows_file, "\n", sep = "")
cat("  ", figure_file, "\n", sep = "")

cat("\nSummary:\n")
cat("  Windows:", nrow(merged_tbl), "\n")
cat("  Candidate windows:", nrow(candidate_windows_tbl), "\n")
cat("  Candidate regions:", nrow(candidate_regions), "\n")
