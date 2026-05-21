####  Libraries
library(ggplot2)
library(maps)
library(ggspatial)
library(dplyr)
library(raster)
library(sp)
library(RColorBrewer) 
library(GGally)
library(scatterplot3d)
library(plotly)
library(fields)
library(viridis)
library(tidyverse)
library(mapdata)
library(patchwork)
library(vroom)
library(grid)     
library(terra)

# Load and filter data
dir_maize <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Phenotype/"
maize_geo <- read.csv(paste0(dir_maize, "taxa_geoloc_pheno.csv")) %>%
  dplyr::filter(
    sp == "Zea mays",
    GEO3 %in% c("Caribbean", "Meso-America", "South America")
  ) %>%
  dplyr::select(2, 6, 7)
names(maize_geo) <- c("Genotypes", "long", "lat")


# Soil N values
data_values <- vroom("/Users/nirwantandukar/Documents/Research/data/Phenotypes/maize_N_NHx.csv")
colnames(data_values)[1] <- "Genotypes"

# Combine
maize_geo <- maize_geo %>%
  left_join(data_values, by = "Genotypes") %>%
  dplyr::select(Genotypes, long, lat, TN_maize)


# #### FOR ALTITUDE
# # # Load the raster
# # ## USED THIS: https://asterweb.jpl.nasa.gov/gdem.asp
# # # DEM: get global elevation (~1 km) and aggregate to ~10 km (~0.083333°)
# library(terra)
# library(geodata)   # install.packages("geodata") if needed
# 
# # 5 arc-min (~0.083333°) global DEM in meters
# dem_5min <- geodata::elevation_global(res = "5", path = tempdir())  # SpatRaster, 1 band
# 
# # If you specifically want 30″ first, then aggregate to 5′:
# # dem_30s  <- geodata::elevation_global(res = "0.5", path = tempdir())
# # dem_5min <- aggregate(dem_30s, fact = 10, fun = mean, na.rm = TRUE)
# 
# # Your lon/lat points (WGS84)
# pts <- vect(maize_geo, geom = c("long","lat"), crs = "EPSG:4326")
# 
# # Extract elevation (meters); bilinear gives smoother values
# elev <- terra::extract(dem_5min, pts, method = "bilinear")
# 
# # Bind back (drop ID column)
# maize_geo$elev_m <- elev[,2]
# 
# # peek
# head(maize_geo[, c("Genotypes","long","lat","elev_m","TN_maize")])




# Convert to spatial points
coordinates(maize_geo) <- ~long + lat
proj4string(maize_geo) <- CRS("+init=epsg:4326")

# Convert back to data frame for ggplot
maize_df <- as.data.frame(maize_geo)
# maize_df <- maize_df[,-4]
# colnames(maize_df)[4] <- "TN_AT"

# Use raw Soil N values (not log-transformed)
maize_df$soil_n <- maize_df$TN_maize
str(maize_df)


# helper: margin() if available, otherwise fall back to grid::unit
marg <- function(top = 0, right = 0, bottom = 0, left = 0, unit = "pt") {
  if (exists("margin", envir = asNamespace("ggplot2"))) {
    ggplot2::margin(top, right, bottom, left)
  } else {
    grid::unit(c(top, right, bottom, left), unit)
  }
}

plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title     = element_text(
      size   = 14,
      face   = "bold",
      hjust  = 0.5,
      margin = margin(b = 10)
    ),
    axis.title.x   = element_text(
      size = 16,
      face = "bold"
    ),
    axis.title.y   = element_text(
      size = 16,
      face = "bold"
    ),
    axis.text.x    = element_text(
      size = 16,
      color = "black"
    ),
    axis.text.y    = element_text(
      size = 16,
      color = "black"
    ),
    axis.line      = element_line(color = "black", linewidth = 0.8),
    axis.ticks     = element_line(color = "black", linewidth = 0.7),
    panel.grid     = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "grey70", linewidth = 0.4),
    legend.direction = "vertical",
    legend.spacing.y = unit(0.2, "cm"),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    plot.margin = margin(15, 15, 15, 15)
  )

## ─────────────────────────────────────────────────────────
##  A) MAIN MAP
## ─────────────────────────────────────────────────────────
main_map <- ggplot() +
  borders(
    "world",
    xlim   = c(-120, -30),
    ylim   = c(-50,  40),
    fill   = "gray95",
    colour = "gray70",
    linewidth = 0.2
  ) +
  geom_point(
    data  = maize_df,
    aes(long, lat, colour = soil_n),
    size  = 0.8,
    alpha = 0.8
  ) +
  coord_fixed(ratio = 1.3, xlim = c(-120, -30), ylim = c(-50, 40)) +
  scale_colour_viridis_c(
    option = "D",  # "D" = viridis, "C" = magma, "B" = cividis
    name   = "Soil N",
    guide  = guide_colourbar(
      direction    = "vertical",
      barwidth     = unit(0.4, "cm"),
      barheight    = unit(4, "cm"),
      ticks.colour = "black",
      title.position = "top",
      title.hjust  = 0.5
    )
  ) +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title       = element_text(size = 16, face = "bold"),
    axis.text        = element_text(size = 16, color = "black"),
    axis.line        = element_line(linewidth = 0.8, colour = "black"),
    axis.ticks       = element_line(linewidth = 0.7, colour = "black"),
    legend.position  = "inside",
    legend.position.inside = c(0.9, 0.2),
    legend.justification = "left",
    legend.title     = element_text(size = 8, face = "bold"),
    legend.text      = element_text(size = 7),
    plot.title       = element_text(face = "bold", size = 11, hjust = 0),
    plot.subtitle    = element_text(size = 9, hjust = 0, margin = marg(bottom = 6)),
    plot.margin      = marg(6, 12, 6, 6)
  )

## ─────────────────────────────────────────────────────────
##  B) HISTOGRAM (larger inset, no grids)
## ─────────────────────────────────────────────────────────

# Calculate quantiles
quantiles <- quantile(maize_df$soil_n, probs = c(0.25, 0.75), na.rm = TRUE)

hist_plot <- ggplot(maize_df, aes(soil_n, fill = after_stat(x))) +
  geom_histogram(bins = 25, color = "black", linewidth = 0.5) # +
#Add vertical lines for quantiles
# geom_vline(xintercept = quantiles,
#            color = "red",
#            linetype = "dashed",
#            size = 0.8,
#            alpha = 0.7)

hist_plot <- hist_plot +
  # Add text annotations
  # annotate("text",
  #          x = quantiles[1],
  #          y = max(ggplot_build(hist_plot)$data[[1]]$count) * 0.95,
  #          label = "Bottom \n 25%",
  #          color = "red",
  #          hjust = 1.1,
  #          size = 2) +
  #  annotate("text",
  #          x = quantiles[2],
  #          y = max(ggplot_build(hist_plot)$data[[1]]$count) * 0.95,
#          label = "Top \n 25%",
#          color = "red",
#          hjust = -0.1,
#          size = 2) +
scale_fill_viridis_c(option = "D", guide = "none") +
  labs(x = "Soil N", y = "Count") +
  theme_minimal(base_size = 10) +
  theme(
    axis.title = element_text(size = 8, face = "bold"),
    axis.text = element_text(size = 7, color = "black"),
    axis.line = element_line(linewidth = 0.7, color = "black"),
    axis.ticks = element_line(linewidth = 0.6, color = "black"),
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  ) 


## ─────────────────────────────────────────────────────────
##  C) COMBINE & SAVE
## ─────────────────────────────────────────────────────────
combined_plot <- main_map +
  inset_element(
    hist_plot,
    left   = 0.001,
    bottom = 0.01,
    right  = 0.52,   # bigger
    top    = 0.44
  )
# 
# 
# ggsave(
#   "SoilN_Map_Landraces_Maize_Romerro.png",
#   plot   = combined_plot,
#   width  = 8,
#   height = 6,
#   dpi    = 300,
#   units = "in",
#   bg     = "white"
# )

getwd()



################################################################################
### LIBRARIES
################################################################################
# Packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(GenomicRanges)
library(ggrepel)
library(viridisLite)
library(tidyverse)
library(vroom)
library(data.table)
library(purrr)
library(qqman)
library(stringr)
library(scales)
library(ggsci)
library(readr)
library(VennDiagram)
library(grid)
library(ggnewscale)  # For multiple color scales in ggplot


# NOTE:
# MAIZE and SORGHUM MANHATTAN.

################################################################################
################################################################################
### MAIZE
################################################################################
################################################################################

################################################################################
### GETTING RAW RESULTS FOR MLM, MLMM, BLINK
################################################################################
mlm_raw <- vroom("/Users/nirwantandukar/Documents/Research/results/GWAS/GAPIT/raw_GWAS_MLM_3PC_N.txt",
                 col_names = TRUE, delim = "\t") %>% dplyr::select(SNP, Chr, Pos, P.value)

BLINK_raw <- vroom("/Users/nirwantandukar/Documents/Research/results/GWAS/GAPIT/raw_GWAS_BLINK_3PC_N.txt",
                   col_names = TRUE, delim = "\t") %>% dplyr::select(SNP, Chr, Pos, P.value)

MLMM_raw <- vroom("/Users/nirwantandukar/Documents/Research/results/GWAS/GAPIT/raw_GWAS_MLMM_3PC_N.txt",
                  col_names = TRUE, delim = "\t") %>% dplyr::select(SNP, Chr, Pos, P.value)

################################################################################
### COMBINING THE FARMCPU RESULTS
################################################################################

# Initialize empty list to hold per-chromosome results
gwas_list <- list()

# Loop through chr1 to chr10
for (chr in 1:10) {
  # Construct the file path
  file_path <- paste0("/Users/nirwantandukar/Documents/Research/results/GWAS/GAPIT/FarmCPU/FarmCPU_TN_3PC_maize_chr", chr, ".rds")
  
  # Read the RDS file
  x <- readRDS(file_path)
  
  # Extract the GWAS results table for TN_maize
  gwas_chr <- x$TN_maize$GWAS
  
  # Append to the list
  gwas_list[[chr]] <- gwas_chr
}

# Combine all into a single data frame
farmcpu_raw <- do.call(rbind, gwas_list)

# Optional: check result
farmcpu_raw <- farmcpu_raw %>%
  dplyr::select(SNP, Chromosome, Position, p.value) %>%
  #rename the columns to match the others
  dplyr::rename(Chr = Chromosome, Pos = Position, P.value = p.value)

farmcpu_raw$log10 <- -log10(farmcpu_raw$P.value)


################################################################################
### COMBINING THE RESULTS
################################################################################

# First, rename columns to standard format for merging
mlm_df     <- mlm_raw     |> dplyr::rename(SNP = SNP, Chromosome = Chr, Position = Pos, MLM     = P.value)
mlmm_df    <- MLMM_raw    |> dplyr::rename(SNP = SNP, Chromosome = Chr, Position = Pos, MLMM    = P.value)
blink_df   <- BLINK_raw   |> dplyr::rename(SNP = SNP, Chromosome = Chr, Position = Pos, BLINK   = P.value)
farmcpu_df <- farmcpu_raw |> dplyr::rename(SNP = SNP, Chromosome = Chr, Position = Pos, FarmCPU = P.value)

# Merge all by SNP, Chromosome, and Position
combined_df <- mlm_df %>%
  full_join(mlmm_df,   by = c("SNP", "Chromosome", "Position")) %>%
  full_join(blink_df,  by = c("SNP", "Chromosome", "Position")) %>%
  full_join(farmcpu_df, by = c("SNP", "Chromosome", "Position"))

# Reorder columns: SNP, Chromosome, Position, then models
combined_df <- combined_df |> dplyr::select(SNP, Chromosome, Position, MLM, MLMM, BLINK, FarmCPU)


################################################################################
##### MANHATTAN PLOT WITH GGPLOT2 (Two threshold lines: 5 and 7)
################################################################################

# 0) Clean obvious NAs (recommended)
combined_df <- combined_df %>%
  filter(!is.na(Chromosome), !is.na(Position))

# 1) For each SNP, find which model has the LOWEST p-value (best signal)
best_per_snp <- combined_df %>%
  rowwise() %>%
  mutate(
    # Get the minimum p-value across models (ignoring NAs)
    min_p = min(c(MLM, MLMM, BLINK, FarmCPU), na.rm = TRUE),
    # Determine which model had that minimum p-value
    best_model = case_when(
      !is.na(MLM) & MLM == min_p ~ "MLM",
      !is.na(MLMM) & MLMM == min_p ~ "MLMM",
      !is.na(BLINK) & BLINK == min_p ~ "BLINK",
      !is.na(FarmCPU) & FarmCPU == min_p ~ "FarmCPU",
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  filter(!is.na(best_model), is.finite(min_p)) %>%
  mutate(
    CHR = as.integer(Chromosome),
    BP = as.numeric(Position),
    P = min_p,
    logp = -log10(min_p)
  ) %>%
  select(SNP, CHR, BP, P, logp, best_model) %>%
  arrange(CHR, BP)

# Check distribution of best models
cat("\n=== Best Model Distribution (all SNPs) ===\n")
print(table(best_per_snp$best_model))

# 2) Chromosome lengths + offsets for cumulative positions
chr_tbl <- best_per_snp %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP, na.rm = TRUE), .groups = "drop") %>%
  arrange(CHR) %>%
  mutate(offset = lag(cumsum(chr_len), default = 0),
         center = offset + chr_len / 2)

# 3) Add cumulative positions
plot_df <- best_per_snp %>%
  left_join(chr_tbl %>% select(CHR, offset), by = "CHR") %>%
  mutate(pos_cum = BP + offset)

# 4) Define thresholds
thr_suggestive <- 5  # -log10(1e-5)
thr_significant <- 7  # -log10(1e-7)

# Check distribution for significant SNPs
cat("\n=== Best Model Distribution (significant SNPs, logp >= 5) ===\n")
print(table(plot_df$best_model[plot_df$logp >= thr_suggestive]))

# 5) Define colors for each GWAS model
model_colors <- c(
  "MLM"     = "#0072B2",  # Blue
  "MLMM"    = "#E69F00",  # Orange
  "BLINK"   = "#009E73",  # Bluish green
  "FarmCPU" = "#D55E00"   # Vermillion
)

# 8) Separate data for background (non-significant) and significant SNPs
bg_df <- plot_df %>% filter(logp < thr_suggestive)
sig_df <- plot_df %>% filter(logp >= thr_suggestive)

# 9) Create Manhattan plot with TWO threshold lines and MODEL-colored SNPs
p_manhattan <- ggplot() +
  # Background points (non-significant, alternating chromosome colors)
  geom_point(data = bg_df, aes(x = pos_cum, y = logp, color = factor(CHR %% 2)),
             alpha = 0.5, size = 1, show.legend = FALSE) +
  scale_color_manual(values = c("0" = "grey65", "1" = "grey40")) +
  
  # New color scale for significant SNPs colored by best model
  ggnewscale::new_scale_color() +
  
  # Significant SNPs colored by their best model
  
  geom_point(data = sig_df, aes(x = pos_cum, y = logp, color = best_model),
             size = 2.5, alpha = 0.85) +
  scale_color_manual(
    values = model_colors,
    name = "Best Model",
    limits = names(model_colors)
  ) +
  
  # Two threshold lines
  geom_hline(yintercept = thr_suggestive, linetype = "dashed",
             color = "black", linewidth = 0.7) +
  geom_hline(yintercept = thr_significant, linetype = "solid",
             color = "black", linewidth = 0.9) +
  
  # Annotations for threshold lines
  annotate("text", x = max(plot_df$pos_cum) * 0.02, y = thr_suggestive + 0.4,
           label = "-log10(p) = 5", hjust = 0, size = 4, color = "black", fontface = "italic") +
  annotate("text", x = max(plot_df$pos_cum) * 0.02, y = thr_significant + 0.4,
           label = "-log10(p) = 7", hjust = 0, size = 4, color = "black", fontface = "italic") +
  
  # Chromosome labels
  scale_x_continuous(breaks = chr_tbl$center, labels = chr_tbl$CHR, expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  
  # Labels
  labs(
    x = "Chromosome",
    y = expression(-log[10](italic(p)))
  ) +
  plot_theme +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.99, 0.99),
    legend.justification = c("right", "top"),
    legend.box = "horizontal",
    legend.background = element_rect(fill = "white", color = "grey70", linewidth = 0.4),
    plot.margin = margin(42, 15, 15, 15)
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))

# Display side-by-side figure
fig1_left_panel <- patchwork::wrap_elements(panel = combined_plot)

fig1_combined <- (fig1_left_panel | p_manhattan) +
  plot_layout(widths = c(1, 2.2)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 18, face = "bold"),
    plot.tag.position = c(0.02, 0.98)
  )

output_file <- "Figs/main/Fig1.png"
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
ggsave(output_file, plot = fig1_combined, width = 16, height = 7, dpi = 300, bg = "white")

if (interactive()) {
  quartz(width = 16, height = 7)
  print(fig1_combined)
}

# ggsave("results/figures/Fig_1_map_manhattan_combined.png",
#        plot = fig1_combined,
#        width = 16, height = 7, dpi = 300, bg = "white")
# 
# ggsave("results/figures/Fig_1_map_manhattan_combined.pdf",
#        plot = fig1_combined,
#        width = 16, height = 7, dpi = 300, bg = "white")
# 





