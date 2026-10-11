#!/usr/bin/env Rscript
################################################################################
## Build the control-phenotype table for the maize landrace panel.
##
##   Rscript make_control_phenotypes.R
##
## Produces one CSV with every trait GWAS_models.R can run, in the layout GAPIT
## wants (column 1 = taxa):
##
##   Taxa, latitude, abs_latitude, longitude, elevation,
##   pH, CEC, aridity, N_0_5, N_5_15
##
## latitude / longitude / elevation / soil N come from files that already exist.
## pH, CEC and aridity are extracted from the same rasters the soil-N phenotype
## came from, at the collection coordinates -- so they are built the same way as
## the trait they are controls for.  Soil pH and CEC are taken at 0-5 cm to match
## the soil-N layer that was scanned.
################################################################################

suppressPackageStartupMessages({ library(data.table) })

gv <- function(k, d) { v <- Sys.getenv(k); if (nzchar(v)) v else d }

PHDIR  <- gv("GN_DIR_PHENO", path.expand("~/Documents/Research/data/Phenotypes"))
RASTER <- gv("GN_DIR_RASTER", path.expand("~/Documents/Research/data/P_prediction/predictors"))
AIRAST <- gv("GN_AI_RASTER", path.expand(
  "~/Documents/Research/data/P_prediction/Global-AI_ET0_v3_annual/ai_v3_yr.tif"))
OUT    <- gv("GN_OUT", path.expand(
  "~/Documents/Research/data/Phenotypes/control_phenotypes_maize_landrace.csv"))

## ---- the pieces that already exist ----------------------------------------
ll <- fread(file.path(PHDIR, "navarro_longlat.csv"), encoding = "UTF-8")
setnames(ll, 1, "Taxa")
el <- fread(file.path(PHDIR, "maize_elevation.csv")); setnames(el, c("Taxa", "elevation"))
ni <- fread(file.path(PHDIR, "maize_nitrogen.csv"));  setnames(ni, 1, "Taxa")
setnames(ni, c("nitrogen_0-5cm_mean_1000_maize", "nitrogen_5-15cm_mean_1000_maize"),
         c("N_0_5", "N_5_15"), skip_absent = TRUE)

d <- Reduce(function(a, b) merge(a, b, by = "Taxa", all = TRUE),
            list(ll[, .(Taxa, longitude = lon, latitude = lat)],
                 el[, .(Taxa, elevation)],
                 ni[, .(Taxa, N_0_5, N_5_15)]))
d[, abs_latitude := abs(latitude)]
cat(sprintf("taxa with coordinates: %d\n", d[!is.na(latitude), .N]))

## ---- pH, CEC, aridity from the rasters -------------------------------------
if (!requireNamespace("terra", quietly = TRUE)) {
  stop("terra is needed to extract pH / CEC / aridity.\n",
       '  install.packages("terra")\n',
       "Everything else in this table does not need it -- rerun once terra is in.")
}
suppressPackageStartupMessages(library(terra))

RAS <- c(pH      = file.path(RASTER, "phh2o_0-5cm_mean_1000.tif"),
         CEC     = file.path(RASTER, "cec_0-5cm_mean_1000.tif"),
         aridity = AIRAST)

pts <- terra::vect(as.data.frame(d[!is.na(latitude) & !is.na(longitude),
                                   .(Taxa, longitude, latitude)]),
                   geom = c("longitude", "latitude"), crs = "EPSG:4326")

for (nm in names(RAS)) {
  f <- RAS[[nm]]
  if (!file.exists(f)) {
    warning("raster not found, column left empty: ", f); d[[nm]] <- NA_real_; next
  }
  r <- terra::rast(f)
  v <- terra::extract(r, pts, method = "bilinear")[, 2]
  tmp <- data.table(Taxa = pts$Taxa, val = as.numeric(v))
  d <- merge(d, tmp, by = "Taxa", all.x = TRUE)
  setnames(d, "val", nm)
  cat(sprintf("%-8s extracted for %d taxa  (range %.3g to %.3g)\n",
              nm, sum(!is.na(d[[nm]])), min(d[[nm]], na.rm = TRUE), max(d[[nm]], na.rm = TRUE)))
}

setcolorder(d, c("Taxa", "latitude", "abs_latitude", "longitude", "elevation",
                 "pH", "CEC", "aridity", "N_0_5", "N_5_15"))
fwrite(d, OUT)
cat("\nwritten: ", OUT, "\n", sep = "")

## ---- how much does each control overlap soil N? -----------------------------
cat("\ncorrelation with soil N 0-5 cm (the trait these are controls for):\n")
for (nm in c("latitude","abs_latitude","longitude","elevation","pH","CEC","aridity","N_5_15")) {
  ok <- complete.cases(d[[nm]], d$N_0_5)
  cat(sprintf("  %-13s r = %+.3f   (n = %d)\n", nm,
              cor(d[[nm]][ok], d$N_0_5[ok]), sum(ok)))
}
