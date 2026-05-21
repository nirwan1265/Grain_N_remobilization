################################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# GWAS
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
################################################################################
### WEBSITES
# Places to download files
# https://www.isric.org/explore/soilgrids/soilgrids-access
# https://geo.nsstc.nasa.gov/SPoRT/modeling/lis/conus3km/geotiff/vsm_percentiles/
# https://nsidc.org/data/spl3sma/versions/3
# https://n5eil01u.ecs.nsidc.org/SMAP/SPL3SMA.003/2015.07.07/
# https://mygeodata.cloud/converter/hdf5-to-geotiff
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9877027/
# https://files.isric.org/public/af250m_nutrient/
# https://esdac.jrc.ec.europa.eu/themes/global-phosphorus
# https://esdac.jrc.ec.europa.eu/resource-type/datasets
# https://daac.ornl.gov/SOILS/guides/Global_Soil_Regolith_Sediment.html
#@ https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=830 - nitrogen deposition
# https://esdac.jrc.ec.europa.eu/content/lucas-2009-topsoil-data#tabs-0-description=1

# Soil Phosphorus retention
# https://files.isric.org/public/other/

################################################################################
### PACKAGES LOAD
################################################################################

# Load packages
library(vroom)
library(bigmemory)
library(biganalytics)
library(compiler)
library(dplyr)
library(qtl2)
library(GAPIT)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# MAIZE
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

################################################################################
### LOAD THE DATA AND PCA
################################################################################

# Read phenotype
myY_N <- read.csv("/Users/nirwantandukar/Documents/Research/data/Phenotypes/maize_N_NHx.csv", header = TRUE)
myY_N <- myY_N[!is.na(myY_N$NHx_maize), ] 
myY_N <- myY_N[,-2]
samples_to_keep <- myY_N$X

# Read PCA and subset to common individuals
pca <- read.csv("/Users/nirwantandukar/Documents/Github/Grain_Nitrogen/results/PCA/PCA_maize_romerro.csv", header = TRUE)
pca <- pca[pca$sample.id %in% samples_to_keep, ]
pca <- pca[, 1:4]  # Use top 3 PCs based on the scree plot


# Loop through chromosomes for GAPIT models
folder_hapmap <- "/Users/nirwantandukar/Documents/Research/data/maize_genotype/RomeroNavarro/hapmap/"
for (chr in 1:10) {
  out_dir=paste0(getwd(),"/MLM/NHx/chr",chr)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  setwd(out_dir)
  hapmap_file <- paste0(folder_hapmap,"chr", chr, ".txt")
  
  myG_N <- read.delim(
    hapmap_file,
    header       = FALSE,   # the HapMap header starts with ‘rs#’
    sep          = "\t",   # HapMap files are TAB‑separated
    quote        = "",     # genotype codes don’t use quotes
    comment.char = "",     # prevents ‘#’ being treated as comment
    check.names  = FALSE,  # keep sample names unchanged
    stringsAsFactors = FALSE
  )
  
  
  myGAPIT <- GAPIT(
    Y=myY_N,
    G=myG_N,
    CV = pca,
    model = c("MLM"),
    file.output = TRUE
  )
  
  # Optional: save results manually if needed
  saveRDS(myGAPIT, file = paste0("GAPIT_MLM_chr", chr, ".rds"))
}