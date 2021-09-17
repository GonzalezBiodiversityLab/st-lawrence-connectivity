# a254
# Bronwyn Rayfield and Jed Lloren, ApexRMS
# Run with R-4.1.1

# Workspace---------------------------------------------------------------------------------------------
# Set environment variable TZ when running on AWS EC2 instance
Sys.setenv(TZ='GMT')
options(stringsAsFactors=FALSE)

# Load packages
library(fasterize)
library(raster)
library(sf)
library(tidyverse)

# Set up directories
# Assumes these directories already exist
rawMapsDir <- "../Data/Spatial"
rawTablesDir <- "../Inputs/RawData/Tables"

# Read in data
# Spatial
landcoverBTSLRaw <- st_read(dsn = file.path(rawMapsDir, "Extrait_Donnees.gdb"), layer="BTSL_SLL_Occ_sol_Land_cover")
siefDataRaw <- st_read(dsn = file.path(rawMapsDir, "Extrait_Donnees.gdb"), layer="SIEF_C08PEEFO")
studyArea <- raster(file.path(rawMapsDir, "study-area.tif"))

# Tabular
landcoverBTSLReclass <- read_csv(file.path(rawTablesDir, "landcoverBTSLReclass.csv"))
forestAgeReclass <- read_csv(file.path(rawTablesDir, "forestAgeReclass.csv"))
depositReclass <- read_csv(file.path(rawTablesDir, "depositReclass.csv"))

# Reformatting data---------------------------------------------------------------------------------------------
# Forest age
# Merge the layer and reclass tables
landcoverMerge <- landcoverBTSLRaw %>%
  left_join(landcoverBTSLReclass)
# Rasterize, crop, and mask to study area
landcoverBTSL <- fasterize(sf = st_cast(landcoverMerge, "MULTIPOLYGON"), 
                       raster = studyArea, # Snap to Outaouais connectivity boundaries
                       field = "Value") %>% 
  mask(., mask=studyArea)

# Forest age
forestAgeMerge <- siefDataRaw %>%
  left_join(forestAgeReclass)

forestAge <- fasterize(sf = st_cast(forestAgeMerge, "MULTIPOLYGON"), 
                       raster = studyArea,
                       field = "Value") %>% 
  mask(., mask=studyArea)
    
# Surficial deposits
depositMerge <- siefDataRaw %>%
  left_join(depositReclass)

surficialDeposits <- fasterize(sf = st_cast(depositMerge, "MULTIPOLYGON"), 
                       raster = studyArea,
                       field = "Value") %>% 
  mask(., mask=studyArea)

# Save outputs---------------------------------------------------------------------------------------------
# Landcover
writeRaster(landcoverBTSL, 
            file.path(rawMapsDir, "BTSL-SLL-Occ-sol-Land-cover.tif"), 
            overwrite = TRUE)
# Forest age
writeRaster(forestAge, 
            file.path(rawMapsDir, "SIEF-C08PEEFO-forest-age.tif"), 
            overwrite = TRUE)
# Surficial deposits
writeRaster(surficialDeposits, 
            file.path(rawMapsDir, "SIEF-C08PEEFO-surficial-deposits.tif"), 
            overwrite = TRUE)