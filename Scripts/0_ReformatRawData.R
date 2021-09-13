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
siefDataRaw <- st_read(dsn = file.path(rawMapsDir, "Extrait_Donnees.gdb"), layer="SIEF_C08PEEFO")
studyArea <- raster(file.path(rawMapsDir, "BTSLOutaouais.tif"))

# Tabular
forestAgeReclass <- read_csv(file.path(rawTablesDir, "forestAgeReclass.csv"))
depositReclass <- read_csv(file.path(rawTablesDir, "depositReclass.csv"))

# Merge the layer and reclass tables according to 1) forest age and 2) surficial deposits
forestAgeMerge <- merge(siefDataRaw, forestAgeReclass, all.x = TRUE, all.y = FALSE)
depositMerge <- merge(siefDataRaw, depositReclass, all.x = TRUE, all.y = FALSE)

# Rasterize, crop, and mask to study area
# Forest age
forestAge <- fasterize(sf = st_cast(forestAgeMerge, "MULTIPOLYGON"), 
                   raster = studyArea, # Snap to evt grid
                   field = "Code") %>% 
  mask(., mask=studyArea)

# Surficial deposits
surficialDeposits <- fasterize(sf = st_cast(depositMerge, "MULTIPOLYGON"), 
                       raster = studyArea, # Snap to evt grid
                       field = "Recode") %>% 
  mask(., mask=studyArea)

# Save outputs---------------------------------------------------------------------------------------------
# Forest age
writeRaster(forestAge, 
            file.path(rawMapsDir, "forestAge.tif"), 
            overwrite = TRUE)

# Surficial deposits
writeRaster(surficialDeposits, 
            file.path(rawMapsDir, "surficialDeposits.tif"), 
            overwrite = TRUE)