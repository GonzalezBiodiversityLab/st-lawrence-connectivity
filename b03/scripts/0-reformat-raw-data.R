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
rawMapsDir <- "../data/spatial"
rawTablesDir <- "../../b01b02/inputs/rawData/tables"

# Read in data
# Spatial
landcoverBTSLRaw <- st_read(dsn = file.path(rawMapsDir, "Extrait_Donnees.gdb"), layer="BTSL_SLL_Occ_sol_Land_cover")
siefDataRaw <- st_read(dsn = file.path(rawMapsDir, "Extrait_Donnees.gdb"), layer="SIEF_C08PEEFO")
studyAreaRaw <- st_read(dsn = file.path(rawMapsDir, "CR_CERQ.gdb"), layer="CR_NIV_01_S")
OutaouaisExtent <- raster(file.path(rawMapsDir, "OutaouaisConnectivityExtent.tif"))

# Tabular
landcoverBTSLReclass <- read_csv(file.path(rawTablesDir, "landcoverBTSLReclass.csv"))
forestAgeReclass <- read_csv(file.path(rawTablesDir, "forestAgeReclass.csv"))
depositReclass <- read_csv(file.path(rawTablesDir, "depositReclass.csv"))

# Reformatting data---------------------------------------------------------------------------------------------
# Forest age
# Merge the layer and reclass tables
landcoverMerge <- landcoverBTSLRaw %>%
  left_join(landcoverBTSLReclass, by = c("CLASSE_GEN" = "CLASSE_GEN")) # Default joins by CLASSE_GEN and CLASSE_DET
# Rasterize, crop, and mask to study area
landcoverBTSL <- fasterize(sf = st_cast(landcoverMerge, "MULTIPOLYGON"), 
                       raster = OutaouaisExtent, # Snap to Outaouais connectivity boundaries
                       field = "Value") %>% 
  mask(., mask=OutaouaisExtent)

# Forest age
forestAgeMerge <- siefDataRaw %>%
  left_join(forestAgeReclass)

forestAge <- fasterize(sf = st_cast(forestAgeMerge, "MULTIPOLYGON"), 
                       raster = OutaouaisExtent,
                       field = "Value") %>% 
  mask(., mask=OutaouaisExtent)
    
# Surficial deposits
depositMerge <- siefDataRaw %>%
  left_join(depositReclass)

surficialDeposits <- fasterize(sf = st_cast(depositMerge, "MULTIPOLYGON"), 
                       raster = OutaouaisExtent,
                       field = "Value") %>% 
  mask(., mask=OutaouaisExtent)

# Study area
studyAreaRaw <- studyAreaRaw %>%
  mutate(ID_Code = c(3, 4, 1, 2))
studyArea <- fasterize(sf = st_cast(studyAreaRaw, "MULTIPOLYGON"), 
                               raster = OutaouaisExtent,
                               field = "ID_Code") %>% 
  mask(., mask=OutaouaisExtent)

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
# Study area
writeRaster(studyArea, 
            file.path(rawMapsDir, "study-area.tif"), 
            overwrite = TRUE)