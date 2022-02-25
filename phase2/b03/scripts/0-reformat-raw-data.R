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
resultsDir <- "b03/model-inputs"
rawMapsDir <- "b03/data/spatial"
rawTablesDir <- "b01b02/data/tables"

# Read in data
# Spatial
landcoverBTSLRaw <- st_read(dsn = file.path(rawMapsDir, "Extrait_Donnees.gdb"), layer="BTSL_SLL_Occ_sol_Land_cover")
siefDataRaw <- st_read(dsn = file.path(rawMapsDir, "Extrait_Donnees.gdb"), layer="SIEF_C08PEEFO")
OutaouaisExtent <- raster(file.path(rawMapsDir, "OutaouaisConnectivityExtent.tif"))
studyArea <- raster(file.path(resultsDir, "b03_studyArea_30m.tif"))
ontarioBoundary <- st_read(dsn = file.path(rawMapsDir, "ON Provincial Border"), layer = "Province")

# Tabular
landcoverBTSLReclass <- read_csv(file.path(rawTablesDir, "landcoverBTSLReclass.csv"))
forestAgeReclass <- read_csv(file.path(rawTablesDir, "forestAgeReclass.csv"))
depositReclass <- read_csv(file.path(rawTablesDir, "depositReclass.csv"))
forestDensityReclass <- read_csv(file.path(rawTablesDir, "speciesForestDensityReclass.csv"))

# Reformatting data---------------------------------------------------------------------------------------------
# Forest age
# Merge the layer and reclass tables
landcoverMerge <- landcoverBTSLRaw %>%
  left_join(landcoverBTSLReclass, by = c("CLASSE_DET" = "CLASSE_DET")) # Default joins by CLASSE_GEN and CLASSE_DET
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
  left_join(depositReclass, by = c("DEP_SUR" = "CODE"))

surficialDeposits <- fasterize(sf = st_cast(depositMerge, "MULTIPOLYGON"), 
                       raster = studyArea,
                       field = "Recode") %>% 
  mask(., mask=studyArea)

# Forest density
forestDensityMerge <- siefDataRaw %>%
  left_join(forestDensityReclass, by = c("CL_DENS" = "DensityName"))

forestDensity <- fasterize(sf = st_cast(forestDensityMerge, "MULTIPOLYGON"), 
                               raster = studyArea,
                               field = "DensityCode") %>% 
  mask(., mask=studyArea)

# Study area
studyAreaRaw <- studyAreaRaw %>%
  mutate(ID_Code = c(3, 4, 1, 2))
studyArea <- fasterize(sf = st_cast(studyAreaRaw, "MULTIPOLYGON"), 
                               raster = OutaouaisExtent,
                               field = "ID_Code") %>% 
  mask(., mask=OutaouaisExtent)

# Ontario boundary
ontarioBoundary <- st_transform(ontarioBoundary, crs = crs(studyArea))
ontarioBoundaryRaster <- fasterize(sf = st_cast(ontarioBoundary, "MULTIPOLYGON"),
                                   raster = studyArea,
                                   field = "OBJECTID") %>% 
  mask(., mask=studyArea)
  
  
  #mask(mask = studyArea)

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
# Forest density
writeRaster(forestDensity, 
            file.path(rawMapsDir, "SIEF-C08PEEFO-forest-density.tif"), 
            overwrite = TRUE)
# Study area
writeRaster(studyArea, 
            file.path(rawMapsDir, "study-area.tif"), 
            overwrite = TRUE)

# Ontario mask
writeRaster(ontarioBoundaryRaster, 
            file.path(rawMapsDir, "ontario-mask.tif"), 
            overwrite = TRUE)
