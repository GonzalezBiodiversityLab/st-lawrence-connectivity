# a254
# Bronwyn Rayfield and Jed Lloren, ApexRMS
# Run with R-1.1

# Workspace---------------------------------------------------------------------------------------------
# Load packages
library(fasterize)
library(raster)
library(sf)
library(tidyverse)

# Set up directories
# Assumes these directories already exist
rawMapsDir <- "../data/spatial"
rawTablesDir <- "../data/tabular"

# Read in data
# Spatial
AP_Raw <- st_read(file.path(rawMapsDir,"AP_REG_S_20210824.shp"))
StudyExtent<-raster(file.path(rawMapsDir,"b03-studyarea.tif"))

# Tabular
reclassTable <- read_csv(file.path(rawTablesDir,"protected-area-reclass.csv"))

# Reformatting data---------------------------------------------------------------------------------------------
# Reclass the study extent to turn 0 values to NA
StudyExtent[StudyExtent < 1] <- NA

# Rasterize the Plaine d'Ottawa shapefile
AP-b03Raw <- fasterize(sf = st_cast(AP_Raw, "POLYGON"), 
                       raster = StudyExtent,
                       field = "OBJECTID") %>% 
  mask(., mask=StudyExtent)

# Reclassify
AP-b03 <- reclassify(AP-b03Raw, reclassTable)

# Save outputs---------------------------------------------------------------------------------------------
# B03 study area
writeRaster(AP-b03, 
            file.path(rawMapsDir, "AP-b03.tif"), 
            overwrite = TRUE)