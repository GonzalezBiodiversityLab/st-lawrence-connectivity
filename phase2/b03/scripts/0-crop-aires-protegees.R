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

# Reformatting data---------------------------------------------------------------------------------------------
# Reclass the study extent to turn 0 values to NA
StudyExtent[StudyExtent < 1] <- NA

# Rasterize the Plaine d'Ottawa shapefile
AP_b03Raw <- fasterize(sf = st_cast(AP_Raw, "MULTIPOLYGON"), 
                       raster = StudyExtent,
                       field = "DESIG_NO") %>% 
  mask(., mask=StudyExtent)

# Filter acquatic environments
AP_b03Raw[AP_b03Raw == 15] <- NA

# Reclassify
reclassTable = matrix(
  c(1, 7, 9, 15, 16, 21, 29, 134, 1 , 1, 1, NA, 1, 1, 1, 1), nrow = 8, ncol = 2,)
AP_b03 <- reclassify(AP_b03Raw, reclassTable)

# Save outputs---------------------------------------------------------------------------------------------
# B03 study area
writeRaster(AP_b03, 
            file.path(rawMapsDir, "AP-b03.tif"), 
            overwrite = TRUE)