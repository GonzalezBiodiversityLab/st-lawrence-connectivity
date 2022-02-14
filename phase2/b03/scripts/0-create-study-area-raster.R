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
rawMapsDir <- "b03/data/spatial"
rawTablesDir <- "b03/data/tabular"

# Read in data
# Spatial
b03Raw <- st_read(file.path(rawMapsDir,"plaine-dottawa.shp"))
studyExtent<-raster(file.path(rawMapsDir,"OutaouaisConnectivityExtent.tif"))

# Tabular
reclassTable <- read_csv(file.path(rawTablesDir,"b03-study-area-reclass.csv"))

# Reformatting data---------------------------------------------------------------------------------------------
# Rasterize the Plaine d'Ottawa shapefile
studyArea <- fasterize(sf = st_cast(b03Raw, "POLYGON"), 
                       raster = studyExtent,
                       field = "FID02") %>% 
  mask(., mask=studyExtent)

# Reclassify
b03_studyarea <- reclassify(studyArea, reclassTable)

# Save outputs---------------------------------------------------------------------------------------------
# B03 study area
writeRaster(b03_studyarea, 
            file.path(rawMapsDir, "b03-studyArea.tif"), 
            overwrite = TRUE)