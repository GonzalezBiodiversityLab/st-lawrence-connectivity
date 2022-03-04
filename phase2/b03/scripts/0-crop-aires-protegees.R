# a254
# Bronwyn Rayfield and Jed Lloren, ApexRMS
# Run with R-1.1

# Workspace ----
# Load packages
library(fasterize)
library(raster)
library(sf)
library(tidyverse)

# Set up directories
# Assumes these directories already exist
rawMapsDir <- "b03/data/spatial"
# rawMaps2Dir <- "../../combined/maps/protected-areas"
rawTablesDir <- "b01b02/data/tables"

# Read in data
# Spatial
protectedAreasRaw <- st_read(file.path(rawMapsDir,"AP_REG_S_20210824.shp"))
# ConsIntRaw <- st_read(dsn = file.path(rawMaps2Dir, "Territoires_interet_BTSL_Juin2019.gdb"), layer = "SitesMulticible")
studyareaRaster <-raster(file.path(rawMapsDir,"b03-studyarea-30m.tif"))

# Reformatting data ----
# Reclass the study extent to turn 0 values to NA
studyareaRaster[studyareaRaster < 1] <- NA

# Rasterize the protected area and conservation interest shapefiles
# Protected areas
protectedAreasRawb03 <- fasterize(sf = st_cast(protectedAreasRaw, "MULTIPOLYGON"), 
                       raster = studyareaRaster,
                       field = "DESIG_NO") %>% 
  mask(., mask=studyareaRaster)

# Conservation interest
# ConsInt <- fasterize(sf = st_cast(ConsIntRaw, "MULTIPOLYGON"), 
#                            raster = studyareaRaster,
#                            field = "SUM") %>% 
#   mask(., mask=studyareaRaster)

# Filter acquatic environments
protectedAreasRawb03[protectedAreasRawb03 == 15] <- NA

# Reclassify
# Protected areas
reclassTable = matrix(c(1, 7, 9, 15, 16, 21, 29, 134, 1 , 1, 1, NA, 1, 1, 1, 1), 
                      nrow = 8, 
                      ncol = 2,)

protectedAreasb03 <- reclassify(protectedAreasRawb03, reclassTable)

# Conservation interest
# ConsInt[ConsInt > 0] <- 1

# Save outputs ----
# B03 protected areas
writeRaster(protectedAreasb03, 
            file.path(rawMapsDir, "b03-protectedAreas.tif"), 
            overwrite = TRUE)
# B03 conservation interest
# writeRaster(ConsInt, 
#             file.path(rawMapsDir, "sites-multicible.tif"), 
#             overwrite = TRUE)