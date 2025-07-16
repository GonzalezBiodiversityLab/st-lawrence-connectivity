# a254
# Sarah Chisholm, ApexRMS
#
# This script mosaics linkage priority rasters across the b01b02 and b03 regions
# for each focal species. Overlapping cells between regions are filled with the
# maximum value. 

## Workspace ----

# Load constants, functions, etc
source("./b03/scripts/0-0-constants.R")

# Parameters
# Template raster
templateRaster <- raster()

# Target crs
targetCrs <- "+proj=lcc +lat_0=44 +lon_0=-68.5 +lat_1=46 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

# Combined b01b02 and b03 extent
targetExtent <- c(-692580, -95520, 116730, 403680) 

# Target resolution
targetResolution <- 30

# Empty raster stack
speciesCorridorRasterStack <- stack()

## Mosaic b01b02 and b03 corridor rasters ----
# Create a template raster at the combined b01b02 and b03 extent with a 30m resolution
extent(templateRaster) <- targetExtent
crs(templateRaster) <- targetCrs
res(templateRaster) <- targetResolution

# Loop over species list to mosaic rasters
for(species in speciesList){

  # Get b03 and b01b02 linkage priority rasters
  b03linkagePriorityRawRaster <- raster(file.path(b03linkagePriorityDir, "v2", str_c("LMproj_", species, "240_linkage_priority.tif")))
  b01b02linkagePriorityRawRaster <- raster(file.path(b01b02linkagePriorityDir, "LinkagePriority_30m", str_c("Corridors_30m_", species, "_linkage_priority.tif")))
  
  # Adjust crs, resolution, and extent
  b03linkagePriorityRaster <- b03linkagePriorityRawRaster %>% 
    projectRaster(to = templateRaster, method = 'ngb')

  b01b02linkagePriorityRaster <- b01b02linkagePriorityRawRaster %>% 
    projectRaster(to = templateRaster, method = 'ngb')

  # Mosaic rasters together and write to disk
  # Maximum value takes precedence in overlapping cells
  mosaic(
    x = b03linkagePriorityRaster,
    y = b01b02linkagePriorityRaster,
    fun = max,
    filename = file.path(b03linkagePriorityDir, str_c(species, "_max_linkage_priority_30m_BTSL.tif")))
  
  # Set NAs to 0
  btslCorridorRaster <- raster(file.path(b03linkagePriorityDir, str_c(species, "_max_linkage_priority_30m_BTSL.tif")))
  btslCorridorRaster[is.na(btslCorridorRaster)] <- 0
  names(btslCorridorRaster) <- species
  
  # Add species corridor rasters to a stack
  speciesCorridorRasterStack <- addLayer(speciesCorridorRasterStack, btslCorridorRaster)
}

# Sum five mosaiced species corridor rasters
allCorridorRaster <- calc(speciesCorridorRasterStack, sum)

# Save summed corridor raster to disk
writeRaster(allCorridorRaster, 
            file.path(b03linkagePriorityDir, "ALL_max_linkage_priority_30m_BTSL.tif"),
            overwrite = TRUE)
