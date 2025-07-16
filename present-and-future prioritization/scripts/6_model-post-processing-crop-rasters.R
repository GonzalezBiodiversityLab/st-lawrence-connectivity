# a254
# Sarah Chisholm, ApexRMS
# Run with R-4.1.2
#
# This script masks model outputs to the BTSL study area and to the BTSL area 
# with a buffer area that extends north of the BTSL. 


# Load constants, functions, etc
source("./b03/scripts/0-constants.R")

# Workspace ----

# Load spatial data
# Primary stratum to be used as clipping mask
primaryStratumRaster <- raster(file.path(b03ProcessedMapsDir, "PrimaryStratum_90m_b03Extent.tif"))

# List of files to crop
files <- list.files(path = b03ResultsMapsDir, 
                    pattern = "tif$", 
                    full.names = TRUE, 
                    recursive = TRUE)

# Create mask rasters ----
# Define reclassification matrices
btslWithBufferReclassMatrix <- matrix(data = c(2, 3, 5, 6, NA, 1, 1, NA), ncol = 2, nrow = 4) 
btslReclassMatrix <- matrix(data = c(2, 3, 5, 6, NA, 1, NA, NA), ncol = 2, nrow = 4) 

# Generate mask rasters of the BTSL with a buffer, and the BTSL only
btslWithBufferMask <- reclassify(x = primaryStratumRaster, 
                               rcl = btslWithBufferReclassMatrix)

btslMask <- reclassify(x = primaryStratumRaster, 
                     rcl = btslReclassMatrix)

# Crop layers to masks ----
# This takes ~ 8 minutes

start_time <- Sys.time() # start the clock
for(file in files){
  
  rawRaster <- raster(file)
  
  # Crop to the BTSL with a buffer
  maskedBtslWithBufferRaster <- mask(rawRaster, btslWithBufferMask)
  
  # Crop to the BTSL only
  maskedBtslRaster <- mask(rawRaster, btslMask)
  
  # Write rasters to disk
  writeRaster(x = maskedBtslWithBufferRaster,
              filename = file,
              format = "GTiff",
              overwrite = TRUE)
  
  writeRaster(x = maskedBtslRaster,
              filename = str_replace(file, ".tif", "_BTSL.tif"),
              format = "GTiff",
              overwrite = TRUE)
}
end_time <- Sys.time() # stop the clock
end_time - start_time # total time to mask layers
