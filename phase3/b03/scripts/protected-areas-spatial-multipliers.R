# a254
# Sarah Chisholm, ApexRMS
#
# This script generates three protected areas spatial multiplier layers: 
# - protected areas + areas of conservation priority
# - sites multicible
# - protected areas + areas of conservation priority + sites multicible

# Load constants, functions, etc
source("./b03/scripts/0-constants.R")

# Workspace ----

# Load spatial data

# Primary stratum
primaryStratumRaster <- raster(file.path(b03ProcessedMapsDir, "PrimaryStratum_90m.tif"))
# b03 extent
b03Extent <- raster(file.path(b03ProcessedMapsDir, "SecondaryStratum_90m_b03Extent.tif"))
# Protected Areas
protectedAreasRaster <- raster(file.path(b03RawMapsDir, "AP-b03.tif"))
# RMN protected areas
rmnShapefile <- st_read(dsn = b03RawMapsDir, layer = "RMN_20210608")
# Atlas priority areas
sitesMulticibleRawRaster <- raster(file.path(b03RawMapsDir, "sites-multicible.tif"))

# Protected areas + RMN protected areas spatial multiplier ----
# Rasterize rmnShapefile
rmnRaster <- fasterize(sf = rmnShapefile,
                       raster = primaryStratumRaster,
                       background = 0)

# Add conservation priotity layers to protected areas layer 
protectedAreasRmnRawRaster <- protectedAreasRaster + rmnRaster
protectedAreasRmnRawRaster[protectedAreasRmnRawRaster == 0] <- -1
protectedAreasRmnRawRaster[protectedAreasRmnRawRaster != -1] <- 0  
protectedAreasRmnRawRaster[protectedAreasRmnRawRaster != 0] <- 1  

# Mask to primaryStratumRaster
protectedAreasRmnRaster <- mask(x = protectedAreasRmnRawRaster,
                                mask = primaryStratumRaster)

# Crop to b03 extent
protectedAreasRmnCroppedRaster <- crop(x = protectedAreasRmnRaster, 
                                       y = b03Extent)

# Atlas priority areas spatial multiplier ----
sitesMulticibleRawRaster[sitesMulticibleRawRaster == 1] <- 0
sitesMulticibleRawRaster[is.na(sitesMulticibleRawRaster)] <- 1

# Mask to primaryStratumRaster
sitesMulticibleRaster <- mask(x = sitesMulticibleRawRaster, 
                              mask = primaryStratumRaster)

# Crop to b03 extent
sitesMulticibleCroppedRaster <- crop(x = sitesMulticibleRaster, 
                                     y = b03Extent)

# Protected areas + RMN protected areas + atlas priority areas spatial multiplier ----
allAreasRaster <- protectedAreasRmnRaster + sitesMulticibleRaster

# Reclassify allAreasRaster
fromValues <- c(0, 1, 2)
toValues <- c(0, 0, 1)

reclassifyMatrix <- cbind(fromValues, toValues)

allAreasReclassifiedRaster <- reclassify(x = allAreasRaster,
                                         rcl = reclassifyMatrix) %>% 
                              mask(mask = primaryStratumRaster)

# Crop to b03 extent
allAreasReclassifiedCroppedRaster <- crop(x = allAreasReclassifiedRaster, 
                                          y = b03Extent)

# Save tifs to disk ----
writeRaster(x = protectedAreasRmnRaster, 
            filename = file.path(b03ProcessedMapsDir, "SpatialMultiplier_ProtectedAreas_b03BufferExtent.tif"),
            overwrite = TRUE)
writeRaster(x = sitesMulticibleRaster, 
            filename = file.path(b03ProcessedMapsDir, "SpatialMultiplier_AtlasPriorityAreas_b03BufferExtent.tif"),
            overwrite = TRUE)
writeRaster(x = allAreasReclassifiedRaster, 
            filename = file.path(b03ProcessedMapsDir, "SpatialMultiplier_ProtectedAreas_AtlasPriorityAreas_b03BufferExtent.tif"),
            overwrite = TRUE)
writeRaster(x = protectedAreasRmnCroppedRaster, 
            filename = file.path(b03ProcessedMapsDir, "SpatialMultiplier_ProtectedAreas_b03Extent.tif"),
            overwrite = TRUE)
writeRaster(x = sitesMulticibleCroppedRaster, 
            filename = file.path(b03ProcessedMapsDir, "SpatialMultiplier_AtlasPriorityAreas_b03Extent.tif"),
            overwrite = TRUE)
writeRaster(x = allAreasReclassifiedCroppedRaster, 
            filename = file.path(b03ProcessedMapsDir, "SpatialMultiplier_AllProtectedAreas_AtlasPriorityAreas_b03Extent.tif"),
            overwrite = TRUE)
