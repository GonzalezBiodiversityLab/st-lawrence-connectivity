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


# Protected Areas
protectedAreasRaster <- raster(file.path(b03RawMapsDir, "AP-b03.tif"))

# Areas of conservation priority
rmnShapefile <- st_read(dsn = b03RawMapsDir, layer = "RMN_20210608")

# Sites multicible
sitesMulticibleRawRaster <- raster(file.path(b03RawMapsDir, "sites-multicible.tif"))

# Protected areas + areas of conservation priority spatial multiplier ----
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

# Sites multicible spatial multiplier ----
sitesMulticibleRawRaster[sitesMulticibleRawRaster == 1] <- 0
sitesMulticibleRawRaster[is.na(sitesMulticibleRawRaster)] <- 1

# Mask to primaryStratumRaster
sitesMulticibleRaster <- mask(x = sitesMulticibleRawRaster, 
                              mask = primaryStratumRaster)

# Protected areas + areas of conservation priority + sites multicible spatial multiplier ----
allAreasRaster <- protectedAreasRmnRaster + sitesMulticibleRaster

# Reclassify allAreasRaster
fromValues <- c(0, 1, 2)
toValues <- c(0, 0, 1)

reclassifyMatrix <- cbind(fromValues, toValues)

allAreasReclassifiedRaster <- reclassify(x = allAreasRaster,
                                         rcl = reclassifyMatrix) %>% 
                              mask(mask = primaryStratumRaster)

# Clip layers to b03 extent ----


# Save tifs to disk ----
writeRaster(x = protectedAreasRmnRaster, 
            filename = file.path(b03ProcessedMapsDir, "SpatialMultiplier_ProtectedAreas_RMN_b03BufferExtent.tif"),
            overwrite = TRUE)
writeRaster(x = sitesMulticibleRaster, 
            filename = file.path(b03ProcessedMapsDir, "SpatialMultiplier_SitesMulticible_b03BufferExtent.tif"),
            overwrite = TRUE)
writeRaster(x = allAreasReclassifiedRaster, 
            filename = file.path(b03ProcessedMapsDir, "SpatialMultiplier_AllProtectedAreas_b03BufferExtent.tif"),
            overwrite = TRUE)
