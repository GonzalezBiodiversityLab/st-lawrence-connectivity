# a254
# Script by Sarah Chisholm

# Load constants, functions, etc
source("./b03/scripts/0-constants.R")

# Workspace ----

# Load spatial data

# Primary stratum
primaryStratumRaster <- raster(file.path(b03ProcessedMapsDir, "PrimaryStratum_90m.tif"))

# b03 extent


# Protected Areas
protectedAreasRaster <- raster(file.path(b03RawMapsDir, "AP-b03.tif"))
rmnShapefile <- st_read(dsn = b03RawMapsDir, layer = "RMN_20210608")

# Rasterize natural area ----
rmnRaster <- fasterize(sf = rmnShapefile,
                       raster = primaryStratumRaster,
                       background = 0)

# Add conservation priotity layers to protected areas layer 
allPAs <- protectedAreasRaster + rmnRaster
allPAs[allPAs == 0] <- 5
allPAs[allPAs != 5] <- 0  
allPAs[allPAs != 0] <- 1  

# Clip to b03 extent


# Save tif to disk
writeRaster(allPAs, file.path(b03ProcessedMapsDir, "SpatialMultiplier_ProtectedAreas_b03BufferExtent.tif"))
