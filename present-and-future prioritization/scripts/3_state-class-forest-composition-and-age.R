# a254
# Sarah Chisholm, ApexRMS
# Run with R-4.1.2
#
# This script imposes forest age and forest composition onto forested state class cells

# Load constants, functions, etc
source("./b03/scripts/0-constants.R")

# Workspace ----
# Load data
# Spatial

# State class raster
stateClassRawRaster <- raster(file.path(b03ProcessedMapsDir, "StateClass2010_90m_coarse_forest.tif"))

# Foret composition raster
forestCompositionRawRaster <- raster(file.path(b03RawMapsDir, "baseline_output_cover_types_0.tif"))
crs(forestCompositionRawRaster) <- crs(stateClassRawRaster)

# Forest age raster
forestAgeRawRaster <- raster(file.path(b03RawMapsDir, "AGE-AVG-0.tif"))
crs(forestAgeRawRaster) <- crs(stateClassRawRaster)

# Create forest mask ---
forestMaskRaster <- stateClassRawRaster
forestMaskRaster[forestMaskRaster == 500] <- 1
forestMaskRaster[forestMaskRaster != 1] <- NA

# Create forest composition mask (0 <- NA)
forestCompositionMaskRaster <- forestCompositionRawRaster
forestCompositionMaskRaster[forestCompositionMaskRaster == 0] <- NA
forestCompositionMaskRaster[!is.na(forestCompositionMaskRaster)] <- 1

# Reclassify forest composition raster ----
# Create reclassification matrix
fromValues <- c(0, 1, 2, 3, 4)
toValues <- c(NA, 1, 2, 2, 3)

reclassMatrix <- cbind(fromValues, toValues)

forestCompositionRaster <- reclassify(x = forestCompositionRawRaster,
                                      rcl = reclassMatrix) 

forestCompositionRaster <- forestCompositionRaster * forestMaskRaster

# Reclassify forest age raster ----
# 1 == 0-30
# 2 == 31-60
# 3 == 61+

# Mask forest age to forest composition mask raster
forestAgeMaskedRaster <- forestAgeRawRaster %>% 
  mask(mask = forestCompositionMaskRaster)

# Create reclassification matrix
fromValues <- c(-1, 30)
toValues <- c(30, 60)
becomesValues <- c(1, 2)

reclassMatrix <- cbind(fromValues, toValues, becomesValues)

forestAgeRaster <- reclassify(x = forestAgeMaskedRaster,
                              rcl = reclassMatrix)

forestAgeRaster[forestAgeRaster >= 61] <- 3

forestAgeRaster <- forestAgeRaster * forestMaskRaster

# Impose forest composition and age onto stateClass raster ----
stateClassRaster <- stateClassRawRaster + (forestCompositionRaster*10) + forestAgeRaster
stateClassRaster[is.na(stateClassRaster)] <- 0

stateClassRaster <- stateClassRawRaster + stateClassRaster

fromValues <- c(100, 400, 410, 500, 700, 800, 810, 1011, 1012, 1013, 1021, 1022, 1023, 1031, 1032, 1033)
toValues <- c(100, 400, 410, 522, 700, 800, 810, 511, 512, 513, 521, 522, 523, 531, 532, 533)
reclassMatrix <- cbind(fromValues, toValues)

stateClassCompositionAgeRaster <- reclassify(x = stateClassRaster,
                                             rcl = reclassMatrix)
# Write tif to disk
writeRaster(stateClassCompositionAgeRaster, 
            file.path(b03ProcessedMapsDir, "StateClass_2010_90m.tif"),
            overwrite = TRUE)


