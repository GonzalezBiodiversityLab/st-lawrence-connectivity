# a254
# Sarah Chisholm, ApexRMS
#
# Crop strata, state class, and protected area spatial multiplier rasters to b03 extent

# Load constants, functions, etc
source("./b03/scripts/0-constants.R")

# Workspace ----
# Load b03 extent raster
b03Extent <- raster(file.path(b03ProcessedMapsDir, "SecondaryStratum_90m_finalExtent.tif"))

# List rasters to crop
strataList <- list(file.path(b03ProcessedMapsDir, "PrimaryStratum_90m.tif"),
                   file.path(b03ProcessedMapsDir, "TertiaryStratum_90m.tif"),
                   file.path(b03ProcessedMapsDir, "StateClass_2010_90m.tif"),
                   file.path(b03ProcessedMapsDir, "SpatialMultiplier_AtlasPriorityAreas_90m.tif"),
                   file.path(b03ProcessedMapsDir, "SpatialMultiplier_ProtectedAreas_AtlasPriorityAreas_90m.tif"),
                   file.path(b03ProcessedMapsDir, "SpatialMultiplier_ProtectedAreas_90m.tif"))

# Crop strata to the final run extent
for(stratum in strataList){
  
  stratumRaster <- raster(stratum) 
  crs(stratumRaster) <- crs(b03Extent)
  
  stratumCroppedRaster <- crop(x = stratumRaster, y = b03Extent) %>% 
                   mask(mask = b03Extent)
  

  writeRaster(stratumCroppedRaster,
              file.path(b03ProcessedMapsDir, str_c(names(stratumRaster), "_finalExtent.tif")),
              overwrite = TRUE)
}