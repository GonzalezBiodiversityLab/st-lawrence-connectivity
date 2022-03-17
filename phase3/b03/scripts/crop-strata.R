# a254
# Sarah Chisholm, ApexRMS
#
# Crop strata and state class raster to b03 extent

# Load constants, functions, etc
source("./b03/scripts/0-constants.R")

# Workspace ----
# Load b03 extent raster
b03Extent <- raster(file.path(b03ProcessedMapsDir, "SecondaryStratum_90m_b03Extent.tif"))

# List rasters to crop
strataList <- list(file.path(b03ProcessedMapsDir, "PrimaryStratum_90m.tif"),
                   file.path(b03ProcessedMapsDir, "TertiaryStratum_90m.tif"),
                   file.path(b03ProcessedMapsDir, "StateClass_2010_90m.tif"))

# Crop strata to focal region
for(stratum in strataList){
  
  stratumRaster <- raster(stratum) %>% 
                   crop(b03Extent)

  writeRaster(stratumRaster,
              file.path(b03ProcessedMapsDir, str_c(names(r), "_b03Extent.tif")))
}

