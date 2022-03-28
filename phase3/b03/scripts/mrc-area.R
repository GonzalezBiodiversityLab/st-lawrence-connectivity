# a254
# Sarah Chisholm, ApexRMS

# The script calculates the area of MRCs pre- and post-cropping down to the final run extent

# Load constants, functions, etc
source("./b03/scripts/0-constants.R")

# Workspace ----
# Load spatial data
secondaryStratumRaster <- raster(file.path(b03ProcessedMapsDir, "SecondaryStratum_90m.tif"))
secondaryStratumCroppedRaster <- raster(file.path(b03ProcessedMapsDir, "SecondaryStratum_90m_finalExtent.tif"))

# Parameters
cellSize <- 90

# Get cell counts of unique values
preCropValueCounts <- freq(secondaryStratumRaster)%>% 
                      as_tibble %>% 
                      rename("Pre-Crop Cell Count" = count)

postCropValueCounts <- freq(secondaryStratumCroppedRaster) %>% 
                       as_tibble %>% 
                       rename("Post-Crop Cell Count" = count)
                       

mrcArea <- preCropValueCounts %>% 
           left_join(postCropValueCounts, by = "value") %>% 
           mutate("Pre-Crop MRC Area (m^2)" = `Pre-Crop Cell Count` * cellSize^2,
                  "Post-Crop MRC Area (m^2)" = `Post-Crop Cell Count` * cellSize^2) %>% 
           rename(Value = value)
