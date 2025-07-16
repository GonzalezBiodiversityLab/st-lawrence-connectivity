# a254
# Sarah Chisholm, ApexRMS

# The script calculates the area of MRCs pre- and post-cropping down to the final run extent

# Load constants, functions, etc
source("./b03/scripts/0-constants.R")

# Workspace ----
# Load spatial data
# Secondary stratum layers
secondaryStratumRaster <- raster(file.path(b03ProcessedMapsDir, "SecondaryStratum_90m.tif"))
secondaryStratumCroppedRaster <- raster(file.path(b03ProcessedMapsDir, "SecondaryStratum_90m_b03Extent.tif"))

# Primary stratum
primaryStratumRaster <- raster(file.path(b03ProcessedMapsDir, "PrimaryStratum_90m.tif"))

# Parameters
cellSize <- 90

# Mask b03 study area from secondary stratum rasters ---
# Set primary stratum crs to match crs of secondary stratum
crs(primaryStratumRaster) <- crs(secondaryStratumRaster)

# Set b03 region (value = 3) to NA
primaryStratumRaster[primaryStratumRaster == 3] <- NA

# Crop masked primary stratum to b03 extent
primaryStratumCroppedRaster <- crop(x = primaryStratumRaster, y = secondaryStratumCroppedRaster)

# Mask secindary stratum layers
secondaryStratumMaskedRaster <- mask(x = secondaryStratumRaster, 
                                     mask = primaryStratumRaster)
secondaryStratumCroppedMaskedRaster <- mask(x = secondaryStratumCroppedRaster, 
                                            mask = primaryStratumCroppedRaster)

# Get cell counts of unique values ----
bufferExtentCounts <- freq(secondaryStratumMaskedRaster)%>% 
                      as_tibble %>% 
                      rename(CellCount = count)

b03ExtentCounts <- freq(secondaryStratumCroppedMaskedRaster) %>% 
                       as_tibble %>% 
                       rename(CropCellCount = count)
                       

mrcArea <- bufferExtentCounts %>% 
           left_join(b03ExtentCounts, by = "value") %>% 
           replace_na(list(CropCellCount = 0)) %>% 
           mutate(Ratio = CropCellCount/CellCount,
                  BufferExtentArea = CellCount * cellSize^2,
                  b03ExtentArea = CropCellCount * cellSize^2) %>% 
           filter(Ratio != 1 & Ratio != 0 & !is.na(value)) %>% 
           select(value, BufferExtentArea, b03ExtentArea, Ratio) %>% 
           rename("MRC ID" = value,
                  "Buffer Extent Area (square meters)" = BufferExtentArea,
                  "b03 Extent Area (square meters)" = b03ExtentArea)

# Save tabular data to disk
write_csv(mrcArea, file.path(b03RawTablesDir, "mrc-area-ratios.csv"))
