# a254
# Bronwyn Rayfield and Jed Lloren, ApexRMS
# Run with R-4.1.1

# Load constants, functions, etc
source("./b03/scripts/0-0-constants.R")

# Workspace ----

# Read in data
# Spatial
landcoverBTSLRaw <- st_read(dsn = file.path(rawMapsDir, "Extrait_Donnees.gdb"), layer="BTSL_SLL_Occ_sol_Land_cover")
siefDataRaw <- st_read(dsn = file.path(rawMapsDir, "Extrait_Donnees.gdb"), layer="SIEF_C08PEEFO")
studyAreaRaw <- st_read(dsn = rawMapsDir, layer = "CR_NIV_01_S")
ontarioBoundary <- st_read(dsn = file.path(rawMapsDir, "ON Provincial Border"), layer = "Province")
outaouaisExtent <- raster(file.path(rawMapsDir, "OutaouaisConnectivityExtent-30m.tif"))

# Tabular
landcoverBTSLReclass <- read_csv(file.path(rawTablesDir, "landcoverBTSLReclass.csv"))
forestAgeReclass <- read_csv(file.path(rawTablesDir, "forestAgeReclass.csv"))
depositReclass <- read_csv(file.path(rawTablesDir, "depositReclass.csv"))
forestDensityReclass <- read_csv(file.path(rawTablesDir, "speciesForestDensityReclass.csv"))

# Reformatting data ----
# Forest age
# Merge the layer and reclass tables
landcoverMerge <- landcoverBTSLRaw %>%
  left_join(landcoverBTSLReclass, by = c("CLASSE_DET" = "CLASSE_DET")) # Default joins by CLASSE_GEN and CLASSE_DET
# Rasterize, crop, and mask to study area
landcoverBTSL <- fasterize(sf = st_cast(landcoverMerge, "MULTIPOLYGON"), 
                       raster = outaouaisExtent, # Snap to Outaouais connectivity boundaries
                       field = "Value") %>% 
  mask(., mask=outaouaisExtent)

# Forest age
forestAgeMerge <- siefDataRaw %>%
  left_join(forestAgeReclass)

forestAge <- fasterize(sf = st_cast(forestAgeMerge, "MULTIPOLYGON"), 
                       raster = outaouaisExtent,
                       field = "Value") %>% 
  mask(., mask=outaouaisExtent)
    
# Surficial deposits
depositMerge <- siefDataRaw %>%
  left_join(depositReclass, by = c("DEP_SUR" = "CODE"))

surficialDeposits <- fasterize(sf = st_cast(depositMerge, "MULTIPOLYGON"), 
                       raster = outaouaisExtent,
                       field = "Recode") %>% 
  mask(., mask=outaouaisExtent)

# Forest density
forestDensityMerge <- siefDataRaw %>%
  left_join(forestDensityReclass, by = c("CL_DENS" = "DensityName"))

forestDensity <- fasterize(sf = st_cast(forestDensityMerge, "MULTIPOLYGON"), 
                               raster = outaouaisExtent,
                               field = "DensityCode") %>% 
  mask(., mask=outaouaisExtent)

# Study area
studyAreaRaw <- studyAreaRaw %>%
  mutate(ID_Code = c(3, 4, 1, 2))
studyArea <- fasterize(sf = st_cast(studyAreaRaw, "MULTIPOLYGON"), 
                               raster = outaouaisExtent,
                               field = "ID_Code") %>% 
  mask(., mask=outaouaisExtent)

# Ontario boundary
ontarioBoundary <- st_transform(ontarioBoundary, crs = crs(outaouaisExtent))
ontarioBoundaryRaster <- fasterize(sf = st_cast(ontarioBoundary, "MULTIPOLYGON"),
                                   raster = outaouaisExtent,
                                   field = "OBJECTID") %>% 
  mask(., mask=outaouaisExtent)

# Save outputs ----
# Landcover
writeRaster(landcoverBTSL, 
            file.path(rawMapsDir, "BTSL-SLL-Occ-sol-Land-cover.tif"), 
            overwrite = TRUE)
# Forest age
writeRaster(forestAge, 
            file.path(rawMapsDir, "SIEF-C08PEEFO-forest-age.tif"), 
            overwrite = TRUE)
# Surficial deposits
writeRaster(surficialDeposits, 
            file.path(rawMapsDir, "SIEF-C08PEEFO-surficial-deposits.tif"), 
            overwrite = TRUE)
# Forest density
writeRaster(forestDensity, 
            file.path(rawMapsDir, "SIEF-C08PEEFO-forest-density.tif"), 
            overwrite = TRUE)
# Study area
writeRaster(studyArea, 
            file.path(rawMapsDir, "b03-studyarea-zones-30m.tif"), 
            overwrite = TRUE)
# Ontario mask
writeRaster(ontarioBoundaryRaster, 
            file.path(rawMapsDir, "ontario-mask.tif"), 
            overwrite = TRUE)
