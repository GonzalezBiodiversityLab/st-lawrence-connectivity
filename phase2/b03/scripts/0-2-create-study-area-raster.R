# a254
# Bronwyn Rayfield and Jed Lloren, ApexRMS
# Run with R-1.1

# Load constants, functions, etc
source("./b03/scripts/0-0-constants.R")

# Workspace ----

# Read in data
# Spatial
plaineDottawa <- st_read(file.path(rawMapsDir,"plaine-dottawa.shp"))
studyExtent <- raster(file.path(rawMapsDir,"OutaouaisConnectivityExtent-30m.tif"))

# Tabular
reclassTable <- read_csv(file.path(rawTablesDir,"b03-study-area-reclass.csv"))

# Reformatting data ----
# Rasterize the Plaine d'Ottawa shapefile
b03Studyarea <- fasterize(sf = st_cast(plaineDottawa, "POLYGON"), 
                       raster = studyExtent,
                       field = "FID02") %>% 
  mask(., mask=studyExtent)

# Reclassify
b03Studyarea <- reclassify(b03Studyarea, reclassTable)

# Save outputs ----
# B03 study area
writeRaster(b03Studyarea, 
            file.path(rawMapsDir, "b03-studyarea-30m.tif"), 
            overwrite = TRUE)