#####################################################################################################
# A254 Outaouais connectivity project                              
# Creating a new study extent to capture the western portion of the St. Lawrence Lowlands ecoregion                          
# 08/2021  
#                                                                                              
# Note: When loading the resulting raster layer use NAD83(CSRS) to WGS 84 (1) if prompted
#
#
# Script created by Jed Lloren with input from Bronwyn Rayfield for ApexRMS                                  
#####################################################################################################

# Packages
library(raster)
library(rgdal)
library(sp)

# Create a new raster named OutaouaisConnectivityExtent
# ncol and nrow were calulated using (xmx-xmn)/90 and (ymx-ymn)/90, respectively
# OutaouaisConnectivityExtent90m <- raster(ncol=3333, nrow=3175, xmn=-691380, xmx=-391410, ymn=117930, ymx=403680)
OutaouaisConnectivityExtent30m <- raster(ncol=9999, nrow=9525, xmn=-691380, xmx=-391410, ymn=117930, ymx=403680)

# Give each of your cells a value
# We'll just be giving the entire raster a value of 1
# values(OutaouaisConnectivityExtent90m) <- rep(1, ncell(OutaouaisConnectivityExtent90m))
values(OutaouaisConnectivityExtent30m) <- rep(1, ncell(OutaouaisConnectivityExtent30m))

# Give your raster the same coordinate reference system as the original raster
# crs(OutaouaisConnectivityExtent90m) <- "+proj=lcc +lat_0=44 +lon_0=-68.5 +lat_1=46 +lat_2=60 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"
crs(OutaouaisConnectivityExtent30m) <- "+proj=lcc +lat_0=44 +lon_0=-68.5 +lat_1=46 +lat_2=60 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"

# Export the resulting raster
# writeRaster(OutaouaisConnectivityExtent90m, "C:/Enter/File/Path/OutaouaisConnectivityExtent.tif", format = "GTiff", overwrite=TRUE)
writeRaster(OutaouaisConnectivityExtent30m, "C:/Enter/File/Path/OutaouaisConnectivityExtent-30m.tif", format = "GTiff", overwrite=TRUE)