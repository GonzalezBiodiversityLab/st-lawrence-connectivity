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

# Create an object with the same as your file
primaryStratum<-raster("primaryStratum.tif")

# Look at your raster metadata
# I will be referring to this raster as "original raster" throughout
primaryStratum

# class      : RasterLayer 
# dimensions : 3175, 4399, 13966825  (nrow, ncol, ncell)
# resolution : 90, 90  (x, y)
# extent     : -491430, -95520, 117930, 403680  (xmin, xmax, ymin, ymax)
# crs        : +proj=lcc +lat_0=44 +lon_0=-68.5 +lat_1=46 +lat_2=60 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs 
# source     : PrimaryStratum.tif 
# names      : PrimaryStratum 
# values     : -32768, 32767  (min, max)

# Create a new raster named OutaouaisConnectivityExtent
# ymx and ymn values were retained from the primaryStratum raster
# xmx (eastern extent) was moved 295.89km west
# xmn (western extent) was moved 199.98km west
# ncol and nrow were calulated using (xmx-xmn)/90 and (ymx-ymn)/90, respectively
# This new raster has an area of 85,716.4275 km^2
OutaouaisConnectivityExtent <- raster(ncol=3333, nrow=3175, xmn=-691380, xmx=-391410, ymn=117930, ymx=403680)

# Check your raster metadata
OutaouaisConnectivityExtent

# class      : RasterLayer 
# dimensions : 3175, 3333, 10582275  (nrow, ncol, ncell)
# resolution : 90, 90  (x, y)
# extent     : -691380, -391410, 117930, 403680  (xmin, xmax, ymin, ymax)
# crs        : NA 

# Give each of your cells a value
# We'll just be giving the entire raster a value of 1
values(OutaouaisConnectivityExtent) <- rep(1, ncell(OutaouaisConnectivityExtent))

# Check your raster metadata to see if it worked
OutaouaisConnectivityExtent

# class      : RasterLayer 
# dimensions : 3175, 3333, 10582275  (nrow, ncol, ncell)
# resolution : 90, 90  (x, y)
# extent     : -691380, -391410, 117930, 403680  (xmin, xmax, ymin, ymax)
# crs        : NA 
# source     : memory
# names      : layer 
# values     : 1, 1  (min, max)

# Plot your raster to see if it works
x11(); plot(OutaouaisConnectivityExtent)

# Give your raster the same coordinate reference system as the original raster
crs(OutaouaisConnectivityExtent) <- crs(primaryStratum)

# Look at your raster metadata to see if it worked
OutaouaisConnectivityExtent

# class      : RasterLayer 
# dimensions : 3175, 3333, 10582275  (nrow, ncol, ncell)
# resolution : 90, 90  (x, y)
# extent     : -691380, -391410, 117930, 403680  (xmin, xmax, ymin, ymax)
# crs        : +proj=lcc +lat_0=44 +lon_0=-68.5 +lat_1=46 +lat_2=60 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs 
# source     : memory
# names      : layer 
# values     : 1, 1  (min, max)

# Compare with the original raster
primaryStratum

# class      : RasterLayer 
# dimensions : 3175, 4399, 13966825  (nrow, ncol, ncell)
# resolution : 90, 90  (x, y)
# extent     : -491430, -95520, 117930, 403680  (xmin, xmax, ymin, ymax)
# crs        : +proj=lcc +lat_0=44 +lon_0=-68.5 +lat_1=46 +lat_2=60 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs 
# source     : PrimaryStratum.tif 
# names      : PrimaryStratum 
# values     : -32768, 32767  (min, max)

# Export the resulting raster
writeRaster(OutaouaisConnectivityExtent, "C:/Enter/File/Path/OutaouaisConnectivityExtent.tif", format = "GTiff", overwrite=TRUE)