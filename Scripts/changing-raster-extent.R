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

# Create a new raster named newRaster
newRaster <- raster(ncol=3333, nrow=3175, xmn=-691380, xmx=-391410, ymn=117930, ymx=403680)

# Check your raster metadata
newRaster

# Give each of your cells a value
values(newRaster) <- rep(1, ncell(newRaster))

# Check your raster metadata to see if it worked
newRaster

# Plot your raster to see if it works
x11(); plot(newRaster)

# Give your raster the same coordinate reference system as the original raster
crs(newRaster) <- crs(primaryStratum)

# Look at your raster metadata to see if it worked
newRaster

# Compare with the original raster
primaryStratum

# Export the resulting raster
writeRaster(newRaster, "C:/Enter/File/Path/newRaster.tif", format = "GTiff", overwrite=TRUE)