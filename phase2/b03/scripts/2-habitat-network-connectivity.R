# a254
# Bronwyn Rayfield, ApexRMS
# Create betweenness and habitat netowrk maps for focal species
#                                                                   
# Inputs:                                                           
#    - Habitat suitability map                                      
#    - Resistance map                                               
#    - Dispersal parameters                                         
# Outputs:                                                          
#    - Habitat network map                                          
#    - Betweenness map  

                
# Load constants, functions, etc
source("./b03/scripts/0-0-constants.R")

# Workspace ----
# Packages
library(igraph)
library(grainscape)

# Settings
options(stringsAsFactors=FALSE, SHAPE_RESTORE_SHX=T, useFancyQuotes = F, digits=10)

# Generate short-range betweenness and link layers ----
for(species in speciesList){
  
  # Read in data
  habitatRaster <- raster(file.path(b03habitatDir, paste0(species, "_habitatPatch_Focal_Filtered_", myResolution, "m.tif")))
  resistanceRaster <- raster(file.path(b03resistanceDir, paste0(species, "_resistance_Focal_", myResolution, "m.tif")))
  
  # Extract mpg ----
  mpg <- MPG(cost = resistanceRaster, patch = habitatRaster)
  
  # Betweenness ----
  btwn <- data.frame(patchId = as.numeric(vertex_attr(mpg$mpg, 'patchId')), btwn = betweenness(mpg$mpg, weights = edge_attr(mpg$mpg, 'lcpPerimWeight'), directed = FALSE))
  #make a look-up table between patch id and betweenness value
  btwn_lookup <- cbind(patchId = btwn$patchId, btwn = 1 / (max(btwn$btwn) - min(btwn$btwn)) * (btwn$btwn - min(btwn$btwn)))
  #replace patch ids with betweeness values
  shortBtwnRaster <- reclassify(mpg$patchId, btwn_lookup)
  
  # Calculate overall EC ----
  dist_mat <- shortest.paths(mpg$mpg, weights = edge_attr(mpg$mpg, 'lcpPerimWeight'))
  
  #Natal
  #matrix of dispersal kernel
  kernel_mat <- exp(-dist_mat / 459)
  rm(dist_mat)
  
  #matrix of product node attributes
  attribute_mat1 <- as.vector(as.numeric(vertex_attr(mpg$mpg, 'patchArea')))
  attribute_mat <- attribute_mat1 %*% t(attribute_mat1)
  
  #matrix of PCvalues
  PC_mat <- kernel_mat * attribute_mat
  rm(attribute_mat)
  
  PCnum <- sum(PC_mat)
  ECNatal <- sqrt(PCnum)
  
  # Convert link rasters to shapefiles ----
  linkPolygon <- rasterToPolygons(mpg$lcpPerimWeight) %>% st_as_sf
  
  # Make binary link layers
  lcpPerimWeightBinary <- mpg$lcpPerimWeight
  lcpPerimWeightBinary[lcpPerimWeightBinary >= 0] <- 1
  lcpPerimWeightBinary[is.na(lcpPerimWeightBinary)] <- 0
  
  # Save outputs ----
  #geotif
  writeRaster(shortBtwnRaster, file.path(b03networkDir, paste0(species, "_ShortBetweenness_", myResolution, "m.tif")), overwrite=TRUE)
  writeRaster(mpg$lcpPerimWeight, file.path(b03networkDir, paste0(species, "_LinkLength_", myResolution, "m_b03.tif")), overwrite=TRUE)
  writeRaster(lcpPerimWeightBinary, file.path(b03networkDir, paste0(species, "_LinkLength_01_", myResolution, "m_b03.tif")), overwrite=TRUE)
  #shapefile
  st_write(linkPolygon, dsn = b03networkDir, layer = paste0(species, "_LinkLength_", myResolution, "_b03"), driver = "ESRI Shapefile", append = FALSE)
}
  
# Generate long-range betweenness layers ----
for(species in speciesList){

  # Read in data
  habitatRaster <- raster(file.path(b03habitatDir, paste0(species, "_habitatPatch_", myResolution, "m.tif")))
  resistanceRaster <- raster(file.path(b03resistanceDir, paste0(species, "_resistance_", myResolution, "m.tif")))
  
  # Extract mpg ----
  mpg <- MPG(cost = resistanceRaster, patch = habitatRaster)
  
  # Betweenness ----
  btwn <- data.frame(patchId = as.numeric(vertex_attr(mpg$mpg, 'patchId')), btwn = betweenness(mpg$mpg, weights = edge_attr(mpg$mpg, 'lcpPerimWeight'), directed = FALSE))
  #make a look-up table between patch id and betweenness value
  btwn_lookup <- cbind(patchId = btwn$patchId, btwn = 1 / (max(btwn$btwn) - min(btwn$btwn)) * (btwn$btwn - min(btwn$btwn)))
  #replace patch ids with betweeness values
  btwnRaster <- reclassify(mpg$patchId, btwn_lookup)
  
  # Mask long-range betweenness rasters to b03 filtered patches
  filteredPatchRaster <- raster(file.path(b03habitatDir, paste0(species, "_habitatPatch_Focal_Filtered_", myResolution, "m.tif")))
  filteredPatchRaster[filteredPatchRaster == 0] <- NA
    
  maskedBtwnRaster <- btwnRaster * filteredPatchRaster
  
  # Calculate overall EC ----
  dist_mat <- shortest.paths(mpg$mpg, weights = edge_attr(mpg$mpg, 'lcpPerimWeight'))
  
  #Natal
  #matrix of dispersal kernel
  kernel_mat <- exp(-dist_mat / 459)
  rm(dist_mat)
  
  #matrix of product node attributes
  attribute_mat1 <- as.vector(as.numeric(vertex_attr(mpg$mpg, 'patchArea')))
  attribute_mat <- attribute_mat1 %*% t(attribute_mat1)
  
  #matrix of PCvalues
  PC_mat <- kernel_mat * attribute_mat
  rm(attribute_mat)
  
  PCnum <- sum(PC_mat)
  ECNatal <- sqrt(PCnum)
  
  # Convert link rasters to shapefiles ----
  linkPolygon <- rasterToPolygons(mpg$lcpPerimWeight) %>% st_as_sf
  
  # Make binary link layers
  lcpPerimWeightBinary <- mpg$lcpPerimWeight
  lcpPerimWeightBinary[lcpPerimWeightBinary >= 0] <- 1
  lcpPerimWeightBinary[is.na(lcpPerimWeightBinary)] <- 0

  # Save outputs ----
  #geotif
  writeRaster(maskedBtwnRaster, file.path(b03networkDir, paste0(species, "_LongBetweenness_", myResolution, "m.tif")), overwrite=TRUE)
  writeRaster(mpg$lcpPerimWeight, file.path(b03networkDir, paste0(species, "_LinkLength_", myResolution, "m.tif")), overwrite=TRUE)
  writeRaster(lcpPerimWeightBinary, file.path(b03networkDir, paste0(species, "_LinkLength_01_", myResolution, "m.tif")), overwrite=TRUE)
  #shapefile
  st_write(linkPolygon, dsn = b03networkDir, layer = paste0(species, "_LinkLength_", myResolution), driver = "ESRI Shapefile", append = FALSE)
}
