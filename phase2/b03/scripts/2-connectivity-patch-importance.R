#####################################################################
# a254      
# Calculate patch importance metrics for focal species
# 03-2022                                       					
#                             
#	  Uses Makurhini R package (see https://connectscape.github.io/Makurhini/)          
# 		- Function: MK_dPCIIC()   
#                                             
#   Inputs (for focal species):
#    -habitat patches
#    -resistance layer
#    -dispersal distance 
#  
#   Outputs:
#    -patch importance and fractions
#                                                                   
# Script by B Rayfield for ApexRMS 									
#####################################################################

# Load constants, functions, etc
source("./b03/scripts/0-0-constants.R")

## Workspace ----

# Packages
library(Makurhini)

# Parameters
species <- "URAM"
dispersalMode <- "Natal" # "Gap"

# Load data
# Tabular
dispersalDistance <- read_csv(file.path(b01b02RawTablesDir, "speciesDispersalParameters.csv"))

## Calculate patch importance for focal species ----

# NB the  MK_dPCIIC() function does not work inside a loop so manually set each species for now
# Loop through species, generating Patch Importance data and saving outputs
# List of species
# speciesList <- dispersalDistance$Species

#for (i in speciesList[c(5)]){
  
  # species <- i
  
  ## Load data for focal species
  # Spatial data
  if(withinBTSL) {
    habitatPatches <- raster(file.path(b03habitatDir, paste0(species, "_habitatPatch_Focal_Filtered_", coarseResolution, "m.tif")))
    resistance <- raster(file.path(b03resistanceDir, paste0(species, "_resistance_Focal_", coarseResolution, "m.tif")))
  } else {
    habitatPatches <- raster(file.path(b03habitatDir, paste0(species, "_habitatPatch_", coarseResolution, "m.tif")))
    resistance <- raster(file.path(b03resistanceDir, paste0(species, "_resistance_", coarseResolution, "m.tif")))
  }
  
  # Species dispersal maximal/or median distance 
  maxdist <- dispersalDistance[[which(dispersalDistance$Species == species), dispersalMode]]
  # Albert distance estimates represent the Median
  prob <- 0.5
  
   # Create raster with unique patch id for each species
  habitatPatchesID <- clump(habitatPatches, directions=8)
  length(unique(habitatPatchesID))
  
  ## PC fractions measurements for focal species in focal Area ----
 
  # PC & fractions
  PCfocal <- MK_dPCIIC(nodes = habitatPatchesID, 
                       distance = list(type = "least-cost", 
                                       resistance = resistance),
                       attribute = NULL,
                       metric = "PC", 
                       probability = prob, 
                       distance_thresholds = maxdist)
  
  # PCfocal is a raster stack : plot individual rasters
  plot(PCfocal)
  
  ## Save single-band geotiffs ----
  if(withinBTSL){
    # dPC
    writeRaster(PCfocal[[2]], 
                file.path(b03patchImportanceDir, 
                          paste0(species, "_PC_", dispersalMode, "_Focal_", coarseResolution, "m.tif")), 
                overwrite=TRUE)  
    #dPCintra
    writeRaster(PCfocal[[3]], 
                file.path(b03patchImportanceDir, 
                          paste0(species, "_PC_Intra_", dispersalMode, "_Focal_", coarseResolution, "m.tif")), 
                overwrite=TRUE)  
    #dPCflux
    writeRaster(PCfocal[[4]], 
                file.path(b03patchImportanceDir, 
                          paste0(species, "_PC_Flux_", dispersalMode, "_Focal_", coarseResolution, "m.tif")), 
                overwrite=TRUE)  
    #dPCconnector
    writeRaster(PCfocal[[5]], 
                file.path(b03patchImportanceDir, 
                          paste0(species, "_PC_Connector_", dispersalMode, "_Focal_", coarseResolution, "m.tif")), 
                overwrite=TRUE)  		
  } else {
    # dPC
    writeRaster(PCfocal[[2]], 
                file.path(b03patchImportanceDir, 
                          paste0(species, "_PC_", dispersalMode, "_", coarseResolution, "m.tif")), 
                overwrite=TRUE)  
    #dPCintra
    writeRaster(PCfocal[[3]], 
                file.path(b03patchImportanceDir, 
                          paste0(species, "_PC_Intra_", dispersalMode, "_", coarseResolution, "m.tif")), 
                overwrite=TRUE)  
    #dPCflux
    writeRaster(PCfocal[[4]], 
                file.path(b03patchImportanceDir, 
                          paste0(species, "_PC_Flux_", dispersalMode, "_", coarseResolution, "m.tif")), 
                overwrite=TRUE)  
    #dPCconnector
    writeRaster(PCfocal[[5]], 
                file.path(b03patchImportanceDir, 
                          paste0(species, "_PC_Connector_", dispersalMode, "_", coarseResolution, "m.tif")), 
                overwrite=TRUE) 
  }
  
  ## Optional - transform outputs using log transform and or rescaling 0 - 1 ----
  
  # Functions
  #dPC
  pc <-PCfocal[[2]] %>%
    calc(., fun=rescaleR)
  
  #dPCintra
  intra <-PCfocal[[3]] %>%
    calc(., fun=rescaleR)
  
  #dPCflux
  flux <- PCfocal[[4]] %>%
    calc(., fun=rescaleR)    
  
  #dPCconnector
  connector <- PCfocal[[5]] %>%
    calc(., fun=function(x){log(x+10^-16)}) %>%
    calc(., fun=rescaleR)    
  
  
  ## Save transformed geotiffs
  if(withinBTSL){
  # dPC
  writeRaster(pc, 
              file.path(b03patchImportanceDir, 
                        paste0(species, "_PC_0-1_", dispersalMode, "_Focal_", coarseResolution, "m.tif")), 
              overwrite=TRUE)  
  #dPCintra
  writeRaster(intra, 
              file.path(b03patchImportanceDir, 
                        paste0(species, "_PC_Intra_0-1_", dispersalMode, "_Focal_", coarseResolution, "m.tif")), 
              overwrite=TRUE)  
  #dPCflux
  writeRaster(flux, 
              file.path(b03patchImportanceDir, 
                        paste0(species, "_PC_Flux_0-1_", dispersalMode, "_Focal_", coarseResolution, "m.tif")), 
              overwrite=TRUE)  
  #dPCconnector
  writeRaster(connector, 
              file.path(b03patchImportanceDir, 
                        paste0(species, "_PC_Connector_log&0-1_", dispersalMode, "_Focal_", coarseResolution, "m.tif")), 
              overwrite=TRUE)
  } else {
    # dPC
    writeRaster(pc, 
                file.path(b03patchImportanceDir, 
                          paste0(species, "_PC_0-1_", dispersalMode, "_", coarseResolution, "m.tif")), 
                overwrite=TRUE)  
    #dPCintra
    writeRaster(intra, 
                file.path(b03patchImportanceDir, 
                          paste0(species, "_PC_Intra_0-1_", dispersalMode, "_", coarseResolution, "m.tif")), 
                overwrite=TRUE)  
    #dPCflux
    writeRaster(flux, 
                file.path(b03patchImportanceDir, 
                          paste0(species, "_PC_Flux_0-1_", dispersalMode, "_", coarseResolution, "m.tif")), 
                overwrite=TRUE)  
    #dPCconnector
    writeRaster(connector, 
                file.path(b03patchImportanceDir, 
                          paste0(species, "_PC_Connector_log&0-1_", dispersalMode, "_", coarseResolution, "m.tif")), 
                overwrite=TRUE)
  }
#}