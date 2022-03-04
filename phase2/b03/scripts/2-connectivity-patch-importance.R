#####################################################################
# a254      
# Calculate patch importance metrics for focal species
# 02-2022                                       					
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

## Workspace ---------------------------------------------------------

# Packages
library(Makurhini)

# Spatial data
studyArea <- focalArea <- st_read(file.path(b03ProcessedMapsDir, "b03-studyarea.shp")) # Study area polygon

# Tabular data
dispersalDistance <- read_csv(file.path(b01b02RawTablesDir, "speciesDispersalParameters.csv"))

## Calculate patch importance for focal species -----------------------------------------

# NB the  MK_dPCIIC() function does not work inside a loop so manually set each species for now
# Loop through species, generating Patch Importance data and saving outputs
# List of species
# speciesList <- dispersalDistance$Species

#for (i in speciesList[c(5)]){
  
#  species <- i
  
  ## Load data for focal species
  # Spatial data
  habitatPatches <- raster(file.path(habitatDir, paste0(species, "_habitatPatch_", myResolution, "m.tif")))
  resistance <- raster(file.path(resistanceDir, paste0(species, "_resistance_", myResolution, "m.tif")))
  
  # Species dispersal maximal/or median distance 
  maxdist <- dispersalDistance[[which(dispersalDistance$Species == species), dispersalMode]]
  # Albert distance estimates represent the Median
  prob <- 0.5
  
  # # NB temporary fix projection of habitat patch and resistance rasters, follow up with Sarah
  # # Reproject habitat and resistance 
  # habitatPatches <- projectRaster(habitatPatchesRaw, crs=crs(studyArea))
  # resistance <- projectRaster(resistanceRaw, crs=crs(studyArea))
  # 
  # # NB temporary fix of habitatPatches layer - check with Sarah
  # habitatPatches[habitatPatches < 1] <- 0
  
  # Crop habitat and resistance to study area
  habitatPatchesFocal <- habitatPatches %>%
    crop(., extent(studyArea), snap="out") %>% # Crop to focal area extent
    mask(., mask=studyArea) %>% # Clip to focal area
    trim(.) # Trim extra white spaces  
  
  resistanceFocal <- resistance %>%
    crop(., extent(studyArea), snap="out") %>% # Crop to focal area extent
    mask(., mask=studyArea) %>% # Clip to focal area
    trim(.) # Trim extra white spaces  
  
  # Create raster with unique patch id for each species
  habitatPatchesIDFocal <- clump(habitatPatchesFocal, directions=8)
  length(unique(habitatPatchesIDFocal))
  
  ## PC fractions measurements for focal species in focal Area ---------------
  
  # PC & fractions
  PCfocal <- MK_dPCIIC(nodes = habitatPatchesIDFocal, 
                       distance = list(type = "least-cost", 
                                       resistance = resistanceFocal),
                       attribute = NULL,
                       metric = "PC", 
                       probability = prob, 
                       distance_thresholds = maxdist)
  
  # PCfocal is a raster stack : plot individual rasters
  plot(PCfocal)
  
  ## Save single-band geotiffs ---------------------------------------------
  # dPC
  writeRaster(PCfocal[[2]], 
              file.path(b03patchImportanceDir, 
                        paste0(species, "_PC_", dispersalMode, "_", myResolution, "m.tif")), 
              overwrite=TRUE)  
  #dPCintra
  writeRaster(PCfocal[[3]], 
              file.path(b03patchImportanceDir, 
                        paste0(species, "_PC_Intra_", dispersalMode, "_", myResolution, "m.tif")), 
              overwrite=TRUE)  
  #dPCflux
  writeRaster(PCfocal[[4]], 
              file.path(b03patchImportanceDir, 
                        paste0(species, "_PC_Flux_", dispersalMode, "_", myResolution, "m.tif")), 
              overwrite=TRUE)  
  #dPCconnector
  writeRaster(PCfocal[[5]], 
              file.path(b03patchImportanceDir, 
                        paste0(species, "_PC_Connector_", dispersalMode, "_", myResolution, "m.tif")), 
              overwrite=TRUE)  		
  
  
  ## Optional - transform outputs using log transform and or rescaling 0 - 1------------
  
  # Functions
  rescaleR <- function(x, new.min = 0, new.max = 1) {
    x.min = suppressWarnings(min(x, na.rm=TRUE))
    x.max = suppressWarnings(max(x, na.rm=TRUE))
    new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))
  }
  
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
  # dPC
  writeRaster(pc, 
              file.path(b03patchImportanceDir, 
                        paste0(species, "_PC_0-1_", dispersalMode, "_", myResolution, "m.tif")), 
              overwrite=TRUE)  
  #dPCintra
  writeRaster(intra, 
              file.path(b03patchImportanceDir, 
                        paste0(species, "_PC_Intra_0-1_", dispersalMode, "_", myResolution, "m.tif")), 
              overwrite=TRUE)  
  #dPCflux
  writeRaster(flux, 
              file.path(b03patchImportanceDir, 
                        paste0(species, "_PC_Flux_0-1_", dispersalMode, "_", myResolution, "m.tif")), 
              overwrite=TRUE)  
  #dPCconnector
  writeRaster(connector, 
              file.path(b03patchImportanceDir, 
                        paste0(species, "_PC_Connector_log&0-1_", dispersalMode, "_", myResolution, "m.tif")), 
              overwrite=TRUE)  		
#}