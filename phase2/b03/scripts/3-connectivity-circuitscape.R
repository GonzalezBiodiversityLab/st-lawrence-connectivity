#####################################################################
# a254      
# Run Circuitscape for focal species
# 02-2022                                       					
#                             
#	  Uses Circuitscape          
#                                             
#   Inputs (for focal species):
#    -resistance layer
#  
#   Outputs:
#    -current density map
#                                                                   
# Script by B Rayfield for ApexRMS 									
#####################################################################

## Workspace ---------------------------------------------------------

# Packages
library(tidyverse)
library(sf)
library(raster)

# Directories
b03Dir <- "b03"
b01b02Dir <- "b01b02"
habitatDir <- file.path(b03Dir, "model-outputs", "spatial", "1.Habitat")
resistanceDir <- file.path(b03Dir, "model-outputs", "spatial", "2.Resistance")
b01b02RawTablesDir <- file.path(b01b02Dir, "inputs", "rawData", "tables")
b03ProcessedMapsDir <- file.path(b03Dir, "model-inputs", "spatial")
b03CircuitscapeDir <- file.path(b03Dir, "model-outputs", "spatial", "5.Circuitscape")

# Spatial data
studyArea <- focalArea <- st_read(file.path(b03ProcessedMapsDir, "b03-studyarea.shp")) # Study area polygon

# Tabular data
dispersalDistance <- read_csv(file.path(b01b02RawTablesDir, "speciesDispersalParameters.csv"))

# Run settings
myResolution <- 30

# Functions
rescaleR <- function(x, new.min = 0, new.max = 1) {
  x.min = suppressWarnings(min(x, na.rm=TRUE))
  x.max = suppressWarnings(max(x, na.rm=TRUE))
  new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))
}

## Run Circuitscape for focal species -----------------------------------------
# Loop through species, generating Circuitscape data and saving outputs
# List of species
speciesList <- dispersalDistance$Species

for (i in speciesList){
  
  species <- i
  
  ## Load data for focal species
  resistance <- raster(file.path(resistanceDir, paste0(species, "_resistance_", myResolution, "m.tif")))
  
  # Create focal regions
  # Create focal region raster for N-S by adding top and bottom rows
  extentNS<-extent(xmin(resistance),xmax(resistance),ymin(resistance)-res(resistance)[1],ymax(resistance)+res(resistance)[1])
  focalRegionNS<-raster(extentNS,nrow=nrow(resistance)+2,ncol=ncol(resistance))
  focalRegionNS[1,]<-1
  focalRegionNS[nrow(focalRegionNS),]<-2
  
  # Create focal region raster for E-W by adding left and right columns
  extentEW<-extent(xmin(resistance)-res(resistance)[2],xmax(resistance)+res(resistance)[2],ymin(resistance),ymax(resistance))
  focalRegionEW<-raster(extentEW,nrow=nrow(resistance),ncol=ncol(resistance)+2)
  focalRegionEW[,1]<-1
  focalRegionEW[,ncol(focalRegionEW)]<-2
  
  # Extend resistance rasters to match NS and EW extents
  resistanceNS<-extend(resistance, extentNS,value=1)
  resistanceNS[is.na(resistanceNS)]<-100
  resistanceEW<-extend(resistance, extentEW,value=1)
  resistanceEW[is.na(resistanceEW)]<-100
  
  # Save focal region and extended resistance rasters to temporary folder as Circuitscape inputs
  # Focal region rasters
  focalRegionNSName = file.path(b03CircuitscapeDir, species, "NSfocalRegion.asc")
  focalRegionEWName = file.path(b03CircuitscapeDir, species, "EWfocalRegion.asc")
  writeRaster(focalRegionNS, focalRegionNSName, overwrite = TRUE)
  writeRaster(focalRegionEW, focalRegionEWName, overwrite = TRUE)
  
  # Extended resistance rasters
  resistanceRasterNSName = file.path(getwd(), b03CircuitscapeDir, species, "NSResistance.asc")
  resistanceRasterEWName = file.path(getwd(), b03CircuitscapeDir, species, "EWResistance.asc")
  writeRaster(resistanceNS, resistanceRasterNSName, overwrite=TRUE)
  writeRaster(resistanceEW, resistanceRasterEWName, overwrite=TRUE)
  
  # Make .ini files and save to output folder
  # Note that these files specify the path for Circuitscape outputs to be written to temp folder
  NS_ini<-c("[circuitscape options]", 
            "data_type = raster",
            "scenario = pairwise",
            "write_cum_cur_map_only = true",
            "write_cur_maps = true",
            paste0("point_file = ", file.path(getwd(), b03CircuitscapeDir, species, "NSfocalRegion.asc")),
            paste0("habitat_file = ", resistanceRasterNSName),
            paste0("output_file = ", file.path(getwd(), b03CircuitscapeDir, species, "NS.out")))
  writeLines(NS_ini,file.path(getwd(), b03CircuitscapeDir, species, "NS.ini"))
  
  EW_ini<-c("[circuitscape options]", 
            "data_type = raster",
            "scenario = pairwise",
            "write_cum_cur_map_only = true",
            "write_cur_maps = true",
            paste0("point_file = ", file.path(getwd(), b03CircuitscapeDir, species, "EWfocalRegion.asc")),
            paste0("habitat_file = ", resistanceRasterEWName),
            paste0("output_file = ", file.path(getwd(), b03CircuitscapeDir, species, "EW.out")))
  writeLines(EW_ini,file.path(getwd(), b03CircuitscapeDir, species, "EW.ini"))
  
  # Make Julia scripts and save to output folder
  NS_run_jl <- file(file.path(getwd(), b03CircuitscapeDir, species, "NSscript.jl"))
  writeLines(c("using Circuitscape", paste0("compute(", "\"", file.path(b03CircuitscapeDir, species, "NS.ini"),"\")")), NS_run_jl)
  close(NS_run_jl)
  
  EW_run_jl <- file(file.path(getwd(), b03CircuitscapeDir, species, "EWscript.jl"))
  writeLines(c("using Circuitscape", paste0("compute(", "\"", file.path(b03CircuitscapeDir, species, "EW.ini"),"\")")), EW_run_jl)
  close(EW_run_jl)
  
  # Make the N-S and E-W Circuitscape run commands and save to output folder
  NS_run <- paste("C:\\Users\\Administrator\\AppData\\Local\\Programs\\Julia-1.7.2\\bin\\julia.exe" , paste0("\"",file.path(getwd(), b03CircuitscapeDir, species, "NSscript.jl"),"\""))
  EW_run <- paste("C:\\Users\\Administrator\\AppData\\Local\\Programs\\Julia-1.7.2\\bin\\julia.exe" , paste0("\"",file.path(getwd(), b03CircuitscapeDir, species, "EWscript.jl"),"\""))
  
  # Run the commands
  system(NS_run)
  system(EW_run)
  
  # Read in NS and EW cum current maps, crop, and combine
  currMapNS <- crop(raster(file.path(b03CircuitscapeDir, species, "NS_cum_curmap.asc")), resistance)
  currMapEW <- crop(raster(file.path(b03CircuitscapeDir, species, "EW_cum_curmap.asc")), resistance)
  currMapOMNI <- currMapEW + currMapNS
  currMapOMNI[is.na(resistance)] <- NA
  
  # Transform outputs using log transform and or rescaling 0 - 1
  currMapOMNI_01 <- currMapOMNI %>%
    calc(., fun=rescaleR)    
  
  currMapOMNI_log_01 <- currMapOMNI %>%
    calc(., fun=function(x){log(x+10^-16)}) %>%
    calc(., fun=rescaleR)    

  # Crop to study area
  currMapOMNIFocal <- currMapOMNI %>%
    crop(., extent(studyArea), snap="out") %>% # Crop to focal area extent
    mask(., mask=studyArea) %>% # Clip to focal area
    trim(.) # Trim extra white spaces  

  currMapOMNI_01Focal <- currMapOMNI_01 %>%
    crop(., extent(studyArea), snap="out") %>% # Crop to focal area extent
    mask(., mask=studyArea) %>% # Clip to focal area
    trim(.) # Trim extra white spaces  

  currMapOMNI_log_01Focal <- currMapOMNI_log_01 %>%
    crop(., extent(studyArea), snap="out") %>% # Crop to focal area extent
    mask(., mask=studyArea) %>% # Clip to focal area
    trim(.) # Trim extra white spaces  
  
  # Write outputs
  writeRaster(currMapOMNI, file.path(b03CircuitscapeDir, species, paste0(species, "_currentdensity_", myResolution, "m.tif")), overwrite=TRUE)  
  writeRaster(currMapOMNI_01, file.path(b03CircuitscapeDir, species, paste0(species, "_currentdensity_01_", myResolution, "m.tif")), overwrite=TRUE)  
  writeRaster(currMapOMNI_log_01, file.path(b03CircuitscapeDir, species, paste0(species, "_currentdensity_log_01_", myResolution, "m.tif")), overwrite=TRUE)  
  writeRaster(currMapOMNIFocal, file.path(b03CircuitscapeDir, species, paste0(species, "_currentdensity_Focal_", myResolution, "m.tif")), overwrite=TRUE)  
  writeRaster(currMapOMNI_01Focal, file.path(b03CircuitscapeDir, species, paste0(species, "_currentdensity_01_Focal_", myResolution, "m.tif")), overwrite=TRUE)  
  writeRaster(currMapOMNI_log_01Focal, file.path(b03CircuitscapeDir, species, paste0(species, "_currentdensity_log_01_Focal_", myResolution, "m.tif")), overwrite=TRUE)  
  
}
