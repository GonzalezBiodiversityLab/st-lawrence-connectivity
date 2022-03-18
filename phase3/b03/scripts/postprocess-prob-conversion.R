library(tidyverse)
library(rsyncrosim)
library(raster)
library(fasterize)
library(sf)
library(scales)


options(tibble.width = Inf, tibble.print_min = Inf)

#### Inputs ####
# Input parameters
resultTag <- c("CON_NC","BAU_NC")
numScenarios <- length(resultTag)
numIterations <- 40
timestep1 = 2010
timestep2 = 2110

# Input files and folders
workingDir <- "C:/Users/bronw/Documents/Apex/Projects/Active/A224_MDDELCC/stsim/ResultsSummary/"

# Primary stratum
spatialInitialConditionsDir <- "C:/Users/bronw/Documents/Apex/Projects/Active/A224_MDDELCC/stsim/SpatialInitialConditions/FullExtent/"
primaryStratumFileName <- "PrimaryStratum.tif"

####MAPPING FUNCTIONS####

#Assemble all the pieces and map them together
#Define the map theme
theme_map <- function(...) {
  theme_minimal() +
    theme(
      #text = element_text(family = "Arial", color = "#22211d"),
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.background = element_rect(fill = "white", color = NA), 
      panel.background = element_rect(fill = "white", color = NA), 
      legend.background = element_rect(fill = "white", color = NA),
      #panel.border = element_rect(colour = "grey", fill=NA, size=1),
      legend.position="bottom",
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(-10,-10,-10,-10),
      plot.margin=grid::unit(c(0,0,0,0), "mm")
    )
}



# STSim
SyncroSimDir <- "C:/Users/bronw/Documents/Apex/SyncroSim/2-2-10/"
modelDir <- c("D:/Apex/Projects/A224/stsim/CON_NC_backup/", 
              "D:/Apex/Projects/A224/stsim/BAU_NC_backup/")
resultScenarioNumber <- c(166, 164)
modelFile <- "BTSL_stconnect.ssim"

#### Create clipping mask
primaryStratum <- raster(paste0(spatialInitialConditionsDir, primaryStratumFileName))
# BTSL map
btslMask <- Which(primaryStratum %in% c(3))
btslMask[btslMask == 0] <- NA

#Shapefile with ecoregions
bt4 <- st_read("C:/Users/bronw/Documents/Apex/Projects/Active/A224_MDDELCC/Data/CERO04_BTSL20201127.shp")
bt6=fasterize(bt4, primaryStratum, field="FID04")
writeRaster(bt6, "C:/Users/bronw/Documents/Apex/Projects/Active/A224_MDDELCC/Data/CERO04_BTSL20201127.tif")

#Saint-Lawrence Lowlands
bt=st_read("C:/Users/bronw/Documents/Apex/Projects/Active/A224_MDDELCC/Data/btsl_90m_polygon.shp")
btp=as(bt, "Spatial")

# State class maps and transitions ---------------------------------------------------
# Loop over scenarios, loop over iterations
for(scn in 1:numScenarios){
  print(paste("Working on scenario", resultTag[scn]))
  
  # Setting up st-sim library, project, scenario
  mySession <- session(SyncroSimDir)
  myLibrary <- ssimLibrary(paste0(modelDir[scn], modelFile), session=mySession)
  myProject <- project(myLibrary, "Definitions")  # Assumes there is only one default project per library
  myScenario <- scenario(myProject, resultScenarioNumber[scn])
  #datasheet(myScenario)
  
  #### Create maps of state class probability ####
  # Get state class map for iteration 1
  print(paste("Working on scenario", resultTag[scn], "iteration 1"))
  stateRas <- datasheetRaster(myScenario, "stsim_OutputSpatialState", iteration = 1, timestep=timestep1)
  # Crop to focal stratum
  stateRas <- stateRas * btslMask
  # Isolate natural areas (forest and wetland)
  natAreas1 <- Which(stateRas %in% c(800, 810, 511, 512, 513, 521, 522, 523, 531, 532, 533))
  
  # Get state class map for iteration 1, timestep 2
  stateRas <- datasheetRaster(myScenario, "stsim_OutputSpatialState", iteration = 1, timestep=timestep2)
  # Crop to focal stratum
  stateRas <- stateRas * btslMask
  # Isolate natural areas (forest and wetland)
  natAreas2 <- Which(stateRas %in% c(800, 810, 511, 512, 513, 521, 522, 523, 531, 532, 533))
  
  # Difference between natural areas in timestep 1 and timestep 2
  natAreas12 <- natAreas1-natAreas2
  # Crop to focal stratum
  natAreas12 <- natAreas12 * btslMask
  
  # Create mask for natural areas in 2010
  natAreasMask <- natAreas1
  natAreasMask[natAreasMask==0] <- NA
  
  # Get state class maps for all other iterations
  for(i in 2:numIterations){
    print(paste("Working on scenario", resultTag[scn], "iteration", i))
    # Read in state class rasters for start and end of time interval
    stateRas <- datasheetRaster(myScenario, "stsim_OutputSpatialState", iteration = i, timestep=timestep1)
    # Crop to focal stratum
    stateRas <- stateRas * btslMask
    # Isolate natural areas (forest and wetland)
    natAreas1 <- Which(stateRas %in% c(800, 810, 511, 512, 513, 521, 522, 523, 531, 532, 533))
    
    # Get state class map for iteration 1, timestep 2
    stateRas <- datasheetRaster(myScenario, "stsim_OutputSpatialState", iteration = i, timestep=timestep2)
    # Crop to focal stratum
    stateRas <- stateRas * btslMask
    # Isolate natural areas (forest and wetland)
    natAreas2 <- Which(stateRas %in% c(800, 810, 511, 512, 513, 521, 522, 523, 531, 532, 533))
    
    # Difference between natural areas in timestep 1 and timestep 2
    natAreas12 <- natAreas12 + (natAreas1 - natAreas2)
    # Crop to focal stratum
    natAreas12 <- natAreas12 * btslMask
  }
  
  # Divide by num iterations to get probability
  # Keep only areas within natural areas
  # Pixel-level summary
  natAreasProbConversion <- (natAreas12 / numIterations)  * natAreasMask
  
  writeRaster(natAreasProbConversion, paste0(workingDir, resultTag[scn], "_", timestep1, "-", timestep2, "_NaturalAreaProbConversion_BTSL.tif"), overwrite=T)

  # Ecoregion level-4 summary
  
  values(natAreasProbConversion)[is.na(values(bt6))]=NA
  mns<-data.frame(zonal(natAreasProbConversion, bt6, fun="mean"))
  sms<-data.frame(zonal(natAreasProbConversion, bt6, fun="sum"))
  ecv<-merge(mns, sms, by="zone")
  names(ecv)<-c("FID04", "mean", "sum")
  #head(ecv)
  
  ecr2=merge(bt4, ecv[c("FID04", "mean")])
  
  write.csv(ecv, paste0(workingDir, resultTag[scn], "_", timestep1, "-", timestep2, "_NaturalAreaProbConversion_BTSL_Ecoregion.csv"), row.names = F)
  st_write(ecr2, paste0(workingDir, resultTag[scn], "_", timestep1, "-", timestep2, "_NaturalAreaProbConversion_BTSL_Ecoregion.shp"))
  
  # Ten quantiles
  qn <- rescale(quantile(
    ecv$mean,
    probs=seq(0, 1, length.out=64)))

  # Min and max vlaues for legend
  mnV <- as.numeric(min(ecv$mean))
  mxV <- as.numeric(max(ecv$mean))

  q=ggplot() +
    geom_polygon(data = btp, aes(x = long+2000, #shadow
                                 y = lat-800,
                                 group=group),
                 color = "grey", size = 1, fill="grey") +
    geom_sf(data = bt, size = 1, fill="white", color="black") + #btsl
    geom_sf(data = ecr2, aes(fill=mean), size = 0.03, color="black") +
    coord_sf() +
    theme_map() +
    labs(x = NULL,
         y = NULL) +
    scale_fill_distiller(palette="Reds",
                         direction = 1,
#                         Name="Mean probability of conversion",
                           # here we use guide_colourbar because it is still a continuous scale
                          guide = guide_colorbar(
                            title = "Mean probability of conversion",
                             direction = "horizontal",
                             barheight = unit(2, units = "mm"),
                             barwidth = unit(50, units = "mm"),
                             draw.ulim = F,
                             title.position = 'top',
                             # some shifting around
                             title.hjust = 0.5,
                             label.hjust = 0.5
                           ))

ggsave(
  paste0(workingDir, resultTag[scn], "_", timestep1, "-", timestep2, "_NaturalAreaProbConversion_BTSL_Ecoregion1.png"), q, height=4, width=11)

}  


