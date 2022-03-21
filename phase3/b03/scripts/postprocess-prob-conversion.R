
# Load constants, functions, etc
source("./b03/scripts/0-constants.R")


options(tibble.width = Inf, tibble.print_min = Inf)


#### Inputs ####
# Input parameters
resultTag <- c("BAU_85")
numScenarios <- length(resultTag)
numIterations <- 10
timestep1 = 2010
timestep2 = 2110

# Input files and folders
fig5Dir <- "../Report/phase3/Figures/Figure5"

# Primary stratum
primaryStratumFileName <- "PrimaryStratum_90m_b03Extent.tif"

# STSim
mySession <- session()
modelDir <- c(b03Libraries)
resultScenarioNumber <- c(211)
modelFile <- "BTSL_stconnect.ssim"

#### Create clipping mask
primaryStratum <- raster(file.path(b03ProcessedMapsDir, primaryStratumFileName))
# BTSL map
btslMask <- Which(primaryStratum %in% c(3))
btslMask[btslMask == 0] <- NA

#Shapefile with ecoregions
bt4 <- st_read(file.path(b03RawMapsDir, "CR_CERQ_NIV_04_S_b03.shp"))
bt6=fasterize(bt4, primaryStratum, field="FID04")
writeRaster(bt6, file.path(b03RawMapsDir, "CR_CERQ_NIV_04_S_b03.tif"), overwrite = TRUE)

#Saint-Lawrence Lowlands
bt=st_read(file.path(b03ProcessedMapsDir, "b03-studyarea.shp"))
btp=as(bt, "Spatial")

# State class maps and transitions ---------------------------------------------------
# Loop over scenarios, loop over iterations
for(scn in 1:numScenarios){
  print(paste("Working on scenario", resultTag[scn]))
  
  # Setting up st-sim library, project, scenario
  myLibrary <- ssimLibrary(name = file.path(modelDir[scn], modelFile), session = mySession, package = "stconnect")
  myProject <- project(ssimObject = myLibrary, project = "Definitions")  # Assumes there is only one default project per library
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
  
  writeRaster(natAreasProbConversion, file.path(b03ResultsMapsDir, paste0(resultTag[scn], "_", timestep1, "-", timestep2, "_NaturalAreaProbConversion_BTSL_b03.tif")), overwrite=T)

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
  file.path(fig5Dir, paste0(resultTag[scn], "_", timestep1, "-", timestep2, "_NaturalAreaProbConversion_BTSL_b03.png")), q, height=4, width=11)

}  


