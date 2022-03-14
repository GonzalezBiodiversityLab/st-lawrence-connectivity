Sys.setenv(TZ='GMT')

library(rgrass7)

options(stringsAsFactors = FALSE)

# Directories
# Directories
projectDir <- "E:/"
gisBase <- "C:/Program Files/GRASS GIS 7.4.3"
gisDbase <- paste0(projectDir, "grass7")
rawTablesDir <- paste0(projectDir, "Inputs/RawData/Tables/")
processedMapsDir <- paste0(projectDir, "Inputs/ProcessedData/Maps/")
habitatDir <- paste0(projectDir, "Outputs/1.Habitat/")
resistanceDir <- paste0(projectDir, "Outputs/2.Resistance/")
networkDir <- paste0(projectDir, "Outputs/3.NetworkConnectivity/")
circuitDir <- paste0(projectDir, "Outputs/4.CircuitConnectivity/")
zonInputDir <- paste0(projectDir, "Outputs/5.HabitatPrioritization/ZonationInputs/")

# Input parameters
# Run this at 30m resolution
myResolution <- 30
speciesList<-c("MAAM", "PLCI", "RASY", "BLBR", "URAM")
conservationCriteriaList<-data.frame(Name=c("habitatSuitability", "betweennessShortBTSL", "dECGapBTSL", "dECNatalBTSL", "currentFlowBTSL"),
                                    FileName=c(paste0("habitatSuitability_01_Focal_", myResolution, "m.tif"), "betweenness_BTSL.tif", "PC_0-1_Gap_240m.tif", "PC_0-1_Natal_240m.tif", paste0("currentdensity_log_01_Focal_", myResolution, ".tif")),
                                    Directory=c(habitatDir, networkDir, networkDir, networkDir, networkDir, circuitDir))
# Set up GRASS mapset for the first time
doGRASSSetup <- F
# Name of BTSL ecological boundary
landcoverName <- paste0("landcover_", myResolution, "m.tif")
protectedAreasName <- paste0("protectedAreas_", myResolution,"m.tif")
studyAreaName <- paste0("studyArea_", myResolution, "m.tif")

###############
# GRASS setup #
###############
if(doGRASSSetup){
  #https://gis.stackexchange.com/questions/183032/create-a-new-grass-database-in-r-with-crs-projection
  # Manually set up empty GRASS database - see GRASSTemplate
  initGRASS(gisBase=gisBase, gisDbase=gisDbase, location=paste0("BTSL_", myResolution, "m"), mapset='PERMANENT', override=TRUE)
  
  execGRASS("g.proj", georef=paste0(processedMapsDir, landcoverName), flags="c")
  
  # Initialize new mapset inheriting projection info
  execGRASS("g.mapset", mapset="HabitatPrioritization", flags="c")
  
  # Load layers into grass database
  for(i in 1:length(speciesList)){
    species<-speciesList[i]
    for(j in 1:nrow(conservationCriteriaList)){
      mapName<-paste0(species, "_", conservationCriteriaList$FileName[j])
      execGRASS("r.in.gdal", input=paste0(conservationCriteriaList$Directory[j], species, "_", conservationCriteriaList$FileName[j]), output=paste0(species, "_", conservationCriteriaList$Name[j]), flags=c("overwrite", "o"))
    }
  }
  execGRASS("r.in.gdal", input=paste0(processedMapsDir, landcoverName), output="landcover", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=paste0(processedMapsDir, studyAreaName), output="studyArea1", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=paste0(processedMapsDir, protectedAreasName), output="protectedAreas", flags=c("overwrite", "o"))
}else{
  initGRASS(gisBase=gisBase, gisDbase=gisDbase, location=paste0("BTSL_", myResolution, "m"), mapset='HabitatPrioritization', override=TRUE)
}

# Set the geographic region
execGRASS('g.region', n='403680', e='-95520', w='-491430', s='117960', res=paste0(myResolution))

# Massage study area and protected areas layers
execGRASS("r.mapcalc", expression='studyArea0=if(studyArea1>0,0)', flags=c('overwrite'))
execGRASS('r.reclass.area', input='protectedAreas', output='protectedAreas150', value=150, mode='greater', flags=c('overwrite', 'd'))

# Create a binary Zonation mask 
# Reclassify the landcover map to idenitfy all natural areas that could be potential habitat
rcl<-read.csv(paste0(rawTablesDir,"speciesLandcoverReclass.csv"),header=TRUE)
attach(rcl)
rcl1<-NULL
for(i in 1:length(speciesList)){
  speciesHabitat<-rcl[rcl[speciesList[i]]>=60,]
  rcl1<-merge(rcl1,speciesHabitat, all.y=T)
}
write.table(c(paste0(rcl1[,'LandcoverCode'],'=1'),'*=NULL'),'rule.txt',sep="",col.names=FALSE,quote=FALSE,row.names=FALSE)
execGRASS('r.reclass',input='landcover',output='Z_mask',rules='rule.txt',flags=c('overwrite'))
execGRASS("r.mapcalc", expression='Z_mask0_BTSL=if(Z_mask==1,0,1)*studyArea1', flags=c('overwrite'))

# Crop habitat suitability and betweenness long to BTSL
for(i in 1:length(speciesList)){
  species<-speciesList[i]
  execGRASS('r.mapcalc',expression=paste0(species, '_habitatSuitabilityBTSL=', species, '_habitatSuitability*studyArea1'), flags=c('overwrite'))
  execGRASS('r.mapcalc',expression=paste0(species, '_betweennessLongBTSL=', species, '_betweennessLong*studyArea1'), flags=c('overwrite'))
}
# Update working names of conervation criteria
conservationCriteriaNames <- conservationCriteriaList$Name 
conservationCriteriaNames[conservationCriteriaNames=="habitatSuitability"] <- "habitatSuitabilityBTSL"
conservationCriteriaNames[conservationCriteriaNames=="betweennessLong"] <- "betweennessLongBTSL"

# Make all areas outside of natural areas NULL using Zonation binary mask
for(i in 1:length(speciesList)){
  species<-speciesList[i]
  for(j in 1:length(conservationCriteriaNames)){
    execGRASS('r.mapcalc',expression=paste0(species, '_', conservationCriteriaNames[j], '_NatAreas=', species, '_', conservationCriteriaNames[j] , '*Z_mask'), flags=c('overwrite'))
  }
}
# Update working names of conervation criteria
conservationCriteriaNames <- paste0(conservationCriteriaNames, "_NatAreas") 

# Rescale habitat suitability 0 - 1
for(i in 1:length(speciesList)){
  species<-speciesList[i]
  execGRASS('r.mapcalc', expression=paste0(species, '_habitatSuitabilityBTSL_NatAreas_01=', species, '_habitatSuitabilityBTSL_NatAreas/100.0'), flags=c('overwrite'))
}  
# Update working names of conervation criteria
conservationCriteriaNames[1] <- paste0(conservationCriteriaNames[1], "_01") 


# Log-transform and rescale connectivity layers 0 - 1
connectivityCriteriaNames <- conservationCriteriaNames[-1]
for(i in 1:length(speciesList)){
  species<-speciesList[i]
  for(j in 1:length(connectivityCriteriaNames)){
    #log transformation
    execGRASS('r.mapcalc',expression=paste0(species, '_', connectivityCriteriaNames[j],'_Log=log(', species, '_', connectivityCriteriaNames[j],',10)'), flags=c('overwrite'))

    #rescale 0 - 1
    statistic <- execGRASS("r.univar", map = paste0(species, '_', connectivityCriteriaNames[j], '_Log'), flags = c("g", "quiet"), intern = TRUE, legacyExec = TRUE)
    mapMax <- format(as.numeric(strsplit(statistic[agrep(statistic, pattern = "max")], "=")[[1]][[2]]), nsmall=2)
    mapMin <- format(as.numeric(strsplit(statistic[agrep(statistic, pattern = "min")], "=")[[1]][[2]]), nsmall=2)
    execGRASS('r.mapcalc', expression=paste0(species, '_', connectivityCriteriaNames[j], '_Log01=1/(',mapMax,'-',mapMin,')*(',species, '_', connectivityCriteriaNames[j], '_Log-', mapMin, ')'), flags=c('overwrite'))
  }
}  
conservationCriteriaNames[2:6] <- paste0(conservationCriteriaNames[2:6], "_Log01") 

# Rescale connectivity layers 0 - 1 without transforming them
connectivityCriteriaNames <- connectivityCriteriaNames[-5]
for(i in 1:length(speciesList)){
  species<-speciesList[i]
  for(j in 1:length(connectivityCriteriaNames)){
    #rescale 0 - 1
    statistic <- execGRASS("r.univar", map = paste0(species, '_', connectivityCriteriaNames[j]), flags = c("g", "quiet"), intern = TRUE, legacyExec = TRUE)
    mapMax <- format(as.numeric(strsplit(statistic[agrep(statistic, pattern = "max")], "=")[[1]][[2]]), nsmall=2)
    mapMin <- format(as.numeric(strsplit(statistic[agrep(statistic, pattern = "min")], "=")[[1]][[2]]), nsmall=2)
    execGRASS('r.mapcalc', expression=paste0(species, '_', connectivityCriteriaNames[j], '_01=1/(',mapMax,'-',mapMin,')*(',species, '_', connectivityCriteriaNames[j], '-', mapMin, ')'), flags=c('overwrite'))
  }
}  

# Calculate max-min of each map
maxMinTable <- data.frame(Species=as.character(), MapName=as.character(), Min=as.numeric(), Max=as.numeric())
k<-0
for(i in 1:length(speciesList)){
  species<-speciesList[i]
  for(j in 1:length(conservationCriteriaNames)){
    k<-k+1
    maxMinTable[k,"Species"]<-species
    maxMinTable[k,"MapName"]<-paste0(species, "_", conservationCriteriaNames[j])
    statistic <- execGRASS("r.univar", map = maxMinTable[k,"MapName"], flags = c("g", "quiet"), intern = TRUE, legacyExec = TRUE)
    maxMinTable[k, "Max"] <- as.numeric(strsplit(statistic[agrep(statistic, pattern = "max")], "=")[[1]][[2]])
    maxMinTable[k, "Min"] <- as.numeric(strsplit(statistic[agrep(statistic, pattern = "min")], "=")[[1]][[2]])
  }
  for(m in 1:length(connectivityCriteriaNames)){
    k<-k+1
    maxMinTable[k,"Species"]<-species
    maxMinTable[k,"MapName"]<-paste0(species, "_", connectivityCriteriaNames[m], '_01')
    statistic <- execGRASS("r.univar", map = maxMinTable[k,"MapName"], flags = c("g", "quiet"), intern = TRUE, legacyExec = TRUE)
    maxMinTable[k, "Max"] <- as.numeric(strsplit(statistic[agrep(statistic, pattern = "max")], "=")[[1]][[2]])
    maxMinTable[k, "Min"] <- as.numeric(strsplit(statistic[agrep(statistic, pattern = "min")], "=")[[1]][[2]])
  }
}
write.csv(maxMinTable, paste0(zonInputDir, "MaxMinTable.csv"), row.names=F)

#Summed species maps
execGRASS('r.mapcalc',expression='all_habitatSuitabilityBTSL_NatAreas=(if(isnull(MAAM_habitatSuitabilityBTSL_NatAreas_01),0,MAAM_habitatSuitabilityBTSL_NatAreas_01)+if(isnull(BLBR_habitatSuitabilityBTSL_NatAreas_01),0,BLBR_habitatSuitabilityBTSL_NatAreas_01)+if(isnull(URAM_habitatSuitabilityBTSL_NatAreas_01),0,URAM_habitatSuitabilityBTSL_NatAreas_01)+if(isnull(PLCI_habitatSuitabilityBTSL_NatAreas_01),0,PLCI_habitatSuitabilityBTSL_NatAreas_01)+if(isnull(RASY_habitatSuitabilityBTSL_NatAreas_01),0,RASY_habitatSuitabilityBTSL_NatAreas_01))*studyArea1', flags=c('overwrite'))
execGRASS('r.mapcalc',expression='all_currentFlowBTSL_NatAreas=(if(isnull(MAAM_currentFlowBTSL_NatAreas_01),0,MAAM_currentFlowBTSL_NatAreas_01)+if(isnull(BLBR_currentFlowBTSL_NatAreas_01),0,BLBR_currentFlowBTSL_NatAreas_01)+if(isnull(URAM_currentFlowBTSL_NatAreas_01),0,URAM_currentFlowBTSL_NatAreas_01)+if(isnull(PLCI_currentFlowBTSL_NatAreas_01),0,PLCI_currentFlowBTSL_NatAreas_01)+if(isnull(RASY_currentFlowBTSL_NatAreas_01),0,RASY_currentFlowBTSL_NatAreas_01))*studyArea1', flags=c('overwrite'))
execGRASS('r.mapcalc',expression='all_betweennessLongBTSL_NatAreas=(if(isnull(MAAM_betweennessLongBTSL_NatAreas_01),0,MAAM_betweennessLongBTSL_NatAreas_01)+if(isnull(BLBR_betweennessLongBTSL_NatAreas_01),0,BLBR_betweennessLongBTSL_NatAreas_01)+if(isnull(URAM_betweennessLongBTSL_NatAreas_01),0,URAM_betweennessLongBTSL_NatAreas_01)+if(isnull(PLCI_betweennessLongBTSL_NatAreas_01),0,PLCI_betweennessLongBTSL_NatAreas_01)+if(isnull(RASY_betweennessLongBTSL_NatAreas_01),0,RASY_betweennessLongBTSL_NatAreas_01))*studyArea1', flags=c('overwrite'))
execGRASS('r.mapcalc',expression='all_betweennessShortBTSL_NatAreas=(if(isnull(MAAM_betweennessShortBTSL_NatAreas_01),0,MAAM_betweennessShortBTSL_NatAreas_01)+if(isnull(BLBR_betweennessShortBTSL_NatAreas_01),0,BLBR_betweennessShortBTSL_NatAreas_01)+if(isnull(URAM_betweennessShortBTSL_NatAreas_01),0,URAM_betweennessShortBTSL_NatAreas_01)+if(isnull(PLCI_betweennessShortBTSL_NatAreas_01),0,PLCI_betweennessShortBTSL_NatAreas_01)+if(isnull(RASY_betweennessShortBTSL_NatAreas_01),0,RASY_betweennessShortBTSL_NatAreas_01))*studyArea1', flags=c('overwrite'))
execGRASS('r.mapcalc',expression='all_dECGapBTSL_NatAreas=(if(isnull(MAAM_dECGapBTSL_NatAreas_01),0,MAAM_dECGapBTSL_NatAreas_01)+if(isnull(BLBR_dECGapBTSL_NatAreas_01),0,BLBR_dECGapBTSL_NatAreas_01)+if(isnull(URAM_dECGapBTSL_NatAreas_01),0,URAM_dECGapBTSL_NatAreas_01)+if(isnull(PLCI_dECGapBTSL_NatAreas_01),0,PLCI_dECGapBTSL_NatAreas_01)+if(isnull(RASY_dECGapBTSL_NatAreas_01),0,RASY_dECGapBTSL_NatAreas_01))*studyArea1', flags=c('overwrite'))
execGRASS('r.mapcalc',expression='all_dECNatalBTSL_NatAreas=(if(isnull(MAAM_dECNatalBTSL_NatAreas_01),0,MAAM_dECNatalBTSL_NatAreas_01)+if(isnull(BLBR_dECNatalBTSL_NatAreas_01),0,BLBR_dECNatalBTSL_NatAreas_01)+if(isnull(URAM_dECNatalBTSL_NatAreas_01),0,URAM_dECNatalBTSL_NatAreas_01)+if(isnull(PLCI_dECNatalBTSL_NatAreas_01),0,PLCI_dECNatalBTSL_NatAreas_01)+if(isnull(RASY_dECNatalBTSL_NatAreas_01),0,RASY_dECNatalBTSL_NatAreas_01))*studyArea1', flags=c('overwrite'))
execGRASS('r.mapcalc',expression='all_connectivityBTSL_NatAreas=all_currentFlowBTSL_NatAreas+all_betweennessLongBTSL_NatAreas+all_betweennessShortBTSL_NatAreas+all_dECGapBTSL_NatAreas+all_dECNatalBTSL_NatAreas', flags=c('overwrite'))
execGRASS('r.mapcalc',expression='all_qualityconnectivityBTSL_NatAreas=all_habitatSuitabilityBTSL_NatAreas+all_currentFlowBTSL_NatAreas+all_betweennessLongBTSL_NatAreas+all_betweennessShortBTSL_NatAreas+all_dECGapBTSL_NatAreas+all_dECNatalBTSL_NatAreas', flags=c('overwrite'))

#Save summed species maps
execGRASS('r.out.gdal',input='all_habitatSuitabilityBTSL_NatAreas',output=paste0(networkDir,'all_habitatSuitability.tif'),format='GTiff',createopt='COMPRESS=LZW',flags=c('overwrite'))
execGRASS('r.out.gdal',input='all_currentFlowBTSL_NatAreas',output=paste0(networkDir,'all_currentFlow.tif'),format='GTiff',createopt='COMPRESS=LZW',flags=c('overwrite'))
execGRASS('r.out.gdal',input='all_betweennessLongBTSL_NatAreas',output=paste0(networkDir,'all_betweennessLong.tif'),format='GTiff',createopt='COMPRESS=LZW',flags=c('overwrite'))
execGRASS('r.out.gdal',input='all_betweennessShortBTSL_NatAreas',output=paste0(networkDir,'all_betweennessShort.tif'),format='GTiff',createopt='COMPRESS=LZW',flags=c('overwrite'))
execGRASS('r.out.gdal',input='all_dECGapBTSL_NatAreas',output=paste0(networkDir,'all_dECGap.tif'),format='GTiff',createopt='COMPRESS=LZW',flags=c('overwrite'))
execGRASS('r.out.gdal',input='all_dECNatalBTSL_NatAreas',output=paste0(networkDir,'all_dECNatal.tif'),format='GTiff',createopt='COMPRESS=LZW',flags=c('overwrite'))
execGRASS('r.out.gdal',input='all_connectivityBTSL_NatAreas',output=paste0(networkDir,'all_connectivity.tif'),format='GTiff',createopt='COMPRESS=LZW',flags=c('overwrite'))
execGRASS('r.out.gdal',input='all_qualityconnectivityBTSL_NatAreas',output=paste0(networkDir,'all_qualityconnectivity.tif'),format='GTiff',createopt='COMPRESS=LZW',flags=c('overwrite'))


# Set geographic region for Zonation analysis (add 1 row of NA pixels all around each map)
execGRASS('g.region', n=paste0('n+',myResolution), e=paste0('e+',myResolution), w=paste0('w-',myResolution), s=paste0('s-',myResolution), flags='c')

# Write out habitat suitability and log-transformed connectivity layers to geotiffs
for(i in 1:length(speciesList)){
  species<-speciesList[i]
  for(j in 1:nrow(conservationCriteriaList)){
    execGRASS('r.patch', input=paste0(species, '_', conservationCriteriaNames[j], ',Z_mask0_BTSL'),output=paste0('Z_', species, '_', conservationCriteriaNames[j]),flags=c('overwrite'))
    execGRASS('r.out.gdal',input=paste0('Z_', species, '_', conservationCriteriaNames[j]),output=paste0(zonInputDir,species, '_', conservationCriteriaNames[j], '_Z.tif'),format='GTiff',flags=c('overwrite'))
  }
}

# Write out non-log-transformed connectivity layers to geotiffs
for(i in 1:length(speciesList)){
  species<-speciesList[i]
  for(j in 1:length(connectivityCriteriaNames)){
    execGRASS('r.patch', input=paste0(species, '_', connectivityCriteriaNames[j], '_01,Z_mask0_BTSL'),output=paste0('Z_', species, '_', connectivityCriteriaNames[j], '_01'),flags=c('overwrite'))
    execGRASS('r.out.gdal',input=paste0('Z_', species, '_', connectivityCriteriaNames[j], '_01'),output=paste0(zonInputDir,species, '_', connectivityCriteriaNames[j], '_01_Z.tif'),format='GTiff',flags=c('overwrite'))
  }
}

# Change the extent of the protected areas raster to match the geographic regions for Zonation analysis
execGRASS('r.mapcalc', expression='Z_protectedAreas=protectedAreas', flags=c('overwrite'))
execGRASS('r.out.gdal',input='Z_protectedAreas',output=paste0(zonInputDir, 'Z_protectedAreas.tif'),format='GTiff',flags=c('overwrite'))

# Reset the geographic region
execGRASS('g.region', n='403680', e='-95520', w='-491430', s='117960', res=paste0(myResolution))
