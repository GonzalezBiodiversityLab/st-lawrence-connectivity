#####################################################################
# a254      
# Prepare Zonation inputs 
# 03-2022                                       					
#                             
#   Inputs (for focal species):
#    -habitat suitability, short and long betweenness,
#    -short and long dEC, current density
#    -protected areas
#  
#   Outputs:
#    -all layers formatted for Zonation input
#                                                                   
# Script by B Rayfield for ApexRMS 									
#####################################################################

# Workspace -------------------

# Load constants, functions, etc
source("./b03/scripts/0-0-constants.R")

options(stringsAsFactors = FALSE)
Sys.setenv(TZ='GMT')

# Input parameters
# List of Zonation feature layers
conservationCriteriaList<-data.frame(Name=c("habitatSuitabilityFocal", "betweennessShortFocal", "betweennessLongFocal", "dECGapFocal", "dECNatalFocal", "currentFlowFocal"),
                                     FileName=c(paste0("habitatSuitability_Focal_", myResolution, "m.tif"), paste0("Betweenness_b03_", myResolution, "m.tif"), paste0("Betweenness_", myResolution, "m.tif"), paste0("PC_Gap_Focal_", coarseResolution, "m.tif"), paste0("PC_Natal_Focal_", coarseResolution, "m.tif"), paste0("currentdensity_Focal_", myResolution, "m.tif")),
                                     Directory=c(b03habitatDir, b03networkDir, b03networkDir, b03patchImportanceDir, b03patchImportanceDir, b03circuitscapeDir))

# Names of additional layers
landcoverName <- paste0("b03-landcover-", myResolution, "m.tif")
protectedAreasName <- paste0("protected-areas-150ha-btsl-focal-", myResolution,"m.tif")
studyAreaName <- paste0("b03-studyarea-", myResolution, "m-trimmed.tif")

# Set up GRASS mapset for the first time
doGRASSSetup <- T

###############
# GRASS setup #
###############
if(doGRASSSetup){
  #https://gis.stackexchange.com/questions/183032/create-a-new-grass-database-in-r-with-crs-projection
  # Manually set up empty GRASS database - see GRASSTemplate
  initGRASS(gisBase=gisBase, gisDbase=gisDbase, location=paste0("BTSL_", myResolution, "m"), mapset='PERMANENT', override=TRUE)
  
  execGRASS("g.proj", georef=file.path(b03ProcessedMapsDir, landcoverName), flags="c")
  
  # Initialize new mapset inheriting projection info
  execGRASS("g.mapset", mapset="Zonation", flags="c")
  
  # Load layers into grass database
  execGRASS("r.in.gdal", input=file.path(b03ProcessedMapsDir, landcoverName), output="landcover", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=file.path(b03ProcessedMapsDir, studyAreaName), output="studyArea1", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=file.path(b03ProcessedMapsDir, protectedAreasName), output="protectedAreas", flags=c("overwrite", "o"))
  for(i in 1:length(speciesList)){
    species<-speciesList[i]
    for(j in 1:nrow(conservationCriteriaList)){
      mapName<-paste0(species, "_", conservationCriteriaList$FileName[j])
      execGRASS("r.in.gdal", input=file.path(conservationCriteriaList$Directory[j], paste0(species, "_", conservationCriteriaList$FileName[j])), output=paste0(species, "_", conservationCriteriaList$Name[j]), flags=c("overwrite", "o"))
    }
  }
  # Manually add BLBR patch importance (480m resolution)
  execGRASS("r.in.gdal", input=file.path(b03patchImportanceDir, "BLBR_PC_Gap_Focal_480m.tif"), output="BLBR_dECGapFocal", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=file.path(b03patchImportanceDir, "BLBR_PC_Natal_Focal_480m.tif"), output="BLBR_dECNatalFocal", flags=c("overwrite", "o"))
}else{
  initGRASS(gisBase=gisBase, gisDbase=gisDbase, location=paste0("BTSL_", myResolution, "m"), mapset='Zonation', override=TRUE)
}

# Set the geographic region
execGRASS('g.region', raster='studyArea1', res=paste0(myResolution))

# Massage study area and protected areas layers
execGRASS("r.mapcalc", expression='studyArea0=if(studyArea1>0,0)', flags=c('overwrite'))
execGRASS('r.reclass.area', input='protectedAreas', output='protectedAreas150', value=150, mode='greater', flags=c('overwrite', 'd'))

# Create a binary Zonation mask 
# Reclassify the landcover map to idenitfy all natural areas that could be potential habitat
rcl<-read.csv(file.path(b01b02RawTablesDir,"speciesLandcoverReclass.csv"),header=TRUE)
attach(rcl)
rcl1<-NULL
for(i in 1:length(speciesList)){
  speciesHabitat<-rcl[rcl[speciesList[i]]>=60,]
  rcl1<-merge(rcl1,speciesHabitat, all.y=T)
}
write.table(c(paste0(rcl1[,'LandcoverCode'],'=1'),'*=NULL'),'rule.txt',sep="",col.names=FALSE,quote=FALSE,row.names=FALSE)
execGRASS('r.reclass',input='landcover',output='Z_mask',rules='rule.txt',flags=c('overwrite'))
execGRASS("r.mapcalc", expression='Z_mask0_BTSL=if(Z_mask==1,0,1)*studyArea1', flags=c('overwrite'))

# Update working names of conservation criteria
conservationCriteriaNames <- conservationCriteriaList$Name 

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
  execGRASS('r.mapcalc', expression=paste0(species, '_habitatSuitabilityFocal_NatAreas_01=', species, '_habitatSuitabilityFocal_NatAreas/100.0'), flags=c('overwrite'))
}  
# Update working names of conservation criteria
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
conservationCriteriaNames[2:length(conservationCriteriaNames)] <- paste0(conservationCriteriaNames[2:length(conservationCriteriaNames)], "_Log01") 

# Rescale connectivity layers 0 - 1 without transforming them
connectivityCriteriaNames <- connectivityCriteriaNames[-length(connectivityCriteriaNames)]
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
write.csv(maxMinTable, file.path(b03zonationDir, "ZonationInputs", "MaxMinTable.csv"), row.names=F)

#Summed species maps
execGRASS('r.mapcalc',expression='all_habitatSuitabilityFocal_NatAreas=(if(isnull(MAAM_habitatSuitabilityFocal_NatAreas_01),0,MAAM_habitatSuitabilityFocal_NatAreas_01)+if(isnull(BLBR_habitatSuitabilityFocal_NatAreas_01),0,BLBR_habitatSuitabilityFocal_NatAreas_01)+if(isnull(URAM_habitatSuitabilityFocal_NatAreas_01),0,URAM_habitatSuitabilityFocal_NatAreas_01)+if(isnull(PLCI_habitatSuitabilityFocal_NatAreas_01),0,PLCI_habitatSuitabilityFocal_NatAreas_01)+if(isnull(RASY_habitatSuitabilityFocal_NatAreas_01),0,RASY_habitatSuitabilityFocal_NatAreas_01))*studyArea1', flags=c('overwrite'))
execGRASS('r.mapcalc',expression='all_currentFlowFocal_NatAreas=(if(isnull(MAAM_currentFlowFocal_NatAreas_Log01),0,MAAM_currentFlowFocal_NatAreas_Log01)+if(isnull(BLBR_currentFlowFocal_NatAreas_Log01),0,BLBR_currentFlowFocal_NatAreas_Log01)+if(isnull(URAM_currentFlowFocal_NatAreas_Log01),0,URAM_currentFlowFocal_NatAreas_Log01)+if(isnull(PLCI_currentFlowFocal_NatAreas_Log01),0,PLCI_currentFlowFocal_NatAreas_Log01)+if(isnull(RASY_currentFlowFocal_NatAreas_Log01),0,RASY_currentFlowFocal_NatAreas_Log01))*studyArea1', flags=c('overwrite'))
execGRASS('r.mapcalc',expression='all_betweennessLongFocal_NatAreas=(if(isnull(MAAM_betweennessLongFocal_NatAreas_01),0,MAAM_betweennessLongFocal_NatAreas_01)+if(isnull(BLBR_betweennessLongFocal_NatAreas_01),0,BLBR_betweennessLongFocal_NatAreas_01)+if(isnull(URAM_betweennessLongFocal_NatAreas_01),0,URAM_betweennessLongFocal_NatAreas_01)+if(isnull(PLCI_betweennessLongFocal_NatAreas_01),0,PLCI_betweennessLongFocal_NatAreas_01)+if(isnull(RASY_betweennessLongFocal_NatAreas_01),0,RASY_betweennessLongFocal_NatAreas_01))*studyArea1', flags=c('overwrite'))
execGRASS('r.mapcalc',expression='all_betweennessShortFocal_NatAreas=(if(isnull(MAAM_betweennessShortFocal_NatAreas_01),0,MAAM_betweennessShortFocal_NatAreas_01)+if(isnull(BLBR_betweennessShortFocal_NatAreas_01),0,BLBR_betweennessShortFocal_NatAreas_01)+if(isnull(URAM_betweennessShortFocal_NatAreas_01),0,URAM_betweennessShortFocal_NatAreas_01)+if(isnull(PLCI_betweennessShortFocal_NatAreas_01),0,PLCI_betweennessShortFocal_NatAreas_01)+if(isnull(RASY_betweennessShortFocal_NatAreas_01),0,RASY_betweennessShortFocal_NatAreas_01))*studyArea1', flags=c('overwrite'))
execGRASS('r.mapcalc',expression='all_dECGapFocal_NatAreas=(if(isnull(MAAM_dECGapFocal_NatAreas_01),0,MAAM_dECGapFocal_NatAreas_01)+if(isnull(BLBR_dECGapFocal_NatAreas_01),0,BLBR_dECGapFocal_NatAreas_01)+if(isnull(URAM_dECGapFocal_NatAreas_01),0,URAM_dECGapFocal_NatAreas_01)+if(isnull(PLCI_dECGapFocal_NatAreas_01),0,PLCI_dECGapFocal_NatAreas_01)+if(isnull(RASY_dECGapFocal_NatAreas_01),0,RASY_dECGapFocal_NatAreas_01))*studyArea1', flags=c('overwrite'))
execGRASS('r.mapcalc',expression='all_dECNatalFocal_NatAreas=(if(isnull(MAAM_dECNatalFocal_NatAreas_01),0,MAAM_dECNatalFocal_NatAreas_01)+if(isnull(BLBR_dECNatalFocal_NatAreas_01),0,BLBR_dECNatalFocal_NatAreas_01)+if(isnull(URAM_dECNatalFocal_NatAreas_01),0,URAM_dECNatalFocal_NatAreas_01)+if(isnull(PLCI_dECNatalFocal_NatAreas_01),0,PLCI_dECNatalFocal_NatAreas_01)+if(isnull(RASY_dECNatalFocal_NatAreas_01),0,RASY_dECNatalFocal_NatAreas_01))*studyArea1', flags=c('overwrite'))
execGRASS('r.mapcalc',expression='all_connectivityFocal_NatAreas=all_currentFlowFocal_NatAreas+all_betweennessShortFocal_NatAreas+all_dECGapFocal_NatAreas+all_dECNatalFocal_NatAreas', flags=c('overwrite'))
execGRASS('r.mapcalc',expression='all_connectivityFocal_NatAreas=all_currentFlowFocal_NatAreas+all_betweennessLongFocal_NatAreas+all_betweennessShortFocal_NatAreas+all_dECGapFocal_NatAreas+all_dECNatalFocal_NatAreas', flags=c('overwrite'))
execGRASS('r.mapcalc',expression='all_qualityconnectivityFocal_NatAreas=all_habitatSuitabilityFocal_NatAreas+all_currentFlowFocal_NatAreas+all_betweennessShortFocal_NatAreas+all_dECGapFocal_NatAreas+all_dECNatalFocal_NatAreas', flags=c('overwrite'))
execGRASS('r.mapcalc',expression='all_qualityconnectivityFocal_NatAreas=all_habitatSuitabilityFocal_NatAreas+all_currentFlowFocal_NatAreas+all_betweennessLongFocal_NatAreas+all_betweennessShortFocal_NatAreas+all_dECGapFocal_NatAreas+all_dECNatalFocal_NatAreas', flags=c('overwrite'))

#Save summed species maps
execGRASS('r.out.gdal',input='all_habitatSuitabilityFocal_NatAreas',output=file.path(b03habitatDir,paste0('ALL_habitatSuitability_NatAreas_01_', myResolution, 'm.tif')),format='GTiff',createopt='COMPRESS=LZW',flags=c('overwrite'))
execGRASS('r.out.gdal',input='all_currentFlowFocal_NatAreas',output=file.path(b03circuitscapeDir, paste0("ALL_currendensity_NatAreas_Log01_", myResolution, "m.tif")),format='GTiff',createopt='COMPRESS=LZW',flags=c('overwrite'))
execGRASS('r.out.gdal',input='all_betweennessLongFocal_NatAreas',output=file.path(b03networkDir,paste0('ALL_betweennessLong_01_', myResolution,'m.tif')),format='GTiff',createopt='COMPRESS=LZW',flags=c('overwrite'))
execGRASS('r.out.gdal',input='all_betweennessShortFocal_NatAreas',output=file.path(b03networkDir,paste0('ALL_betweennessShort_01_', myResolution, 'm.tif')),format='GTiff',createopt='COMPRESS=LZW',flags=c('overwrite'))
execGRASS('r.out.gdal',input='all_dECGapFocal_NatAreas',output=file.path(b03patchImportanceDir,'ALL_dECGap_01_240m.tif'),format='GTiff',createopt='COMPRESS=LZW',flags=c('overwrite'))
execGRASS('r.out.gdal',input='all_dECNatalFocal_NatAreas',output=file.path(b03patchImportanceDir,'ALL_dECNatal_01_240m.tif'),format='GTiff',createopt='COMPRESS=LZW',flags=c('overwrite'))
execGRASS('r.out.gdal',input='all_connectivityFocal_NatAreas',output=file.path(b03networkDir,'ALL_connectivity.tif'),format='GTiff',createopt='COMPRESS=LZW',flags=c('overwrite'))
execGRASS('r.out.gdal',input='all_qualityconnectivityFocal_NatAreas',output=file.path(b03networkDir,'ALL_qualityConnectivity.tif'),format='GTiff',createopt='COMPRESS=LZW',flags=c('overwrite'))


# Set geographic region for Zonation analysis (add 1 row of NA pixels all around each map)
execGRASS('g.region', n=paste0('n+',myResolution), e=paste0('e+',myResolution), w=paste0('w-',myResolution), s=paste0('s-',myResolution), flags='c')

# Write out habitat suitability and log-transformed connectivity layers to geotiffs
for(i in 1:length(speciesList)){
  species<-speciesList[i]
  for(j in 1:nrow(conservationCriteriaList)){
    execGRASS('r.patch', input=paste0(species, '_', conservationCriteriaNames[j], ',Z_mask0_BTSL'),output=paste0('Z_', species, '_', conservationCriteriaNames[j]),flags=c('overwrite'))
    execGRASS('r.out.gdal',input=paste0('Z_', species, '_', conservationCriteriaNames[j]),output=file.path(b03zonationDir, "ZonationInputs", paste0(species, '_', conservationCriteriaNames[j], '_Z.tif')),format='GTiff',flags=c('overwrite'))
  }
}

# Write out non-log-transformed connectivity layers to geotiffs
for(i in 1:length(speciesList)){
  species<-speciesList[i]
  for(j in 1:length(connectivityCriteriaNames)){
    execGRASS('r.patch', input=paste0(species, '_', connectivityCriteriaNames[j], '_01,Z_mask0_BTSL'),output=paste0('Z_', species, '_', connectivityCriteriaNames[j], '_01'),flags=c('overwrite'))
    execGRASS('r.out.gdal',input=paste0('Z_', species, '_', connectivityCriteriaNames[j], '_01'),output=file.path(b03zonationDir, "ZonationInputs", paste0(species, '_', connectivityCriteriaNames[j], '_01_Z.tif')),format='GTiff',flags=c('overwrite'))
  }
}

# Change the extent of the protected areas raster to match the geographic regions for Zonation analysis
execGRASS('r.mapcalc', expression='Z_protectedAreas=protectedAreas', flags=c('overwrite'))
execGRASS('r.out.gdal',input='Z_protectedAreas',output=file.path(b03zonationDir, "ZonationInputs", 'Z_protectedAreas.tif'),format='GTiff',flags=c('overwrite'))

# Reset the geographic region
execGRASS('g.region', raster='StudyArea1', res=paste0(myResolution))
