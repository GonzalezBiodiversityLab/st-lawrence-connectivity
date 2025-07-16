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
conservationCriteriaList<-data.frame(Name=c("habitatSuitability", "betweenness", "betweennessBTSL", "dECGap", "dECNatal", "currentFlow", "currentFlowBTSL"),
                                    FileName=c(paste0("habitatSuitability_", myResolution, "m.tif"), "betweenness_01.tif", "betweenness_BTSL_01.tif", "dECGap_BTSL_01.tif", "dECNatal_BTSL_01.tif", "30_FULL_curmap_01.tif", "30_FULL_curmap_BTSL_01.tif"),
                                    Directory=c(habitatDir, networkDir, networkDir, networkDir, networkDir, circuitDir, circuitDir))
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
  execGRASS("g.mapset", mapset="HabitatNetworks", flags="c")
  
  # Load layers into grass database
  for(i in 1:length(speciesList)){
    species<-speciesList[i]
    #Habitat patches
    execGRASS("r.in.gdal", input=paste0(habitatDir, species, '_habitatPatch_30m.tif'), output=paste0(species, "_habitatPatch"), flags=c("overwrite", "o"))
    #Resistance
    execGRASS("r.in.gdal", input=paste0(resistanceDir, species, '_resistance_30m.tif'), output=paste0(species, "_resistance"), flags=c("overwrite", "o"))
    #Links
    execGRASS("r.in.gdal", input=paste0(networkDir, species, '_linkWeight.asc'), output=paste0(species, "_linkWeight"), flags=c("overwrite", "o"))
    #Connectivity analysis results
    for(j in 1:nrow(conservationCriteriaList)){
      execGRASS("r.in.gdal", input=paste0(conservationCriteriaList$Directory[j], species, "_", conservationCriteriaList$FileName[j]), output=paste0(species, "_", conservationCriteriaList$Name[j]), flags=c("overwrite", "o"))
    }
  }
  execGRASS("r.in.gdal", input=paste0(processedMapsDir, landcoverName), output="landcover", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=paste0(processedMapsDir, studyAreaName), output="studyArea1", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=paste0(processedMapsDir, protectedAreasName), output="protectedAreas0", flags=c("overwrite", "o"))
}else{
  initGRASS(gisBase=gisBase, gisDbase=gisDbase, location=paste0("BTSL_", myResolution, "m"), mapset='HabitatNetworks', override=TRUE)
}

# Set the geographic region
#execGRASS('g.region', n='403680', e='-95520', w='-491430', s='124830')
execGRASS('g.region', n='403680', e='-95520', w='-491430', s='117960', res=paste0(myResolution))

for(i in 1:length(speciesList)){
  species<-speciesList[i]
#  execGRASS('r.mapcalc',expression=paste0(species, '_habitat=if(', species, '_habitatSuitability>=60,1,0)'), flags=c('overwrite'))
#  execGRASS('r.mapcalc',expression=paste0(species, '_links=if(isnull(', species, '_linkWeight),0,1)'), flags=c('overwrite'))
#  execGRASS('r.mapcalc',expression=paste0(species, '_betweennessClipped=if(isnull(', species, '_betweenness),0,',species,'_betweenness)*studyArea1'), flags=c('overwrite'))
#  execGRASS('r.mapcalc',expression=paste0(species, '_betweennessBTSL01=if(isnull(', species, '_betweennessBTSL),0,',species,'_betweennessBTSL)*studyArea1'), flags=c('overwrite'))
  execGRASS('r.mapcalc',expression=paste0(species, '_dECGap01=if(isnull(', species, '_dECGap),0,',species,'_dECGap)*studyArea1'), flags=c('overwrite'))
  execGRASS('r.mapcalc',expression=paste0(species, '_dECNatal01=if(isnull(', species, '_dECNatal),0,',species,'_dECNatal)*studyArea1'), flags=c('overwrite'))

#  execGRASS('r.mapcalc',expression=paste0(species, '_habitatSuitabilityBTSL=', species, '_habitatSuitability*studyArea1'), flags=c('overwrite'))
#  execGRASS('r.mapcalc',expression=paste0(species, '_habitatPatchBTSL=', species, '_habitatPatch*studyArea1'), flags=c('overwrite'))
  execGRASS('r.mapcalc',expression=paste0(species, '_resistanceBTSL=', species, '_resistance*studyArea1'), flags=c('overwrite'))
  
  # #Code to rescale betweenness maps to 0-1 after clipping them to BTSL *** NOT WORKING ***
  # execGRASS('r.mapcalc',expression=paste0(species, '_betweennessClipped=', species, '_betweenness*studyArea1*1000000'), flags=c('overwrite'))
  # statistic <- execGRASS("r.univar", map = paste0(species, '_betweennessClipped'), flags = c("g", "quiet"), intern = TRUE, legacyExec = TRUE)
  # mapMax <- floor(as.numeric(strsplit(statistic[agrep(statistic, pattern = "max")], "=")[[1]][[2]]))
  # mapMin <- as.numeric(strsplit(statistic[agrep(statistic, pattern = "min")], "=")[[1]][[2]])
  # execGRASS('r.rescale',input=paste0(species, '_betweennessClipped'), from=c(mapMin, mapMax), output=paste0(species, '_betweennessClipped01'), to=c(0,1000000), flags=c('overwrite'))
}

#Summed species maps
execGRASS('r.mapcalc',expression='all_habitatSum=MAAM_habitat+BLBR_habitat+URAM_habitat+PLCI_habitat+RASY_habitat', flags=c('overwrite'))
execGRASS('r.mapcalc',expression='all_habitatPatch=MAAM_habitatPatch+BLBR_habitatPatch+URAM_habitatPatch+PLCI_habitatPatch+RASY_habitatPatch', flags=c('overwrite'))
execGRASS('r.mapcalc',expression='all_habitatBinary=if(all_habitat>0,1,0)', flags=c('overwrite'))
execGRASS('r.mapcalc',expression='all_habitatSuitability=MAAM_habitatSuitability*MAAM_habitat+BLBR_habitatSuitability*BLBR_habitat+URAM_habitatSuitability*URAM_habitat+PLCI_habitatSuitability*PLCI_habitat+RASY_habitatSuitability*RASY_habitat', flags=c('overwrite'))
execGRASS('r.stats',input='all_habitatSuitability', flags=c('c'))
execGRASS('r.mapcalc',expression='all_links=MAAM_links+BLBR_links+URAM_links+PLCI_links+RASY_links', flags=c('overwrite'))
execGRASS('r.mapcalc',expression='all_currentFlowBTSL=MAAM_currentFlowBTSL+BLBR_currentFlowBTSL+URAM_currentFlowBTSL+PLCI_currentFlowBTSL+RASY_currentFlowBTSL', flags=c('overwrite'))
execGRASS('r.mapcalc',expression='all_betweenness=MAAM_betweennessClipped+BLBR_betweennessClipped+URAM_betweennessClipped+PLCI_betweennessClipped+RASY_betweennessClipped', flags=c('overwrite'))
#execGRASS('r.mapcalc',expression='all_betweenness01=MAAM_betweennessClipped01+BLBR_betweennessClipped01+URAM_betweennessClipped01+PLCI_betweennessClipped01+RASY_betweennessClipped01', flags=c('overwrite'))
execGRASS('r.mapcalc',expression='all_betweennessBTSL=MAAM_betweennessBTSL01+BLBR_betweennessBTSL01+URAM_betweennessBTSL01+PLCI_betweennessBTSL01+RASY_betweennessBTSL01', flags=c('overwrite'))
execGRASS('r.mapcalc',expression='all_dECGap=BLBR_dECGap01+MAAM_dECGap01+PLCI_dECGap01+RASY_dECGap01+URAM_dECGap01', flags=c('overwrite')) #+BLBR_dECGap01+PLCI_dECGap01+RASY_dECGap01
execGRASS('r.mapcalc',expression='all_dECNatal=BLBR_dECNatal01+MAAM_dECNatal01+PLCI_dECNatal01+RASY_dECNatal01+URAM_dECNatal01', flags=c('overwrite')) #+BLBR_dECNatal01+PLCI_dECNatal01+RASY_dECNatal01

execGRASS('r.mapcalc',expression='landcoverBTSL=landcover*studyArea1', flags=c('overwrite'))

#Sum all connectivity maps across all species
execGRASS('r.mapcalc',expression='all_connectivity=all_betweenness+all_betweennessBTSL+all_dECGap+all_dECNatal+all_currentFlowBTSL', flags=c('overwrite'))
execGRASS('r.mapcalc',expression='all_habitat=MAAM_habitat_binary_30m+BLBR_habitat_binary_30m+URAM_habitat_binary_30m+PLCI_habitat_binary_30m+RASY_habitat_binary_30m', flags=c('overwrite'))
write.table(c('0=NULL','1=1','2=1','3=1','4=1','5=1'),'rule.txt',sep="",col.names=FALSE,quote=FALSE,row.names=FALSE)
execGRASS('r.reclass',input='all_habitat',output='all_habitat_NA',rules='rule.txt',flags=c('overwrite'))


#Save to geotiffs
execGRASS('r.out.gdal',input='all_habitatPatch',output=paste0(networkDir,'all_habitatPatch.tif'),format='GTiff',createopt='COMPRESS=LZW',flags=c('overwrite'))

execGRASS('r.out.gdal',input='all_habitat_suitability_30m',output=paste0(conOutputDir,'all_habitat_suitability_30m.tif'),format='GTiff',createopt='COMPRESS=LZW',flags=c('overwrite'))
execGRASS('r.out.gdal',input='all_habitat_suitability',output=paste0(conOutputDir,'all_habitat_suitability.tif'),format='GTiff',createopt='COMPRESS=LZW',flags=c('overwrite'))
execGRASS('r.out.gdal',input='all_betweenness',output=paste0(conOutputDir,'all_betweenness.tif'),format='GTiff',createopt='COMPRESS=LZW',flags=c('overwrite'))
execGRASS('r.out.gdal',input='all_dECGap',output=paste0(conOutputDir,'all_dECGap.tif'),format='GTiff',createopt='COMPRESS=LZW',flags=c('overwrite'))
execGRASS('r.out.gdal',input='all_dECNatal',output=paste0(conOutputDir,'all_dECNatal.tif'),format='GTiff',createopt='COMPRESS=LZW',flags=c('overwrite'))
execGRASS('r.out.gdal',input='all_connectivity_all_species',output=paste0(conOutputDir,'all_connectivity_all_species.tif'),format='GTiff',createopt='COMPRESS=LZW',flags=c('overwrite'))




execGRASS("r.mapcalc", expression='studyArea0=if(studyArea1>0,0)', flags=c('overwrite'))

#Protected Areas
write.table(c('1=1','0=NULL'), 'rule.txt', sep="", col.names=FALSE, quote=FALSE, row.names=FALSE)
execGRASS('r.reclass', input='protectedAreas01', output='protectedAreasNA1', rules='rule.txt', flags=c('overwrite'))
execGRASS('r.clump', input='protectedAreasNA1', output='protectedAreasId', flags=c('overwrite', 'd'))
execGRASS('r.reclass.area', input='protectedAreasId', output='protectedAreasId150', value=150, mode='greater', flags=c('overwrite', 'd'))
execGRASS('r.to.vect.area', input='protectedAreasId', output='protectedAreasId_vect', flags=c('overwrite'))
execGRASS('r.to.vect.area', input='protectedAreasId150', output='protectedAreasId150_vect', flags=c('overwrite'))




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

# Set geographic region for Zonation analysis (add 1 row of NA pixels all around each map)
execGRASS('g.region', n=paste0('n+',myResolution), e=paste0('e+',myResolution), w=paste0('w-',myResolution), s=paste0('s-',myResolution), flags='c')

# Make all areas outside of natural areas NULL and write out geotiffs
for(i in 1:length(speciesList)){
  species<-speciesList[i]
  for(j in 1:nrow(conservationCriteriaList)){
    execGRASS('r.mapcalc',expression=paste0('Z_',species, '_', conservationCriteriaList$Name[j], '=Z_mask*studyArea1*', species, '_', conservationCriteriaList$Name[j]), flags=c('overwrite'))
    execGRASS('r.out.gdal',input=paste0('Z_',species, '_', conservationCriteriaList$Name[j]),output=paste0(zonInputDir,species, '_', conservationCriteriaList$Name[j], '_Z.tif'),format='GTiff',flags=c('overwrite'))
  }
}

# Change the extent of the protected areas raster to match the geographic regions for Zonation analysis
execGRASS('r.mapcalc', expression='Z_protectedAreas=protectedAreas', flags=c('overwrite'))
execGRASS('r.out.gdal',input='Z_protectedAreas',output=paste0(zonInputDir, 'Z_protectedAreas.tif'),format='GTiff',flags=c('overwrite'))
