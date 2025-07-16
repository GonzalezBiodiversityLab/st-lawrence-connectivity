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
linkageDir <- paste0(projectDir, "Outputs/6.LinkagePrioritization/LinkagePriority_30m/")

# Input parameters
# Run this at 30m resolution
myResolution <- 30
speciesList<-c("MAAM", "PLCI", "RASY", "BLBR", "URAM")
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
  execGRASS("g.mapset", mapset="LinkagePrioritization", flags="c")
  
  # Load layers into grass database
  for(i in 1:length(speciesList)){
    species<-speciesList[i]
    execGRASS("r.in.gdal", input=paste0(linkageDir, 'Corridors_30m_', species, '_linkage_priority.tif'), output=paste0(species, "_linkages"), flags=c("overwrite", "o"))
  }
  execGRASS("r.in.gdal", input=paste0(processedMapsDir, landcoverName), output="landcover", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=paste0(processedMapsDir, studyAreaName), output="studyArea1", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=paste0(processedMapsDir, protectedAreasName), output="protectedAreas", flags=c("overwrite", "o"))
}else{
  initGRASS(gisBase=gisBase, gisDbase=gisDbase, location=paste0("BTSL_", myResolution, "m"), mapset='LinkagePrioritization', override=TRUE)
}

# Set the geographic region
execGRASS('g.region', n='403680', e='-95520', w='-491430', s='117960', res=paste0(myResolution))

# Massage study area and protected areas layers
execGRASS("r.mapcalc", expression='studyArea0=if(studyArea1>0,0)', flags=c('overwrite'))
execGRASS('r.reclass.area', input='protectedAreas', output='protectedAreas150', value=150, mode='greater', flags=c('overwrite', 'd'))

#Summed species maps
execGRASS('r.mapcalc',expression='all_linkages=if(isnull(MAAM_linkages),0,MAAM_linkages)+if(isnull(BLBR_linkages),0,BLBR_linkages)+if(isnull(URAM_linkages),0,URAM_linkages)+if(isnull(PLCI_linkages),0,PLCI_linkages)+if(isnull(RASY_linkages),0,RASY_linkages)', flags=c('overwrite'))
execGRASS('r.out.gdal',input='all_linkages',output=paste0(linkageDir,'Corridors_30m_AllSpecies_linkage_priority.tif'),format='GTiff',createopt='COMPRESS=LZW',flags=c('overwrite'))
execGRASS('r.out.gdal',input='all_linkages',output=paste0(linkageDir,'all_linkages.tif'),format='GTiff',createopt='COMPRESS=LZW',flags=c('overwrite'))
