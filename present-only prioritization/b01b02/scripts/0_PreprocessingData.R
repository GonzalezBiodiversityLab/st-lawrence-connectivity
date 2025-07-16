library(rgrass7)

Sys.setenv(TZ='GMT')
options(stringsAsFactors=FALSE)

# Input parameters
# Resolution in meters
# Run this once for 10m and once for 30m
myResolution <- 30
# Set up GRASS mapset for the first time
doGRASSSetup <- F

# Directories
projectDir <- "E:/"
gisBase <- "C:/Program Files/GRASS GIS 7.4.3"
gisDbase <- paste0(projectDir, "grass7")
rawMapsDir <- paste0(projectDir, "Inputs/RawData/Maps/")
rawTablesDir <- paste0(projectDir, "Inputs/RawData/Tables/")
processedMapsDir <- paste0(projectDir, "Inputs/ProcessedData/Maps/")

# Raw data filenames
landcoverBTSLPolyName <- "PASL_Merge_QL"
landcoverBufferName <- "Occ_sol_2014_recl_FED_10m_aout2017.tif"
landcoverQuebecName <- "utilisation_territoire.tif"
landcoverBTSLName <- "Ras_FED_aout2017_calss_det.tif"
roadName <- "BDTQ_Roads_10m.tif"
ageName <- "Age_FullExtent.tif"
densityName <- "Densite_FullExtent.tif"
depositName <- "RCL_DEPOT.tif"
drainageName <- "RCL_DRAINAGE.tif"
studyAreaName <- "BTSL_StudyArea_10m.tif"
#waterName <- "Water_StLaurent.tif"

###############
# GRASS setup #
###############
if(doGRASSSetup){
  #https://gis.stackexchange.com/questions/183032/create-a-new-grass-database-in-r-with-crs-projection
  # Manually set up empty GRASS database - see GRASSTemplate
  initGRASS(gisBase=gisBase, gisDbase=gisDbase, location=paste0('BTSL_', myResolution, 'm'), mapset='PERMANENT', override=TRUE)
  
  execGRASS("g.proj", georef=paste0(rawMapsDir, landcoverBTSLName), flags="c")
  
  # Initialize new mapset inheriting projection info
  execGRASS("g.mapset", mapset = "RawData", flags="c")

  # Load layers into grass database
  execGRASS("v.in.ogr", input=rawMapsDir, layer=landcoverBTSLPolyName, output="rawDataLandcoverBTSLPoly", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=paste0(rawMapsDir, landcoverBufferName), output="rawDataLandcoverBuffer", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=paste0(rawMapsDir, landcoverQuebecName), output="rawDataLandcoverQuebec", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=paste0(rawMapsDir, landcoverBTSLName), output="rawDataLandcoverBTSL", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=paste0(rawMapsDir, ageName), output="rawDataForestAge", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=paste0(rawMapsDir, densityName), output="rawDataForestDensity", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=paste0(rawMapsDir, depositName), output="rawDataDeposit", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=paste0(rawMapsDir, drainageName), output="rawDataDrainage", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=paste0(rawMapsDir, studyAreaName), output="rawDataStudyArea", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=paste0(rawMapsDir, roadName), output="rawDataRoad", flags=c("overwrite", "o"))
  #execGRASS("r.in.gdal", input=paste0(rawMapsDir, waterName), output="rawDataWaterStLaurent", flags=c("overwrite", "o"))
}else{
  initGRASS(gisBase=gisBase, gisDbase=gisDbase, location=paste0('BTSL_', myResolution, 'm'), mapset="RawData", override=TRUE)
}

# Set the geographic region
#execGRASS('g.region', n='403680', e='-95520', w='-491430', s='124830')
execGRASS('g.region', n='403680', e='-95520', w='-491430', s='117960', res=paste0(myResolution))


execGRASS('r.resamp.stats', input="landcoverBTSL1", output=paste0("landcoverBTSL1_", myResolution, "m"), method="mode", flags=c("overwrite"))
execGRASS('r.resamp.stats', input="landcoverBTSL1", output=paste0("landcoverBTSL1_", myResolution, "m"), method="mode", flags=c("overwrite"))
execGRASS('r.resamp.stats', input="landcoverBTSL1", output=paste0("landcoverBTSL1_", myResolution, "m"), method="mode", flags=c("overwrite"))

execGRASS('r.out.gdal', input='landcover',output=paste0(processedMapsDir, 'landcover_', myResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))
execGRASS('r.out.gdal', input='forestAge',output=paste0(processedMapsDir, 'forestAge_', myResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))
execGRASS('r.out.gdal', input='rawDataForestDensity', output=paste0(processedMapsDir, 'forestDensity_', myResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW',flags=c('overwrite'))
execGRASS('r.out.gdal', input='rawDataDrainage', output=paste0(processedMapsDir, 'drainage_', myResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))
execGRASS('r.out.gdal', input='deposit', output=paste0(processedMapsDir, 'deposit_', myResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))
execGRASS('r.out.gdal', input='studyArea1', output=paste0(processedMapsDir, 'studyArea_', myResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))

###############################################
# Study Area, Age, Density, Drainage, Deposit #
###############################################
# NB No need to reclassify forest density or drainage
# Study area
# Fix NA values in study area layer (255=NULL)
write.table(c('1=1','255=NULL'), 'rule.txt', sep="", col.names=FALSE, quote=FALSE, row.names=FALSE)
execGRASS('r.reclass', input='rawDataStudyArea', output='studyArea1', rules='rule.txt', flags=c('overwrite'))

# Forest age - reclassify
rcl<-read.csv(paste0(rawTablesDir, "forestAgeReclass.csv"), header=TRUE)[,c('Value','Code')]
write.table(paste0(rcl[,'Value'], '=', rcl[,'Code']), 'rule.txt', sep="", col.names=FALSE, quote=FALSE, row.names=FALSE)
execGRASS('r.reclass', input='rawDataForestAge', output='forestAge', rules='rule.txt', flags=c('overwrite'))

# Surficial deposit - reclassify
rcl<-read.csv(paste0(rawTablesDir, "depositReclass.csv"),header=TRUE)[,c('Value','Recode')]
write.table(paste0(rcl[,'Value'], '=', rcl[,'Recode']), 'rule.txt', sep="", col.names=FALSE, quote=FALSE, row.names=FALSE)
execGRASS('r.reclass', input='rawDataDeposit', output='deposit', rules='rule.txt', flags=c('overwrite'))


##################
# Landcover BTSL #
##################
# Reclass tables
speciesLandcoverReclass <- read.csv(paste0(rawTablesDir, "speciesLandcoverReclass.csv"), header=TRUE)
landcoverReclass <- read.csv(paste0(rawTablesDir, "landcoverBTSLReclass.csv"), header=TRUE)
roadReclassBDTQ <- read.csv(paste0(rawTablesDir, "roadReclassBDTQ.csv"), header=TRUE)

# Reclassify landcover in BTSL to match Albert et al. classes
write.table(paste0(landcoverReclass[,'Value'],'=', landcoverReclass[,'Code']), 'rule.txt', sep="", col.names=FALSE, quote=FALSE, row.names=FALSE)
execGRASS('r.reclass', input='rawDataLandcoverBTSL', output='landcoverBTSL1', rules='rule.txt', flags=c('overwrite'))
execGRASS('r.resamp.stats', input="landcoverBTSL1", output=paste0("landcoverBTSL1_", myResolution, "m"), method="mode", flags=c("overwrite"))

# Fix broken linear features by extracting them, rasterizing them at higher resolution, and then reimposing them on the landcover raster
# Extract linear agriculture (i.e. "Milieu agricole non cultivé")
# Note that characters with accents are not retained in the attribute table therefore use "Milieu agricole non cultivÃ©"
agLinearCode <- speciesLandcoverReclass$LandcoverCode[speciesLandcoverReclass$LandcoverName == "AgricultureLinearElements"]
agLinearDetailedClasses <- landcoverReclass$SQL_CLASSE_DET[landcoverReclass$Code == agLinearCode]
whereString = paste0("CLASSE_DET IN ('",paste(agLinearDetailedClasses, collapse="','"),"')")
execGRASS("v.extract", input="rawDataLandcoverBTSLPoly", layer='1', type="area", where=whereString, output="landcoverAgLinearPoly", flags=c("overwrite"))
execGRASS('g.region', res='2')
execGRASS("v.to.rast", input="landcoverAgLinearPoly", use="val", value=agLinearCode, output="landcoverAgLinearRas_2m", flags=c("overwrite"))
execGRASS('g.region', res=paste0(myResolution))
execGRASS('r.resamp.stats', input="landcoverAgLinearRas_2m", output=paste0("landcoverAgLinearRas_", myResolution, "m"), method="maximum", flags=c("overwrite"))
# Impose linear agriculture on landcover map
execGRASS('r.patch', input=paste0("landcoverAgLinearRas_", myResolution, "m,landcoverBTSL1_", myResolution, "m"), output=paste0("landcoverBTSL2_", myResolution, "m"), flags=c("overwrite"))

# Extract roads from landcover map
minorRoadCode <- speciesLandcoverReclass$LandcoverCode[speciesLandcoverReclass$LandcoverName == "MinorRoads"]
majorRoadCode <- speciesLandcoverReclass$LandcoverCode[speciesLandcoverReclass$LandcoverName == "MajorRoads"]
# Minor
roadMinorDetailedClasses <- landcoverReclass$SQL_CLASSE_DET[landcoverReclass$Code == minorRoadCode]
whereString = paste0("CLASSE_DET IN ('",paste(roadMinorDetailedClasses, collapse="','"),"')")
execGRASS("v.extract", input="rawDataLandcoverBTSLPoly", layer='1', type="area", where=whereString, output="landcoverMinorRoadPASLPoly", flags=c("overwrite"))
execGRASS('g.region', res='2')
execGRASS("v.to.rast", input="landcoverMinorRoadPASLPoly", use="val", value=minorRoadCode, output="landcoverMinorRoadPASLRas_2m", flags=c("overwrite"))
execGRASS('g.region', res=paste0(myResolution))
execGRASS('r.resamp.stats', input="landcoverMinorRoadPASLRas_2m", output=paste0("landcoverMinorRoadPASLRas_", myResolution, "m"), method="maximum", flags=c("overwrite"))
# Major
roadMajorDetailedClasses <- landcoverReclass$SQL_CLASSE_DET[landcoverReclass$Code == majorRoadCode]
whereString = paste0("CLASSE_DET IN ('",paste(roadMajorDetailedClasses, collapse="','"),"')")
execGRASS("v.extract", input="rawDataLandcoverBTSLPoly", layer='1', type="area", where=whereString, output="landcoverMajorRoadPASLPoly", flags=c("overwrite"))
execGRASS('g.region', res='2')
execGRASS("v.to.rast", input="landcoverMajorRoadPASLPoly", use="val", value=majorRoadCode, output="landcoverMajorRoadPASLRas_2m", flags=c("overwrite"))
execGRASS('g.region', res=paste0(myResolution))
execGRASS('r.resamp.stats', input="landcoverMajorRoadPASLRas_2m", output=paste0("landcoverMajorRoadPASLRas_", myResolution, "m"), method="maximum", flags=c("overwrite"))
# All roads PASL
execGRASS('r.mapcalc', expression=paste0("landcoverRoadPASL1 = if(isnull(landcoverMajorRoadPASLRas_", myResolution, "m),0,landcoverMajorRoadPASLRas_", myResolution, "m) + if(isnull(landcoverMinorRoadPASLRas_", myResolution, "m),0,landcoverMinorRoadPASLRas_", myResolution, "m)"), flags=c("overwrite"))
write.table(c('0=0', paste0(minorRoadCode, '=', minorRoadCode), paste0(majorRoadCode, '=', majorRoadCode), paste0(minorRoadCode + majorRoadCode, '=', majorRoadCode)), 'rule.txt', sep="", col.names=FALSE, quote=FALSE, row.names=FALSE)
execGRASS('r.reclass', input='landcoverRoadPASL1', output='landcoverRoadPASL', rules='rule.txt', flags=c('overwrite'))

# Add roads from BDTQ database
# Rescale raw BDTQ road data to myResolution
execGRASS('r.resamp.stats', input="rawDataRoad", output=paste0("rawDataRoad_", myResolution, "m"), method="maximum", flags=c("overwrite"))
rcl<-read.csv(paste0(rawTablesDir, "roadReclassBDTQ.csv"), header=TRUE)[,c('Value','Code')]
write.table(paste0(rcl[,'Value'], '=', rcl[,'Code']), 'rule.txt', sep="", col.names=FALSE, quote=FALSE, row.names=FALSE)
execGRASS('r.reclass', input=paste0("rawDataRoad_", myResolution, "m"), output='landcoverRoadBDTQ', rules='rule.txt', flags=c('overwrite'))

# Combine all roads into a sinlge layer
# NB: assume the PASL roads are correct and then add the BDTQ roads to them
execGRASS('r.mapcalc', expression='landcoverRoad = if(landcoverRoadPASL>0,landcoverRoadPASL,landcoverRoadBDTQ)', flags=c("overwrite"))
execGRASS('r.null', map='landcoverRoad', setnull='0')

# Impose roads on BTSL landcover map
execGRASS('r.patch', input=paste0("landcoverRoad,landcoverBTSL2_", myResolution, "m"), output='landcoverBTSL', flags=c('overwrite'))

####################
# Landcover Buffer #
####################
# Reclassify landcover in Quebec to match buffer classes
rcl<-read.csv(paste0(rawTablesDir, "landcoverQuebecReclass.csv"),header=TRUE)[,c('Code','Recode')]
write.table(paste0(rcl[,'Code'],'=',rcl[,'Recode']),'rule.txt',sep="",col.names=FALSE,quote=FALSE,row.names=FALSE)
execGRASS('r.reclass',input='rawDataLandcoverQuebec',output='landcoverQuebecClasseGen' ,rules='rule.txt',flags=c('overwrite'))

# Fill out buffer with Quebec landcover
execGRASS('r.patch', input="rawDataLandcoverBuffer,landcoverQuebecClasseGen", output='landcoverBufferFilled', flags=c('overwrite'))

# Reclassify landcover in buffer to match Albert et al. classes
rcl<-read.csv(paste0(rawTablesDir, "landcoverBufferReclass.csv"),header=TRUE)[,c('Value','Code')]
write.table(paste0(rcl[,'Value'],'=',rcl[,'Code']),'rule.txt',sep="",col.names=FALSE,quote=FALSE,row.names=FALSE)
execGRASS('r.reclass',input='landcoverBufferFilled',output='landcoverBuffer' ,rules='rule.txt',flags=c('overwrite'))

#####################################
# Combine BTSL and Buffer landcover #
#####################################
execGRASS('r.patch', input='landcoverBTSL,landcoverBuffer', output='landcover', flags=c('overwrite'))

############################
# Save outputs to geotiffs #
############################
execGRASS('r.out.gdal', input='landcover',output=paste0(processedMapsDir, 'landcover_', myResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))
execGRASS('r.out.gdal', input='forestAge',output=paste0(processedMapsDir, 'forestAge_', myResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))
execGRASS('r.out.gdal', input='rawDataForestDensity', output=paste0(processedMapsDir, 'forestDensity_', myResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW',flags=c('overwrite'))
execGRASS('r.out.gdal', input='rawDataDrainage', output=paste0(processedMapsDir, 'drainage_', myResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))
execGRASS('r.out.gdal', input='deposit', output=paste0(processedMapsDir, 'deposit_', myResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))
execGRASS('r.out.gdal', input='studyArea1', output=paste0(processedMapsDir, 'studyArea_', myResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))
