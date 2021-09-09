# a254
# Bronwyn Rayfield and Jed Lloren, ApexRMS
# Run with R-4.1.1                                
# N.B. Spatial layers were originally received as geodatabase and manually converted to shapefiles before reading into R

# JL
# Workspace---------------------------------------------------------------------------------------------
# Set environment variable TZ when running on AWS EC2 instance
Sys.setenv(TZ='GMT')
options(stringsAsFactors=FALSE)

# Load packages
library(rgrass7)

# Input parameters
# Resolution in meters
# Run this once for 10m and once for 30m
myResolution <- 30
# Set up GRASS mapset for the first time
doGRASSSetup <- F

# JL
# Set up directories
# Assumes these directories already exist
gisBase <- "C:/Program Files/GRASS GIS 7.8"
gisDbase <- "../grass7"
rawMapsDir <- "../Data/Spatial"
rawTablesDir <- "../Inputs/RawData/Tables"
processedMapsDir <- "../Inputs/ProcessedData/Maps"

# Raw data filenames
studyareaName <- "OutaouaisConnectivityExtent.tif"
protectedareaName <- "AP_REG_S_20210824.shp"
landcoverBTSLPolyName <-  "BTSL_SLL_Occ_sol_Land_cover.shp"
roadName <- "AQ_Routes_l_20210701.shp"
siefName <- "SIEF_C08PEEFO.shp"
#landcoverBufferName <- "utilisation_territoire_2016"

###############
# GRASS setup #
###############
if(doGRASSSetup){
  #https://gis.stackexchange.com/questions/183032/create-a-new-grass-database-in-r-with-crs-projection
  # Manually set up empty GRASS database - see GRASSTemplate
  initGRASS(gisBase=gisBase, gisDbase=gisDbase, location=paste0('BTSL_', myResolution, 'm'), mapset='PERMANENT', override=TRUE)
  
  execGRASS("g.proj", georef=file.path(rawMapsDir, landcoverBTSLPolyName), flags="c")
  
  # Initialize new mapset inheriting projection info
  execGRASS("g.mapset", mapset = "RawData", flags="c")

  # Load layers into grass database
  execGRASS("v.in.ogr", input=rawMapsDir, layer=siefName, output="rawDataSiefName", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=file.path(rawMapsDir, landcoverBTSLPolyName), output="rawDataLandcoverBTSL", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=paste0(rawMapsDir, ageName), output="rawDataForestAge", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=paste0(rawMapsDir, densityName), output="rawDataForestDensity", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=paste0(rawMapsDir, depositName), output="rawDataDeposit", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=paste0(rawMapsDir, drainageName), output="rawDataDrainage", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=paste0(rawMapsDir, studyAreaName), output="rawDataStudyArea", flags=c("overwrite", "o"))
  execGRASS("v.in.ogr", input=rawMapsDir, layer=roadName, output="rawDataRoad", flags=c("overwrite", "o"))
  #execGRASS("r.in.gdal", input=paste0(rawMapsDir, waterName), output="rawDataWaterStLaurent", flags=c("overwrite", "o"))
}else{
  initGRASS(gisBase=gisBase, gisDbase=gisDbase, location=paste0('BTSL_', myResolution, 'm'), mapset="RawData", override=TRUE)
}

# Set the geographic region
#execGRASS('g.region', n='403680', e='-95520', w='-491430', s='124830')
execGRASS('g.region', n='403680', e='-95520', w='-491430', s='117960')


###############################################
# Study Area, Age, Density, Drainage, Deposit #
###############################################
# NB No need to reclassify forest density or drainage
# Study area
# Fix NA values in study area layer (255=NULL)
write.table(c('1=1','255=NULL'), 'rule.txt', sep="", col.names=FALSE, quote=FALSE, row.names=FALSE)
execGRASS('r.reclass', input='rawDataStudyArea', output='studyArea1', rules='rule.txt', flags=c('overwrite'))

# Forest age - reclassify
rcl<-read.csv(paste0(rawTablesDir, "ageReclass.csv"), header=TRUE)[,c('Value','Code')]
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
execGRASS('r.reclass', input='rawDataLandcoverBTSL', output='landcoverBTSL1_', rules='rule.txt', flags=c('overwrite'))

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
execGRASS('r.patch', input="landcoverAgLinearRas,landcoverBTSL1", output='landcoverBTSL2', flags=c("overwrite"))

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
# Minor
roadMinorDetailedClasses <- roadReclassBDTQ$SQL_DESCRIPTIO[roadReclassBDTQ$Code == minorRoadCode]
whereString = paste0("DESCRIPTIO IN ('",paste(roadMinorDetailedClasses, collapse="','"),"')")
execGRASS("v.extract", input="rawDataRoad", type="line", where=whereString, output="landcoverMinorRoadBDTQ", flags=c("overwrite"))
execGRASS('g.region', res='2')
execGRASS("v.to.rast", input="landcoverMinorRoadBDTQ", use="val", value=minorRoadCode, output="landcoverMinorRoadBDTQRas_2m", flags=c("overwrite"))
execGRASS('g.region', res=paste0(myResolution))
execGRASS('r.resamp.stats', input="landcoverMinorRoadBDTQRas_2m", output=paste0("landcoverMinorRoadBDTQRas_", myResolution, "m"), method="maximum", flags=c("overwrite"))
# Major
roadMajorDetailedClasses <- roadReclassBDTQ$SQL_DESCRIPTIO[roadReclassBDTQ$Code == majorRoadCode]
whereString = paste0("DESCRIPTIO IN ('",paste(roadMajorDetailedClasses, collapse="','"),"')")
execGRASS("v.extract", input="rawDataRoad", type="line", where=whereString, output="landcoverMajorRoadBDTQ", flags=c("overwrite"))
execGRASS('g.region', res='2')
execGRASS("v.to.rast", input="landcoverMajorRoadBDTQ", use="val", value=majorRoadCode, output="landcoverMajorRoadBDTQRas_2m", flags=c("overwrite"))
execGRASS('g.region', res=paste0(myResolution))
execGRASS('r.resamp.stats', input="landcoverMajorRoadBDTQRas_2m", output=paste0("landcoverMajorRoadBDTQRas_", myResolution, "m"), method="maximum", flags=c("overwrite"))
# All roads BDTQ
execGRASS('r.mapcalc', expression=paste0("landcoverRoadBDTQ1 = if(isnull(landcoverMajorRoadBDTQRas_", myResolution, "m),0,landcoverMajorRoadBDTQRas_", myResolution, "m) + if(isnull(landcoverMinorRoadBDTQRas_", myResolution, "m),0,landcoverMinorRoadBDTQRas_", myResolution, "m)"), flags=c("overwrite"))
write.table(c('0=0', paste0(minorRoadCode, '=', minorRoadCode), paste0(majorRoadCode, '=', majorRoadCode), paste0(minorRoadCode + majorRoadCode, '=', majorRoadCode)), 'rule.txt', sep="", col.names=FALSE, quote=FALSE, row.names=FALSE)
execGRASS('r.reclass', input='landcoverRoadBDTQ1', output='landcoverRoadBDTQ', rules='rule.txt', flags=c('overwrite'))

# Combine all roads into a sinlge layer
# NB: assume the PASL roads are correct and then add the BDTQ roads to them
execGRASS('r.mapcalc', expression='landcoverRoad = if(landcoverRoadPASL>0,landcoverRoadPASL,landcoverRoadBDTQ)', flags=c("overwrite"))
execGRASS('r.null', map='landcoverRoad', setnull='0')

# Impose roads on BTSL landcover map
execGRASS('r.patch', input='landcoverRoad,landcoverBTSL2', output='landcoverBTSL', flags=c('overwrite'))

####################
# Landcover Buffer #
####################
# Reclassify landcover in buffer to match Albert et al. classes
rcl<-read.csv(paste0(rawTablesDir, "landcoverBufferReclass.csv"),header=TRUE)[,c('Value','Code')]
write.table(paste0(rcl[,'Value'],'=',rcl[,'Code']),'rule.txt',sep="",col.names=FALSE,quote=FALSE,row.names=FALSE)
execGRASS('r.reclass',input='rawDataLandcoverBuffer',output='landcoverBuffer' ,rules='rule.txt',flags=c('overwrite'))

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
