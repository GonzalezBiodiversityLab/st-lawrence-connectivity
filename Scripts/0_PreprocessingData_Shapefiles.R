# a254
# Bronwyn Rayfield and Jed Lloren, ApexRMS
# Run with R-4.1.1                                
# NB: Spatial layers were originally received as geodatabase and manually converted to shapefiles in QGIS before running code on line 32
# NB: OutaouaisConnectivityExtent.tif was created previously to determine the Phase 4 study extent

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

#JL
# Raw data filenames
countyName <- "MRC_s_2021_02.shp"
densityName <- "OutaouaisConnectivityExtent.tif"
forestageName <- "SIEF-C08PEEFO-forest-age.tif"
landcoverBTSLPolyName <-  "BTSL_SLL_Occ_sol_Land_cover.shp"
landcoverBTSLName <-  "BTSL-SLL-Occ-sol-Land-cover.tif"
landcoverBufferName <- "utilisation_territoire_2018.tif"
privatelandName <- "RMN_20210608.shp"
protectedareaName <- "AP_REG_S_20210824.shp"
roadName <- "AQ_Routes_l_20210701.shp"
studyareaName <- "study-area.tif"
surficialdepositName <- "SIEF-C08PEEFO-surficial-deposits.tif"
# No layer for drainage

# GRASS setup---------------------------------------------------------------------------------------------
if(doGRASSSetup){
  #https://gis.stackexchange.com/questions/183032/create-a-new-grass-database-in-r-with-crs-projection
  # Manually set up empty GRASS database - see GRASSTemplate
  initGRASS(gisBase=gisBase, gisDbase=gisDbase, location=paste0('BTSL_', myResolution, 'm'), mapset='PERMANENT', override=TRUE)
  
  execGRASS("g.proj", georef=file.path(rawMapsDir, studyareaName), flags="c")
  
  # Initialize new mapset inheriting projection info
  execGRASS("g.mapset", mapset = "RawData", flags="c")
  
  # JL
  # Load layers into grass database
  execGRASS("v.in.ogr", input=file.path(rawMapsDir, roadName), output="rawDataRoad", flags=c("overwrite", "o"))
  execGRASS("v.in.ogr", input=file.path(rawMapsDir, protectedareaName), output="rawDataProtectedArea", flags=c("overwrite", "o"))
  execGRASS("v.in.ogr", input=file.path(rawMapsDir, countyName), output="rawDataCounty", flags=c("overwrite", "o"))
  execGRASS("v.in.ogr", input=file.path(rawMapsDir, landcoverBTSLPolyName), output="rawDataLandcoverBTSLPoly", flags=c("overwrite", "o"))
  execGRASS("v.in.ogr", input=file.path(rawMapsDir, privatelandName), output="rawDataPrivateLand", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=file.path(rawMapsDir, densityName), output="rawDataExtent", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=file.path(rawMapsDir, forestageName), output="rawDataForestAge", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=file.path(rawMapsDir, landcoverBTSLName), output="rawDataLandcoverBTSL", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=file.path(rawMapsDir, landcoverBufferName), output="rawDataLandcoverBuffer", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=file.path(rawMapsDir, studyareaName), output="rawDataStudyArea", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=file.path(rawMapsDir, surficialdepositName), output="rawDataSurficialDeposits", flags=c("overwrite", "o"))
}else{
  initGRASS(gisBase=gisBase, gisDbase=gisDbase, location=paste0('BTSL_', myResolution, 'm'), mapset="RawData", override=TRUE)
}

# Set the geographic region
execGRASS('g.region', n='403680', e='-391410', w='-691380', s='117930')

# check your geographic location
execGRASS("g.region", flags = "p") # only included in the working file, to be removed later

# Reclassify layers---------------------------------------------------------------------------------------------
# NB No need to reclassify forest density or drainage
# Study area
# Fix NA values in study area layer: i.e. 2, 3, and 4 = NULL
rclStudyArea<-read.csv(file.path(rawTablesDir, "studyareaReclass.csv"), header=TRUE)[,c('Value','Code')]
write.table(rclStudyArea, file="../Inputs/RawData/Tables/studyareaRules.txt", sep="=", col.names=FALSE, quote=FALSE, row.names=FALSE)
execGRASS('r.reclass', input='rawDataStudyArea', output='studyArea', rules=file.path(rawTablesDir, "studyareaRules.txt"), flags=c('overwrite'))

# Forest age - reclassify
rclAge<-read.csv(file.path(rawTablesDir, "forestAgeReclass.csv"), header=TRUE)[,c('Value','Code')]
write.table(rclAge, file="../Inputs/RawData/Tables/forestAgeRules.txt", sep="=", col.names=FALSE, quote=FALSE, row.names=FALSE)
execGRASS('r.reclass', input='rawDataForestAge', output='forestAge', rules=file.path(rawTablesDir, "forestAgeRules.txt"), flags=c('overwrite'))

# Surficial deposit - reclassify
rclDeposit<-read.csv(file.path(rawTablesDir, "depositReclass.csv"),header=TRUE)[,c('Value','Recode')]
write.table(rclDeposit, file="../Inputs/RawData/Tables/depositRules.txt", sep="=", col.names=FALSE, quote=FALSE, row.names=FALSE)
execGRASS('r.reclass', input='rawDataSurficialDeposits', output='deposit', rules=file.path(rawTablesDir, "depositRules.txt"), flags=c('overwrite'))

# Landcover BTSL---------------------------------------------------------------------------------------------
# Reclass tables
speciesLandcoverReclass <- read.csv(file.path(rawTablesDir, "speciesLandcoverReclass.csv"), header=TRUE)
landcoverReclassRaw <- read.csv(file.path(rawTablesDir, "landcoverBTSLReclass.csv"), header=TRUE)
roadReclassBDTQ <- read.csv(file.path(rawTablesDir, "roadReclassBDTQ.csv"), header=TRUE)

# Reclassify landcover in BTSL to match Albert et al. classes
landcoverReclass<-landcoverReclassRaw[, c("Value", "Code")]
write.table(landcoverReclass, file="../Inputs/RawData/Tables/landReclassRule.txt", sep="=", col.names=FALSE, quote=FALSE, row.names=FALSE)
execGRASS('r.reclass', input='rawDataLandcoverBTSL', output='landcoverBTSLReclass', rules=file.path(rawTablesDir, "landReclassRule.txt"), flags=c('overwrite'))

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
execGRASS('r.patch', input=paste0("landcoverAgLinearRas_", myResolution, "m,landcoverBTSLReclass"), output=paste0("landcoverBTSL", myResolution, "m"), flags=c("overwrite"))

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

# Add roads from AQ database
# Minor
roadMinorDetailedClasses <- roadReclassAQ$SQL_DESCRIPTIO[roadReclassAQ$Code == minorRoadCode]
whereString = paste0("DESCRIPTIO IN ('",paste(roadMinorDetailedClasses, collapse="','"),"')")
execGRASS("v.extract", input="rawDataRoad", type="line", where=whereString, output="landcoverMinorRoadAQ", flags=c("overwrite"))
execGRASS('g.region', res='2')
execGRASS("v.to.rast", input="landcoverMinorRoadAQ", use="val", value=minorRoadCode, output="landcoverMinorRoadAQRas_2m", flags=c("overwrite"))
execGRASS('g.region', res=paste0(myResolution))
execGRASS('r.resamp.stats', input="landcoverMinorRoadAQRas_2m", output=paste0("landcoverMinorRoadAQRas_", myResolution, "m"), method="maximum", flags=c("overwrite"))
# Major
roadMajorDetailedClasses <- roadReclassAQ$SQL_DESCRIPTIO[roadReclassAQ$Code == majorRoadCode]
whereString = paste0("DESCRIPTIO IN ('",paste(roadMajorDetailedClasses, collapse="','"),"')")
execGRASS("v.extract", input="rawDataRoad", type="line", where=whereString, output="landcoverMajorRoadAQ", flags=c("overwrite"))
execGRASS('g.region', res='2')
execGRASS("v.to.rast", input="landcoverMajorRoadAQ", use="val", value=majorRoadCode, output="landcoverMajorRoadAQRas_2m", flags=c("overwrite"))
execGRASS('g.region', res=paste0(myResolution))
execGRASS('r.resamp.stats', input="landcoverMajorRoadAQRas_2m", output=paste0("landcoverMajorRoadAQRas_", myResolution, "m"), method="maximum", flags=c("overwrite"))
# All roads AQ
execGRASS('r.mapcalc', expression=paste0("landcoverRoadAQ1 = if(isnull(landcoverMajorRoadAQRas_", myResolution, "m),0,landcoverMajorRoadAQRas_", myResolution, "m) + if(isnull(landcoverMinorRoadAQRas_", myResolution, "m),0,landcoverMinorRoadAQRas_", myResolution, "m)"), flags=c("overwrite"))
write.table(c('0=0', paste0(minorRoadCode, '=', minorRoadCode), paste0(majorRoadCode, '=', majorRoadCode), paste0(minorRoadCode + majorRoadCode, '=', majorRoadCode)), 'rule.txt', sep="", col.names=FALSE, quote=FALSE, row.names=FALSE)
execGRASS('r.reclass', input='landcoverRoadAQ1', output='landcoverRoadAQ', rules='rule.txt', flags=c('overwrite'))

# Combine all roads into a single layer
# NB: assume the PASL roads are correct and then add the BDTQ roads to them
execGRASS('r.mapcalc', expression='landcoverRoad = if(landcoverRoadPASL>0,landcoverRoadPASL,landcoverRoadBDTQ)', flags=c("overwrite"))
execGRASS('r.null', map='landcoverRoad', setnull='0')

# Impose roads on BTSL landcover map
execGRASS('r.patch', input='landcoverRoad,landcoverBTSL2', output='landcoverBTSL', flags=c('overwrite'))

####################
# Landcover Buffer #
####################
# Reclassify landcover in buffer to match Albert et al. classes
rclBuffer<-read.csv(file.path(rawTablesDir, "landcoverBufferReclass.csv"),header=TRUE)[,c('Value','Code')]
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
