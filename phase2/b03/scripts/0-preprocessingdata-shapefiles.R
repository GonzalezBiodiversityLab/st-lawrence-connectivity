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
b03Dir <- "b03"
b01b02Dir <- "b01b02"
gisDbase <- file.path(b03Dir, "grass7")
b03RawMapsDir <- file.path(b03Dir, "data", "spatial")
b01b02RawMapsDir <- file.path(b01b02Dir, "data", "spatial")
b03RawTablesDir <- file.path(b03Dir, "data", "tabular")
b01b02RawTablesDir <- file.path(b01b02Dir, "data", "tables")
b03ProcessedMapsDir <- file.path(b03Dir, "model-outputs")

#JL
# Raw data filenames
countyName <- "MRC_s_2021_02.shp"
extentName <- "OutaouaisConnectivityExtent.tif"
forestageName <- "SIEF-C08PEEFO-forest-age.tif"
forestDensityName <- "SIEF-C08PEEFO-forest-density.tif"
landcoverBTSLPolyName <- "BTSL_SLL_Occ_sol_Land_cover.shp"
landcoverBTSLName <- "BTSL-SLL-Occ-sol-Land-cover.tif"
landcoverBufferName <- "utilisation_territoire_2018.tif"
landcoverBufferFillName <- "Occ_sol_2014_recl_FED_10m_aout2017.tif"
privatelandName <- "RMN_20210608.shp"
protectedareaName <- "AP_REG_S_20210824.shp"
roadName <- "AQ_Routes_l_20210701.shp"
studyareaName <- "b03-studyArea.tif"
surficialdepositName <- "SIEF-C08PEEFO-surficial-deposits.tif"
drainageName <- "RCL_DRAINAGE.tif"

# GRASS setup---------------------------------------------------------------------------------------------
if(doGRASSSetup){
  #https://gis.stackexchange.com/questions/183032/create-a-new-grass-database-in-r-with-crs-projection
  # Manually set up empty GRASS database - see GRASSTemplate
  initGRASS(gisBase=gisBase, gisDbase=gisDbase, location=paste0('BTSL_', myResolution, 'm'), mapset='PERMANENT', override=TRUE)
  
  execGRASS("g.proj", georef=file.path(b03RawMapsDir, studyareaName), flags="c")
  
  # Initialize new mapset inheriting projection info
  execGRASS("g.mapset", mapset = "RawData", flags="c")
  
  # JL
  # Load layers into grass database
  execGRASS("v.in.ogr", input=file.path(b03RawMapsDir, roadName), output="rawDataRoad", flags=c("overwrite", "o"))
  execGRASS("v.in.ogr", input=file.path(b03RawMapsDir, protectedareaName), output="rawDataProtectedArea", flags=c("overwrite", "o"))
  execGRASS("v.in.ogr", input=file.path(b03RawMapsDir, countyName), output="rawDataCounty", flags=c("overwrite", "o"))
  execGRASS("v.in.ogr", input=file.path(b03RawMapsDir, landcoverBTSLPolyName), output="rawDataLandcoverBTSLPoly", flags=c("overwrite", "o"))
  execGRASS("v.in.ogr", input=file.path(b03RawMapsDir, privatelandName), output="rawDataPrivateLand", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=file.path(b03RawMapsDir, extentName), output="rawDataExtent", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=file.path(b03RawMapsDir, forestageName), output="rawDataForestAge", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=file.path(b03RawMapsDir, forestDensityName), output="rawDataForestDensity", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=file.path(b03RawMapsDir, landcoverBTSLName), output="rawDataLandcoverBTSL", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=file.path(b03RawMapsDir, landcoverBufferName), output="rawDataLandcoverBuffer", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=file.path(b01b02RawMapsDir, landcoverBufferFillName), output="rawDataLandcoverBufferFill", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=file.path(b03RawMapsDir, studyareaName), output="rawDataStudyArea", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=file.path(b03RawMapsDir, surficialdepositName), output="rawDataSurficialDeposits", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=file.path(b03RawMapsDir, drainageName), output="rawDataDrainage", flags=c("overwrite", "o"))
}else{
  initGRASS(gisBase=gisBase, gisDbase=gisDbase, location=paste0('BTSL_', myResolution, 'm'), mapset="RawData", override=TRUE)
}

# Set the geographic region
execGRASS('g.region', n='403680', e='-391410', w='-691380', s='117930') # This data should already be set when running line 53

# check your geographic location
execGRASS("g.region", flags = "p") # only included in the working file, to be removed later

# Reclassify layers---------------------------------------------------------------------------------------------
# NB No need to reclassify forest density or drainage
# Study area
# Fix NA values in study area layer: i.e. 2, 3, and 4 = NULL
# See 0-create-study-area-raster.R
# rclStudyArea <- read.csv(file.path("b01b02", "data", "tables", "studyareaReclass.csv"), header=TRUE)[,c('Value','Code')]
# write.table(rclStudyArea, file="b01b02/data/tables/studyareaRules.txt", sep="=", col.names=FALSE, quote=FALSE, row.names=FALSE)
# execGRASS('r.reclass', input='rawDataStudyArea', output='b03-studyArea-test', rules=file.path("b01b02", "data", "tables", "studyareaRules.txt"), flags=c('overwrite'))

# Forest age - reclassify
rclAge<-read.csv(file.path(b01b02RawTablesDir, "forestAgeReclass.csv"), header=TRUE)[,c('Value','Code')]
write.table(rclAge, file=file.path(b01b02RawTablesDir, "forestAgeRules.txt"), sep="=", col.names=FALSE, quote=FALSE, row.names=FALSE)
execGRASS('r.reclass', input='rawDataForestAge', output='forestAge', rules=file.path(b01b02RawTablesDir, "forestAgeRules.txt"), flags=c('overwrite'))

# Surficial deposit - reclassify
rclDeposit<-read.csv(file.path(b01b02RawTablesDir, "depositReclass.csv"),header=TRUE)[,c('Value','Recode')]
write.table(rclDeposit, file=file.path(b01b02RawTablesDir, "depositRules.txt"), sep="=", col.names=FALSE, quote=FALSE, row.names=FALSE)
execGRASS('r.reclass', input='rawDataSurficialDeposits', output='surficialDeposit', rules=file.path(b01b02RawTablesDir, "depositRules.txt"), flags=c('overwrite'))

# Landcover BTSL---------------------------------------------------------------------------------------------
# Reclass tables
speciesLandcoverReclass <- read.csv(file.path(b01b02RawTablesDir, "speciesLandcoverReclass.csv"), header=TRUE)
landcoverReclassRaw <- read.csv(file.path(b01b02RawTablesDir, "landcoverBTSLReclass.csv"), header=TRUE)
roadReclassBDTQ <- read.csv(file.path(b01b02RawTablesDir, "roadReclassBDTQ.csv"), header=TRUE)
roadReclassAQ <- read.csv(file.path(b03RawTablesDir, "roadReclassAQ.csv"), header=TRUE)

# Reclassify landcover in BTSL to match Albert et al. classes
landcoverReclass <- landcoverReclassRaw[, c("Value", "Code")]
write.table(landcoverReclass, file=file.path(b01b02RawTablesDir, "landReclassRule.txt"), sep="=", col.names=FALSE, quote=FALSE, row.names=FALSE)
execGRASS('r.reclass', input='rawDataLandcoverBTSL', output='landcoverBTSLReclass', rules=file.path(b01b02RawTablesDir, "landReclassRule.txt"), flags=c('overwrite'))

# Fix broken linear features by extracting them, rasterizing them at higher resolution, and then reimposing them on the landcover raster
# Extract linear agriculture (i.e. "Milieu agricole non cultiv?")
agLinearCode <- speciesLandcoverReclass$LandcoverCode[speciesLandcoverReclass$LandcoverName == "AgricultureLinearElements"]
agLinearDetailedClasses <- landcoverReclassRaw$SQL_CLASSE_DET[landcoverReclassRaw$Code == agLinearCode]
whereString <- paste0("CLASSE_DET IN ('",paste(agLinearDetailedClasses, collapse="','"),"')")
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
roadMinorDetailedClasses <- landcoverReclassRaw$SQL_CLASSE_DET[landcoverReclassRaw$Code == minorRoadCode]
whereString = paste0("CLASSE_DET IN ('",paste(roadMinorDetailedClasses, collapse="','"),"')")
execGRASS("v.extract", input="rawDataLandcoverBTSLPoly", layer='1', type="area", where=whereString, output="landcoverMinorRoadPASLPoly", flags=c("overwrite"))
execGRASS('g.region', res='2')
execGRASS("v.to.rast", input="landcoverMinorRoadPASLPoly", use="val", value=minorRoadCode, output="landcoverMinorRoadPASLRas_2m", flags=c("overwrite"))
execGRASS('g.region', res=paste0(myResolution))
execGRASS('r.resamp.stats', input="landcoverMinorRoadPASLRas_2m", output=paste0("landcoverMinorRoadPASLRas_", myResolution, "m"), method="maximum", flags=c("overwrite"))
# Major
roadMajorDetailedClasses <- landcoverReclassRaw$SQL_CLASSE_DET[landcoverReclassRaw$Code == majorRoadCode]
whereString = paste0("CLASSE_DET IN ('",paste(roadMajorDetailedClasses, collapse="','"),"')")
execGRASS("v.extract", input="rawDataLandcoverBTSLPoly", layer='1', type="area", where=whereString, output="landcoverMajorRoadPASLPoly", flags=c("overwrite"))
execGRASS('g.region', res='2')
execGRASS("v.to.rast", input="landcoverMajorRoadPASLPoly", use="val", value=majorRoadCode, output="landcoverMajorRoadPASLRas_2m", flags=c("overwrite"))
execGRASS('g.region', res=paste0(myResolution))
execGRASS('r.resamp.stats', input="landcoverMajorRoadPASLRas_2m", output=paste0("landcoverMajorRoadPASLRas_", myResolution, "m"), method="maximum", flags=c("overwrite"))
# All roads PASL
execGRASS('r.mapcalc', expression=paste0("landcoverRoadPASL1 = if(isnull(landcoverMajorRoadPASLRas_", myResolution, "m),0,landcoverMajorRoadPASLRas_", myResolution, "m) + if(isnull(landcoverMinorRoadPASLRas_", myResolution, "m),0,landcoverMinorRoadPASLRas_", myResolution, "m)"), flags=c("overwrite"))
write.table(c('0=0', paste0(minorRoadCode, '=', minorRoadCode), paste0(majorRoadCode, '=', majorRoadCode), paste0(minorRoadCode + majorRoadCode, '=', majorRoadCode)), paste0(b03RawTablesDir, '/roadPASLrule.txt'), sep="", col.names=FALSE, quote=FALSE, row.names=FALSE)
execGRASS('r.reclass', input='landcoverRoadPASL1', output='landcoverRoadPASL', rules=paste0(b03RawTablesDir, '/roadPASLrule.txt'), flags=c('overwrite'))

# # Add roads from AQ database
# # Minor
roadMinorDetailedClasses <- roadReclassAQ$ClsRte_SQL[roadReclassAQ$Code == minorRoadCode]
whereString = paste0("ClsRte IN ('",paste(roadMinorDetailedClasses, collapse="','"),"')")
execGRASS("v.extract", input="rawDataRoad", type="line", where=whereString, output="landcoverMinorRoadAQ", flags=c("overwrite"))
execGRASS('g.region', res='2')
execGRASS("v.to.rast", input="landcoverMinorRoadAQ", use="val", value=minorRoadCode, output="landcoverMinorRoadAQRas_2m", flags=c("overwrite"))
execGRASS('g.region', res=paste0(myResolution))
execGRASS('r.resamp.stats', input="landcoverMinorRoadAQRas_2m", output=paste0("landcoverMinorRoadAQRas_", myResolution, "m"), method="maximum", flags=c("overwrite"))
# Major
roadMajorDetailedClasses <- roadReclassAQ$ClsRte_SQL[roadReclassAQ$Code == majorRoadCode]
whereString = paste0("ClsRte IN ('",paste(roadMajorDetailedClasses, collapse="','"),"')")
execGRASS("v.extract", input="rawDataRoad", type="line", where=whereString, output="landcoverMajorRoadAQ", flags=c("overwrite"))
execGRASS('g.region', res='2')
execGRASS("v.to.rast", input="landcoverMajorRoadAQ", use="val", value=majorRoadCode, output="landcoverMajorRoadAQRas_2m", flags=c("overwrite"))
execGRASS('g.region', res=paste0(myResolution))
execGRASS('r.resamp.stats', input="landcoverMajorRoadAQRas_2m", output=paste0("landcoverMajorRoadAQRas_", myResolution, "m"), method="maximum", flags=c("overwrite"))
# All roads AQ
execGRASS('r.mapcalc', expression=paste0("landcoverRoadAQ1 = if(isnull(landcoverMajorRoadAQRas_", myResolution, "m),0,landcoverMajorRoadAQRas_", myResolution, "m) + if(isnull(landcoverMinorRoadAQRas_", myResolution, "m),0,landcoverMinorRoadAQRas_", myResolution, "m)"), flags=c("overwrite"))
write.table(c('0=0', paste0(minorRoadCode, '=', minorRoadCode), paste0(majorRoadCode, '=', majorRoadCode), paste0(minorRoadCode + majorRoadCode, '=', majorRoadCode)), paste0(b03RawTablesDir, '/rule.txt'), sep="", col.names=FALSE, quote=FALSE, row.names=FALSE)
execGRASS('r.reclass', input='landcoverRoadAQ1', output='landcoverRoadAQ', rules=paste0(b03RawTablesDir, '/rule.txt'), flags=c('overwrite'))

# Combine all roads into a single layer
# NB: assume the PASL roads are correct and then add the AQ roads to them
execGRASS('r.mapcalc', expression='landcoverRoad = if(landcoverRoadPASL>0,landcoverRoadPASL,landcoverRoadAQ)', flags=c("overwrite"))
execGRASS('r.null', map='landcoverRoad', setnull='0')

# Impose roads on BTSL landcover map
execGRASS('r.patch', input=paste0('landcoverRoad,landcoverBTSL', myResolution, 'm'), output='landcoverBTSL', flags=c('overwrite'))

####################
# Landcover Buffer #
####################
# Reclassify landcover in Quebec to match buffer classes
rcl<-read.csv(file.path(b01b02RawTablesDir, "landcoverQuebecReclass.csv"),header=TRUE)[,c('Code','Recode')]
write.table(paste0(rcl[,'Code'],'=',rcl[,'Recode']),file.path(b03RawTablesDir,'landCoverQuebecReclassRule.txt'),sep="",col.names=FALSE,quote=FALSE,row.names=FALSE)
execGRASS('r.reclass',input='rawDataLandcoverBuffer',output='landcoverBufferClasseGen',rules=file.path(b03RawTablesDir,'landCoverQuebecReclassRule.txt'),flags=c('overwrite'))

# Fill out buffer with Quebec landcover
execGRASS('r.patch', input="rawDataLandcoverBufferFill,landcoverBufferClasseGen", output='landcoverBufferFilled', flags=c('overwrite'))

# Reclassify landcover in buffer to match Albert et al. classes
rcl<-read.csv(file.path(b01b02RawTablesDir, "landcoverBufferReclass.csv"),header=TRUE)[,c('Value','Code')]
write.table(paste0(rcl[,'Value'],'=',rcl[,'Code']),file.path(b03RawTablesDir, 'landcoverBufferReclassRule.txt'),sep="",col.names=FALSE,quote=FALSE,row.names=FALSE)
execGRASS('r.reclass',input='landcoverBufferFilled',output='landcoverBuffer' ,rules=file.path(b03RawTablesDir, 'landcoverBufferReclassRule.txt'),flags=c('overwrite'))

#####################################
# Combine BTSL and Buffer landcover #
#####################################
execGRASS('r.patch', input='landcoverBTSL,landcoverBuffer', output='landcover', flags=c('overwrite'))

############################
# Save outputs to geotiffs #
############################
execGRASS('r.out.gdal', input='landcover',output=paste0(b03ProcessedMapsDir, '/b03_landcover_', myResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))
execGRASS('r.out.gdal', input='forestAge',output=paste0(b03ProcessedMapsDir, '/b03_forestAge_', myResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))
execGRASS('r.out.gdal', input='rawDataForestDensity', output=paste0(b03ProcessedMapsDir, '/b03_forestDensity_', myResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW',flags=c('overwrite'))
execGRASS('r.out.gdal', input='rawDataDrainage', output=paste0(b03ProcessedMapsDir, '/b03_drainage_', myResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))
execGRASS('r.out.gdal', input='surficialDeposit', output=paste0(b03ProcessedMapsDir, '/b03_deposit_', myResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))
execGRASS('r.out.gdal', input='rawDataStudyArea', output=paste0(b03ProcessedMapsDir, '/b03_studyArea_', myResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))
