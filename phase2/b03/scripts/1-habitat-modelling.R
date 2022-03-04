# a254
# Bronwyn Rayfield and Sarah Chisholm, ApexRMS
# Run with R-4.1.2  

# Load constants, functions, etc
source("./b03/scripts/0-0-constants.R")

# Workspace ----

# Input parameters
# Data filenames
landcoverName <- paste0("b03-landcover-", myResolution, "m.tif")
ageName <- paste0("b03-forest-age-", myResolution, "m.tif")
densityName <- paste0("b03-forest-density-", myResolution, "m.tif")
depositName <- paste0("b03-deposit-", myResolution, "m.tif")
drainageName <- paste0("b03-drainage-", myResolution, "m.tif")
studyAreaName <- paste0("b03-studyarea-", myResolution, "m.tif")

###############
# GRASS setup #
###############
if(doGRASSSetup){
  #https://gis.stackexchange.com/questions/183032/create-a-new-grass-database-in-r-with-crs-projection
  # Manually set up empty GRASS database - see GRASSTemplate
  initGRASS(gisBase=gisBase, gisDbase=gisDbase, location=paste0("BTSL_", myResolution, "m"), mapset='PERMANENT', override=TRUE)

  execGRASS("g.proj", georef=file.path(b03ProcessedMapsDir, landcoverName), flags="c")
  
  # Initialize new mapset inheriting projection info
  execGRASS("g.mapset", mapset="HabitatModelling", flags="c")
  
  # Load layers into grass database
  execGRASS("r.in.gdal", input=file.path(b03ProcessedMapsDir, landcoverName), output="landcover", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=file.path(b03ProcessedMapsDir, ageName), output="forestAge", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=file.path(b03ProcessedMapsDir, densityName), output="forestDensity", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=file.path(b03ProcessedMapsDir, depositName), output="deposit", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=file.path(b03ProcessedMapsDir, drainageName), output="drainage", flags=c("overwrite", "o"))
  execGRASS("r.in.gdal", input=file.path(b03ProcessedMapsDir, studyAreaName), output="studyArea1", flags=c("overwrite", "o"))
}else{
  initGRASS(gisBase=gisBase, gisDbase=gisDbase, location=paste0("BTSL_", myResolution, "m"), mapset='HabitatModelling', override=TRUE)
}

# Set the geographic region
execGRASS('g.region', n='403680', e='-391410', w='-691380', s='117930') 

#########################################################
# Landscape analyses to prepare for habitat suitability #
#########################################################
# Create study area layer with 0's instead of 1's
write.table('1=0','rule.txt',sep="",col.names=FALSE,quote=FALSE,row.names=FALSE)
execGRASS('r.reclass',input='studyArea1',output='studyArea0',rules='rule.txt',flags=c('overwrite'))

# Read in species parameters
speciesLandcoverReclass <- read.csv(file.path(b01b02RawTablesDir, "speciesLandcoverReclass.csv"), header=TRUE)
speciesHabitatPatchParams <- read.csv(file.path(b01b02RawTablesDir, "speciesHabitatPatchParameters.csv"), header=TRUE)
speciesResistanceReclass <- read.csv(file.path(b01b02RawTablesDir, "speciesResistanceReclass.csv"), header=TRUE)

# Isolate fallow lands and reclss as code = 99
fallowCode <- speciesLandcoverReclass$LandcoverCode[speciesLandcoverReclass$LandcoverName %in% c("Fallow", "FallowLinearElements")]
write.table(c(paste0(fallowCode, '=99'),'*=NULL'),'rule.txt',sep="",col.names=FALSE,quote=FALSE,row.names=FALSE)
execGRASS('r.reclass',input='landcover',output='fallow',rules='rule.txt',flags=c('overwrite'))
# Add fallow to density map
execGRASS('r.patch',input='fallow,forestDensity',output='forestDensityFallow',flags=c('overwrite'))
# Add fallow to age map
execGRASS('r.patch',input='fallow,forestAge',output='forestAgeFallow',flags=c('overwrite'))

# Isolate major roads
majorRoadCode <- speciesLandcoverReclass$LandcoverCode[speciesLandcoverReclass$LandcoverName == "MajorRoads"]
write.table(c(paste0(majorRoadCode, '=1'),'*=NULL'),'rule.txt',sep="",col.names=FALSE,quote=FALSE,row.names=FALSE)
execGRASS('r.reclass',input='landcover',output='majorRoad',rules='rule.txt',flags=c('overwrite'))
# Calculate distance from major roads
execGRASS('r.grow.distance',input='majorRoad',distance='distanceFromMajorRoad',flags=c('overwrite'))

# Isolate minor roads    
minorRoadCode <- speciesLandcoverReclass$LandcoverCode[speciesLandcoverReclass$LandcoverName == "MinorRoads"]
write.table(c(paste0(minorRoadCode, '=1'),'*=NULL'),'rule.txt',sep="",col.names=FALSE,quote=FALSE,row.names=FALSE)
execGRASS('r.reclass',input='landcover',output='minorRoad',rules='rule.txt',flags=c('overwrite'))
# Calculate distance from minor roads
execGRASS('r.grow.distance',input='minorRoad',distance='distanceFromMinorRoad',flags=c('overwrite'))

# Isolate urban
urbanCode <- speciesLandcoverReclass$LandcoverCode[speciesLandcoverReclass$LandcoverName == "Urban"]
write.table(c(paste0(urbanCode, '=1'),'*=NULL'),'rule.txt',sep="",col.names=FALSE,quote=FALSE,row.names=FALSE)
execGRASS('r.reclass',input='landcover',output='urban',rules='rule.txt',flags=c('overwrite'))
# Calculate distance from urban
execGRASS('r.grow.distance',input='urban',distance='distanceFromUrban',flags=c('overwrite'))

# Isolate forest
forestCode <- speciesLandcoverReclass$LandcoverCode[speciesLandcoverReclass$LandcoverName %in% c("ForestDeciduous", "ForestMixed", "ForestConiferous")]
write.table(c(paste0(forestCode, '=1'),'*=NULL'),'rule.txt',sep="",col.names=FALSE,quote=FALSE,row.names=FALSE)
execGRASS('r.reclass',input='landcover',output='forest',rules='rule.txt',flags=c('overwrite'))
# Calculate distance from forest
execGRASS('r.grow.distance',input='forest',distance='distanceFromForest',flags=c('overwrite'))

# Isolate forest and fallow
forestFallowCode <- speciesLandcoverReclass$LandcoverCode[speciesLandcoverReclass$LandcoverName %in% c("ForestDeciduous", "ForestMixed", "ForestConiferous", "Fallow")]
write.table(c(paste0(forestFallowCode, '=1'),'*=NULL'),'rule.txt',sep="",col.names=FALSE,quote=FALSE,row.names=FALSE)
execGRASS('r.reclass',input='landcover',output='forestFallow',rules='rule.txt',flags=c('overwrite'))

# Isolate non-forest
execGRASS('r.mapcalc', expression='nonforest1=if(isnull(forest),1,0)', flags=c('overwrite'))
write.table(c('1=1','*=NULL'), 'rule.txt', sep="", col.names=FALSE, quote=FALSE, row.names=FALSE)
execGRASS('r.reclass', input='nonforest1', output='nonforest', rules='rule.txt', flags=c('overwrite'))
# Calculate distance from nonforest (aka forest edge)
execGRASS('r.grow.distance',input='nonforest',distance='distanceFromForestEdge',flags=c('overwrite'))

# Isolate wetlands
wetlandCode <- speciesLandcoverReclass$LandcoverCode[speciesLandcoverReclass$LandcoverName %in% c("WetlandsOpen", "WetlandsTreed")]
write.table(c(paste0(wetlandCode, '=1'),'*=NULL'),'rule.txt', sep="", col.names=FALSE, quote=FALSE, row.names=FALSE)
execGRASS('r.reclass', input='landcover', output='wetlands', rules='rule.txt', flags=c('overwrite'))


#######################
# Habitat suitability #
#######################
for(i in 1:length(speciesList)){
  species<-speciesList[i]
  
  # Landcover
  # Reclassify landcover based on suitability
  write.table(paste0(speciesLandcoverReclass[,'LandcoverCode'],'=',speciesLandcoverReclass[,species]),'rule.txt',sep="",col.names=FALSE,quote=FALSE,row.names=FALSE)
  execGRASS('r.reclass',input='landcover',output=paste0(species, '_landcover'),rules='rule.txt',flags=c('overwrite'))
  
  # Forest density
  # Reclassify forest density based on suitability
  rcl<-read.csv(paste0(b01b02RawTablesDir, "/speciesForestDensityReclass.csv"),header=TRUE)[,c('DensityCode',species)]
  write.table(paste0(rcl[,'DensityCode'],'=',rcl[,species]),'rule.txt',sep="",col.names=FALSE,quote=FALSE,row.names=FALSE)
  execGRASS('r.reclass',input='forestDensityFallow',output=paste0(species,'_forestDensity'),rules='rule.txt',flags=c('overwrite'))
  
  # Forest age
  # Reclassify forest age based on suitability
  rcl<-read.csv(paste0(b01b02RawTablesDir, "/speciesForestAgeReclass.csv"),header=TRUE)[,c('AgeCode',species)]
  write.table(paste0(rcl[,'AgeCode'],'=',rcl[,species]),'rule.txt',sep="",col.names=FALSE,quote=FALSE,row.names=FALSE)
  execGRASS('r.reclass',input='forestAgeFallow',output=paste0(species,'_forestAge'),rules='rule.txt',flags=c('overwrite'))
  
  # Drainage
  # Reclassify drainage based on suitability
  rcl<-read.csv(paste0(b01b02RawTablesDir, "/speciesDrainageReclass.csv"),header=TRUE)[,c('DrainageCode',species)]
  write.table(paste0(rcl[,'DrainageCode'],'=',rcl[,species]),'rule.txt',sep="",col.names=FALSE,quote=FALSE,row.names=FALSE)
  execGRASS('r.reclass',input='drainage',output=paste0(species,'_drainage'),rules='rule.txt',flags=c('overwrite'))
  
  # Deposit
  # Reclassify deposit based on suitability
  rcl<-read.csv(paste0(b01b02RawTablesDir, "/speciesDepositReclass.csv"),header=TRUE)[,c('DepositCode',species)]
  write.table(paste0(rcl[,'DepositCode'],'=',rcl[,species]),'rule.txt',sep="",col.names=FALSE,quote=FALSE,row.names=FALSE)
  execGRASS('r.reclass',input='deposit',output=paste0(species,'_deposit'),rules='rule.txt',flags=c('overwrite'))
  
  # Roads
  # Reclassify distance from major roads based on species suitability
  rcl<-read.csv(paste0(b01b02RawTablesDir,"/speciesDistanceToMajorRoadReclass.csv"),header=TRUE)[,c('From','To',species)]
  write.table(paste0(rcl[,'From'],' thru ',rcl[,'To'],'=',rcl[,species]),'rule.txt',sep="",col.names=FALSE,quote=FALSE,row.names=FALSE)
  execGRASS('r.reclass',input='distanceFromMajorRoad',output=paste0(species,'_distanceFromMajorRoad'),rules='rule.txt',flags=c('overwrite'))
  execGRASS('r.mapcalc',expression=paste0(species,'_distanceFromMajorRoad=',species,'_distanceFromMajorRoad'),flags=c('overwrite'))
  
  # Reclassify distance from minor roads based on species suitability
  rcl<-read.csv(paste0(b01b02RawTablesDir,"/speciesDistanceToMinorRoadReclass.csv"),header=TRUE)[,c('From','To',species)]
  write.table(paste0(rcl[,'From'],' thru ',rcl[,'To'],'=',rcl[,species]),'rule.txt',sep="",col.names=FALSE,quote=FALSE,row.names=FALSE)
  execGRASS('r.reclass',input='distanceFromMinorRoad',output=paste0(species,'_distanceFromMinorRoad'),rules='rule.txt',flags=c('overwrite'))
  execGRASS('r.mapcalc',expression=paste0(species,'_distanceFromMinorRoad=',species,'_distanceFromMinorRoad'),flags=c('overwrite'))
  
  # Urban
  # Reclassify distance from urban based on species suitability
  rcl<-read.csv(paste0(b01b02RawTablesDir,"/speciesDistanceToUrbanReclass.csv"),header=TRUE)[,c('From','To',species)]
  write.table(paste0(rcl[,'From'],' thru ',rcl[,'To'],'=',rcl[,species]),'rule.txt',sep="",col.names=FALSE,quote=FALSE,row.names=FALSE)
  execGRASS('r.reclass',input='distanceFromUrban',output=paste0(species,'_distanceFromUrban'),rules='rule.txt',flags=c('overwrite'))
    
  # Forest
  # Reclassify distance from forest based on species suitability
  rcl<-read.csv(paste0(b01b02RawTablesDir,"/speciesDistanceToForestReclass.csv"),header=TRUE)[,c('From','To',species)]
  write.table(paste0(rcl[,'From'],' thru ',rcl[,'To'],'=',rcl[,species]),'rule.txt',sep="",col.names=FALSE,quote=FALSE,row.names=FALSE)
  execGRASS('r.reclass',input='distanceFromForest',output=paste0(species,'_distanceFromForest'),rules='rule.txt',flags=c('overwrite'))

  # Forest edge
  # Reclassify distance from forest edge (nonforest) based on species suitability
  rcl<-read.csv(paste0(b01b02RawTablesDir, "/speciesDistanceToForestEdgeReclass.csv"), header=TRUE)[, c('From', 'To', species)]
  write.table(paste0(rcl[,'From'], ' thru ', rcl[,'To'], '=', rcl[,species]), 'rule.txt', sep="", col.names=FALSE, quote=FALSE, row.names=FALSE)
  execGRASS('r.reclass',input='distanceFromForestEdge',output=paste0(species,'_distanceFromForestEdge'),rules='rule.txt',flags=c('overwrite'))

  # Combine landcover, forest age, forest density, drainage, deposit type, distance to minor and major roads, distance to urban
  execGRASS('r.mapcalc', expression=paste0(species, '_habitatSuitability', '=', species, '_landcover*if(isnull(', species, '_forestDensity),100,',species, '_forestDensity)/100*if(isnull(', species, '_forestAge),100,', species, '_forestAge)/100*if(isnull(', species, '_drainage),100,', species, '_drainage)/100*if(isnull(', species, '_deposit),100,', species, '_deposit)/100*', species, '_distanceFromMinorRoad/100*', species,'_distanceFromMajorRoad/100*', species, '_distanceFromUrban/100*', species, '_distanceFromForest/100*', species, '_distanceFromForestEdge/100'), flags=c('overwrite'))

  if(species == 'RASY'){
    #Identify potential suitable habitat (habitat suitability > 60)
    execGRASS('r.mapcalc', expression=paste0(species, '_habitatSuitabilityPotential', '=' ,species, '_habitatSuitability'), flags=c('overwrite'))
    write.table(c('0 thru 59=NULL', '60 thru 100=1'), 'rule.txt', sep="", col.names=FALSE, quote=FALSE, row.names=FALSE)
    execGRASS('r.reclass', input=paste0(species, '_habitatSuitability'), output=paste0(species, '_habitatPotential'), rules='rule.txt', flags=c('overwrite'))
    #Keep only terrestrial habitat within potential suitable habitat
    execGRASS('r.mapcalc', expression=paste0(species, '_habitatTerrestrialPotential','=if(isnull(forestFallow),0,forestFallow)*', species, '_habitatPotential'), flags=c('overwrite'))
    
    #Distance from terrestrial habitat
    execGRASS('r.grow.distance', input=paste0(species,'_habitatTerrestrialPotential'), distance=paste0(species, '_distanceFromHabitatTerrestrial1'), flags=c('overwrite'))
    #Reclassify distance from terrestrial habitat based on species suitability
    write.table(c('0 thru 200=100','200 thru 600=70','600 thru 200000=59','*=NULL'), 'rule.txt', sep="", col.names=FALSE, quote=FALSE, row.names=FALSE)
    execGRASS('r.reclass', input=paste0(species,'_distanceFromHabitatTerrestrial1'), output=paste0(species, '_distanceFromHabitatTerrestrial'), rules='rule.txt', flags=c('overwrite'))
    
    #Identify aquatic habitat (i.e. wetlands within disperal distance from terrestrial habitat)
    execGRASS('r.mapcalc', expression=paste0(species, '_wetlands=wetlands', '*', species, '_distanceFromHabitatTerrestrial'), flags=c('overwrite'))
    write.table(c('0 thru 59=NULL','60 thru 100=1'), 'rule.txt', sep="", col.names=FALSE, quote=FALSE, row.names=FALSE)
    execGRASS('r.reclass', input=paste0(species, '_wetlands'), output=paste0(species, '_habitatAquatic'), rules='rule.txt', flags=c('overwrite'))
    
    #Distance from aquatic habitat
    execGRASS('r.grow.distance', input=paste0(species, '_habitatAquatic'), distance=paste0(species, '_distanceFromHabitatAquatic1'), flags=c('overwrite'))
    #Reclassify distance from aquatic habitat based on species suitability
    write.table(c('0 thru 200=100', '200 thru 600=70', '600 thru 200000=59', '*=NULL'), 'rule.txt', sep="", col.names=FALSE, quote=FALSE, row.names=FALSE)
    execGRASS('r.reclass', input=paste0(species, '_distanceFromHabitatAquatic1'), output=paste0(species, '_distanceFromHabitatAquatic'), rules='rule.txt', flags=c('overwrite'))
    
    # #Identify terrestrial habitat within dispersal distance from aquatic habitat    
    # execGRASS('r.mapcalc',expression=paste0(species,'_habitat_terrestrial_close_to_wet=',species,'_habitat_terrestrial_potential*',species,'_distance_from_habitat_aquatic_reclass/100'),flags=c('overwrite'))
    # write.table(c('0 thru 59=NULL','60 thru 100=1'),'rule.txt',sep="",col.names=FALSE,quote=FALSE,row.names=FALSE)
    # execGRASS('r.reclass',input=paste0(species,'_habitat_terrestrial_close_to_wet'),output=paste0(species,'_habitat_terrestrial'),rules='rule.txt',flags=c('overwrite'))
    
    # Combine habitat suitability and distance to aquatic habitat
    execGRASS('r.mapcalc', expression=paste0(species, '_habitatSuitability=', species, '_habitatSuitabilityPotential*', species, '_distanceFromHabitatAquatic/100'), flags=c('overwrite'))
    
    # #Combine aquatic and terrestrial habitats 
    # execGRASS('r.patch',input=paste0(species,'_habitat_aquatic,',species,'_habitat_terrestrial'),output=paste0(species,'_habitat'),flags=c('overwrite'))
  }
  
  # Identify suitable habitat (habitat suitability > 60)
  suitabilityThreshold <- speciesHabitatPatchParams$SuitabilityThreshold[speciesHabitatPatchParams$Species==species]
  write.table(c(paste0('0 thru ', suitabilityThreshold-1, '=NULL'), paste0(suitabilityThreshold, ' thru 100=1')),'rule.txt',sep="",col.names=FALSE,quote=FALSE,row.names=FALSE)
  execGRASS('r.reclass',input=paste0(species,'_habitatSuitability'),output=paste0(species,'_habitat'),rules='rule.txt',flags=c('overwrite'))
  # Replace NAs with 0's to make a binary habitat layer
  execGRASS('r.mapcalc', expression=paste0(species, '_habitatBinary=if(isnull(', species, '_habitat),','0,',species, '_habitat)'), flags=c('overwrite'))
  
  
  #############################
  # Habitat patch delineation #
  #############################
  if(species %in% speciesFillGaps){
    #Keep a copy of habitat map with gaps
    execGRASS('r.mapcalc', expression=paste0(species, '_habitatWithGaps', '=', species, '_habitat'), flags=c('overwrite'))
    
    # Distance from habitat patches
    execGRASS('r.grow.distance', input=paste0(species, '_habitat'), distance=paste0(species, '_habitatDistance'), flags=c('overwrite'))
    
    # High-pass filter: moving window to identify gaps in forest based on gap-smoothing parameter
    gap <- as.numeric(subset(speciesHabitatPatchParams, Species==species, select=Gap))
    # Moving window size is the diameter of a circle (number of pixels). It must be an odd number.
    neighborhoodSize <- gap / myResolution + 1
    execGRASS('r.neighbors', input=paste0(species, '_habitatDistance'), output=paste0(species, '_habitatSmooth'), size=neighborhoodSize, method='maximum', flags=c('overwrite','c'))
    
    # Reclassify pixels in gaps to create bridges across gaps
    #Note that when myResolution=30 use this but for myResolution=10 we used the commented line
    write.table(c(paste0('0 thru ', gap/2, '=1'), paste0(gap/2+1, ' thru 20000000=NULL')), 'rule.txt', sep="", col.names=FALSE, quote=FALSE, row.names=FALSE)
    # write.table(c(paste0('0 thru ', gap/2-1, '=1'), paste0(gap/2, ' thru 20000000=NULL')), 'rule.txt', sep="", col.names=FALSE, quote=FALSE, row.names=FALSE)
    execGRASS('r.reclass', input=paste0(species, '_habitatSmooth'), output=paste0(species, '_habitatSmoothReclass'), rules='rule.txt', flags=c('overwrite'))
    execGRASS('r.mapcalc', expression=paste0(species, '_habitatSmoothReclass=', species, '_habitatSmoothReclass'), flags=c('overwrite'))
    
    # Combine suitable patches and bridges across gaps
    execGRASS('r.patch', input=paste0(species,'_habitatWithGaps,', species, '_habitatSmoothReclass'), output=paste0(species, '_habitat'), flags=c('overwrite'))
  }
  
  if(species == 'RASY'){
    # Assign patch ids to habitat layer
    execGRASS('r.clump', input=paste0(species, '_habitat'), output=paste0(species, '_habitatClump'), flags=c('overwrite'))
    # Calculate how much wetland is in each habitat patch
    execGRASS('r.stats.zonal', base=paste0(species, '_habitatClump'), cover='wetlands', method='count', output=paste0(species, '_habitatWetlandsCount'), flags=c('overwrite'))
    # Suitable habitat must have at least 3 pixels of wetland within its borders
    write.table(c('0 thru 2=NULL','3 thru 1000000=1'),'rule.txt', sep="", col.names=FALSE, quote=FALSE, row.names=FALSE)
    execGRASS('r.reclass', input=paste0(species, '_habitatWetlandsCount'), output=paste0(species, '_habitatWetlandInside'), rules='rule.txt', flags=c('overwrite'))
    
    # Identify patches below minimum area threshold
    minPatchSize<-as.numeric(subset(speciesHabitatPatchParams, Species==species, select=MinPatchSize))
    execGRASS('r.reclass.area', input=paste0(species, '_habitatWetlandInside'), output=paste0(species, '_habitatSmall'), value=minPatchSize-0.01, mode='lesser', flags=c('overwrite', 'd'))
    
    # Identify patches above minimum area threshold
    execGRASS('r.reclass.area', input=paste0(species, '_habitatWetlandInside'), output=paste0(species, '_habitatLarge'), value=minPatchSize, mode='greater', flags=c('overwrite', 'd'))
  }
  else{
    # Identify patches below minimum area threshold
    minPatchSize<-as.numeric(subset(speciesHabitatPatchParams, Species==species, select=MinPatchSize))
    execGRASS('r.reclass.area', input=paste0(species,'_habitat'), output=paste0(species,'_habitatSmall'), value=minPatchSize-0.01, mode='lesser', flags=c('overwrite','d'))
    
    # Identify patches above minimum area threshold
    execGRASS('r.reclass.area', input=paste0(species,'_habitat'), output=paste0(species,'_habitatLarge'), value=minPatchSize, mode='greater', flags=c('overwrite','d'))
  }
  
  # Replace NAs with 0's to make a binary habitat layer
  execGRASS('r.mapcalc', expression=paste0(species, '_habitatLargeBinary=if(isnull(', species, '_habitatLarge),','0,',species, '_habitatLarge)'), flags=c('overwrite'))
  
  
  ########################
  # Landscape resistance #
  ########################
  # Reclass: patches above minimum patch=900; patches below minimum patch size=910
  habitatTooSmallCode <- speciesResistanceReclass$LandcoverCode[speciesResistanceReclass$LandcoverName=="HabitatTooSmall"]
  write.table(c(paste0('1=',habitatTooSmallCode),'*=NULL'), 'rule.txt', sep="", col.names=FALSE, quote=FALSE, row.names=FALSE)
  execGRASS('r.reclass', input=paste0(species,'_habitatSmall'), output=paste0(species,'_habitatSmallReclass'), rules='rule.txt', flags=c('overwrite'))
  habitatCode <- speciesResistanceReclass$LandcoverCode[speciesResistanceReclass$LandcoverName=="Habitat"]
  write.table(c(paste0('1=',habitatCode),'*=NULL'), 'rule.txt', sep="", col.names=FALSE, quote=FALSE, row.names=FALSE)
  execGRASS('r.reclass', input=paste0(species,'_habitatLarge'), output=paste0(species,'_habitatLargeReclass'), rules='rule.txt', flags=c('overwrite'))
  
  # Impose suitable habitat on landcover map
  execGRASS('r.patch', input=paste0(species, '_habitatLargeReclass,', species, '_habitatSmallReclass'), output=paste0(species,'_habitatReclass'), flags=c('overwrite'))
  execGRASS('r.patch', input=paste0(species, '_habitatReclass,landcover'), output=paste0(species,'_landcoverResistance'), flags=c('overwrite'))
  
  # Reclassify landcover with resistance values
  rcl<-read.csv(paste0(b01b02RawTablesDir,"/speciesResistanceReclass.csv"),header=TRUE)[,c('LandcoverCode',species)]
  write.table(paste0(rcl[,'LandcoverCode'],'=',rcl[,species]), 'rule.txt', sep="", col.names=FALSE, quote=FALSE, row.names=FALSE)
  execGRASS('r.reclass', input=paste0(species, '_landcoverResistance'), output=paste0(species,'_resistance'), rules='rule.txt', flags=c('overwrite'))
  
  ###########################################
  # Save outputs as geotiff for full region #
  ###########################################
  execGRASS('r.out.gdal', input=paste0(species, '_habitatSuitability'), output=paste0(b03habitatDir, "/", species, '_habitatSuitability_', myResolution, 'm.tif'), format='GTiff',createopt='COMPRESS=LZW', flags=c('overwrite'))
  execGRASS('r.out.gdal', input=paste0(species, '_habitatLargeBinary'), output=paste0(b03habitatDir, "/", species, '_habitatPatch_', myResolution, 'm.tif'), format='GTiff',createopt='COMPRESS=LZW', flags=c('overwrite'))
  execGRASS('r.out.gdal', input=paste0(species, '_resistance'), output=paste0(b03resistanceDir, "/", species, '_resistance_', myResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))
  
  ################################################################################
  # Save habitat patch and resistance map at coarser resolutions for full region #
  ################################################################################
  # 240m
  # Set the geographic region
  newResolution<-240
  execGRASS('g.region', n='403770', e='-391380', w='-691380', s='117930', res=paste0(newResolution))
  
  # Save maps as geotiffs
  execGRASS('r.out.gdal', input=paste0(species,'_habitatLargeBinary'), output=paste0(b03habitatDir, "/", species, '_habitatPatch_', newResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))
  execGRASS('r.out.gdal', input=paste0(species,'_resistance'), output=paste0(b03resistanceDir, '/', species,'_resistance_', newResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))
  
  # 480m
  # Set the geographic region
  newResolution<-480
  execGRASS('g.region', n='404010', e='-391380', w='-691380', s='117930', res=paste0(newResolution))
  
  # Save maps as geotiffs
  execGRASS('r.out.gdal', input=paste0(species,'_habitatLargeBinary'), output=paste0(b03habitatDir, "/", species,'_habitatPatch_', newResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))
  execGRASS('r.out.gdal', input=paste0(species,'_resistance'), output=paste0(b03resistanceDir, "/", species,'_resistance_', newResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))
  
  # Reset the geographic region
  execGRASS('g.region', n='403680', e='-391410', w='-691380', s='117930', res=paste0(myResolution))
  
  
  ###############################
  # Rescale habitat suitability #
  ###############################
  # call habitatSuitability layers from b03habitatDir and apply rescale function
  suitabilityMap_01 <- raster(paste0(b03habitatDir, "/", species, "_habitatSuitability_", myResolution, 'm.tif')) %>%
    calc(., fun=rescaleR)
  
  # Save maps as geotiffs
  writeRaster(suitabilityMap_01, paste0(b03habitatDir, "/", species, "_habitatSuitability_01_", myResolution, 'm.tif'), overwrite=TRUE)   
  
  # Add maps to mapset (to clip to focal region next)
  execGRASS("r.in.gdal", input=paste0(b03habitatDir, "/", species, "_habitatSuitability_01_", myResolution, 'm.tif'), output=paste0(species, "_habitatSuitability_01"), flags=c("overwrite", "o"))
  
  #############################
  # Clip maps to focal region # 
  #############################
  # 30m
  # Apply mask of focal region
  execGRASS('r.mask', raster="studyArea0",flags=c("overwrite"))
  # Clip maps to mask
  execGRASS('r.mapcalc', expression=paste0(species, "_habitatSuitability_Focal_", myResolution, "m=", species, "_habitatSuitability * 1"), flags=c("overwrite"))
  execGRASS('r.mapcalc', expression=paste0(species, "_habitatSuitability_01_Focal_", myResolution, "m=", species, "_habitatSuitability_01 * 1"), flags=c("overwrite"))
  execGRASS('r.mapcalc', expression=paste0(species, "_habitatPatch_Focal_", myResolution, "m=", species, "_habitatLargeBinary * 1"), flags=c("overwrite"))
  execGRASS('r.mapcalc', expression=paste0(species, "_resistance_Focal_", myResolution, "m=", species, "_resistance * 1"), flags=c("overwrite"))
  # Save clipped maps as geotiffs
  execGRASS('r.out.gdal', input=paste0(species, "_habitatSuitability_Focal_", myResolution, "m"), output=paste0(b03habitatDir, '/', species, '_habitatSuitability_Focal_', myResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))
  execGRASS('r.out.gdal', input=paste0(species, "_habitatSuitability_01_Focal_", myResolution, "m"), output=paste0(b03habitatDir, '/', species, '_habitatSuitability_01_Focal_', myResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))
  execGRASS('r.out.gdal', input=paste0(species, "_habitatPatch_Focal_", myResolution, "m"), output=paste0(b03habitatDir, '/', species, '_habitatPatch_Focal_', myResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))
  execGRASS('r.out.gdal', input=paste0(species, "_resistance_Focal_", myResolution, "m"), output=paste0(b03resistanceDir, '/', species, '_resistance_Focal_', myResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))
  # Remove mask from mapset
  execGRASS('r.mask', raster="studyArea0",flags=c("r"))
  
  # 240m
  # Set the geographic region
  newResolution<-240
  execGRASS('g.region', n='403770', e='-391380', w='-691380', s='117930', res=paste0(newResolution))
  # Apply mask of focal region
  execGRASS('r.mask', raster="studyArea0",flags=c("overwrite"))
  # Clip maps to mask
  execGRASS('r.mapcalc', expression=paste0(species, "_habitatPatch_Focal_", newResolution, "m=", species, "_habitatLargeBinary * 1"), flags=c("overwrite"))
  execGRASS('r.mapcalc', expression=paste0(species, "_resistance_Focal_", newResolution, "m=", species, "_resistance * 1"), flags=c("overwrite"))
  # Save clipped maps as geotiffs
  execGRASS('r.out.gdal', input=paste0(species, "_habitatPatch_Focal_", newResolution, "m"), output=paste0(b03habitatDir, '/', species, '_habitatPatch_Focal_', newResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))
  execGRASS('r.out.gdal', input=paste0(species, "_resistance_Focal_", newResolution, "m"), output=paste0(b03resistanceDir, '/', species, '_resistance_Focal_', newResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))
  # Remove mask from mapset
  execGRASS('r.mask', raster="studyArea0",flags=c("r"))
  
  # 480m
  # Set the geographic region
  newResolution<-480
  execGRASS('g.region', n='404010', e='-391380', w='-691380', s='117930', res=paste0(newResolution))
  # Apply mask of focal region
  execGRASS('r.mask', raster="studyArea0",flags=c("overwrite"))
  # Clip maps to mask
  execGRASS('r.mapcalc', expression=paste0(species, "_habitatPatch_Focal_", newResolution, "m=", species, "_habitatLargeBinary * 1"), flags=c("overwrite"))
  execGRASS('r.mapcalc', expression=paste0(species, "_resistance_Focal_", newResolution, "m=", species, "_resistance * 1"), flags=c("overwrite"))
  # Save clipped maps as geotiffs
  execGRASS('r.out.gdal', input=paste0(species, "_habitatPatch_Focal_", newResolution, "m"), output=paste0(b03habitatDir, '/', species, '_habitatPatch_Focal_', newResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))
  execGRASS('r.out.gdal', input=paste0(species, "_resistance_Focal_", newResolution, "m"), output=paste0(b03resistanceDir, '/', species, '_resistance_Focal_', newResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))
  # Remove mask from mapset
  execGRASS('r.mask', raster="studyArea0",flags=c("r"))
  
  # Reset the geographic region
  execGRASS('g.region', n='403680', e='-391410', w='-691380', s='117930', res=paste0(myResolution))
  
  ########################################
  # Filter patch sizes within focal area #
  ########################################
  #  - apply patch filter to size of wetland patches for RASY?
  # if(species == 'RASY'){
  #   # Assign patch ids to habitat layer
  #   execGRASS('r.clump', input=paste0(species, '_habitat'), output=paste0(species, '_habitatClump'), flags=c('overwrite'))
  #   # Calculate how much wetland is in each habitat patch
  #   execGRASS('r.stats.zonal', base=paste0(species, '_habitatClump'), cover='wetlands', method='count', output=paste0(species, '_habitatWetlandsCount'), flags=c('overwrite'))
  #   # Suitable habitat must have at least 3 pixels of wetland within its borders
  #   write.table(c('0 thru 2=NULL','3 thru 1000000=1'),'rule.txt', sep="", col.names=FALSE, quote=FALSE, row.names=FALSE)
  #   execGRASS('r.reclass', input=paste0(species, '_habitatWetlandsCount'), output=paste0(species, '_habitatWetlandInside'), rules='rule.txt', flags=c('overwrite'))
  #   
  #   # Identify patches below minimum area threshold
  #   minPatchSize<-as.numeric(subset(speciesHabitatPatchParams, Species==species, select=MinPatchSize))
  #   execGRASS('r.reclass.area', input=paste0(species, '_habitatWetlandInside'), output=paste0(species, '_habitatSmall'), value=minPatchSize-0.01, mode='lesser', flags=c('overwrite', 'd'))
  #   
  #   # Identify patches above minimum area threshold
  #   execGRASS('r.reclass.area', input=paste0(species, '_habitatWetlandInside'), output=paste0(species, '_habitatLarge'), value=minPatchSize, mode='greater', flags=c('overwrite', 'd'))
  # }
  # else{
  
  # Identify patches below minimum area threshold
  minPatchSize<-as.numeric(subset(speciesHabitatPatchParams, Species==species, select=MinPatchSize))
  execGRASS('r.reclass.area', input=paste0(species, "_habitatPatch_Focal_", myResolution, "m"), output=paste0(species, "_habitatSmall_Focal_", myResolution, "m"), value=minPatchSize-0.01, mode='lesser', flags=c('overwrite','d'))
  # Identify patches above minimum area threshold
  execGRASS('r.reclass.area', input=paste0(species, "_habitatPatch_Focal_", myResolution, "m"), output=paste0(species, "_habitatLarge_Focal_", myResolution, "m"), value=minPatchSize, mode='greater', flags=c('overwrite','d'))

  # Replace NAs with 0's to make a binary habitat layer
  execGRASS('r.mapcalc', expression=paste0(species, '_habitatPatch_Focal_Filtered_', myResolution, 'm=if(isnull(', species, '_habitatLarge_Focal_', myResolution,'m),0,',species, '_habitatLarge_Focal_', myResolution, 'm)'), flags=c('overwrite'))
  
  # Save maps as geotiffs
  # Apply mask of focal region
  execGRASS('r.mask', raster="studyArea0",flags=c("overwrite"))
  execGRASS('r.out.gdal', input=paste0(species, "_habitatPatch_Focal_Filtered_", myResolution, "m"), output=paste0(b03habitatDir, '/', species, '_habitatPatch_Focal_Filtered_', myResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))
  # Remove mask from mapset
  execGRASS('r.mask', raster="studyArea0",flags=c("r"))
  
  # Save filtered habitat patch at coarser resolutions 
  # 240m
  # Set the geographic region
  newResolution<-240
  execGRASS('g.region', n='403770', e='-391380', w='-691380', s='117930', res=paste0(newResolution))
  # Apply mask of focal region
  execGRASS('r.mask', raster="studyArea0",flags=c("overwrite"))
  execGRASS('r.out.gdal', input=paste0(species, "_habitatPatch_Focal_Filtered_", myResolution, "m"), output=paste0(b03habitatDir, '/', species, '_habitatPatch_Focal_Filtered_', newResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))
  # Remove mask from mapset
  execGRASS('r.mask', raster="studyArea0",flags=c("r"))
  
  # 480m
  # Set the geographic region
  newResolution<-480
  execGRASS('g.region', n='404010', e='-391380', w='-691380', s='117930', res=paste0(newResolution))
  # Apply mask of focal region
  execGRASS('r.mask', raster="studyArea0",flags=c("overwrite"))
  execGRASS('r.out.gdal', input=paste0(species, "_habitatPatch_Focal_Filtered_", myResolution, "m"), output=paste0(b03habitatDir, '/', species, '_habitatPatch_Focal_Filtered_', newResolution, 'm.tif'), format='GTiff', createopt='COMPRESS=LZW', flags=c('overwrite'))
  # Remove mask from mapset
  execGRASS('r.mask', raster="studyArea0",flags=c("r"))
  
  # Reset the geographic region
  execGRASS('g.region', n='403680', e='-391410', w='-691380', s='117930', res=paste0(myResolution))
}

#############################################
# Generate habitat patch summary statistics #
#############################################
# Get the total area of the study region
studyareaRaster <- raster(file.path(b03ProcessedMapsDir, studyAreaName))

focalCount <- freq(studyareaRaster) %>% 
              as_tibble %>% 
              filter(value == 1) %>% 
              pull(count)

focalArea <- focalCount * 0.09

# Calculate patch statistics and summarize in a table
patchSummaryStatistics <- map_dfr(speciesList,
~{ 
  # Get habitat patch and suitability maps
  habitatPatchRaster <- raster(paste0(b03habitatDir, '/', .x, '_habitatPatch_Focal_Filtered_', myResolution, 'm.tif'))
  habitatSuitabilityRaster <- raster(paste0(b03habitatDir, '/', .x, "_habitatSuitability_Focal_", myResolution, 'm.tif'))
  
  # Calculate total patch area
  patchCount <- freq(habitatPatchRaster)
  patchArea <- patchCount[[2,2]] * 0.09
  
  # Calculate percent of suitable habitat
  percentSuitableHabitat <- (patchArea / focalArea) * 100
  
  # Calculate the number of unique patches
  patchIds <- clump(habitatPatchRaster)
  numPatches <- length(unique(patchIds))
  
  # Calculate the mean patch area
  meanPatchArea <- patchArea / numPatches
  
  # Calculate the mean habitat quality in patches
  meanHabitatQuality <- zonal(x = habitatSuitabilityRaster,
                              z = habitatPatchRaster,
                              fun = 'mean') %>%
                        as_tibble %>% 
                        filter(zone == 1) %>% 
                        pull(mean)
  
  output <- tibble(
    Species = .x,
    'Suitable Habitat (%)' = percentSuitableHabitat,
    'Total Habitat Patch Area (Ha)' = patchArea,
    'Number of Habitat Patches' = numPatches,
    'Mean Habitat Patch Area (Ha)' = meanPatchArea,
    'Mean Habitat Quality' = meanHabitatQuality
  )
})

# Save tabular summary to disk
write_csv(patchSummaryStatistics, file.path(b03Dir, "model-outputs", "tabular", "patch-summary-statistics-focal-30m.csv"))
