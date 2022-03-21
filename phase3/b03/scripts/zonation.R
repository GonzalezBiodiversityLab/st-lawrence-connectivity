
# Load constants, functions, etc
source("./b03/scripts/0-constants.R")

options(stringsAsFactors = FALSE)
Sys.setenv(TZ='GMT')


# Input parameters
numIterations <- 10

# List of Zonation feature layers
conservationCriteriaList<-data.frame(Weight=c(1, 1, 1),
                                     FileName=c("HabitatSuitability",
                                                "HabitatPatch",
                                                "OMNI_cum_curmap"),
                                     Directory=rep(file.path(b03zonationDir, "ZonationInputs"), 3))

# Names of additional layers
protectedAreasName <- file.path(b03ProcessedMapsDir, "SpatialMultiplier_ProtectedAreas_b03Extent.tif")
natAreasMaskFilename <-file.path(b03ProcessedMapsDir, "b03-natural-areas-2010.tif")

# Zonation parameters
zonationExeFileName <- "C:\\Users\\Administrator\\AppData\\Local\\zonation 4.0.0rc1_compact\\bin\\zig4.exe"
# Removal rule
# Reference table of removal rules
# Code                      Name
# 1     Basic core-area Zonation
# 2    Additive benefit function
# 3        Target based planning
# 4 Generalized benefit function
# 5                       Random 
removalRuleCode <- 2
warpFactor <- 1000
edgeRemoval <- TRUE  # 1=TRUE, 0=FALSE
numberEdgePointsToAdd <- 0
# Zonation file paths
zonationTemplateSettingsFileName <- file.path(b03ProcessedTabularDir, "zonation-settings-template.dat")
zonationSettingsFileName <- file.path(b03zonationDir, "ZonationInputs", "RunSetting.dat")

# Create Zonation settings file --------------------------
natAreasMask <- raster(natAreasMaskFilename)

protectedAreas <-raster(protectedAreasName)
protectedAreasNatAreas <- protectedAreas * natAreasMask
writeRaster(protectedAreasNatAreas, file.path(b03zonationDir, "ZonationInputs", "Z_protectedAreas.tif"))

# Read in Zonation settings template
# NB, this file is used for the ALL species analysis and for the individual species analyses
zonationSet <- readLines(zonationTemplateSettingsFileName)
# Overwrite template with user-input values
zonationSet[2] <- paste("removal rule =", removalRuleCode)
zonationSet[3] <- paste("warp factor =", warpFactor)
zonationSet[4] <- paste("edge removal =", edgeRemoval)
zonationSet[5] <- paste("add edge points =", numberEdgePointsToAdd)
zonationSet[15] <- paste("mask file =", file.path(b03zonationDir, "ZonationInputs", "Z_protectedAreas.tif"))
writeLines(zonationSet, zonationSettingsFileName)

# Run for all species combined -----------------------
# Create directory if it doesn't already exist
dir.create(file.path(b03zonationDir, "ALL"), showWarnings = FALSE)


# Zonation file paths
zonationBiodiversityFeaturesFileName <- file.path(b03zonationDir, "ALL", "BiodiversityFeatureList.spp")
zonationOutputFileName <- file.path(b03zonationDir, "ALL", "ALL_Zonation.txt")
zonationRunFileName <- file.path(b03zonationDir, "ALL", "ZonationRun.bat")

# Make Zonation input file of biodiversity features
# Create difference map for habitat suitability and current density for each species
for(i in 1:length(speciesList)){
    species<-speciesList[i]
    habSuit2010 <- raster(file.path(b03Libraries, "BTSL_stconnect.ssim.output", "Scenario-218", "stconnect_HSOutputHabitatSuitability", paste0("HabitatSuitability.", species, ".it1.ts2010.tif")))
    habSuit2110 <- raster(file.path(b03Libraries, "BTSL_stconnect.ssim.output", "Scenario-218", "stconnect_HSOutputHabitatSuitability", paste0("HabitatSuitability.", species, ".it1.ts2110.tif")))
    habSuit2010NatAreas <- habSuit2010  * natAreasMask
    habSuit2110NatAreas <- habSuit2110  * natAreasMask    
    diffHabSuit <- abs(habSuit2110NatAreas - habSuit2010NatAreas)
#    diffHabSuit01 <-rescaleR(diffHabSuit, 0, 1)

    currentDensity2010 <- raster(file.path(b03Libraries, "BTSL_stconnect.ssim.output", "Scenario-220", "stconnect_CCOutputCumulativeCurrent", paste0("OMNI_cum_curmap.", species, ".it1.ts2010.tif")))
    currentDensity2110 <- raster(file.path(b03Libraries, "BTSL_stconnect.ssim.output", "Scenario-220", "stconnect_CCOutputCumulativeCurrent", paste0("OMNI_cum_curmap.", species, ".it1.ts2110.tif")))
    currentDensity2010NatAreas <- currentDensity2010  * natAreasMask
    currentDensity2110NatAreas <- currentDensity2110  * natAreasMask    
    diffCurrentDensity <- abs(currentDensity2110NatAreas - currentDensity2010NatAreas)

    writeRaster(diffHabSuit, file.path(b03zonationDir, "ZonationInputs", paste0(species, "_HabitatSuitabilityChange2010-2110.tif")), overwrite=TRUE)    
    writeRaster(diffCurrentDensity, file.path(b03zonationDir, "ZonationInputs", paste0(species, "_CurrentDensityChange2010-2110.tif")), overwrite=TRUE)    
}



zonationSpp <- tibble(X1 = double(),X2 = integer(),X3=integer(),X4=integer(),X5=integer(),X6=character())
for(i in 1:length(speciesList)){
  species<-speciesList[i]
  for(iteration in 1:numIterations){
    habSuit2010 <- raster(file.path(b03Libraries, "BTSL_stconnect.ssim.output", "Scenario-218", "stconnect_HSOutputHabitatSuitability", paste0("HabitatSuitability.", species, ".it", iteration, ".ts2010.tif")))
    habSuit2010NatAreas <- habSuit2010  * natAreasMask
    writeRaster(habSuit2010NatAreas, file.path(b03zonationDir, "ZonationInputs", paste0("HabitatSuitability.", species, ".it", iteration, ".ts2010.tif")), overwrite = TRUE)    
    
    habPatch2010 <- raster(file.path(b03Libraries, "BTSL_stconnect.ssim.output", "Scenario-218", "stconnect_HSOutputHabitatPatch", paste0("HabitatPatch.", species, ".it", iteration, ".ts2010.tif")))
    habPatch2010NatAreas <- habPatch2010  * natAreasMask
    writeRaster(habPatch2010NatAreas, file.path(b03zonationDir, "ZonationInputs", paste0("HabitatPatch.", species, ".it", iteration, ".ts2010.tif")), overwrite = TRUE)    
    
    currentDensity2010 <- raster(file.path(b03Libraries, "BTSL_stconnect.ssim.output", "Scenario-220", "stconnect_CCOutputCumulativeCurrent", paste0("OMNI_cum_curmap.", species, ".it", iteration, ".ts2010.tif")))
    currentDensity2010NatAreas <- currentDensity2010  * natAreasMask
    writeRaster(currentDensity2010NatAreas, file.path(b03zonationDir, "ZonationInputs", paste0("OMNI_cum_curmap.", species, ".it", iteration, ".ts2010.tif")), overwrite = TRUE)    
    
    conservationCriteriaFilePaths <- file.path(conservationCriteriaList$Directory, paste0(conservationCriteriaList$FileName, ".", species, ".it", iteration, ".ts2010.tif"))
    zonationSpp <- bind_rows(zonationSpp, data.frame(X1=conservationCriteriaList$Weight, X2=0, X3=1, X4=1, X5=1, X6=conservationCriteriaFilePaths))
  }
  zonationSpp <- bind_rows(zonationSpp, data.frame(X1=0.5, X2=0, X3=1, X4=1, X5=1, X6=file.path(b03zonationDir, "ZonationInputs", paste0(species, "_HabitatSuitabilityChange2010-2110.tif"))))
  zonationSpp <- bind_rows(zonationSpp, data.frame(X1=0.5, X2=0, X3=1, X4=1, X5=1, X6=file.path(b03zonationDir, "ZonationInputs", paste0(species, "_CurrentDensityChange2010-2110.tif"))))
}

write.table(zonationSpp, zonationBiodiversityFeaturesFileName, row.names = FALSE, col.names=FALSE)

# Make Zonation run file
zonationRun <- paste0("-r ", "\"", zonationSettingsFileName, "\" ", "\"", zonationBiodiversityFeaturesFileName, "\" ", "\"", zonationOutputFileName, "\" ", "0.0 0 1.0 0")
writeLines(zonationRun, zonationRunFileName)

# Run Zonation
system(paste0("\"", zonationExeFileName, "\" ", zonationRun))

