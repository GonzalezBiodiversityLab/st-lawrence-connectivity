#####################################################################
# a254      
# Run Zonation for all focal species and for individual focal species 
# 03-2022                                       					
#     
#   Uses Zonation
#
#   Inputs (for focal species):
#    -habitat suitability, short and long betweenness,
#    -short and long dEC, current density
#    -protected areas
#  
#   Outputs:
#    -Zonation outputs including priority natural areas
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
conservationCriteriaList<-data.frame(Weight=c(1, 1, 0.5, 0.5, 1),
                                     FileName=c("habitatSuitabilityFocal_NatAreas_01_Z.tif",
                                                "betweennessShortFocal_NatAreas_Log01_Z.tif",
                                                "dECGapFocal_NatAreas_01_Z.tif",
                                                "dECNatalFocal_NatAreas_01_Z.tif",
                                                "currentFlowFocal_NatAreas_Log01_Z.tif"),
                                     Directory=rep(file.path(b03zonationDir, "ZonationInputs"),5))

# Names of additional layers
protectedAreasName <- file.path(b03zonationDir, "Z_protectedAreas.tif")

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
# Read in Zonation settings template
# NB, this file is used for the ALL species analysis and for the individual species analyses
zonationSet <- readLines(zonationTemplateSettingsFileName)
# Overwrite template with user-input values
zonationSet[2] <- paste("removal rule =", removalRuleCode)
zonationSet[3] <- paste("warp factor =", warpFactor)
zonationSet[4] <- paste("edge removal =", edgeRemoval)
zonationSet[5] <- paste("add edge points =", numberEdgePointsToAdd)
zonationSet[15] <- paste("mask file =", protectedAreasName)
writeLines(zonationSet, zonationSettingsFileName)

# Run for all species combined -----------------------
# Create directory if it doesn't already exist
dir.create(file.path(b03zonationDir, "ALL"), showWarnings = FALSE)

# Zonation file paths
zonationBiodiversityFeaturesFileName <- file.path(b03zonationDir, "ALL", "BiodiversityFeatureList.spp")
zonationOutputFileName <- file.path(b03zonationDir, "ALL", "ALL_Zonation.txt")
zonationRunFileName <- file.path(b03zonationDir, "ALL", "ZonationRun.bat")

# Make Zonation input file of biodiversity features
zonationSpp <- tibble(X1 = double(),X2 = integer(),X3=integer(),X4=integer(),X5=integer(),X6=character())
for(i in 1:length(speciesList)){
  species<-speciesList[i]
  conservationCriteriaFilePaths <- file.path(conservationCriteriaList$Directory, paste(species, conservationCriteriaList$FileName, sep="_"))
  zonationSpp <- bind_rows(zonationSpp, data.frame(X1=conservationCriteriaList$Weight, X2=0, X3=1, X4=1, X5=1, X6=conservationCriteriaFilePaths))
}  
write.table(zonationSpp, zonationBiodiversityFeaturesFileName, row.names = FALSE, col.names=FALSE)

# Make Zonation run file
zonationRun <- paste0("-r ", "\"", zonationSettingsFileName, "\" ", "\"", zonationBiodiversityFeaturesFileName, "\" ", "\"", zonationOutputFileName, "\" ", "0.0 0 1.0 0")
writeLines(zonationRun, zonationRunFileName)

# Run Zonation
system(paste0("\"", zonationExeFileName, "\" ", zonationRun))


# Run for each species individually --------------------
# Note that we do not need to change the Zonation settings file
for(i in 1:length(speciesList)){
  species<-speciesList[i]
  
  # Create directories if they don't already exist
  dir.create(file.path(b03zonationDir, species), showWarnings = FALSE)
  
  # Zonation file paths 
  zonationBiodiversityFeaturesFileName <- file.path(b03zonationDir, species, "BiodiversityFeatureList.spp")
  zonationOutputFileName <- file.path(b03zonationDir, species, paste(species, "Zonation.txt", sep="_"))
  zonationRunFileName <- file.path(b03zonationDir, species, "ZonationRun.bat")
  
  # Make Zonation input file of biodiversity features
  conservationCriteriaFilePaths <- file.path(conservationCriteriaList$Directory, paste(species, conservationCriteriaList$FileName, sep="_"))
  zonationSpp <- tibble(X1=conservationCriteriaList$Weight, X2=0, X3=1, X4=1, X5=1, X6=conservationCriteriaFilePaths)
  write.table(zonationSpp, zonationBiodiversityFeaturesFileName, row.names = FALSE, col.names=FALSE)
  
  # Make Zonation run file
  zonationRun <- paste0("-r ", "\"", zonationSettingsFileName, "\" ", "\"", zonationBiodiversityFeaturesFileName, "\" ", "\"", zonationOutputFileName, "\" ", "0.0 0 1.0 0")
  writeLines(zonationRun, zonationRunFileName)
  
  # Run Zonation
  system(paste0("\"", zonationExeFileName, "\" ", zonationRun))
}
