# a254
# Bronwyn Rayfield, ApexRMS
#
# This script generates transition multiplier values for three climate values 
# based off of forest age/composisiont output from Landis-II runs

## Workspace ----
options(tibble.width = Inf, tibble.print_min = Inf)

# Input files and folders
SyncroSimDir <- "C:/Program Files/SyncroSim"
modelFile <- "BTSL_stconnect-Fix-11-3-2022/BTSL_stconnect.ssim"

# Save Datasheets to Library
# Set to T to save datasheets back to library
saveDatasheets <- T  

# Setting up st-sim library and project
mySession <- session()
myLibrary <- ssimLibrary(file.path(b03ProcessedTabularDir, modelFile), session=mySession)
myProject <- project(myLibrary, "Definitions")  # Assumes there is only one default project per library

## Read in state classes and create forest transitions ----
# Read in stateClasses
stateClasses <- datasheet(myProject, "stsim_StateClass")

# Create forest transition classes - named as Class:Subclass
srcClassesForest <- stateClasses[stateClasses$StateLabelXID == "Foret",]
srcClassesForest$ID <- as.numeric(srcClassesForest$ID)

# Setup a dataframe with all combinations of to/from forest state classes
transitionTypesForest <- srcClassesForest[rep(seq_len(nrow(srcClassesForest)), each = nrow(srcClassesForest)), ]
for (i in 1:nrow(srcClassesForest)) {
  for (j in 1:nrow(srcClassesForest)) {
    row <- ((i-1)*nrow(srcClassesForest)) + j
    # print(row)
    transitionTypesForest$destClass[row] <- srcClassesForest[j,]$Name
    transitionTypesForest$destID[row] <- srcClassesForest[j,]$ID
  }
}

# Remove the rows with same source and destination Class
transitionTypesForest <- transitionTypesForest[transitionTypesForest$Name != transitionTypesForest$destClass,]

# Add columns with Name and ID for Transition Types. Note that TransitionTypeID = (srcID * 1000) + destID 
transitionTypesForest$Transition <- paste0(transitionTypesForest$Name, "->", transitionTypesForest$destClass)
transitionTypesForest$TransitionID <- transitionTypesForest$ID * 1000 + transitionTypesForest$destID


## Transition Multipliers - Baseline Climate ----
# Scenario-specific inputs
scenarioName <- "Transition Multipliers: Climate Baseline" # Name of st-sim scenario for transition multipliers
landisName <- "landis_baseline_age_transitions.csv" #Name of Landis csv from Guillaume

# Read in 20-year transition multiplier probabilities from Guillaume
transMultRaw <- read_csv(file.path(b03RawTablesDir, landisName))

# Merge Landis output file with stsim state classes
transMultMergedstsim <- mutate(transMultRaw, srcClassID = (500 + from_cover*10 + age_before), destClassID = (500 + to_cover * 10 + age_after)) %>%
  left_join(stateClasses, by=c("srcClassID"="ID")) %>%
  left_join(stateClasses, by=c("destClassID"="ID")) %>%
  filter(Name.x != Name.y) %>%
  mutate(TransitionGroupID = paste0(Name.x, "->", Name.y, " [Type]"),
         Year = from_time + initYear, 
         Timestep = Year, 
         TertiaryStratumID = land_type, 
         Amount20year = probability, 
         probabilityModified = ifelse(probability == 1, 0.9999, probability),
         numTimesteps = to_time - from_time,
         Amount1year = 1 - (1-probabilityModified)^(1/numTimesteps)) %>%
  select(names(transMultRaw), Timestep, TertiaryStratumID, TransitionGroupID, srcClassID, destClassID, Amount20year, probabilityModified, numTimesteps, Amount1year) 

write_csv(transMultMergedstsim, file.path(b03ProcessedTabularDir, str_c("mergedSTSimDefinitions_", landisName)))

# Expand Landis output file to include all transitions
# Transitions that do not occur need to have a multiplier value = 0
transMult <- transMultMergedstsim %>%
  select(Timestep, TertiaryStratumID, TransitionGroupID, Amount=Amount1year) 
# Make transitionTypesForest$TransitionGroupID into a factor so that it can be expanded
transMult$TransitionGroupID <- factor(transMult$TransitionGroupID, levels=paste(transitionTypesForest$Transition, "[Type]"))

transMultExpanded <- transMult %>% 
  expand(Timestep, TertiaryStratumID, TransitionGroupID)

# Merge expanded transitions to get all transition multipliers
transMultAll <- transMult %>% 
  right_join(transMultExpanded) %>% 
  mutate(Amount=replace_na(Amount,0))

# Save to st-sim library
myScenario <- scenario(myProject, scenarioName)
if (saveDatasheets) saveDatasheet(myScenario, transMultAll, "stsim_TransitionMultiplierValue")

# Error saving datasheet to library:
# Cannot translate ID value for column 'TertiaryStratumID': 1
# Warning messages:
#   1: In if ((class(data) != "list") | (class(data[[1]]) != "data.frame")) { :
#       the condition has length > 1 and only the first element will be used
#     2: In if (class(cDat) != "data.frame") { :
#         the condition has length > 1 and only the first element will be used

# Save csv to disk
write_csv(transMultAll, file.path(b03ProcessedTabularDir, "stsim-transition-mutiplier-values-baseline.csv"))

## Transition Multipliers - RCP 4.5 Climate ----
# Scenario-specific inputs
scenarioName <- "Transition Multipliers: Climate RCP 4.5" # Name of st-sim scenario for transition multipliers
landisName <- "landis_rcp45_age_transitions.csv" #Name of Landis csv from Guillaume

# Read in 20-year transition multiplier probabilities from Guillaume
transMultRaw <- read_csv(file.path(b03RawTablesDir, landisName))

# Merge Landis output file with stsim state classes
transMultMergedstsim <- mutate(transMultRaw, srcClassID = (500 + from_cover*10 + age_before), destClassID = (500 + to_cover * 10 + age_after)) %>%
  left_join(stateClasses, by=c("srcClassID"="ID")) %>%
  left_join(stateClasses, by=c("destClassID"="ID")) %>%
  filter(Name.x != Name.y) %>%
  mutate(TransitionGroupID = paste0(Name.x, "->", Name.y, " [Type]"),
         Year = from_time + initYear, 
         Timestep = Year, 
         TertiaryStratumID = land_type, 
         Amount20year = probability, 
         probabilityModified = ifelse(probability == 1, 0.9999, probability),
         numTimesteps = to_time - from_time,
         Amount1year = 1 - (1-probabilityModified)^(1/numTimesteps)) %>%
  select(names(transMultRaw), Timestep, TertiaryStratumID, TransitionGroupID, srcClassID, destClassID, Amount20year, probabilityModified, numTimesteps, Amount1year) 

write_csv(transMultMergedstsim, file.path(b03ProcessedTabularDir, str_c("mergedSTSimDefinitions_", landisName)))

# Expand Landis output file to include all transitions
# Transitions that do not occur need to have a multiplier value = 0
transMult <- transMultMergedstsim %>%
  select(Timestep, TertiaryStratumID, TransitionGroupID, Amount=Amount1year) 
# Make transitionTypesForest$TransitionGroupID into a factor so that it can be expanded
transMult$TransitionGroupID <- factor(transMult$TransitionGroupID, levels=paste(transitionTypesForest$Transition, "[Type]"))

transMultExpanded <- transMult %>% 
  expand(Timestep, TertiaryStratumID, TransitionGroupID)

# Merge expanded transitions to get all transition multipliers
transMultAll <- transMult %>% 
  right_join(transMultExpanded) %>% 
  mutate(Amount=replace_na(Amount,0))

# Save to st-sim library
myScenario <- scenario(myProject, scenarioName)
if (saveDatasheets) saveDatasheet(myScenario, transMultAll, "stsim_TransitionMultiplierValue")

# Save csv to disk
write_csv(transMultAll, file.path(b03ProcessedTabularDir, "stsim-transition-mutiplier-values-rcp45.csv"))


## Transition Multipliers - RCP 8.5 Climate ----
# Scenario-specific inputs
scenarioName <- "Transition Multipliers: Climate RCP 8.5" # Name of st-sim scenario for transition multipliers
landisName <- "landis_rcp85_age_transitions.csv" #Name of Landis csv from Guillaume

# Read in transition multiplier probabilities from Guillaume
transMultRaw <- read_csv(file.path(b03RawTablesDir, landisName))

# Merge Landis output file with stsim state classes
transMultMergedstsim <- mutate(transMultRaw, srcClassID = (500 + from_cover*10 + age_before), destClassID = (500 + to_cover * 10 + age_after)) %>%
  left_join(stateClasses, by=c("srcClassID"="ID")) %>%
  left_join(stateClasses, by=c("destClassID"="ID")) %>%
  filter(Name.x != Name.y) %>%
  mutate(TransitionGroupID = paste0(Name.x, "->", Name.y, " [Type]"),
         Year = from_time + initYear, 
         Timestep = Year, 
         TertiaryStratumID = land_type, 
         Amount20year = probability, 
         probabilityModified = ifelse(probability == 1, 0.9999, probability),
         numTimesteps = to_time - from_time,
         Amount1year = 1 - (1-probabilityModified)^(1/numTimesteps)) %>%
  select(names(transMultRaw), Timestep, TertiaryStratumID, TransitionGroupID, srcClassID, destClassID, Amount20year, probabilityModified, numTimesteps, Amount1year) 

write_csv(transMultMergedstsim, file.path(b03ProcessedTabularDir, str_c("mergedSTSimDefinitions_", landisName)))

# Expand Landis output file to include all transitions
# Transitions that do not occur need to have a multiplier value = 0
transMult <- transMultMergedstsim %>%
  select(Timestep, TertiaryStratumID, TransitionGroupID, Amount=Amount1year) 
# Make transitionTypesForest$TransitionGroupID into a factor so that it can be expanded
transMult$TransitionGroupID <- factor(transMult$TransitionGroupID, levels=paste(transitionTypesForest$Transition, "[Type]"))

transMultExpanded <- transMult %>% 
  expand(Timestep, TertiaryStratumID, TransitionGroupID)

# Merge expanded transitions to get all transition multipliers
transMultAll <- transMult %>% 
  right_join(transMultExpanded) %>% 
  mutate(Amount=replace_na(Amount,0))

# Save to st-sim library
myScenario <- scenario(myProject, scenarioName)
if (saveDatasheets) saveDatasheet(myScenario, transMultAll, "stsim_TransitionMultiplierValue")

# Save csv to disk
write_csv(transMultAll, file.path(b03ProcessedTabularDir, "stsim-transition-mutiplier-values-rcp85.csv"))
