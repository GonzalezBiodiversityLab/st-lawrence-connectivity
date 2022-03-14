# a254
# Sarah Chisholm and Bronwyn Rayfield, ApexRMS
# Run with R-4.1.2
#
# This script is used to define constants and load
# necessary functions in a reproducible way across scripts

# Workspace ----
# Set environment variable TZ when running on AWS EC2 instance
Sys.setenv(TZ='UTC')

# Load packages
library(raster)
library(tidyverse)


# Parameters ----
## Global parameters ----


## Transition Multiplier parameters ----
initYear <- 2010

# Directories ----
## Core directories ----
# gisBase <- "C:/Program Files/GRASS GIS 7.8"
b03Dir <- "b03"
b01b02Dir <- "b01b02"
# gisDbase <- file.path(b03Dir, "grass7")
# RscriptDir<-paste0(b03Dir, "RScripts/")

# b03 
b03RawMapsDir <- file.path(b03Dir, "data", "spatial")
b03RawTablesDir <- file.path(b03Dir, "data", "tabular")
b03ProcessedMapsDir <- file.path(b03Dir, "model-inputs", "spatial")
b03ProcessedTabularDir <- file.path(b03Dir, "model-inputs", "tabular")
b03Libraries <- file.path(b03Dir, "libraries")
# b02 
b01b02RawMapsDir <- file.path(b01b02Dir, "inputs", "rawData", "maps")
b01b02RawTablesDir <- file.path(b01b02Dir, "inputs", "rawData", "tables")

## Composite directories ----
# b03habitatDir <- file.path(b03Dir, "model-outputs", "spatial", "1.Habitat")
# b03resistanceDir <- file.path(b03Dir, "model-outputs", "spatial", "2.Resistance")
# b03networkDir <- file.path(b03Dir, "model-outputs", "spatial", "3.NetworkConnectivity")
# b03patchImportanceDir <- file.path(b03Dir, "model-outputs", "spatial", "4.PatchImportance")
# b03circuitscapeDir <- file.path(b03Dir, "model-outputs", "spatial", "5.Circuitscape")

# Functions ----
