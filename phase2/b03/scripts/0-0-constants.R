# a254
# Sarah Chisholm and Bronwyn Rayfield, ApexRMS
# Run with R-4.1.2
#
# This script is used to define constants and load
# necessary functions in a reproducible way across scripts

# Workspace ----
# Set environment variable TZ when running on AWS EC2 instance
Sys.setenv(TZ='UTC')
options(stringsAsFactors = FALSE)

# Load packages
library(raster)
library(tidyverse)
library(rgrass7)
library(rgdal)
library(sf)
library(fasterize)

# Parameters ----
## Global parameters ----

# Resolution in meters
myResolution <- 30

# Coarse resolution in meters
coarseResolution <- 240

# Set up GRASS mapset for the first time
doGRASSSetup <- F

# List of species codes
speciesList <- c("MAAM", "BLBR", "URAM", "PLCI", "RASY")

# Species that should have gaps within patches filled
speciesFillGaps <- c('MAAM','URAM')

## Habitat network analysis parameters ----

# Set clipping threshold
# Patches with less than clipAreaThreshold proportion of their area within the 
# ecological boundary will be removed
clipAreaThreshold <- 0.8

## Protected areas parameters ----

# Minimum percentage of the protected area that must fall in water for the area 
# to be considered "aquatic" (and therefore removed from consideration)
waterThreshold <- 50 

# Protected areas must be at least this size (in ha) to be retained within the buffer
bufferThreshold <- 900 

# Protected areas must be at least this size (in ha) to be retained within the study area
studyAreaThreshold <- 150 

# Directories ----
## Core directories ----
gisBase <- "C:/Program Files/GRASS GIS 7.8"
b03Dir <- "b03"
b01b02Dir <- "b01b02"
gisDbase <- file.path(b03Dir, "grass7")
RscriptDir<-paste0(b03Dir, "RScripts/")

# b03 
b03RawMapsDir <- file.path(b03Dir, "data", "spatial")
b03RawTablesDir <- file.path(b03Dir, "data", "tabular")
b03ProcessedMapsDir <- file.path(b03Dir, "model-inputs", "spatial")
b03ProcessedTabularDir <- file.path(b03Dir, "model-inputs", "tabular")
# b02 
b01b02RawMapsDir <- file.path(b01b02Dir, "inputs", "rawData", "maps")
b01b02RawTablesDir <- file.path(b01b02Dir, "inputs", "rawData", "tables")

## Composite directories ----
b03habitatDir <- file.path(b03Dir, "model-outputs", "spatial", "1.Habitat")
b03resistanceDir <- file.path(b03Dir, "model-outputs", "spatial", "2.Resistance")
b03networkDir <- file.path(b03Dir, "model-outputs", "spatial", "3.NetworkConnectivity")
b03patchImportanceDir <- file.path(b03Dir, "model-outputs", "spatial", "4.PatchImportance")
b03circuitscapeDir <- file.path(b03Dir, "model-outputs", "spatial", "5.Circuitscape")
b03zonationDir <-file.path(b03Dir, "model-outputs", "spatial", "6.Zonation")

# Functions ----
# Rescale a raster layer
rescaleR <- function(x, new.min = 0, new.max = 1) {
  x.min = suppressWarnings(min(x, na.rm=TRUE))
  x.max = suppressWarnings(max(x, na.rm=TRUE))
  new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))
}

# Get GRASS vector attribute table
v.get.att <- function(vector_name, sep){
  # Get attributes
  att <- execGRASS("v.db.select", map=vector_name, separator=sep, intern=T)
  
  # Format as dataframe
  tc <- textConnection(att)
  df <- read.table(tc, header = TRUE, sep=sep)
  close(tc)
  
  # Return resulting dataframe
  return(df)
}
