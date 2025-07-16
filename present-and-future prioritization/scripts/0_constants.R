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
library(sf)
library(tidyverse)
library(fasterize)
library(tidyverse)
library(rsyncrosim)
library(scales)

# Parameters ----
## Global parameters ----

# List of species codes
speciesList <- c("MAAM", "BLBR", "URAM", "PLCI", "RANA")


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
b03ResultsMapsDir <- file.path(b03Dir, "model-outputs", "spatial")
b03ResultsTabularDir <- file.path(b03Dir, "model-outputs", "tabular")
b03zonationDir <- file.path(b03Dir, "model-outputs", "zonation")

# Functions ----
# Rescale a raster layer
rescaleR <- function(x, new.min = 0, new.max = 1) {
  x.min = suppressWarnings(min(x, na.rm=TRUE))
  x.max = suppressWarnings(max(x, na.rm=TRUE))
  new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min))
}

# Mapping functions
#Assemble all the pieces and map them together
#Define the map theme
theme_map <- function(...) {
  theme_minimal() +
    theme(
      #text = element_text(family = "Arial", color = "#22211d"),
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.background = element_rect(fill = "white", color = NA), 
      panel.background = element_rect(fill = "white", color = NA), 
      legend.background = element_rect(fill = "white", color = NA),
      #panel.border = element_rect(colour = "grey", fill=NA, size=1),
      legend.position="bottom",
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(-10,-10,-10,-10),
      plot.margin=grid::unit(c(0,0,0,0), "mm")
    )
}
