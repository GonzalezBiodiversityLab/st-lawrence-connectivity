#####################################################################
# a254      
# Species performance curves from Zonation 
# 03-2022                                       					
#     
#Inputs:
#   -Zonation curves
#   -Zonation feature info
#  
#   Outputs:
#    -Plot of species performance curves: Proportion landscape lost vs. Conservation Criteria retained
#                                                                   
# Script by B Rayfield for ApexRMS 									
#####################################################################

# Workspace -------------------

library(reshape2)
library(ggplot2)
library(wesanderson)

# Load constants, functions, etc
source("./b03/scripts/0-0-constants.R")

options(stringsAsFactors = FALSE)
Sys.setenv(TZ='GMT')

# Zonation output directory
zonationOutputDirName <- file.path(b03zonationDir, "ALL")

# Zonation performance curves -----------------------
criteria.curves<-read.table(file.path(zonationOutputDirName, "ALL_Zonation.curves_BR.txt"), header=FALSE)
criteria.names<-read.table(file.path(zonationOutputDirName, "ALL_Zonation.features_info_BR.txt"), header=FALSE)
criteria.names<-as.character(criteria.names[,ncol(criteria.names)])

input.names<-c()
for(i in 1: length(criteria.names)){
  input.names[i]<-strsplit(criteria.names[i], "ZonationInputs/")[[1]][2]
}

criteria.curves.names<-c("Prop_landscape_lost",	"cost_needed_for_top_fraction",	"min_prop_rem",	"ave_prop_rem",	"W_prop_rem",	"ext-1", "ext-2", input.names)
names(criteria.curves)<-criteria.curves.names
attach(criteria.curves)

MAAM<-criteria.curves[,8:(8+5), drop=FALSE]
BLBR<-criteria.curves[,14:(14+5), drop=FALSE]
URAM<-criteria.curves[,20:(20+5), drop=FALSE]
PLCI<-criteria.curves[,26:(26+5), drop=FALSE]
RASY<-criteria.curves[,32:(32+5), drop=FALSE]

criteria.curves$BLBR<-rowMeans(BLBR)
criteria.curves$MAAM<-rowMeans(MAAM)
criteria.curves$PLCI<-rowMeans(PLCI)
criteria.curves$RASY<-rowMeans(RASY)
criteria.curves$URAM<-rowMeans(URAM)


# Plot ---------------------
criteria.curves<-criteria.curves[,-(2:(ncol(criteria.curves)-5))]
df <- melt(criteria.curves ,  id.vars = 'Prop_landscape_lost', variable.name = 'species')
x11(); ggplot(df, aes(Prop_landscape_lost,value)) + geom_line(aes(colour = species), size=1.5) + scale_color_manual(values=wes_palette(n=5, name="Darjeeling1")) + xlab("Proportion de milieux naturels perdus") + ylab("Proportion des critères de conservation retenus")
ggsave(file.path(b03zonationDir, "speciesCurves_b03.png"), width = 16, height = 13, units = "cm")

