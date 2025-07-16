################################################################################
## 05
## Composite Figures
################################################################################
#date: 20201127
#author: Kyle T. Martins

#Script is for plotting composite figures with output from the previous scripts

#Set working directory to the source file location 

#Link references: 
#https://www.datanovia.com/en/lessons/combine-multiple-ggplots-into-a-figure/
#https://www.rdocumentation.org/packages/ggpubr/versions/0.4.0/topics/ggarrange

for(i in 1:2){setwd("..")}

library(sf)
library(reshape2)
library(stringr)
library(ggplot2)
library(viridis)
library(scales)
library(cowplot)
library(magick)
library(ggimage)
library(ggthemes)
library(ggdark)
library(RColorBrewer)
library(scales)
library(ggpubr)

fld="./Results/ConnectivityAnalysis/Circuit_Trends_Ecoregion/"
fl1=list.files(fld)[1]
fl2=list.files(fld)[2]

img1=image_read(paste0(fld, fl1))
img1=image_ggplot(img1)
img2=image_read(paste0(fld, fl2))
img2=image_ggplot(img2)

gcomp=ggarrange(img1, img2, nrow=1, ncol=2)

ggsave("./Results/ConnectivityAnalysis/Composite_trial/Trial1.png", gcomp)
