################################################################################
## 03
## Ecoregion Connectivity Analysis
################################################################################
#date: 20201127
#author: Kyle T. Martins

#Set working directory to the source file location 

for(i in 1:2){setwd("..")}

library(raster)
library(sf)
library(fasterize)
library(reshape2)
library(stringr)
library(ggplot2)
library(viridis)
library(scales)
library(cowplot)
library(magick)
library(ggimage)

#File organization directory
fod=read.csv(paste0("./Data/GeospatialData/5.ConnectivityAnalysis/",
                    "ConnectivityFileDirectory.csv"))

#Data with mean values per ecoregion
ecv=read.csv(paste0("./Data/GeospatialData/5.ConnectivityAnalysis/",
                    "ConnectivityTabularSummary.csv"))
head(ecv)

#Shapefile with ecoregions
ecr=
st_read(paste0(
  "G:/My Drive/4 - Eco2urbDataBase/Eco2urbProjectRepo/2020/MELCC/",
  "Phase-III/Data/GeospatialData/1.ReferenceData/5.CadresEcologiques/",
  "CERO04_BTSL20201127.shp"))

#Saint-Lawrence Lowlands
bt=st_read(paste0("G:/My Drive/",
            "./4 - Eco2urbDataBase/Eco2urbProjectRepo/2020/MELCC/Phase-III/",
            "./Data/GeospatialData/1.ReferenceData/3.BTSLDelimitation/",
            "btsl_90m_polygon.shp"))

btp=as(bt, "Spatial")

head(ecv)

fll=sort(unique(ecv$fls2))

#Specify the path for exporting the data
path2="/Results/ConnectivityAnalysis/"

####RANGE ANALYSIS####

#Look through the distribution of mean difference values across runs and decide
#on thresholds for negative, neutral and positive changes

#Extract the difference data
dif=ecv[ecv$type=="Difference",]
head(dif)
hist(dif$mean, xlab="Mean difference")
hist(dif$mean, xlab="Mean difference") #skewed towards the negative
qs=quantile(round(dif$mean,4))
#Neutral change : between -0.0012 and 0.0001
#Negative change : below -0.0012
#Positive change : above 0.0001

ecv$meanRounded=round(ecv$mean, 4)
ecv$Change=""
ecv[ecv$type=="Difference",]$Change=
  ifelse(ecv[ecv$type=="Difference",]$meanRounded<=qs[4] & 
                    ecv[ecv$type=="Difference",]$meanRounded>=qs[2], 
                  "Neutral", 
                  ifelse(ecv[ecv$type=="Difference",]$meanRounded<qs[2], 
                         "Negative", "Positive" ))
sort(unique(ecv$mean))
summary(ecv$mean<=0.0001)

boxplot(mean~Change, data=ecv[ecv$type=="Difference",])
table(ecv[ecv$type=="Difference",]$Change)

ecv$Change=as.factor(ecv$Change)
ecv$Change=factor(ecv$Change, levels=c("Negative", "Neutral", "Positive"))

####MAPPING FUNCTIONS####

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

####MAPPING PARAMETERS####

#Get the range of connectivity and difference values to assign the same range 
#to all the map legends

sort(ecv$mean, decreasing=TRUE)[1:10]
head(ecv)

#See the range of values per species across scenarios and years
aggregate(mean~sp, data=ecv[ecv$type=="Yearly",], max)
aggregate(mean~sp, data=ecv[ecv$type=="Yearly",], min)

aggregate(mean~sp, data=ecv[ecv$type=="Difference",], max)
aggregate(mean~sp, data=ecv[ecv$type=="Difference",], min)

mnmx=data.frame(type=c("Yearly", "Difference"), min="", max="",
                stringsAsFactors=FALSE)

mnmx[mnmx$type=="Yearly",c("min", "max")]=c(
  round(0,5), round(0.033,5))

mnmx[mnmx$type=="Difference",c("min", "max")]=c(-0.011,0.0037)

qnyr <- rescale(quantile(
  ecv[ecv$fls2=="curmap_BAU45_Resistance.BLBR.ts2010.tif",]$mean, 
  probs=seq(0, 1, length.out=64)))

qndf <- rescale(quantile(
  ecv[ecv$fls2=="curmap_BAU45_Resistance.BLBR.ts2110-2010.tif",]$mean, 
  probs=seq(0, 1, length.out=64)))


####LOOP THROUGH MAPS####

nfls=length(fll)

for(i in 1:nfls){
ecv2=ecv[ecv$fls2==fll[i],]
ecr2=merge(ecr, ecv2[c("FID04", "mean", "Change")])

if(sum(ecr2$mean)!=0){
nm=ifelse(unique(ecv2$type)=="Yearly",
          paste0("Current Density ", unique(ecv2$strt)), 
          "Current Difference 2010-2110")

col=ifelse(unique(ecv2$type)=="Yearly", "plasma", "viridis")

fld=ifelse(unique(ecv2$type)=="Yearly", "Circuit_Yearly_Ecoregion", 
           "Circuit_Difference_Ecoregion")

sp=as.character(unique(ecv2$sp))

an_file=paste0(
  "./Data/GeospatialData/5.ConnectivityAnalysis/Circuit_animals_png/",
  sp, ".png")

if(unique(ecv2$type)=="Yearly"){
  qn<-qnyr
  mnV<- as.numeric(mnmx[mnmx$type=="Yearly",]$min)
  mxV<-as.numeric(mnmx[mnmx$type=="Yearly",]$max)
} else {
  qn<-qndf
  mnV<- as.numeric(mnmx[mnmx$type=="Difference",]$min)
  mxV<-as.numeric(mnmx[mnmx$type=="Difference",]$max)
}

q=ggplot() + 
  geom_polygon(data = btp, aes(x = long+2000, #shadow
                               y = lat-800,
                               group=group),
               color = "grey", size = 1, fill="grey")+
  geom_sf(data = bt, size = 1, fill="white", color="black")+ #btsl
  geom_sf(data = ecr2, aes(fill=mean), size = 0.03, color="black")+
  #geom_sf(data=riv, fill=colst[7], size=0.2 )
  #geom_raster(data=bm4, aes(x=x,y=y, fill=bm4[,1])) + #raster layer
  coord_sf()+
  theme_map()+
  labs(x = NULL, 
       y = NULL)+
  scale_fill_viridis( #viridis color schale for index values
    values = qn,
    option = col, 
    #direction = -1,
    #specifies the legend title
    name = nm,
    #breaks=c(0, 0.15,0.5,0.85),
    #labels=labelst, #specifies the labels used in the legend
    limits=c(mnV, mxV),
    # here we use guide_colourbar because it is still a continuous scale
    guide = guide_colorbar(
      direction = "horizontal",
      barheight = unit(2, units = "mm"),
      barwidth = unit(50, units = "mm"),
      draw.ulim = F,
      title.position = 'top',
      # some shifting around
      title.hjust = 0.5,
      label.hjust = 0.5
    ))
#q

q<-ggdraw(q) + 
  draw_image(an_file, x = 0.70, y = 1.25, hjust = 1, vjust = 1, 
             scale=0.2)

ggsave(
paste0(".", path2, fld,"/",
       gsub("[.]tif","" , fll[i]), ".png"), q, height=6, width=6)

if(unique(ecv2$type)!="Yearly"){
q=ggplot() + 
  geom_polygon(data = btp, aes(x = long+2000, #shadow
                               y = lat-800,
                               group=group),
               color = "grey", size = 1, fill="grey")+
  geom_sf(data = bt, size = 1, fill="white", color="black")+ #btsl
  geom_sf(data = ecr2, aes(fill=Change), size = 0.03, color="black")+
  #geom_sf(data=riv, fill=colst[7], size=0.2 )
  #geom_raster(data=bm4, aes(x=x,y=y, fill=bm4[,1])) + #raster layer
  coord_sf()+
  theme_map()+
  labs(x = NULL, 
       y = NULL)+
  scale_fill_manual( #viridis color scale for index values
    values = c("#87a5ffff", "#fff787ff", "#ff5555ff"),
    #direction = -1,
    #specifies the legend title
    name = nm
    #breaks=c(0, 0.15,0.5,0.85),
    #labels=labelst, #specifies the labels used in the legend
    
    # here we use guide_colourbar because it is still a continuous scale
  )
#q

q<-ggdraw(q) + 
  draw_image(an_file, x = 0.70, y = 1.25, hjust = 1, vjust = 1, 
             scale=0.2)

ggsave(
  paste0(".", path2, fld,"_Binned", "/",
         gsub("[.]tif","" , fll[i]), ".png"), q, height=6, width=6)

}
}
}
