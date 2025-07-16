################################################################################
## 02
## Mapping Circuit Connectivity
################################################################################
#date: 20201127
#author: Kyle T. Martins

#Script is for plotting results from the connectivity analysis


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
dim(fod)

#Data with mean, max and min values per ecoregion
ecv=read.csv(paste0("./Data/GeospatialData/5.ConnectivityAnalysis/",
                    "ConnectivityTabularSummary.csv"))[,-1]

#Delimitation of the BTSL
bt=st_read(paste0("G:/My Drive/",
              "./4 - Eco2urbDataBase/Eco2urbProjectRepo/2020/MELCC/Phase-III/",
            "./Data/GeospatialData/1.ReferenceData/3.BTSLDelimitation/",
              "btsl_90m_polygon.shp"))

#Delimitation of the ecoregions
ecor=st_read(paste0(
           "G:/My Drive/4 - Eco2urbDataBase/Eco2urbProjectRepo/2020/MELCC/",
           "Phase-III/Data/GeospatialData/1.ReferenceData/5.CadresEcologiques/",
           "CERO04_BTSL20201127.shp"))

plot(st_geometry(bt))
plot(st_geometry(ecor), add=TRUE)

btp=as(bt, "Spatial")

#Yearly data
flsyr=list.files(paste0(
  "./Data/GeospatialData/5.ConnectivityAnalysis/Circuit_Yearly/"))
flsyr=flsyr[grepl("tif", flsyr)]

#Difference data
flsdf=list.files(paste0(
  "./Data/GeospatialData/5.ConnectivityAnalysis/Circuit_Differences/"))
flsdf=flsdf[grepl("tif", flsdf)]

fld=c("Circuit_Yearly", "Circuit_Differences")

l1=as.list(fld)
l1[[2]]=flsdf; l1[[1]]=flsyr
names(l1)=fld
l1=melt(l1)
head(l1)
names(l1)=c("file", "fld")
nfls=dim(l1)[1]

l1=merge(l1, fod, by.x="file", by.y="fls2")
head(l1)

#Specify the path for exporting the data
path2="/Results/ConnectivityAnalysis/"

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
aggregate(max~sp, data=ecv[ecv$type=="Yearly",], max)
aggregate(min~sp, data=ecv[ecv$type=="Yearly",], min)

aggregate(max~sp, data=ecv[ecv$type=="Difference",], max)
aggregate(min~sp, data=ecv[ecv$type=="Difference",], min)


mnmx=data.frame(type=c("Yearly", "Difference"), min="", max="",
                stringsAsFactors=FALSE)

mnmx[mnmx$type=="Yearly",c("min", "max")]=c(
  round(0,5), round(0.22,5))

mnmx[mnmx$type=="Difference",c("min", "max")]=c(-0.18,0.11)

#Take a random set of raster files to set the quantile distribution of the 
#color gradients for the set of files

flsl=which(l1$file=="curmap_BAU45_Resistance.BLBR.ts2010.tif")
ra=raster(paste0(
  "./Data/GeospatialData/5.ConnectivityAnalysis/",
  l1[flsl,2],"/", l1[flsl,1]))
bm4 <- as(ra, "SpatialPixelsDataFrame")
bm4 <- as.data.frame(bm4); 
qnyr <- rescale(quantile(bm4[,1], probs=seq(0, 1, 
                                          length.out=64)))

flsl=which(l1$file=="curmap_BAU45_Resistance.BLBR.ts2110-2010.tif")
ra=raster(paste0(
  "./Data/GeospatialData/5.ConnectivityAnalysis/",
  l1[flsl,2],"/", l1[flsl,1]))
bm4 <- as(ra, "SpatialPixelsDataFrame")
bm4 <- as.data.frame(bm4); 
qndf <- rescale(quantile(bm4[,1], probs=seq(0, 1, 
                                            length.out=64)))

####MAPPING LOOP####

#left off at i = 164
for(i in 164:nfls){
ra=raster(paste0(
"./Data/GeospatialData/5.ConnectivityAnalysis/",
l1[i,2],"/", l1[i,1]))

bm4 <- as(ra, "SpatialPixelsDataFrame")
bm4 <- as.data.frame(bm4); 
head(bm4)

if(sum(bm4[,1])!=0){
col=ifelse(l1$fld[i]=="Circuit_Yearly", "plasma", "viridis")

nm=ifelse(l1$fld[i]=="Circuit_Yearly",
          paste0("Current Density ", l1$strt[i]), 
          "Current Difference 2010-2110")

if(l1$fld[i]=="Circuit_Yearly"){
  qn<-qnyr
  mnV<- as.numeric(mnmx[mnmx$type=="Yearly",]$min)
  mxV<-as.numeric(mnmx[mnmx$type=="Yearly",]$max)
  } else {
  qn<-qndf
  mnV<- as.numeric(mnmx[mnmx$type=="Difference",]$min)
  mxV<-as.numeric(mnmx[mnmx$type=="Difference",]$max)
  }

sp=as.character(l1$sp[i])

an_file=paste0(
  "./Data/GeospatialData/5.ConnectivityAnalysis/Circuit_animals_png/",
  sp, ".png")

q=ggplot() + 
  geom_polygon(data = btp, aes(x = long+2000, #shadow
                                y = lat-800,
                                group=group),
               color = "grey", size = 1, fill="grey")+
  geom_sf(data = bt, size = 1, fill="white", color="black")+ #btsl
  #geom_sf(data=riv, fill=colst[7], size=0.2 )
  geom_raster(data=bm4, aes(x=x,y=y, fill=bm4[,1])) + #raster layer
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


png(paste0(".", path2, as.character(l1$fld[i]),"/",
           gsub("[.]tif", "", l1$file[i]), ".png"), #Export the data
    res=300, width=6, height=6, units="in") #specifying high resolution images
print(ggdraw(q) + 
        draw_image(an_file, x = 0.70, y = 1.25, hjust = 1, vjust = 1, 
                   scale=0.2))
dev.off()
}
print((round(i/nfls,4))*100)
}

#End of script#