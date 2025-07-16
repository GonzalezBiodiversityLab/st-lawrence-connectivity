################################################################################
## 04
## Graphical Summaries
################################################################################
#date: 20201127
#author: Kyle T. Martins

#Script is for plotting graphical summaries of connectivity results

#Set working directory to the source file location 

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
library(stringr)

#File organization directory
fod=read.csv(paste0("./Data/GeospatialData/5.ConnectivityAnalysis/",
                    "ConnectivityFileDirectory.csv"))

#Data with mean values per ecoregion
ecv=read.csv(paste0("./Data/GeospatialData/5.ConnectivityAnalysis/",
                    "ConnectivityTabularSummary.csv"))

#Shapefile with ecoregions
ecr=
  st_read(paste0(
    "G:/My Drive/4 - Eco2urbDataBase/Eco2urbProjectRepo/2020/MELCC/",
    "Phase-III/Data/GeospatialData/1.ReferenceData/5.CadresEcologiques/",
    "CERO04_BTSL20201127.shp"))

#Raster with the ecoregions
rec=raster(paste0(
  "Data/GeospatialData/1.ReferenceData/5.CadresEcologiques/",
  "CERO04_BTSL20201127.tif"
))
#plot(rec)

#Shapefile of level 2 ecoregions
elv2=st_read(paste0("Data/GeospatialData/1.ReferenceData/5.CadresEcologiques/",
        "CERO02_BTSL20201203.shp"))

#Get the ecoregion names for level 3 
econm=read.csv(paste0("./Data/GeospatialData/5.ConnectivityAnalysis/",
                "ConnectivityLV3EcoregionNames.csv"))

#Get the unique list of file names
fll=unique(ecv$fls2)

#Group them at the level of SCENARIO X SPECIES 
fll2=str_split_fixed(as.character(fll), "[.]", 3)
fll2=unique(paste(fll2[,1], fll2[,2], sep="."))

#Specify the path for exporting the data
path2="/Results/ConnectivityAnalysis/"

#Specify the folder where should be exported
fld="Circuit_Trends_Ecoregion"

####DATA FORMATTING####

#Merge in the ecoregion names
head(ecv)
summary(ecv$FID03 %in% econm$FID03)
summary(duplicated(econm$FID03))
ecv=merge(ecv, econm, by="FID03")
head(ecv)

#Get the areas of the ecoregions
head(ecr)
ecr$areaHa=as.numeric(round(st_area(ecr)*0.0001,2))

#Merge in the areas to the ecv table (summary table with stats per ecoregion)
head(ecv)
ecv=merge(ecv, data.frame(ecr)[c("FID04", "areaHa")], by="FID04", all.x=TRUE)
head(ecv)

#Calculate the number of cells per ecoregion (level 4) and merge that to ecv
f1=freq(rec)
f1=f1[!is.na(f1[,1]),]
f1=data.frame(f1)
names(f1)=c("FID04", "count")

ecv=merge(ecv, f1, by="FID04")
head(ecv)

#Quick check 
plot(c(ecv$sum/ecv$mean)~ecv$areaHa, xlab="Area Ecoregion", 
     ylab="Sum / Mean Value") #odd that there's an outlier...and that the 
#relationship isn't 1 : 1 
abline(a=0, b=1)

#OK here we get a perfect 1 : 1 relationship, which is what we wanted
plot(c(ecv$sum/ecv$mean)~ecv$count, xlab="No. Cells Ecoregion", 
     ylab="Sum / Mean Value") #odd that there's an outlier...and that the 
#relationship isn't 1 : 1 
abline(a=0, b=1)

#Calcualte the mean value per ecoregion at level 3 from level 4 values
l3m=aggregate(cbind(sum, count)~FID03+fls2, data=ecv, sum)
head(l3m)
l3m$mean=round(l3m$sum/l3m$count,6)
head(l3m)

#To this df, add the information specific to the years, the species, etc. to 
#help with data management
head(ecv)
scinfo=unique(ecv[c("fls2", "scen", "sp", 
                    "strt", "end", "type")])
summary(duplicated(scinfo$fls2))
l3m=merge(l3m, scinfo, by="fls2")
head(l3m)

#To this df, add the level 3 ecoregion names
l3m=merge(l3m, econm, by="FID03")
head(l3m)

#Define the level 2 ecoregion codes to level 3 and 4 codes
rec #raster file with ecoregions level 4
elv2 #shapefile with ecoregions level 2
head(elv2)
rl2=fasterize(elv2, rec, field="FID02")
plot(rl2)
ct24=crosstab(rl2, rec)
ct24=melt(ct24)
names(ct24)=c("FID02", "FID04", "count")
ct24=ct24[ct24$count!=0,]
head(ct24)
summary(duplicated(ct24$FID04)) #check that no duplicates...
#Conduct the merges
lk1=unique(data.frame(ecr)[c("FID04", "FID03")])
summary(duplicated(lk1$FID04))
lk1=merge(lk1, ct24[c("FID04", "FID02")], by="FID04", all.x=TRUE)
lk1[is.na(lk1$FID02),] #checked in QGIS and will be assigning this ecoregion 
#level 4 the value of LV2 64 (southern portion)
lk1[is.na(lk1$FID02),]$FID02=64
#Export this reference file, since quite useful, incidentally
write.csv(lk1, 
          paste0("Data/GeospatialData/1.ReferenceData/5.CadresEcologiques/",
                 "CadreEcoloCodeCorrespondanceBTSL.csv"), row.names=FALSE)

l3m=merge(l3m, unique(lk1[c("FID03", "FID02")]), by="FID03")

#Add columns to distinguish the scenario types
l3m$scen=as.character(l3m$scen)
l3m$lusc=substr(l3m$scen,1,nchar(l3m$scen)-2)
l3m$clsc=str_sub(l3m$scen, start= -2)
head(l3m)

#Reasert the ordering of rows in the spreadsheet
class(l3m$strt)
l3m=l3m[with(l3m, order(type, scen, sp, strt, FID03)),]
#View(l3m)

#Use l3m mean values when plotting the trends in the ecological data


####GRAPHICAL SUMMARIES####

####ECOREGION GRAPHS####

for(i in 1:length(fll2)){

#Isolate the yearly data for a given species and scenario
ecv2=l3m[grepl(fll2[i], l3m$fls2) & l3m$type=="Yearly",]

#NOT NEEDED 20201203- Aggregate the data at the level of the ecoregion (level 3)
# ecv3=aggregate(mean~FID03+strt+NOM_ENS_PH, ecv2, mean)
# head(ecv3) #NOT NEEDED 20201203

#Calculate the cumulative percent difference
l1=split(ecv2, ecv2$FID03)

for(w in 1:length(names(l1))){
xx=l1[[w]]
xx$percD=round((xx$mean-xx$mean[1])/xx$mean[1],4)
l1[[w]]=xx
}
ecv3=do.call(rbind, l1)
#View(ecv3)

#Select the areas of interest to show up on the graphs
sel=c(230, 9666, 137, 88, 82, 184)

#Graping parametrs
nms=unique(ecv3[ecv3$FID03 %in% sel,]$NOM_ENS_PH)
nms=paste0("NIV03_", unique(ecv3[ecv3$FID03 %in% sel,]$FID03))

colors=brewer.pal(length(sel), "Paired")

sp=as.character(unique(ecv2$sp))

an_file=paste0(
  "./Data/GeospatialData/5.ConnectivityAnalysis/Circuit_animals_png/",
  sp, ".png")

q<-ggplot()+
  #geom_point()+
  geom_line(data=ecv3[!ecv3$FID03 %in% sel,], 
            aes(x=strt, y=percD, group=FID03), color="grey", 
            alpha=0.5, size=1)+
  geom_line(data=ecv3[ecv3$FID03 %in% sel,], 
            aes(x=strt, y=percD, colour=as.factor(FID03), 
                           group=FID03), size=1.5)+
  theme_gdocs()+
  theme(plot.title = element_text(family = "Fira Sans Condensed"),
        #plot.background = element_rect(fill = "grey10"),
        panel.background = element_blank(),
        #panel.grid.major = element_line(color = "grey30", size = 0.2),
        #panel.grid.minor = element_line(color = "grey30", size = 0.2),
        legend.background = element_blank(),
        axis.ticks = element_blank(),
        legend.key = element_blank(),
        legend.position = "top",
        legend.direction="horizontal",
        legend.text=element_text(size=10))+
  xlab("Year")+
  ylab("Percent Change in Mean Connectivity")+
  scale_y_continuous(labels = scales::percent, limits=c(-0.3, 0.15))+
  scale_colour_manual(labels=as.character(nms), values=colors, 
                      name="Ecoregion")

q<-ggdraw(q) + 
  draw_image(an_file, x = 0.80, y = 0.8,  hjust = 1, vjust = 1, 
             scale=0.2)

ggsave(
  paste0(".", path2, fld, "/",
         gsub("[.]tif","" , fll[i]), ".png"), q, height=6, width=6)

print(round((i/length(fll2))*100,2))  
}


####SCENARIO x ECOREGION GRAPHS####

#Here's we're plotting the mean change in current density for each ecoregion, 
#the lines represent the results for individual scenarios


combs=unique(l3m[l3m$type=="Yearly", c("FID03", "sp")])
combs$filename=paste0(combs$sp,"-", "LV03_", combs$FID03, ".png")

fld="Circuit_Trends_EcoregionXScenario"

for(i in 1:dim(combs)[1]){
  
  #head(l3m)
  tik=l3m[l3m$FID03==combs[i,1] & l3m$sp==combs[i,2] 
          & l3m$type=="Yearly", ]
  
  #Calculate the cumulative percent difference
  l1=split(tik, tik$scen)
  
  for(w in 1:length(names(l1))){
    xx=l1[[w]]
    xx$percD=round((xx$mean-xx$mean[1])/xx$mean[1],4)
    l1[[w]]=xx
  }
  tik=do.call(rbind, l1)
  
  colors=brewer.pal(length(unique(tik$scen)), "Paired")
  
  sp=as.character(unique(tik$sp))
  
  an_file=paste0(
    "./Data/GeospatialData/5.ConnectivityAnalysis/Circuit_animals_png/",
    sp, ".png")
  
  q<-ggplot()+
    geom_line(data=tik, 
              aes(x=strt, y=percD, colour=as.factor(scen), 
                  group=scen), size=1.5)+
    theme_gdocs()+
    theme(plot.title = element_text(family = "Fira Sans Condensed"),
          #plot.background = element_rect(fill = "grey10"),
          panel.background = element_blank(),
          #panel.grid.major = element_line(color = "grey30", size = 0.2),
          #panel.grid.minor = element_line(color = "grey30", size = 0.2),
          legend.background = element_blank(),
          axis.ticks = element_blank(),
          legend.key = element_blank(),
          legend.position = "top",
          legend.direction="horizontal",
          legend.text=element_text(size=10))+
    xlab("Year")+
    ylab("Percent Change in Mean Connectivity")+
    scale_y_continuous(labels = scales::percent, limits=c(-0.3, 0.15))+
    scale_colour_manual(values=colors, name="Scenario")
  
  q<-ggdraw(q) + 
    draw_image(an_file, x = 0.80, y = 0.8,  hjust = 1, vjust = 1, 
               scale=0.2)
  
  ggsave(
    paste0(".", path2, fld, "/",
           combs[i,3]), q, height=6, width=6)
  
  print(round((i/dim(combs)[1])*100,2))  
  
}


####BTSL X SCENARIO x ECOREGION GRAPHS####

#These graphs regroup all the ecoregion level 3 data for each species in a 
#single graph, there are two separate pannels for each of the two main level
#2 divisions of the BTSL

spl=as.character(unique(l3m$sp))

fld="Circuit_Trends_BTSLEcoregionXScenario"

lv2.labs <- c("Middle Saint-Laurent Plain", "Upper Saint-Laurent Plain")
names(lv2.labs) <- c("51", "64")

flnms=paste0(spl, "_BTSL_Scenario_Ecoregions.png")

for(i in 1:length(flnms)){
  
  #head(l3m)
  tik=l3m[l3m$sp==spl[i] & l3m$type=="Yearly", ]
  head(tik); dim(tik)
  
  #Define the variable to organize the calculation of cumulative differences
  tik$FID03scen=paste0(tik$FID03, "-", tik$scen)
  length(unique(tik$FID03scen))
  
  #Calculate the cumulative percent difference
  l1=split(tik, tik$FID03scen)
  l1[[1]]
  
  for(w in 1:length(names(l1))){
    xx=l1[[w]]
    xx$percD=round((xx$mean-xx$mean[1])/xx$mean[1],4)
    l1[[w]]=xx
  }
  tik=do.call(rbind, l1)
  
  colors=brewer.pal(length(unique(tik$scen)), "Paired")
  
  sp=as.character(unique(tik$sp))
  
  an_file=paste0(
    "./Data/GeospatialData/5.ConnectivityAnalysis/Circuit_animals_png/",
    sp, ".png")
  
  q<-ggplot()+
    geom_line(data=tik, 
              aes(x=strt, y=percD, colour=as.factor(scen), 
                  group=FID03scen), size=1.5)+
    facet_wrap(.~FID02, labeller = labeller(FID02 = lv2.labs))+
    theme_gdocs()+
    theme(plot.title = element_text(family = "Fira Sans Condensed"),
          #plot.background = element_rect(fill = "grey10"),
          panel.background = element_blank(),
          #panel.grid.major = element_line(color = "grey30", size = 0.2),
          #panel.grid.minor = element_line(color = "grey30", size = 0.2),
          legend.background = element_blank(),
          axis.ticks = element_blank(),
          legend.key = element_blank(),
          legend.position = "top",
          legend.direction="horizontal",
          legend.text=element_text(size=10))+
    xlab("Year")+
    ylab("Percent Change in Mean Connectivity")+
    scale_y_continuous(labels = scales::percent, limits=c(-0.3, 0.15))+
    scale_colour_manual(values=colors, name="Scenario")
  
  q<-ggdraw(q) + 
    draw_image(an_file, x = 0.80, y = 0.8,  hjust = 1, vjust = 1, 
               scale=0.2)
  
  ggsave(
    paste0(".", path2, fld, "/",
           flnms[i]), q, height=6, width=8)
  
  print(round((i/length(flnms))*100,2))  
  
}

####BTSL X SCENARIO X ECOREGION GRAPHS - GRIDDED####

spl=as.character(unique(l3m$sp))

fld="Circuit_Trends_BTSLEcoregionXScenario_Grid"

flnms=paste0(spl, "_BTSL_Scenario_Ecoregions_Grid.png")

head(l3m)

for(i in 1:length(flnms)){
  
  #head(l3m)
  tik=l3m[l3m$sp==spl[i] & l3m$type=="Yearly", ]
  head(tik); dim(tik)
  
  #Define the variable to organize the calculation of cumulative differences
  tik$FID03scen=paste0(tik$FID03, "-", tik$scen)
  length(unique(tik$FID03scen))
  
  #Calculate the cumulative percent difference
  l1=split(tik, tik$FID03scen)
  l1[[1]]
  
  for(w in 1:length(names(l1))){
    xx=l1[[w]]
    xx$percD=round((xx$mean-xx$mean[1])/xx$mean[1],4)
    l1[[w]]=xx
  }
  tik=do.call(rbind, l1)
  
  #Select the areas of interest to show up on the graphs
  sel=c(230, 9666, 137, 88, 82, 184)
  
  #nms=unique(tik[tik$FID03 %in% sel,]$NOM_ENS_PH)
  nms=paste0("NIV03_", unique(tik[tik$FID03 %in% sel,]$FID03))
  
  colors=brewer.pal(length(sel), "Paired")
  
  sp=as.character(unique(tik$sp))
  
  an_file=paste0(
    "./Data/GeospatialData/5.ConnectivityAnalysis/Circuit_animals_png/",
    sp, ".png")
  
  q<-ggplot()+
    geom_line(data=tik[!tik$FID03 %in% sel,], 
              aes(x=strt, y=percD, group=FID03), color="grey", 
              alpha=0.5, size=1)+
    geom_line(data=tik[tik$FID03 %in% sel,], 
              aes(x=strt, y=percD, colour=as.factor(FID03), 
                  group=FID03), size=1.5)+
    facet_grid(clsc~lusc)+
    theme_gdocs()+
    theme(plot.title = element_text(family = "Fira Sans Condensed"),
          #plot.background = element_rect(fill = "grey10"),
          panel.background = element_blank(),
          #panel.grid.major = element_line(color = "grey30", size = 0.2),
          #panel.grid.minor = element_line(color = "grey30", size = 0.2),
          legend.background = element_blank(),
          axis.ticks = element_blank(),
          legend.key = element_blank(),
          legend.position = "top",
          legend.direction="horizontal",
          legend.text=element_text(size=10))+
    xlab("Year")+
    ylab("Percent Change in Mean Connectivity")+
    scale_y_continuous(labels = scales::percent, limits=c(-0.3, 0.15))+
    scale_colour_manual(values=colors, name="Ecoregion", 
                        labels=as.character(nms))
  
  q<-ggdraw(q) + 
    draw_image(an_file, x = 0.61, y = 1.42,  hjust = 1, vjust = 1, 
               scale=0.1)
  
  ggsave(
    paste0(".", path2, fld, "/",
           flnms[i]), q, height=6, width=8)
  
  print(round((i/length(flnms))*100,2))  
  
}


#End of script#