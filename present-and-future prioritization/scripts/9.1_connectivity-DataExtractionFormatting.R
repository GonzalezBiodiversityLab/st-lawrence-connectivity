################################################################################
## 01
## Organizing data from the connectivity analysis
################################################################################
#date: 20201126
#author: Kyle T. Martins

#Set working directory to the source file location 

for(i in 1:9){setwd("..")}

library(raster)
library(sf)
library(fasterize)
library(reshape2)
library(stringr)
library(ggplot2)

####DOWNLOAD THE DATA####

#Load the updated delimitation of the BTSL (sent 20201126)
bt=st_read(paste0("G:/My Drive/",
  "./4 - Eco2urbDataBase/Eco2urbProjectRepo/2020/MELCC/Phase-III/",
  "./Data/GeospatialData/1.ReferenceData/3.BTSLDelimitation/",
                  "btsl_90m_polygon.shp"))
areaBTSL=as.numeric(st_area(bt))*0.0001

#Load the delimitations of the ecoregions at the two levels of interest
ecor4= #Level 4, data sent from the MELCC
st_read(paste0(
  "G:/My Drive/4 - Eco2urbDataBase/Eco2urbProjectRepo/2020/MELCC/Phase-III/",
  "Data/GeospatialData/1.ReferenceData/5.CadresEcologiques/CER04.shp"))
  
ecor3= #Level 3
  st_read(paste0("G:/My Drive/",
                 "./4 - Eco2urbDataBase/Eco2urbDataRepo/GeospatialData/",
                 "Conservation/Quebec/",
                 "BTSL/CadresEcologiques/", 
        "CR_NIV_03_S.shp"))

plot(st_geometry(bt))
plot(st_geometry(ecor), add=TRUE)

head(ecor3)
econm=unique(data.frame(ecor3)[c("FID03", "NOM_ENS_PH")])
write.csv(econm, 
      paste0("G:/My Drive/4 - Eco2urbDataBase/Eco2urbProjectRepo/2020/MELCC/",
      "Phase-III/Data/GeospatialData/5.ConnectivityAnalysis/",
      "ConnectivityLV3EcoregionNames.csv"), row.names=FALSE)

#Load an example of the raster data to get the right projection
rex=raster(paste0("C:/Users/eco2urb/Documents/Downloads/",
              "curmap_CON85_Resistance.MAAM.ts2030.tif"))

bt=st_transform(bt, projection(rex))
ecor4=st_transform(ecor4, projection(rex))
ecor3=st_transform(ecor3, projection(rex))


####PREP THE ECOREGION LAYER####

#Intersect the ecoregions with the BTSL delimitation to get the right delim
bt=st_buffer(bt, 0)
ecor4=st_cast(ecor4, "MULTIPOLYGON")
ecor4=st_buffer(ecor4, 0)
bt3=st_intersection(ecor3, st_geometry(bt))
bt3=st_cast(bt3, "MULTIPOLYGON")
bt4=st_intersection(ecor4, st_geometry(bt))
bt4=st_cast(bt4, "MULTIPOLYGON")

#Merge the ecoregion level 3 codes into the shapefile for level 4
head(bt4)
head(bt3)
#Some editing had to be done because of sliver polygons the two levels were 
#not completely nested
bt5=st_intersection(bt3, bt4[c("FID04", "geometry")])
bt5$areaHa=as.numeric(st_area(bt5))*0.0001
bt5=data.frame(bt5)[c("areaHa", "FID03", "FID04")]
head(bt5)
bt5=aggregate(areaHa~FID03+FID04, bt5, sum)
bt5=aggregate(areaHa~FID03+FID04, bt5, max)
bt5=bt5[bt5$areaHa>1,]
bt5[duplicated(bt5$FID04),]
bt5[bt5$FID04 %in% c("320", "527"),]
sort(bt5$areaHa)
bt5=bt5[bt5$areaHa>10,]
bt5[duplicated(bt5$FID04),] #No more duplicates
#Merge the level 3 codes into the level 4 polygon
bt4=merge(bt4, bt5[c("FID03", "FID04")], by="FID04")
bt4$FID03=as.factor(bt4$FID03)
head(bt4)
plot(bt4[c("FID03")])

#Refresh the area and length calculations
bt4$SHAPE_Area=st_area(bt4)
bt4$SHAPE_Leng=st_length(bt4)

#Some of the classifications at level 3 need to be corrected further
#see if there are still any remaining slivver polygons
ars=aggregate(SHAPE_Area~FID03, bt4, sum)
ars$percent=round(ars$SHAPE_Area/sum(ars$SHAPE_Area)*100,2)
ars[order(ars$percent),]
#161 should be 359
bt4[bt4$FID03==161,]$FID03=359

#Some of the classifications at level 4 need to be corrected further
#see if there are still any remaining slivver polygons
ars=aggregate(SHAPE_Area~FID04, bt4, sum)
ars$percent=round(ars$SHAPE_Area/sum(ars$SHAPE_Area)*100,2)
ars[order(ars$percent),]
#319 and 1617 should be 391
bt4[bt4$FID04==319,]$FID04=391
bt4[bt4$FID04==1617,]$FID04=391

#Make sure that the surface area for the ecoregion delim hasn't changed from the
#original delim of the BTSL
sum(st_area(bt4))/sum(st_area(bt)) #very similar but slightly different

length(unique(bt4$FID04)) #109 different categories at level 4
length(unique(bt4$FID03)) #18 different levels for level 3

st_write(bt4, 
         paste0(
           "G:/My Drive/4 - Eco2urbDataBase/Eco2urbProjectRepo/2020/MELCC/",
           "Phase-III/Data/GeospatialData/1.ReferenceData/5.CadresEcologiques/",
           "CERO04_BTSL20201127.shp"), delete_layer = TRUE)

#Fasterize the ecoregion file
head(bt4)
summary(duplicated(bt4$FID04))
bt6=fasterize(bt4, rex, field="FID04")
writeRaster(bt6, paste0(
          "G:/My Drive/4 - Eco2urbDataBase/Eco2urbProjectRepo/2020/MELCC/",
          "Phase-III/Data/GeospatialData/1.ReferenceData/5.CadresEcologiques/",
          "CERO04_BTSL20201127.tif"))
plot(bt6)
summary(is.na(values(bt6)))
3512717*90*90/sum(st_area(bt4)) #compare the total areas 
3512717*90*90/sum(st_area(bt)) #compare the total areas

####SORT THROUGH THE CONNECTIVITY DATA####

fls=list.files(paste0("C:/Users/eco2urb/Documents/Downloads"))
fls2=fls#[grepl("-", fls)]

spl1=str_split_fixed(fls2,"[.]",  4)
scen=gsub("curmap_|_Resistance", "", spl1[,1])
end=as.numeric(gsub("ts", "", str_split_fixed(spl1[,3], "-", 2)[,1]))
strt=as.numeric(gsub("ts", "", str_split_fixed(spl1[,3], "-", 2)[,2]))
sp=spl1[,2]

scdf=data.frame(scen, sp, strt, end, fls2)
head(scdf)

#Quick correction to the start and end for the yearly data
scdf[is.na(scdf$strt),]$strt=scdf[is.na(scdf$strt),]$end

#Flag whether yearly data or difference data
scdf$type=ifelse(grepl("-",scdf$fls2), "Difference", "Yearly")
head(scdf)

#Simplifying a bit and only looking at comparisons if between the first 
#and last timestep
scdf=scdf[
scdf$type=="Difference" & (scdf$strt=="2010" & scdf$end=="2110") |
  scdf$type=="Yearly",]

#Make sure that have the same number of replicates per scenario and per species
table(scdf$scen)
table(scdf$sp)

#Create a list for the files
l1=split(as.character(scdf$fls2), as.character(scdf$fls2))

#Define the number of iterations in the loop
nmfls=length(names(l1))


#Some of the files had problems with them, so they need to be re-processed
flt=c("curmap_BAU85_Resistance.BLBR.ts2050.tif",
"curmap_BAU85_Resistance.MAAM.ts2050.tif",
"curmap_BAU85_Resistance.MAAM.ts2110-2010.tif",
"curmap_BAU85_Resistance.MAAM.ts2110.tif",
"curmap_BAU85_Resistance.PLCI.ts2050.tif",
"curmap_BAU85_Resistance.PLCI.ts2090.tif",
"curmap_BAU85_Resistance.URAM.ts2070.tif")

fltn=which(scdf$fls2 %in% flt)

for(i in fltn){
  #for(i in 1:nmfls){
fld=ifelse(scdf$type[i]=="Difference","Circuit_Differences", "Circuit_Yearly")
flnm=as.character(scdf$fls2[i])
r1=raster(paste0("C:/Users/eco2urb/Documents/Downloads/", 
              flnm))
#plot(r1)
#plot(bt6, add=TRUE)
values(r1)[is.na(values(bt6))]=NA
mns=data.frame(zonal(r1, bt6, fun="mean"))
sms=data.frame(zonal(r1, bt6, fun="sum"))
mxs=data.frame(zonal(r1, bt6, fun="max"))
mis=data.frame(zonal(r1, bt6, fun="min"))
mns=merge(mns, sms, by="zone")
mns=merge(mns, mxs, by="zone")
mns=merge(mns, mis, by="zone")
head(mns)
l1[[i]]=mns
 writeRaster(r1, 
 paste0("./My Drive/4 - Eco2urbDataBase/Eco2urbProjectRepo/2020/MELCC/",
              "Phase-III/Data/GeospatialData/5.ConnectivityAnalysis/",
              fld,"/", flnm), overwrite=TRUE)
print(round(i/nmfls,4)*100)
} #start 4:14 , end : 5:12 more or less 1 hour


#Format and export the tabular data
for(i in 1:nmfls){
l1[[i]]$file=names(l1)[i]
}

l2=do.call(rbind, l1)
rownames(l2)=1:dim(l2)[1]
names(l2)=c("FID04", "mean", "sum", "max","min", "fls2")
head(l2)
head(scdf)
l3=merge(l2, scdf, by="fls2")
bt4lk=data.frame(bt4)[c("FID04", "FID03")]
bt4lk=unique(bt4lk)
summary(duplicated(bt4lk$FID04))
l4=merge(l3, bt4lk, by="FID04")
head(l4)

write.csv(l4, paste0(
  "./My Drive/4 - Eco2urbDataBase/Eco2urbProjectRepo/2020/MELCC/",
  "Phase-III/Data/GeospatialData/5.ConnectivityAnalysis/",
  "ConnectivityTabularSummary.csv"
), row.names=FALSE)

#Export the table summarizing the different runs
write.csv(scdf, paste0(
  "./My Drive/4 - Eco2urbDataBase/Eco2urbProjectRepo/2020/MELCC/",
  "Phase-III/Data/GeospatialData/5.ConnectivityAnalysis/",
  "ConnectivityFileDirectory.csv"), row.names=FALSE
)


#End of scipt#