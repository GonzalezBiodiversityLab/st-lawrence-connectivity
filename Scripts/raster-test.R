library(raster)

install.packages("rgdal")
primaryStratum<-raster("primaryStratum.tif")
primaryStratum
phase3Extent <- extent(primaryStratum)
phase3Extent
phase3Extent$xmin
phase3Extent[[xmin]]
newExtent <- extent(-491430-90, -95520 ,  117930 ,  403680 )
newExtent
x <- raster(ncol=36, nrow=18, xmn=-1000, xmx=1000, ymn=-100, ymx=900)
x
valuex(s) <- 1:ncell(x)
values(x) <- 1:ncell(x)
x11(); plot(x)
values(x) <- rep(1, ncell(x))
x11(); plot(x)
crs(x) <- crs(primaryStratum)
x
