Sys.setenv(TZ='GMT')

library(raster)
#install.packages("devtools")
#library("devtools")
#install_github("achubaty/grainscape")
library(grainscape)
library(foreign)
library(igraph)
library(plyr)

options(stringsAsFactors = FALSE)

# Directories
projectDir <- "c:/gitprojects/a254/phase2/b03/"
gisBase <- "C:/Program Files/GRASS GIS 7.8"
gisDbase <- paste0(projectDir, "grass7")
rawTablesDir <- paste0(projectDir, "data/tabular/")
b01b02RawTablesDir <- paste0(projectDir, "../b01b02/data/tabular/")
processedMapsDir <- paste0(projectDir, "model-outputs/")
habitatDir <- paste0(projectDir, "model-outputs/1.Habitat/")
resistanceDir <- paste0(projectDir, "model-outputs/2.Resistance/")
networkDir <- paste0(projectDir, "model-outputs/3.NetworkConnectivity/")
RscriptDir<-paste0(projectDir,"RScripts/")

# Input parameters
# speciesList<-c("MAAM", "PLCI", "RASY", "BLBR", "URAM")
# Run on MAAM for now
speciesList<-c("MAAM")
# Run this at 30m resolution
myResolution <- 30
#read in table of focal species dispersal parameters
DISP<-read.csv(paste0(b01b02RawTablesDir,"/speciesDispersalParameters.csv"),header=T)
# Set clipping threshold
# Patches with less than clipAreaThreshold proportion of their area within the ecological boundary will be removed
clipAreaThreshold <- 0.8
# Load in BTSL ecological boundary
ecobound <- raster(paste0(processedMapsDir, "b03-study-area-30m.tif"))
# Read in scripts for connectivity analyses
source(paste0(RscriptDir,"ECstandard.R"))
source(paste0(RscriptDir,"ECNodeImportance_standard.R"))

# Set up output tables
mpgAnalysisTime<-data.frame(Species=character(), NumPixels=integer(), NumNodes=integer(), NumLinks=integer(), 
                            NumNodesClip=integer(), NumLinksClip=integer(), extractSeconds=double(), clipSeconds=double(), 
                            summarySeconds=double(), btwnSeconds=double(), summaryBTSLSeconds=double(), btwnBTSLSeconds=double(), 
                            dECSeconds=double(), stringsAsFactors = FALSE)

#table of summary statistics about the mpg (habitat network) for the full extent
netsum<-data.frame(Species=character(), numNodes=integer(), numLinks=integer(), totArea=double(), meanArea=double(),
                   totQuality=double(), meanQuality=double(), totAreaQuality=double(), meanAreaQuality=double(),
                   totCost=double(), meanCost=double(), overallECgap=double(), overallECnatal=double())
#                   modularity=double(), numCommunities=integer(), meanNodesperComm=double(), meanAreaperComm=double(), meanAreaQualityperComm=double())

#table of summary statistics about the mpg (habitat network) within the BTSL
netsumBTSL<-data.frame(Species=character(), numNodes=integer(), numLinks=integer(), totArea=double(), meanArea=double(),
                       totQuality=double(), meanQuality=double(), totAreaQuality=double(), meanAreaQuality=double(),
                       totCost=double(), meanCost=double(), overallECgap=double(), overallECnatal=double())
#                   modularity=double(), numCommunities=integer(), meanNodesperComm=double(), meanAreaperComm=double(), meanAreaQualityperComm=double())

############################
# Extract habitat networks #
############################
for(i in 1:length(speciesList)){
  species<-speciesList[i]
  
  PatchMap<-raster(paste0(habitatDir, "/", species, '_habitatPatch_', myResolution, 'm.tif'))
  CostMap<-raster(paste0(resistanceDir, "/", species, '_resistance_', myResolution, 'm.tif'))
  habitatquality<-raster(paste0(habitatDir, "/", species, '_habitatSuitability_', myResolution, 'm.tif'))
  
  ptm <- proc.time()
  mpg<-MPG(cost=CostMap,patch=PatchMap)
  elapsed_time<-proc.time() - ptm
  
  #populate table tracking time for each species
  mpgAnalysisTime[i,"Species"]<-species
  mpgAnalysisTime[i,"NumNodes"]<-nrow(graphdf(mpg)[[1]]$v)
  mpgAnalysisTime[i,"NumLinks"]<-nrow(graphdf(mpg)[[1]]$e)
  mpgAnalysisTime[i,"NumPixels"]<-ncell(PatchMap)
  mpgAnalysisTime[i,"extractSeconds"]<-elapsed_time[3]
  
  #save mpg
  write.csv(graphdf(mpg)[[1]]$v,paste0(networkDir,species,"_patchStatsMPG.csv"),row.names=FALSE)
  write.csv(graphdf(mpg)[[1]]$e,paste0(networkDir,species,"_linkStatsMPG.csv"),row.names=FALSE)
  writeRaster(mpg$patchId,paste0(networkDir,species,"_patchId.tif"), overwrite=TRUE)
  writeRaster(mpg$lcpLinkId,paste0(networkDir,species,"_linkId.tif"), overwrite=TRUE)
  writeRaster(mpg$lcpPerimWeight,paste0(networkDir,species,"_linkWeight.tif"), overwrite=TRUE)
  writeRaster(mpg$patchId,paste0(networkDir,species,"_patchId.asc"), overwrite=TRUE)
  writeRaster(mpg$lcpLinkId,paste0(networkDir,species,"_linkId.asc"), overwrite=TRUE)
  writeRaster(mpg$lcpPerimWeight,paste0(networkDir,species,"_linkWeight.asc"), overwrite=TRUE)
  #save.image(paste0(networkDir,species,"MPG.RData"))
}  
write.csv(mpgAnalysisTime,paste0(networkDir,"AllSpeciesmpgAnalysisTime.csv"),row.names=FALSE)
  
  ######################################
  # Clip habitat networks to ecoregion #
  ######################################
  # If this section is being run after the previous sections then load in the necessary files from disk
 for(i in 1:length(speciesList)){
     species <- speciesList[i]
     original_patchid<-raster(paste0(networkDir, species, "_patchId.asc"))
     linkid<-raster(paste0(networkDir, species, "_linkId.asc"))
     linkStats<-read.csv(paste0(networkDir, species, "_linkStatsMPG.csv"), header=TRUE)
     patchStats<-read.csv(paste0(networkDir, species, "_patchStatsMPG.csv"), header=TRUE)
     
    ptm <- proc.time()
    
     #Set patchid, linkid, patch and link stats for full extent
#     original_patchid<-mpg$patchId
#     linkid<-mpg$lcpLinkId
#     patchStats<-graphdf(mpg)[[1]]$v
#     linkStats<-graphdf(mpg)[[1]]$e
    
    #clip patchid based on ecological boundary
    ecol_patchid<-original_patchid*ecobound
    #which nodes are not in ecological region?
    patchidNOTecol_list<-setdiff(unique(original_patchid), unique(ecol_patchid))
    
    #area of each patch inside ecological region
    count_in<-freq(ecol_patchid)
    #area of each patch outside ecological region
    replacelist<-matrix(c(NA,1,1,NA),ncol=2)
    region_OUTecol<-raster::reclassify(ecobound,replacelist)
    patchid_OUTecol<-original_patchid*region_OUTecol
    count_out<-freq(patchid_OUTecol)
    
    #which patches straddle the boundary?
    overlap_patches_in<-count_in[count_in[,1] %in% count_out[,1],]
    overlap_patches_out<-count_out[count_out[,1] %in% count_in[,1],]
    
    #calculate their percernt area inside the ecological region
    overlap_patches<-data.frame(id=overlap_patches_in[,1],incount=overlap_patches_in[,2],outcount=overlap_patches_out[,2])
    overlap_patches$percentin<-overlap_patches$incount/(overlap_patches$incount+overlap_patches$outcount)
    
    #identify patches to remove if they have less than clipAreaThreshold proportion of their area inside the ecological region
    remove_patches_overlap<-overlap_patches[overlap_patches$percentin<clipAreaThreshold,]
    keep_patches_overlap<-overlap_patches[overlap_patches$percentin>=clipAreaThreshold,]
    remove_patches_notoverlap<-setdiff(patchidNOTecol_list, keep_patches_overlap$id)
    remove_patches<-union(remove_patches_notoverlap, remove_patches_overlap$id)
    
    #Links to delete
    node1_rowids<-which(linkStats$e1 %in% remove_patches)
    node2_rowids<-which(linkStats$e2 %in% remove_patches)
    rowids_del<-union(node1_rowids, node2_rowids)
    #draw clipped link layer
    replacelist<-cbind(-linkStats[rowids_del,"linkId"],NA)
    linkidmap_ecol<-raster::reclassify(linkid,replacelist)
    linkStats_ecol<-linkStats[-rowids_del,]

    #draw clipped patch layer with patches on the boundary removes if their area inside ecol region is < 0.8
    #also removes patches that were only connected to other patches in the ecoregion indirectly via paths outside the ecoregion
    replacelist<-cbind(c(0,setdiff(unique(original_patchid), unique(c(linkStats_ecol$e1, linkStats_ecol$e2)))), NA)
    patchidmap_ecol<-raster::reclassify(original_patchid, replacelist)    

    #Define graph object using data.frame
    landscape.graph<-graph.data.frame(data.frame(linkStats_ecol[,c("e2","e1")]),directed=FALSE) #,vertices=patchquality)
    #find biggest cluster
    C<-clusters(landscape.graph)
    if(C$no > 1){
      bigclus<-which(C$csize==max(C$csize))
      bigClusNodes<-as.numeric(V(landscape.graph)$name[which(C$membership==bigclus)])
      
      #Define graph object in which all the nodes belong to a single cluster
      landscape.graph.NotInBigClus<-delete.vertices(landscape.graph,which(C$membership==bigclus))
      linksNotInBigClus<-matrix(as.numeric(get.edgelist(landscape.graph.NotInBigClus)), ecount(landscape.graph.NotInBigClus), 2)
      
      links_to_delete<-c(linkStats_ecol$linkId[which(linkStats_ecol$e1 %in% linksNotInBigClus[,1] & linkStats_ecol$e2 %in% linksNotInBigClus[,2])],
                             linkStats_ecol$linkId[which(linkStats_ecol$e2 %in% linksNotInBigClus[,1] & linkStats_ecol$e1 %in% linksNotInBigClus[,2])])
      linkStats_ecol<-linkStats_ecol[-which(linkStats_ecol$linkId %in% links_to_delete),]
      
      #draw clipped link layer
      replacelist<-cbind(-links_to_delete,NA)
      linkidmap_ecol<-raster::reclassify(linkidmap_ecol,replacelist)
      
      #draw clipped patchid layer
      #remove patches that are not in bigCluster and patches that are connected to the big cluster via patches outside the ecoregion
#      peninsula_patches<-setdiff(unique(patchidmap_ecol),unique(c(linkStats_ecol$e1, linkStats_ecol$e2)))
#      replacelist<-cbind(c(peninsula_patches, as.numeric(V(landscape.graph)$name[which(C$membership!=bigclus)])), NA)
      replacelist<-cbind(as.numeric(V(landscape.graph)$name[which(C$membership!=bigclus)]), NA)
      patchidmap_ecol<-raster::reclassify(patchidmap_ecol,replacelist)
    }
    
    patchStats_ecol<-patchStats[which(patchStats$patchId%in%unique(c(linkStats_ecol$e1, linkStats_ecol$e2))),]
    elapsed_time<-proc.time() - ptm
    
    #populate table tracking time for each species
    mpgAnalysisTime[i,"Species"]<-species
    mpgAnalysisTime[i,"NumNodesClip"]<-nrow(patchStats_ecol)
    mpgAnalysisTime[i,"NumLinksClip"]<-nrow(linkStats_ecol)
    mpgAnalysisTime[i,"clipSeconds"]<-elapsed_time[3]

    #Write outputs of clipping
    writeRaster(patchidmap_ecol, filename=paste0(networkDir, species, "_patchId_BTSL.tif"), overwrite=TRUE)
    writeRaster(linkidmap_ecol, filename=paste0(networkDir, species, "_linkId_BTSL.tif"), overwrite=TRUE)
    write.csv(patchStats_ecol, paste0(networkDir, species, "_patchStatsMPG_BTSL.csv"), row.names=FALSE)
    write.csv(linkStats_ecol, paste0(networkDir, species, "_linkStatsMPG_BTSL.csv"), row.names=FALSE)
}
write.csv(mpgAnalysisTime,paste0(networkDir,"AllSpeciesmpgAnalysisTime.csv"),row.names=FALSE)
    
    ############################
    # Analyze habitat networks #
    ############################
    # If this section is being run after the previous sections then set up output tables and load in the necessary files from disk
for(i in 4:5){#1:length(speciesList)){
    species <- speciesList[i]

  #species-specific estimates of gap-crossing and natal median dispersal distances
  d50GAP<-DISP$Gap[DISP$Species == species]
  d50NATAL<-DISP$Natal[DISP$Species == species]
  
  #dispersal distance coefficients
  coefficientGAP<-log(0.5)/d50GAP
  coefficientNATAL<-log(0.5)/d50NATAL
  
  # Full extent network calculations
  # nodes<-patchStats
  # links<-linkStats
  # patchId<-original_patchid
  nodes<-read.csv(paste0(networkDir,species,"_patchStatsMPG.csv"),header=TRUE)
  links<-read.csv(paste0(networkDir,species,"_linkStatsMPG.csv"),header=TRUE)
  patchId<-raster(paste0(networkDir,species,"_patchId.asc"))

  #produce patch-level summary of habitat quality
  habitatquality<-raster(paste0(habitatDir, species, "_habitatSuitability_30m.tif"))
  nodeQuality<-data.frame(zonal(habitatquality,patchId,fun='mean'))
  
  #add in node quality
  nodes<-merge(nodes, nodeQuality, by.x="patchId", by.y="zone")
  nodes<-data.frame(name=nodes$patchId, area=nodes$patchArea, quality=nodes$mean, areaquality=(nodes$patchArea*nodes$mean/100))
  
  #define graph object using data.frame
  landscape.graph<-graph.data.frame(links, directed=FALSE, vertices=nodes)
  
  # Network Summary Statistics
  ptm <- proc.time()
  netsum[i, "Species"]<-species
  netsum[i, "numNodes"]<-vcount(landscape.graph)
  netsum[i, "numLinks"]<-ecount(landscape.graph)
  netsum[i, "totArea"]<-sum(as.numeric(V(landscape.graph)$area))
  netsum[i, "meanArea"]<-mean(as.numeric(V(landscape.graph)$area))
  netsum[i, "totQuality"]<-sum(as.numeric(V(landscape.graph)$quality))
  netsum[i, "meanQuality"]<-mean(as.numeric(V(landscape.graph)$quality))
  netsum[i, "totAreaQuality"]<-sum(as.numeric(V(landscape.graph)$areaquality))
  netsum[i, "meanAreaQuality"]<-mean(as.numeric(V(landscape.graph)$areaquality))
  netsum[i, "totCost"]<-sum(as.numeric(E(landscape.graph)$lcpPerimWeight))
  netsum[i, "meanCost"]<-mean(as.numeric(E(landscape.graph)$lcpPerimWeight))
  #netsum[i, "overallECgap"]<-overall.indices.standard(landscape.graph, coefficientGAP,'lcpPerimWeight','areaquality')
  #netsum[i, "overallECnatal"]<-overall.indices.standard(landscape.graph, coefficientNATAL,'lcpPerimWeight','areaquality')
  
  # ### edge.betweenness.community
  # ebc <- edge.betweenness.community(landscape.graph, modularity=TRUE, membership=TRUE)
  #
  # netsum[i, "modularity"]<-max(ebc$modularity)
  # netsum[i, "numCommunities"]<-length(ebc)
  # netsum[i, "meanNodesperComm"]<-mean(sizes(ebc))
  #
  # V(landscape.graph)$comm<-membership(ebc)
  # Vmatrix<-data.frame(V(landscape.graph)$patchId, V(landscape.graph)$area, V(landscape.graph)$areaquality, V(landscape.graph)$comm)
  # names(Vmatrix)<-c("v1", "v2", "v3", "v4")
  # netsum[i, "meanAreaperComm"]<-mean(ddply(Vmatrix, "v4", function(df) sumarea=sum(df$v2))[,2])
  # netsum[i, "meanAreaQualityperComm"]<-mean(ddply(Vmatrix, "v4", function(df) sumarea=sum(df$v3))[,2])
  elapsed_time<-proc.time() - ptm
  #populate table tracking time for each species
  mpgAnalysisTime[i,"summarySeconds"]<-elapsed_time[3]
  
  # Betweenness
  #tabular betweenness output
  ptm <- proc.time()
  btwn<-data.frame(patchId=as.numeric(V(landscape.graph)$name), btwn=betweenness(landscape.graph, weights=E(landscape.graph)$lcpPerimWeight, directed=FALSE))
  
  #raster betweenness output raw
  #replace patch ids with betweeness values
  btwnMap<-reclassify(patchId, btwn)
  
  #raster betweenness output 0 - 1
  #make a look-up table between patch id and betweenness value
  btwn_lookup<-cbind(patchId=btwn$patchId, btwn=1/(max(btwn$btwn)-min(btwn$btwn))*(btwn$btwn-min(btwn$btwn)))
  #replace patch ids with betweeness values 0 - 1
  btwnMap01<-reclassify(patchId, btwn_lookup)
  elapsed_time<-proc.time() - ptm
  
  #populate table tracking time for each species
  mpgAnalysisTime[i,"btwnSeconds"]<-elapsed_time[3]
  
  #Save betweenness outputs
  btwnName<-paste0(species, "_betweenness.csv")
  btwnMapName<-paste0(species, "_betweenness.tif")
  btwnMapName01<-paste0(species, "_betweenness_01.tif")
  writeRaster(btwnMap, filename=paste0(networkDir, btwnMapName), overwrite=TRUE)
  writeRaster(btwnMap01, filename=paste0(networkDir, btwnMapName01), overwrite=TRUE)
  write.csv(btwn, paste0(networkDir, btwnName), row.names=F)
  
  # BTSL extent
  # nodes<-patchStats_ecol
  # links<-linkStats_ecol
  # patchId<-patchidmap_ecol
  nodes<-read.csv(paste0(networkDir,species,"_patchStatsMPG_BTSL.csv"),header=TRUE)
  links<-read.csv(paste0(networkDir,species,"_linkStatsMPG_BTSL.csv"),header=TRUE)
  patchId<-raster(paste0(networkDir,species,"_patchId_BTSL.tif"))
  
  #produce patch-level summary of habitat quality
  habitatquality<-raster(paste0(habitatDir, species, "_habitatSuitability_30m.tif"))
  nodeQuality<-data.frame(zonal(habitatquality,patchId,fun='mean'))
  
  #add in node quality
  nodes<-merge(nodes, nodeQuality, by.x="patchId", by.y="zone")
  nodes<-data.frame(name=nodes$patchId, area=nodes$patchArea, quality=nodes$mean, areaquality=(nodes$patchArea*nodes$mean/100))
  
  #define graph object using data.frame
  landscape.graph.clipped<-graph.data.frame(links, directed=FALSE, vertices=nodes)
  
  
  # Network Summary Statistics
  ptm <- proc.time()
  netsumBTSL[i, "Species"]<-species
  netsumBTSL[i, "numNodes"]<-vcount(landscape.graph.clipped)
  netsumBTSL[i, "numLinks"]<-ecount(landscape.graph.clipped)
  netsumBTSL[i, "totArea"]<-sum(as.numeric(V(landscape.graph.clipped)$area))
  netsumBTSL[i, "meanArea"]<-mean(as.numeric(V(landscape.graph.clipped)$area))
  netsumBTSL[i, "totQuality"]<-sum(as.numeric(V(landscape.graph.clipped)$quality))
  netsumBTSL[i, "meanQuality"]<-mean(as.numeric(V(landscape.graph.clipped)$quality))
  netsumBTSL[i, "totAreaQuality"]<-sum(as.numeric(V(landscape.graph.clipped)$areaquality))
  netsumBTSL[i, "meanAreaQuality"]<-mean(as.numeric(V(landscape.graph.clipped)$areaquality))
  netsumBTSL[i, "totCost"]<-sum(as.numeric(E(landscape.graph.clipped)$lcpPerimWeight))
  netsumBTSL[i, "meanCost"]<-mean(as.numeric(E(landscape.graph.clipped)$lcpPerimWeight))
  netsumBTSL[i, "overallECgap"]<-overall.indices.standard(landscape.graph.clipped, coefficientGAP,'lcpPerimWeight','areaquality')
  netsumBTSL[i, "overallECnatal"]<-overall.indices.standard(landscape.graph.clipped, coefficientNATAL,'lcpPerimWeight','areaquality')
  
  # ### edge.betweenness.community
  # ebc <- edge.betweenness.community(landscape.graph.clipped, modularity=TRUE, membership=TRUE)
  #
  # netsumBTSL[i, "modularity"]<-max(ebc$modularity)
  # netsumBTSL[i, "numCommunities"]<-length(ebc)
  # netsumBTSL[i, "meanNodesperComm"]<-mean(sizes(ebc))
  #
  # V(landscape.graph.clipped)$comm<-membership(ebc)
  # Vmatrix<-data.frame(V(landscape.graph.clipped)$patchId, V(landscape.graph.clipped)$area, V(landscape.graph.clipped)$areaquality, V(landscape.graph.clipped)$comm)
  # names(Vmatrix)<-c("v1", "v2", "v3", "v4")
  # netsumBTSL[i, "meanAreaperComm"]<-mean(ddply(Vmatrix, "v4", function(df) sumarea=sum(df$v2))[,2])
  # netsumBTSL[i, "meanAreaQualityperComm"]<-mean(ddply(Vmatrix, "v4", function(df) sumarea=sum(df$v3))[,2])
  elapsed_time<-proc.time() - ptm
  #populate table tracking time for each species
  mpgAnalysisTime[i,"summaryBTSLSeconds"]<-elapsed_time[3]
  
  # Betweenness
  
  #tabular betweenness output
  ptm <- proc.time()
  btwn<-data.frame(patchId=as.numeric(V(landscape.graph.clipped)$name), btwn=betweenness(landscape.graph.clipped, weights=E(landscape.graph.clipped)$lcpPerimWeight, directed=FALSE))
  
  #raster betweenness output raw
  #replace patch ids with betweeness values
  btwnMap<-reclassify(patchId, btwn)

  #raster betweenness output 0 - 1
  #make a look-up table between patch id and betweenness value
  btwn_lookup<-cbind(patchId=btwn$patchId, btwn=1/(max(btwn$btwn)-min(btwn$btwn))*(btwn$btwn-min(btwn$btwn)))
  #replace patch ids with betweeness values 0 - 1
  btwnMap01<-reclassify(patchId, btwn_lookup)
  elapsed_time<-proc.time() - ptm
  
  #populate table tracking time for each species
  mpgAnalysisTime[i,"btwnBTSLSeconds"]<-elapsed_time[3]
  
  #Save betweennness outputs
  btwnName<-paste0(species, "_betweenness_BTSL.csv")
  btwnMapName<-paste0(species, "_betweenness_BTSL.tif")
  btwnMapName01<-paste0(species, "_betweenness_BTSL_01.tif")
  writeRaster(btwnMap, filename=paste0(networkDir, btwnMapName), overwrite=TRUE)
  writeRaster(btwnMap01, filename=paste0(networkDir, btwnMapName01), overwrite=TRUE)
  write.csv(btwn, paste0(networkDir, btwnName), row.names=F)
  
  # dEC #
#   if(species %in% c("MAAM","URAM")){
# 
#     #get a list of patches in the ecological region
#     #patchList<-unique(patchId*ecobound)
# 
#     startNode<-1
#     endNode<-startNode+netsumBTSL[i,"numNodes"]-1
# 
#     inputparams<-data.frame(#Node=as.numeric(patchList),
#       Node=as.numeric(V(landscape.graph.clipped)$name[startNode:endNode]),
#       coefficientGAP=as.numeric(coefficientGAP),
#       coefficientNATAL=as.numeric(coefficientNATAL),
#       linkWeight='lcpPerimWeight',
#       nodeWeight='areaquality',
#       EC_Gap=netsumBTSL[i,'overallECgap'],
#       EC_Natal=netsumBTSL[i,'overallECnatal'])
#     ptm <- proc.time()
#     dEC<-data.frame(t(apply(inputparams,1,function(x,y){node.importance(x,y)},y=landscape.graph.clipped)))
# 
#     #Raster dEC output
#     #dEC gap
#     #make a look-up table raw
#     dEC_lookup<-cbind(patchId=dEC$Node, dEC=as.numeric(dEC$dEC_GAP))
#     #replace patch ids with dEC values 0 - 1
#     dEC_GAPMap<-reclassify(patchId, dEC_lookup)
#     #make a look-up table 0 - 1
#     dEC_lookup<-cbind(patchId=dEC$Node, dEC=1/(max(as.numeric(dEC$dEC_GAP))-min(as.numeric(dEC$dEC_GAP)))*(as.numeric(dEC$dEC_GAP)-min(as.numeric(dEC$dEC_GAP))))
#     #replace patch ids with dEC values 0 - 1
#     dEC_GAPMap01<-reclassify(patchId, dEC_lookup)
# 
#     #dEC natal
#     #make a look-up table raw
#     dEC_lookup<-cbind(patchId=dEC$Node, dEC=as.numeric(dEC$dEC_NATAL))
#     #replace patch ids with dEC values
#     dEC_NATALMap<-reclassify(patchId, dEC_lookup)
#     #make a look-up table 0 - 1
#     dEC_lookup<-cbind(patchId=dEC$Node, dEC=1/(max(as.numeric(dEC$dEC_NATAL))-min(as.numeric(dEC$dEC_NATAL)))*(as.numeric(dEC$dEC_NATAL)-min(as.numeric(dEC$dEC_NATAL))))
#     #replace patch ids with dEC values
#     dEC_NATALMap01<-reclassify(patchId, dEC_lookup)
# 
#     write.csv(dEC, paste0(networkDir, species, "_dEC.csv"), row.names=FALSE)
#     writeRaster(dEC_GAPMap, filename=paste0(networkDir, species, "_dECGap_BTSL.tif"), overwrite=TRUE)
#     writeRaster(dEC_NATALMap, filename=paste0(networkDir, species, "_dECNatal_BTSL.tif"), overwrite=TRUE)
#     writeRaster(dEC_GAPMap01, filename=paste0(networkDir, species, "_dECGap_BTSL_01.tif"), overwrite=TRUE)
#     writeRaster(dEC_NATALMap01, filename=paste0(networkDir, species, "_dECNatal_BTSL_01.tif"), overwrite=TRUE)
#     elapsed_time<-proc.time() - ptm
# 
#     #populate table tracking time for each species
#     mpgAnalysisTime[i,"dECSeconds"]<-elapsed_time[3]
#   }
# }

mpgAnalysisTime$extractMinutes<-mpgAnalysisTime$extractSeconds/60
mpgAnalysisTime$clipMinutes<-mpgAnalysisTime$clipSeconds/60
mpgAnalysisTime$summaryMinutes<-mpgAnalysisTime$summarySeconds/60
mpgAnalysisTime$btwnMinutes<-mpgAnalysisTime$btwnSeconds/60
mpgAnalysisTime$dECMinutes<-mpgAnalysisTime$dECSeconds/60

write.csv(mpgAnalysisTime,paste0(networkDir,"AllSpeciesMPGAnalysisTime.csv"),row.names=FALSE)
write.csv(netsum,paste0(networkDir,"NetworkSummaryBR1.csv"),row.names=FALSE)
write.csv(netsumBTSL,paste0(networkDir,"NetworkSummary_BTSLBR1.csv"),row.names=FALSE)
