library(ade4)
library(ape)
library(gbm)
library(sp)
library(ecospat)
library(rgeos)
library(raster)
library(rworldmap)
library(parallel)
library(reshape)
library(ggplot2)
library(rgdal)

setwd("H:/ty_jm/")
DataSpecies1 <- read.csv("sord_csv/100_1.csv")
head(DataSpecies1)
myRespName1 <- 'Calamagrostis.purpurea' 
myResp1 <- as.numeric(DataSpecies1[,myRespName1])
myRespXY1<- DataSpecies1[,c("X_WGS84","Y_WGS84")]

DataSpecies2 <- read.csv("sord_csv/100_2.csv")
head(DataSpecies2)
myRespName2 <- 'Calamagrostis.purpurea' 
myResp2 <- as.numeric(DataSpecies2[,myRespName2])
myRespXY2 <- DataSpecies2[,c("X_WGS84","Y_WGS84")]


myExpl = stack(raster( "bioclimate/current/bio1.tif.tif"),
               raster( "bioclimate/current/bio2.tif.tif"),
               raster( "bioclimate/current/bio3.tif.tif"),
               raster( "bioclimate/current/bio4.tif.tif"),
               
               raster( "bioclimate/current/bio8.tif.tif"),
               
               raster( "bioclimate/current/bio10.tif.tif"),
               
               raster( "bioclimate/current/bio12.tif.tif"),
               
               raster( "bioclimate/current/bio15.tif.tif"))
  


xy.sp1<-subset(myRespXY1)#Bromus_erectus
xy.sp2<-subset(myRespXY2) #Daucus_carota

env.sp1<-extract(myExpl,xy.sp1)
env.sp2<-extract(myExpl,xy.sp2)
env.bkg<-na.exclude(values(myExpl))
#################################### PCA-ENVIRONMENT ##################################
pca.cal <- dudi.pca(env.bkg, center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)

pca.cal$eig[1:4]
round(pca.cal$eig/sum(pca.cal$eig)*100,2)[1:2]
# predict the scores on the axes
scores.bkg <- pca.cal$li #scores for background climate
scores.sp1 <- suprow(pca.cal,env.sp1)$lisup #scores for sp1£¬
scores.sp2 <- suprow(pca.cal,env.sp2)$lisup #scores for sp2
# calculation of occurence density (niche z)
z1 <- ecospat.grid.clim.dyn(scores.bkg, scores.bkg, scores.sp1,R=100)
z2 <- ecospat.grid.clim.dyn(scores.bkg, scores.bkg, scores.sp2,R=100)
plot(z1$z.uncor)
points(scores.sp1)
plot(z2$z.uncor)
points(scores.sp2)




ecospat.niche.overlap(z1,z2 ,cor = TRUE)
#################################### stability S in space ##################################
geozS<-ecospat.niche.dynIndexProjGeo(z1,z2,myExpl,index="stability")
plot(geozS,main="Stability")
points(xy.sp1,col="red")
points(xy.sp2,col="blue")

#####Run a niche equivalency test based on two species occurrence density grids.####
ecospat.niche.equivalency.test (z1, z2, rep=100, alternative="greater", ncores = 1)


#####Run a niche similarity test based on two species occurrence density grids.#####
ecospat.niche.similarity.test (z1, z2, rep=100, alternative = "greater",
                               rand.type = 1, ncores= 1)

######Calculate niche expansion, stability and unfilling########
ecospat.niche.dyn.index (z1, z2, intersection=NA)

###Niche Categories and Species Density###
ecospat.plot.niche.dyn(z1, z2, quant=0.25, interest=2,
                       title= "Niche Overlap", name.axis1="PC1",
                       name.axis2="PC2")

ecospat.shift.centroids(scores.sp1,scores.sp2,env.sp1,env.sp2,col = "red")
