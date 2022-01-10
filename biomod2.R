
library(sp)
library(raster)
library(parallel)
library(reshape)
library(ggplot2)
library(biomod2)
library(rgdal)


setwd("E:/ty_jm")
DataSpecies <- read.csv("sord_csv/100_rarefied_points.csv")
head(DataSpecies)
myRespName <-'Calamagrostis.purpurea' 
myResp <- as.numeric(DataSpecies[,myRespName])
myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]



myExpl = stack(raster( "bioclimate/current/bio1.tif.tif"),
               raster( "bioclimate/current/bio2.tif.tif"),
               raster( "bioclimate/current/bio3.tif.tif"),
               raster( "bioclimate/current/bio4.tif.tif"),
              
               raster( "bioclimate/current/bio8.tif.tif"),
             
               raster( "bioclimate/current/bio10.tif.tif"),
               
               raster( "bioclimate/current/bio12.tif.tif"),
              
               raster( "bioclimate/current/bio15.tif.tif"),
              
               raster( "bioclimate/current/aspect.tif.tif"),
            
               raster( "bioclimate/current/slope.tif.tif"))

setwd("E:/ty_jm/80/100/BC26")

myBiomodData <- BIOMOD_FormatingData(resp.var = myResp, 
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName,
                                     PA.nb.rep = 3,  
                                     PA.nb.absences = 2638,
                                     PA.strategy = 'random')

write.csv(myBiomodData@coord,file = "./PA_BC26.csv")

myBiomodData

plot(myBiomodData)
myBiomodOption <- BIOMOD_ModelingOptions()



myBiomodModelOut <- BIOMOD_Modeling( 
                      myBiomodData,   
                      models = c("GLM", "GBM",  "CTA", "ANN",  "FDA", "RF",
                                 "MAXENT.Phillips"),
                      models.options = myBiomodOption, 
                      NbRunEval=5,  
                      DataSplit=70,
                      Prevalence=0.5,
                      VarImport=3, 
                      models.eval.meth = c('TSS','ROC','KAPPA'),
                      SaveObj = TRUE, 
                      rescal.all.models = TRUE,
                      do.full.models = FALSE,
                      modeling.id = paste(myRespName,"FirstModeling",sep=""))



myBiomodModelOut
myBiomodModelEval <- get_evaluations(myBiomodModelOut)
dimnames(myBiomodModelEval)
write.csv(myBiomodModelEval,file = "./myBiomodModelEval_BC26.csv")
KAPPA=myBiomodModelEval["KAPPA","Testing.data",,,]
write.csv(KAPPA,"./KAPPA_BC26.CSV")
ROC=myBiomodModelEval["ROC","Testing.data",,,]
write.csv(ROC,"./ROC_BC26.CSV")
TSS=myBiomodModelEval["TSS","Testing.data",,,]
write.csv(TSS,"./TSS_BC26.CSV")
varimportan=get_variables_importance(myBiomodModelOut)
write.csv(varimportan,"./varimportan_BC26.csv")


myBiomodEM <- BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut, 
                                       chosen.models = 'all', 
                                       em.by='all', 
                                       eval.metric = c('TSS','ROC'), 
                                       eval.metric.quality.threshold = c(0.6,0.8),
                                       prob.mean = T, 
                                       prob.cv = T, 
                                       prob.ci = T,
                                       prob.ci.alpha = 0.05, 
                                       prob.median = T, 
                                       committee.averaging = T, 
                                       prob.mean.weight = T,
                                       prob.mean.weight.decay = 'proportional' )

myBiomodEM




EM_evaluations=get_evaluations(myBiomodEM)
write.csv(EM_evaluations,"/.EM_evaluations_BC26.csv")
myBiomodProj <- BIOMOD_Projection( modeling.output = myBiomodModelOut, 
                                   new.env = myExpl, 
                                   proj.name = 'current', 
                                   selected.models = 'all', 
                                   binary.meth = 'TSS', 
                                   compress = 'xz', 
                                   build.clamping.mask = FALSE,
                                   output.format='.img')



myBiomodProj
setwd("E:/ty_jm/80/100/BC26")
list.files("E:/ty_jm/80/100/BC26/proj_curent")
plot(myBiomodProj)
myCurrentProj <- get_predictions(myBiomodProj)
myCurrentProj


myBiomodEF <- BIOMOD_EnsembleForecasting(projection.output = myBiomodProj,
                                         EM.output = myBiomodEM)


myBiomodEF
plot(myBiomodEF)




setwd("E:/ty_jm")

myExplFuture = stack(raster( "bioclimate/future/BC26/bio1.tif.tif"),
                     raster( "bioclimate/future/BC26/bio2.tif.tif"),
                     raster( "bioclimate/future/BC26/bio3.tif.tif"),
                     raster( "bioclimate/future/BC26/bio4.tif.tif"),
                    
                     raster( "bioclimate/future/BC26/bio8.tif.tif"),
                  
                     raster( "bioclimate/future/BC26/bio10.tif.tif"),
                    
                     raster( "bioclimate/future/BC26/bio12.tif.tif"),
                   
                     raster( "bioclimate/future/BC26/bio15.tif.tif"),
                   
                     raster( "bioclimate/future/BC26/aspect.tif.tif"),
                   
                     raster( "bioclimate/future/BC26/slope.tif.tif"))
setwd("E:/ty_jm/80/100/BC26")

myBiomodProjFuture <- BIOMOD_Projection( 
                      modeling.output = myBiomodModelOut,
                      new.env = myExplFuture,
                      proj.name = 'future', 
                      selected.models = 'all', 
                      binary.meth = 'TSS', 
                      compress = 'xz', 
                      build.clamping.mask = FALSE,
                      clamping.mask = T, 
                      output.format = '.img')





plot(myBiomodProjFuture)
myBiomodEF_future <- BIOMOD_EnsembleForecasting(projection.output = myBiomodProjFuture,
                                         EM.output = myBiomodEM)


myBiomodEF_future
plot(myBiomodEF_future)

setwd("E:/ty_jm/80/100/BC26/Calamagrostis.purpurea/proj_current")
raster_<- stack('proj_current_Calamagrostis.purpurea_ensemble.grd')
writeRaster(raster_, file="proj_current.asc", format="ascii", overwrite=TRUE, bylayer=TRUE, suffix=names(raster_))

setwd("E:/ty_jm/80/100/BC26/Calamagrostis.purpurea/proj_future")
raster_<- stack('proj_future_Calamagrostis.purpurea_ensemble.grd')
writeRaster(raster_, file="proj_future.asc", format="ascii", overwrite=TRUE, bylayer=TRUE, suffix=names(raster_))

setwd("E:/ty_jm/80/100/BC26")