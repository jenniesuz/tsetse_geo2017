#*****************Figures 2 and 3***********************
library("geostatsp")
library("INLABMA") 
library("INLA")
library("rgdal")
library("maptools")
library("sp")
library("rgeos")
library("proj4")
library("graphics")
library("spdep")
library("raster")
#**********count data******************
gp.feb10<- readShapePoints("gp.feb10.shp", CRS("+proj=utm +zone=36 +north 
                                               +datum=WGS84 +units=m +no_defs")) 
gs.feb10<- readShapePoints("gs.feb10.shp", CRS("+proj=utm +zone=36 +north 
                                                +datum=WGS84 +units=m +no_defs"))
#***********landsat variables and elevation***************
elev <- raster("env_var_elev_mean500m.tif") # create list of covariates
b7 <- raster("env_var_B7_mean500m.tif")
ndvi <- raster("env_var_ndvi_mean500m.tif")
temp <- raster("env_var_temp_mean500m.tif")
covlist <- list(b7=b7
                ,ndvi=ndvi
                ,elev=elev
                ,temp=temp)
#******with spatially structured errors*********
#*******using GLGM*************
gp.feb10.fit <- glgm(formula=log10(men_cnt+1) ~ ndvi+b7+temp+elev,data=gp.feb10,grid=ndvi,covariates=covlist,family="gaussian",
                     priorCI=list(sd=c(0.01,5),range=c(500,5000),sdNugget=c(0.01,5)),
                     control.compute=list(dic=T,cpo=T))

sum(log(gp.feb10.fit$inla$cpo[[1]])) # -25.5
#write.csv(gp.feb10.fit$inla$cpo[[1]],"model_gpCPOscores.csv")


cposcores <- read.csv("mode_gpCPOscores.csv",header=T)
plot(gp.feb10$trap_id,log(cposcores[,2]),bty="n")


gp.feb10.fit$inla$.args$data
#write.csv(gp.feb10.fit$inla$.args$data,file="gp.feb10.fit.inla.args.data.csv")

gp.feb10.fit$inla$dic
gp.feb10.fit$parameters$summary # DIC 39.9

dat <- gp.feb10.fit$raster$predict.mean
predict.0.025 <- gp.feb10.fit$raster$predict.0.025quant
predict.0975 <- gp.feb10.fit$raster$predict.0.975quant

coefs <- gp.feb10.fit$parameters$summary

writeRaster(gp.feb10.fit$raster,"model_gp10.tif",overwrite=T)
writeRaster(gp.feb10.fit$raster$predict.mean,"model_logpredvals_gp10.tif",overwrite=T)
writeRaster(predict.0.025,"model_logpred0.25_gp10.tif",overwrite=T)
writeRaster(predict.0975,"model_logpred0.975_gp10.tif",overwrite=T)

gpraster <- stack("model_gp10.tif")
writeRaster(10^gpraster[[12]],"model_gp10.sd.tif",overwrite=T)

gsraster <- stack("model_gs10.tif")
writeRaster(10^gsraster[[12]],"model_gs10.sd.tif",overwrite=T)


predict.0.025 <- raster("model_logpred0.25_gp10.tif")
predict.0975 <- raster("model_logpred0.975_gp10.tif")

temp <- predict.0975 - predict.0.025
reclassify(temp,cbind(5,14,NA),"model_uncertainty_gp10.tif",overwrite=T)
temp <- raster("model_uncertainty_gp10.tif")
writeRaster(temp,"model_uncertainty_gp10.tif",overwrite=T)

#*************With unstructured errors***************************
gp.feb10.fitb <- inla(log10(men_cnt+1) ~ ndvi+elev+b7+temp+f(rep(1,60),model="iid"), data=gp.feb10.fit$inla$.args$data, family='gaussian'
                      ,control.compute=list(dic=T,cpo=T))
# DIC 66.88
sum(log(gp.feb10.fitb$cpo[[1]])) # -37.54

temp <- raster("model_logpred0.975_gp10.tif")
writeRaster(10^temp,"model_pred0.975_gp10.tif",overwrite=T)

# [1] "random.ID"          "random.mean"        "random.sd"          "random.0.025quant"  "random.0.5quant"    "random.0.975quant"  "random.mode"       
# [8] "random.kld"         "random.exp"    


# #******************G. swynnertoni**************************
#******with spatially structured errors*********
#*******using GLGM*************
gs.feb10.fit <- glgm(formula=log10(men_cnt+1) ~ ndvi+elev+b7,data=gs.feb10,grid=ndvi,covariates=covlist,family="gaussian",
                     priorCI=list(sd=c(0.01,5),range=c(500,5000),sdNugget=c(0.01,5)),control.compute=list(dic=T,cpo=T))


sum(log(gs.feb10.fit$inla$cpo[[1]])) # - 8.4

gs.feb10.fit$inla$.args$data
#write.csv(gs.feb10.fit$inla$.args$data,file="gs.feb10.fit.inla.args.data.csv")

summary(gs.feb10.fit$inla) # DIC 13.37
gs.feb10.fit$parameters$summary
#write.csv(gs.feb10.fit$parameters$summary,file="gs.feb10.fit.parameterssummary.csv")

dat <- gs.feb10.fit$raster$predict.mean
predict.0.025 <- gs.feb10.fit$raster$predict.0.025quant
predict.0975 <- gs.feb10.fit$raster$predict.0.975quant

writeRaster(10^(predict.0975-predict.0.025),"model_uncertainty_gs10.tif",overwrite=T)


coefs <- gs.feb10.fit$parameters$summary
writeRaster(gs.feb10.fit$raster,"model_gs10.tif",overwrite=T)
writeRaster(gs.feb10.fit$raster$predict.mean,"model_logpredvals_gs10.tif",overwrite=T)
writeRaster(predict.0.025,"model_logpred0.25_gs10.tif",overwrite=T)
writeRaster(predict.0975,"model_logpred0.975_gs10.tif",overwrite=T)

writeRaster(gs.feb10.fit$raster$random.mean,"model_lograndommean_gs10.tif",overwrite=T)
writeRaster(gs.feb10.fit$raster$random.sd,"model_lograndomsd_gs10.tif",overwrite=T)
writeRaster(gs.feb10.fit$raster$predict.sd,"model_logpredictsd_gs10.tif",overwrite=T)


writeRaster(10^gs.feb10.fit$raster$predict.mean,"model_predvals_gs10.tif",overwrite=T)
#write.csv(coefs,"model_coefslog_gs10.csv")

predict.0.025 <- raster("model_logpred0.25_gs10.tif")
predict.0975 <- raster("model_logpred0.975_gs10.tif")

temp <- predict.0975 - predict.0.025
writeRaster(temp,"model_uncertainty_gs10.tif",overwrite=T)

temp <- raster("model_logpred0.975_gs10.tif")
writeRaster(10^temp,"model_pred0.975_gs10.tif",overwrite=T)

#******with unstructured errors*********
gs.feb10.fitb <- inla(log10(men_cnt+1) ~ ndvi+elev+b7+f(rep(1,60),model="iid")
                      , data=gs.feb10.fit$inla$.args$data, family='gaussian',control.compute=list(dic=T,cpo=T))
# DIC 45.88
sum(log(gs.feb10.fitb$cpo[[1]])) # -23.4

