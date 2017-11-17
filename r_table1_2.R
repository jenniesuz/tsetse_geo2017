


#**********Tables 1 and 2********************
#*******Regression using landsat data********

#********************************************
library("rgdal")       # required packages
library("maptools")
library("sp")
library("rgeos")
library("proj4")
library("rms") 
library("graphics")
library("raster")
library("INLABMA", lib.loc="~/R/win-library/3.1")   # this line may been to be amended depending on where package installed
library("INLA", lib.loc="~/R/win-library/3.1")      # this line may need to be amended depending on where package installed
#**********************************************


#***************Function for summarising raster data around each trap*****
calc.in.buffer <- function(file="mosaic_band7.tif"                        # raster file
                           ,buffer.val=500                                # value around point to average over
                           ,summary.func                                  # for this purpose will be average
                           ,point.locations=trap_loc){                    # point locations
  
  r <- raster(file)
  
  extract.values <- extract(r,point.locations,buffer=buffer.val)
  
  summarise.values <- as.numeric(lapply(extract.values,summary.func,na.rm=T)) # sum/ mean etc
  
  return(cbind.data.frame(point.locations,summarise.values))
}
#***********************************************************************



#**********************Data******************************

#******************2010 data from inside Serengeti National Park*******************
gp.feb10<- readShapePoints("gp.feb10.shp", CRS("+proj=utm +zone=36 +north 
                                               +datum=WGS84 +units=m +no_defs"))
gs.feb10<- readShapePoints("gs.feb10.shp", CRS("+proj=utm +zone=36 +north 
                                                +datum=WGS84 +units=m +no_defs"))
#***********************************************************************************

#*******************Calculate average values of band 7, ndvi and elevation around each trap**************************************
elev <- calc.in.buffer(file="env_var_elev_clip.tif",summary.func=mean,point.locations=coordinates(gs.feb10),buffer.val=500)
ndvi <- calc.in.buffer(file="env_var_ndvi.tif",summary.func=mean,point.locations=coordinates(gs.feb10),buffer.val=500)
b7 <- calc.in.buffer(file="env_var_B7_feb15new.tif",summary.func=mean,point.locations=coordinates(gs.feb10),buffer.val=500)
temp <- calc.in.buffer(file="env_var_temp_100m.tif",summary.func=mean,point.locations=coordinates(gs.feb10),buffer.val=500)

temp <- temp$summarise.values
elev <- elev$summarise.values
b7 <- b7$summarise.values
ndvi <- ndvi$summarise.values
#****************************************************************************************************************************



#*********************Model comparison*******************************

#*************G pallidipes****************************
dat <- list(men_cnt=gp.feb10$men_cnt,temp=temp,elev=elev,b7=b7,ndvi=ndvi,trap=seq(from=1,to=length(gp.feb10$men_cnt),1))

b7.fit <- inla(log10(men_cnt+1) ~  b7, data=dat
               , family="gaussian"
               ,control.compute=list(dic=T,cpo=T)) #DIC 97.63 
sum(log(b7.fit$cpo[[1]])) # -48.86919

ndvi.fit <- inla(log10(men_cnt+1) ~ ndvi
                 , data=dat, family="gaussian"
                 ,control.compute=list(dic=T,cpo=T)) #DIC 119.6 
sum(log(ndvi.fit$cpo[[1]])) # -59.58283

temp.fit <- inla(log10(men_cnt+1) ~ temp, data=dat, 
                 family="gaussian",control.compute=list(dic=T,cpo=T)) #DIC 114.57 
sum(log(temp.fit$cpo[[1]])) #-57.80628

elev.fit <- inla(log10(men_cnt+1) ~ elev, data=dat,
                 family="gaussian",control.compute=list(dic=T,cpo=T))
sum(log(elev.fit$cpo[[1]])) 

b7elev.fit <- inla(log10(men_cnt+1) ~ b7+elev, data=dat
                   , family='gaussian'
                   ,control.compute=list(dic=T,cpo=T)) #DIC 95.37 
sum(log(b7elev.fit$cpo[[1]])) # -47.80202

b7elevtemp.fit <- inla(log10(men_cnt+1) ~ b7+temp+elev
                       , data=dat, family='gaussian'
                       ,control.compute=list(dic=T,cpo=T)) # 94 
sum(log(b7elevtemp.fit$cpo[[1]])) #-47.18462

b7elevndvi.fit <- inla(log10(men_cnt+1) ~ b7+ndvi+elev
                       , data=dat, family='gaussian'
                       ,control.compute=list(dic=T,cpo=T)) # 94
sum(log(b7elevndvi.fit$cpo[[1]])) #-33.07

b7elevtempndvi.fit <- inla(log10(men_cnt+1) ~ b7+temp+elev+ndvi, data=dat
                           , family='gaussian',control.compute=list(dic=T,cpo=T)) 
sum(log(b7elevtempndvi.fit$cpo[[1]])) # DIC 64 -32.68554
#*****************************************************************************


#*************G swynnertoni****************************
dat <- list(men_cnt=gs.feb10$men_cnt,temp=temp,elev=elev,b7=b7,ndvi=ndvi)

b7.fit <- inla(log10(men_cnt+1) ~ b7, data=dat
               , family="gaussian"
               ,control.compute=list(dic=T,cpo=T)) #DIC 81.66 
sum(log(b7.fit$cpo[[1]])) # DIC 64  -41.10231


ndvi.fit <- inla(log10(men_cnt+1) ~ ndvi
                 , data=dat, family="gaussian"
                 ,control.compute=list(dic=T,cpo=T)) #DIC 75.3 
sum(log(ndvi.fit$cpo[[1]])) # -37.97368

temp.fit <- inla(log10(men_cnt+1) ~ temp, data=dat, 
                 family="gaussian",control.compute=list(dic=T,cpo=T)) #DIC 83.4 
sum(log(temp.fit$cpo[[1]])) #-42.0198


elev.fit <- inla(log10(men_cnt+1) ~ elev, data=dat, 
                 family="gaussian",control.compute=list(dic=T,cpo=T)) #DIC 83.4 
sum(log(elev.fit$cpo[[1]])) #-42.0198

ndvielev.fit <- inla(log10(men_cnt+1) ~ ndvi+elev, data=dat
                   , family='gaussian'
                   ,control.compute=list(dic=T,cpo=T)) #DIC 47.6 
sum(log(ndvielev.fit$cpo[[1]])) # -29.30074

ndvielevb7.fit <- inla(log10(men_cnt+1) ~ ndvi+elev+b7
                       , data=dat, family='gaussian'
                       ,control.compute=list(dic=T,cpo=T)) # 47 
sum(log(ndvielevb7.fit$cpo[[1]])) #-24.19505

ndvielevtemp.fit <- inla(log10(men_cnt+1) ~ ndvi+elev+temp
                       , data=dat, family='gaussian'
                       ,control.compute=list(dic=T,cpo=T)) # 47 
sum(log(ndvielevtemp.fit$cpo[[1]])) # -27.31647

b7elevtempndvi.fit <- inla(log10(men_cnt+1) ~ b7+temp+elev+ndvi, data=dat
                           , family='gaussian',control.compute=list(dic=T,cpo=T)) 
sum(log(b7elevtempndvi.fit$cpo[[1]])) # -27.31647 # DIC 48 mean 




# #***************applying cloud mask to raster layers**********
# r <- raster("env_var_temp_100m.tif")
# cloud <- readShapePoly("cloud_mask.shp",proj4string=CRS("+proj=utm +zone=36 +north 
#                                                +datum=WGS84 +units=m +no_defs"))
# mask(r,cloud,filename="env_var_temp100m.tif",inverse=T)
# r <- raster("env_var_B7_feb15new.tif")
# reclassify(r,matrix(c(0,NA),ncol=2),filename="env_var_temp_100m.tif",overwrite=T)
# r <- raster("env_var_temp_100m.tif")
# reclassify(r,matrix(c(0,0.12,1,0.12,3,0),ncol=3,byrow=T),filename="env_var_b7rec.tif",overwrite=T)
# 
# reg.points <- readShapePoints("500mregpoints.shp",proj4string=CRS("+proj=utm +zone=36 +north 
#                                                +datum=WGS84 +units=m +no_defs"))
# 
# ndvi_grid <- calc.in.buffer(file="env_var_ndvi.tif",summary.func=mean,point.locations=reg.points,buffer.val=500)
# temp_grid <- calc.in.buffer(file="env_var_temp_100m.tif",summary.func=mean,point.locations=reg.points,buffer.val=500)
# 
# #reg.points$ndvi <- ndvi_grid$summarise.values
# reg.points$temp <- temp_grid$summarise.values
# 
# ##write.csv(data.frame(reg.points$b7,reg.points$elev),file="serreg.points2.csv")
# 
# reg.points2 <- SpatialPointsDataFrame(coordinates(reg.points),
#                                    data=data.frame(reg.points),proj4string=CRS("+proj=utm +zone=36 +north 
#                                      +datum=WGS84 +units=m +no_defs"))
# 
# drv <- "ESRI Shapefile"
# writeOGR(reg.points2[,1:4],layer="regpointstemp",driver=drv,dsn=".")

