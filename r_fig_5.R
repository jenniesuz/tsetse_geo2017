library("rgdal")
library("maptools")
library("sp")
library("rgeos")
library("proj4")
library("raster")

source("r_2015_trapdist.R") # calculate distance of traps from protected area boundary
#*******************Attach 2015 predicted values to 2015 observed values and distance to boundary************************
#************************Data*************************************
gp.feb15 <- readShapePoints("gp.feb15.shp",CRS("+proj=utm +zone=36 +north 
                            +datum=WGS84 +units=m +no_defs"))
gs.feb15 <- readShapePoints("gs.feb15.shp",CRS("+proj=utm +zone=36 +north 
                                               +datum=WGS84 +units=m +no_defs"))



gp.feb10<- readShapePoints("gp.feb10.shp", CRS("+proj=utm +zone=36 +north 
                                               +datum=WGS84 +units=m +no_defs")) 
gs.feb10<- readShapePoints("gs.feb10.shp", CRS("+proj=utm +zone=36 +north 
                                               +datum=WGS84 +units=m +no_defs"))

#*******************2010 residuals***************************
gpmod <- calc.in.buffer(file="model_predvals_gp10.tif"    # Glossina pallidipes
                        ,summary.func=mean,point.locations=coordinates(gp.feb10)
                        ,buffer.val=500)
gpmod <- gpmod$summarise.values


gsmod <- calc.in.buffer(file="model_predvals_gs10.tif"  # Glossina swynnertoni
                        ,summary.func=mean,point.locations=coordinates(gs.feb10),buffer.val=500)
gsmod <- gsmod$summarise.values

#*******************************************************
tiff("S2Fig.tiff", height = 4, width = 5, units = 'in', compression="lzw", res=300)

par(mfrow=c(1,2),mar=c(4,4,1,0.5),cex=0.7)                    
plot(log10(gp.feb10$men_cnt)
     ,log10(gp.feb10$men_cnt) - log10(gpmod),bty="n"
     ,xlab=expression(paste(Log[10]," observed"))
     ,ylab= expression(paste(Log[10]," observed - predicted"))
     ,ylim=c(-2,0.5)
     ,xlim=c(-0.5,2.5)
     ,main="a")
plot(log10(gs.feb10$men_cnt)
     ,log10(gs.feb10$men_cnt) - log10(gsmod),bty="n"
     ,xlab=expression(paste(Log[10]," observed"))
     ,ylab= expression(paste(Log[10]," observed - predicted"))
     ,ylim=c(-2,0.5)
     ,xlim=c(-0.5,2.5)
     ,main="b")
dev.off()
#**************Fig ******************************

#*********************Model predicted values********************************************
gpmod <- calc.in.buffer(file="model_predvals_gp10.tif"    # Glossina pallidipes
                        ,summary.func=mean,point.locations=coordinates(gp.feb15)
                        ,buffer.val=500)
gpmod <- gpmod$summarise.values
#******************95% credible intervals***********************
gplower <- calc.in.buffer(file="model_pred0.25_gp10.tif"    # Glossina pallidipes
                          ,summary.func=mean,point.locations=coordinates(gp.feb15)
                          ,buffer.val=500)
gplower<- gplower$summarise.values

gpupper <- calc.in.buffer(file="model_pred0.975_gp10.tif"    # Glossina pallidipes
                          ,summary.func=mean,point.locations=coordinates(gp.feb15)
                          ,buffer.val=500)
gpupper<- gpupper$summarise.values
#*****gs*******
gsmod <- calc.in.buffer(file="model_predvals_gs10.tif"  # Glossina swynnertoni
                        ,summary.func=mean,point.locations=coordinates(gs.feb15),buffer.val=500)
gsmod <- gsmod$summarise.values

gslower <- calc.in.buffer(file="model_pred0.25_gs10.tif"    # Glossina pallidipes
                          ,summary.func=mean,point.locations=coordinates(gp.feb15)
                          ,buffer.val=500)
gslower<- gslower$summarise.values

gsupper <- calc.in.buffer(file="model_pred0.975_gs10.tif"    # Glossina pallidipes
                          ,summary.func=mean,point.locations=coordinates(gp.feb15)
                          ,buffer.val=500)
gsupper<- gsupper$summarise.values

#******************Attach trap distances and predicted values********************************
gp.feb15@data$trap.dist <- pallid@data$trap.dist
gp.feb15@data$pa <- pallid@data$pa
gp.feb15@data$pred <- gpmod
gp.feb15@data$upper <- gpupper
gp.feb15@data$lower <- gplower

gs.feb15@data$pred <- gsmod # attach the predicted values to the observed dataset.
gs.feb15@data$trap.dist <- swyn@data$trap.dist
gs.feb15@data$pa <- swyn@data$pa
gs.feb15@data$lower <- gslower
gs.feb15@data$upper <- gsupper
#***********************PLOTS BY DISTANCE*************************************
plotmoddist <- function(y=gp.feb15$men_cnt,u=gp.feb15$sem_u,l=gp.feb15$sem_l, label="a"){
  #plot(gp.feb15$trap.dist,y/max(y)*100+1,bty="n"
  plot(gp.feb15$trap.dist,y+1,bty="n"
       ,type="n"
       ,log="y"
       ,yaxt="n"
       ,ylim=c(1,1000)
       ,xlab=' '
       ,ylab=' '
       ,main=label)

  segments(gp.feb15$trap.dist,u+1,y1=l+1,
           col="darkgrey")

  points(gp.feb15$trap.dist[1:18],y[1:18]+1,col="#e41a1c",pch=0)
  points(gp.feb15$trap.dist[19:36],y[19:36]+1,col="#377eb8",pch=1)
  points(gp.feb15$trap.dist[37:54],y[37:54]+1,col="#4daf4a",pch=2)
  points(gp.feb15$trap.dist[55:72],y[55:72]+1,col="#984ea3",pch=5)
  
  abline(v=0,lty=2)
  axis(2, at = yticks+1, labels = yticks, col.axis="black", las=2)
}
#**************Fig 4******************************
tiff("Fig5.tiff", height = 8, width = 7.5, units = 'in', compression="lzw", res=300)

par(mfrow=c(3,2),mar=c(3,3,1,0.5),oma=c(2,2,0,2),cex=0.7)                    
yticks=c(0,2,5,10,50,100,200,500,1000)
plotmoddist()
mtext("Mean count per trap per day",outer=T,line=0,side=2,cex=0.8
     ,at=.85 )
plotmoddist(y=gs.feb15$men_cnt,u=gs.feb15$sem_u
            ,l=gs.feb15$sem_l,label="b")
legend("topright",title="Transect",legend=c("GGR","IGRS","IGRN","SNP")
       ,col=c("#e41a1c","#377eb8","#4daf4a","#984ea3"),bty="n"
       ,pch=c(0,1,2,5),cex=1.3)

plotmoddist(y=gp.feb15@data$pred,u=gp.feb15@data$upper
            ,l=gp.feb15@data$lower,label="c")
plotmoddist(y=gs.feb15@data$pred,u=gs.feb15@data$upper
            ,l=gs.feb15@data$lower,label="d")

gp.resids <- log10(gp.feb15@data$men_cnt+1)-log10(gp.feb15@data$pred+1)
plot(gp.feb15$trap.dist,gp.resids
     ,bty="n"
     ,type="n"
     ,xlab=' '
     ,ylab=' '
     ,ylim=c(-3,1)
     ,main="e")
points(gp.feb15$trap.dist[1:18],gp.resids[1:18],col="#e41a1c",pch=0)
points(gp.feb15$trap.dist[19:36],gp.resids[19:36],col="#377eb8",pch=1)
points(gp.feb15$trap.dist[37:54],gp.resids[37:54],col="#4daf4a",pch=2)
points(gp.feb15$trap.dist[55:72],gp.resids[55:72],col="#984ea3",pch=5)
mtext(expression(paste(Log[10]," observed - predicted"))
      ,outer=T,line=0,side=2,at=.18,cex=0.8)
abline(v=0,lty=2)

gs.resids <- log10(gs.feb15@data$men_cnt+1)-log10(gs.feb15@data$pred+1)
plot(gp.feb15$trap.dist,gs.resids
     ,bty="n"
     ,type="n"
     ,xlab=' '
     ,ylab=' '
     ,ylim=c(-3,1)
     ,main="f")
points(gp.feb15$trap.dist[1:18],gs.resids[1:18],col="#e41a1c",pch=0)
points(gp.feb15$trap.dist[19:36],gs.resids[19:36],col="#377eb8",pch=1)
points(gp.feb15$trap.dist[37:54],gs.resids[37:54],col="#4daf4a",pch=2)
points(gp.feb15$trap.dist[55:72],gs.resids[55:72],col="#984ea3",pch=5)
abline(v=0,lty=2)



mtext("Predicted abundance",outer=T,line=0,side=2,padj=0,cex=0.8)
mtext("Distance (m) from protected area boundary"
      ,outer=T,line=0,side=1,padj=0,cex=0.8)
dev.off()



