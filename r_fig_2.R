


#******FIGURE 2: Numbers of tsetse caught per trap per day by habitat type***********
#**************2010 data from >10 km inside the Serengeti National Park*************

#********************************************************************
library("rgdal")      # required packages
library("maptools")
library("sp")
library("rgeos")
library("proj4")
library("plyr")
library("raster")
#********************************************************************



#*********Function for calculating standard error of the mean********
sem <- function(x) {  x <- x[! x %in% NA]
                      se.mean <- sd(x)/sqrt(length(x)) 
                      return (se.mean)
}
#*******************************************************************



#***Read in shapefiles containing trap locations and mean counts**********
gp.feb10 <- readShapePoints("gp.feb10.shp"                             # G. pallidipes 2010
                            , CRS("+proj=utm +zone=36 +north 
                                   +datum=WGS84 +units=m +no_defs"))
                            

gs.feb10 <- readShapePoints("gs.feb10.shp"                             # G. swynnertoni 2010
                            , CRS("+proj=utm +zone=36 +north 
                                   +datum=WGS84 +units=m +no_defs"))
                           

gp.feb10$habitat <- factor(gp.feb10$habitat,c("grassland"              # Convert habitat column to factor
                                             ,"savannah"
                                             ,"open_wood"
                                             ,"dense_wood"))
                          

gs.feb10$habitat <- factor(gs.feb10$habitat,c("grassland"
                                             ,"savannah"
                                             ,"open_wood"
                                             ,"dense_wood"))

gp.feb10$habitat[gp.feb10$habitat %in% NA] <- "dense_wood"   # change 'riparian to dense woodland as per Reed et al
gs.feb10$habitat[gs.feb10$habitat %in% NA] <- "dense_wood"
#**************************************************************************



#********************Mean counts by habitat********************************
# G. pallidipes
gplog10meanhabit <- ddply(gp.feb10@data                      # log10 mean 
                          ,.(habitat)
                          ,summarize,mean=mean(log10(men_cnt))) 

gplog10semhabit <- ddply(gp.feb10@data                       # log10 sem
                         ,.(habitat),summarize
                         ,sem=sem(log10(men_cnt)))


gpmeanhabit <- 10^gplog10meanhabit$mean-1                          # mean

gpsemhabitu <- 10^(gplog10meanhabit$mean + gplog10semhabit$sem)-1  # sem
gpsemhabitl <- 10^(gplog10meanhabit$mean - gplog10semhabit$sem)-1

10^(gplog10meanhabit$mean + gplog10semhabit$sem)-1 - 10^( gplog10meanhabit$mean)-1

# G. swynnertoni

gslog10meanhabit <- ddply(gs.feb10@data
                          ,.(habitat)
                          ,summarize,mean=mean(log10(men_cnt)))

gslog10semhabit <- ddply(gs.feb10@data
                         ,.(habitat)
                         ,summarize,sem=sem(log10(men_cnt)))


gsmeanhabit <- 10^gslog10meanhabit$mean-1                               # mean

gssemhabitu <- 10^(gslog10meanhabit$mean + gslog10semhabit$sem)-1       # sem
gssemhabitl <- 10^(gslog10meanhabit$mean - gslog10semhabit$sem)-1


10^(gslog10meanhabit$mean + gslog10semhabit$sem)-1 - 10^( gslog10meanhabit$mean)-1

#********************************************************************



#**************Plot************************************
yticks=c(0,2,5,10,25,50,100,300)

tiff("Fig2.tiff", height =3, width = 4, units = 'in', compression="lzw", res=400)
par(mfcol=c(1,2),mar=c(3,3,1,0.5),oma=c(2,2,0,2),cex=0.4)

boxplot(gp.feb10$men_cnt ~ gp.feb10@data$habitat,frame=F,names=c( " "," "," "," ")
        ,main="a",log="y",ylim=c(0.2,500))

mtext("Grassland       Savannah        Open wood         Dense wood"
      ,side=1,line=1,cex=0.3,at=2.5)

mtext("Habitat type",side=1,line=3,cex=0.4) 

boxplot(gs.feb10$men_cnt ~ gp.feb10@data$habitat,frame=F,names=c( " "," "," "," ")
        ,main="b",log="y",ylim=c(0.2,500))

mtext("Mean count per trap per day",outer=T,line=0,side=2,padj=0,cex=0.4)

mtext("Grassland       Savannah        Open wood         Dense wood"
      ,side=1,line=1,cex=0.3,at=2.5)

mtext("Habitat type",side=1,line=3,cex=0.4)  

dev.off()
#****************************************************************************



# ***********Effect of habitat on trap counts************************************
gphabitat <- aov(log10(gp.feb10$men_cnt + 1) ~ gp.feb10$habitat) # effect of habitat on count per trap per day
gshabitat <- aov(log10(gs.feb10$men_cnt + 1) ~ gp.feb10$habitat)
gphabtuk <- TukeyHSD(gphabitat)
gshabtuk <- TukeyHSD(gshabitat)
#*****************************************************************************



