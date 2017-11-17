library("rgdal")
library("maptools")
library("sp")
library("rgeos")
library("proj4")
library("raster")

#******************Calculate distance from protected area boundary************
pallid <- readShapePoints("gp.feb15.shp",CRS("+proj=utm +zone=36 +north  
                                             +datum=WGS84 +units=m +no_defs"))                 # read in summary data - converted to shapefiles via count_shapfiles.R
swyn <- readShapePoints("gs.feb15.shp",CRS("+proj=utm +zone=36 +north 
                                           +datum=WGS84 +units=m +no_defs"))
pas <- readShapePoints("trapsinvsout.shp", CRS("+proj=utm +zone=36 +north 
                                               +datum=WGS84 +units=m +no_defs"))                 # this contains information whether the traps are inside or outside the protected area
pallid$pa <- pas$pa                                                           # attach the column stating whether inside or outside
swyn$pa <- pas$pa
boundary <- readShapeLines("bound_linenew.shp", CRS("+proj=utm +zone=36 +north 
                                                 +datum=WGS84 +units=m +no_defs"))                 # shapefile containing protected area boundary line
trap.dist <- numeric(length(pallid[,1]))                                      # calculate distance of each trap location from the boundary
for(i in 1:length(trap.dist)){
  trap.dist[i] <- gDistance(pallid[i,],boundary)
}
pallid$trap.dist <- trap.dist                                                 # append trap distance from boundary to shapefile dataframe
swyn$trap.dist <- trap.dist
pallid$trap.dist[pallid$pa %in% 1] <- -pallid$trap.dist[pallid$pa %in% 1]     # use the column 'pa' to say whether the trap is inside or outside the protected area
swyn$trap.dist[swyn$pa %in% 1] <- -swyn$trap.dist[swyn$pa %in% 1]
#***************************************************************************
write.csv(swyn@data,"swyn_feb15.csv")
