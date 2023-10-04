#Set active directory to be the same as where script file resides
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(terra)

focalmax <- function(x, r){
  #establish aggregating factor when radius is too large
  fc = floor(r/res(x)[1]/9+1)
  x1 = x
  if(fc > 1){
    x1 <- aggregate(x, fact = fc, fun = 'max',  na.rm=TRUE)
  }
  #create a focal weights matrix of the appropriate size for given radius
  fm <- focalMat(x1, d=r, type = 'circle')
  #exclude outer portion of circle and ensure max/min values are only mutlplied by 1
  fm.na <- ifelse(fm > 0, 1, NA)
  x1.max <- focal(x1, fm.na, fun='max', na.rm=T)
  #restore resolution in result
  if(fc > 1){
    x1.max <- project(x1.max, x)
  }
  return(x1.max)} 

focalmin <- function(x, r){
  #establish aggregating factor when radius is too large
  fc = floor(r/res(x)[1]/9+1)
  x1 = x
  if(fc > 1){
    x1 <- aggregate(x, fact = fc, fun = 'min',  na.rm=TRUE)
  }
  #create a focal weights matrix of the appropriate size for given radius
  fm <- focalMat(x1, d=r, type = 'circle')
  #exclude outer portion of circle and ensure max/min values are only mutlplied by 1
  fm.na <- ifelse(fm > 0, 1, NA)
  x1.min <- focal(x1, fm.na, fun='min', na.rm=T)
  #restore resolution in result
  if(fc > 1){
    x1.min <- project(x1.min, x)
  }
  return(x1.min)} 

focalmed <- function(x, r){
  #establish aggregating factor when radius is too large
  fc = floor(r/res(x)[1]/9+1)
  x1 = x
  if(fc > 1){
    x1 <- aggregate(x, fact = fc, fun = 'mean',  na.rm=TRUE)
  }
  #create a focal weights matrix of the appropriate size for given radius
  fm <- focalMat(x1, d=r, type = 'circle')
  #exclude outer portion of circle and ensure max/min values are only multiplied by 1
  fm.na <- ifelse(fm > 0, 1, NA)
  x1.mean <- focal(x1, fm.na, fun='median', na.rm=T)
  #restore resolution in result
  if(fc > 1){
    x1.mean <- project(x1.mean, x)
  }
  return(x1.mean)} 


hillpos <- function(xmax, xmin, xmed){#relative slope position
  x.pos <- (dm - xmed)/(xmax - xmed+0.5)
  x.neg <- (dm - xmed)/(xmed - xmin+0.5)
  x.pos <- ifel(x.pos > 0, x.pos,0)
  x.neg <- ifel(x.neg < 0, x.neg,0)
  p <- ((x.pos+x.neg)+1)/2
  return(p)
}
#compound slope position
comphillpos = function(max1,max2,max3,min1,min2,min3,med1,med2,med3){
  x.pos1 <- hillpos(max1,min1,med1)
  x.pos2 <- hillpos(max2,min2,med2)
  x.pos3 <- hillpos(max3,min3,med3)
  x.pos.r1 <- focalmax(x.pos2, 100) - focalmin(x.pos2, 100)
  x.pos.1 <- x.pos1*x.pos.r1 + x.pos2*(x.pos.r1*-1+1)
  x.pos.r2 <- focalmax(x.pos3, 500) - focalmin(x.pos3, 500)
  x.pos <- x.pos.1*x.pos.r2 + x.pos3*(x.pos.r2*-1+1)
  return(x.pos)
}
meanhillpos = function(max1,max2,max3,min1,min2,min3,med1,med2,med3){
  x.pos1 <- hillpos(max1,min1,med1)
  x.pos2 <- hillpos(max2,min2,med2)
  x.pos3 <- hillpos(max3,min3,med3)
  
  x.pos <- (x.pos1+x.pos2+x.pos3)/3
  return(x.pos)
}


dem <- rast('dem.tif')
# dm <- rast('D:/GIS/gsmnp/US_DEM/us_dem/w001001.adf')
#get smaller extent
extnt <- c(xmin= -250000, xmax=000000, ymin=200000, ymax=500000)#lmich
# extnt <- c(xmin= 100000, xmax=300000, ymin=50000, ymax=200000)#nohio
dm <- terra::crop(dem, extnt)

#focal analysis of two different neighborhood radius

xmax = focalmax(dm, 500)

xmin = focalmin(dm, 500)

xmax2 = focalmax(dm, 2500)

xmin2 = focalmin(dm, 2500)

x.med <- focalmed(dm, 500)

x.med2 <- focalmed(dm, 2500)

x.slope <-  tan(terrain(dm, v="slope", unit='radians'))*100

x.slope.med <- focalmed(x.slope, 500)

#local slope position
x.pos <- hillpos(xmax,xmin,x.med)
#mesoscale slope position
x.pos2 <- hillpos(xmax2,xmin2,x.med2)

rel <- xmax-xmin
rel2 <- xmax2-xmin2

#create hillshade 
xslope <- terrain(dm, v='slope', unit='radians')
xaspect <- terrain(dm, v='aspect', unit='radians')
x.shade <- shade(xslope, xaspect)

# #previous analysis keying geomorphic positions: 1.rise, 2.talf, 3.dip, 4.interfluve, 5.crest, 6.toeslope, 7.tread, 8.riser, 9.shoulder, 10.backslope, 11.footslope.
# plot(x.slope.med>1)
# 
# geomorph.flats <- ifel(
#   x.slope.med < 1 & x.slope < 2, 
#   ifel(dm - x.med > 0.5, 1,ifel(dm - x.med > -0.5,2,3)),
#   NA
# )
# 
# plot(x.slope.med < 5)
# 
# geomorph.slopes <- ifel(is.na(geomorph.flats),
#                         ifel(x.rslope %in% 0,
#                              ifel(x.pos2 >= 0.67,
#                                   ifel(x.slope.med < 5,4,5),ifel(x.pos2 < 0.33,6,7)),
#                              NA),
#                         geomorph.flats)
# 
# geomorph <- ifel(is.na(geomorph.slopes),
#                         ifel(x.slope.med < 5,8,
#                              ifel(x.pos >= 0.67, 9, ifel(x.pos >= 0.33, 10,11))),
#                  geomorph.slopes)
# 
# 
# plot(geomorph)
# writeRaster(geomorph, 'geomorph.tif', overwrite=T)
# writeRaster(x.shade, 'x.shade.tif', overwrite=T)

##########alternative geomorphic analysis
#relative position giving more weight to local position where the contrast in mesoscale position is highest.
x.pos.r <- focalmax(x.pos2, 500) - focalmin(x.pos2, 500)
x.pos4 <- x.pos*x.pos.r + x.pos2*(x.pos.r*-1+1)
x.pos.class <- ifel(x.pos4 >= 0.67, 3,ifel(x.pos4 >= 0.33, 2,1))

#relative slope areas, wherein between a slope of 2 to 10 percent the status as sloping or not is relative to its neighborhood
x.rslope <- ifel((x.slope > x.slope.med & x.slope >= 2) | x.slope >= 10, 1,0)

#slope landscape, whether most of the area is flat, intermediate, or hilly
x.slope.med.class <- ifel(
  x.slope.med < 0.5 & x.slope < 2,1, ifel(x.slope.med < 2 & x.slope < 10, 2,3))

geomorph <- x.pos.class + x.rslope*10 + x.slope.med.class*100

writeRaster(geomorph, 'geomorph3.tif', overwrite=T)
writeRaster(x.shade, 'x.shade3.tif', overwrite=T)

classnumbers = c(101,102,103,201,202,203,211,212,213,301,302,303,311,312,313)
classnames = c('dip','flat','rise',
               'valley', 'plain terrace', 'interfluve',
               'riserfoot', 'riserback', 'risershoulder',
               'toeslope', 'hill terrace', 'summit',
               'footslope', 'backslope', 'shoulder'
               )
geomorphclasses <- data.frame(class=classnumbers,name= classnames)
