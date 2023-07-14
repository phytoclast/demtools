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
  x.pos <- (dm - xmed)/(xmax - xmed)
  x.neg <- (dm - xmed)/(xmed - xmin)
  x.pos <- ifel(x.pos > 0, x.pos,0)
  x.neg <- ifel(x.neg < 0, x.neg,0)
  p <- ((x.pos+x.neg)+1)/2
  return(p)
}


dem <- rast('dem.tif')


extnt <- c(xmin= -250000, xmax=000000, ymin=200000, ymax=500000)#lmich
# extnt <- c(xmin= 100000, xmax=300000, ymin=50000, ymax=200000)#nohio
dm <- terra::crop(dem, extnt)
plot(dm)

xmax = focalmax(dm, 500)

xmin = focalmin(dm, 500)

xmax2 = focalmax(dm, 2500)

xmin2 = focalmin(dm, 2500)

x.med <- focalmed(dm, 500)

x.med2 <- focalmed(dm, 2500)

x.slope <-  tan(terrain(dm, v="slope", unit='radians'))*100

x.slope.med <- focalmed(x.slope, 500)


x.pos <- hillpos(xmax,xmin,x.med)
x.pos2 <- hillpos(xmax2,xmin2,x.med2)

rel <- xmax-xmin
rel2 <- xmax2-xmin2

xslope <- terrain(dm, v='slope', unit='radians')
xaspect <- terrain(dm, v='aspect', unit='radians')
x.shade <- shade(xslope, xaspect)

# plot(x.slope.med>1)
# 
# geomorph.flats <- ifel(
#   x.slope.med < 1 & x.slope.med < 2, 
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

##########alt
x.pos4 <- x.pos*x.pos.r + x.pos2*(x.pos.r*-1+1)
x.rslope <- ifel((x.slope > x.slope.med & x.slope >= 2) | x.slope >= 10, 1,0)

# writeRaster(x.pos4, 'x.pos4.tif', overwrite=T)
# x.pos3 <- ifel(x.pos2 >= 0.67 | x.pos >= 0.67, 3, ifel(x.pos2 < 0.33 | x.pos < 0.33,1,2))

x.slope.med.class <- ifel(
  x.slope.med < 0.5 & x.slope < 2,1, ifel(x.slope.med < 2 & x.slope < 10, 2,3))

x.pos.class <- ifel(x.pos4 >= 0.67, 3,ifel(x.pos4 >= 0.33, 2,1))

geomorph <- x.pos.class + x.rslope*10 + x.slope.med.class*100

plot(geomorph)
writeRaster(geomorph, 'geomorph2.tif', overwrite=T)
writeRaster(x.shade, 'x.shade2.tif', overwrite=T)


x.pos.r <- focalmax(x.pos2, 500) - focalmin(x.pos2, 500)

# writeRaster(x.pos.r, 'x.pos.r.tif')



writeRaster(x.pos, 'x.pos.tif', overwrite=T)
writeRaster(x.pos2, 'x.pos2.tif', overwrite=T)