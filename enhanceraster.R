library(terra)
library(dplyr)


#climate grid (mainly temperature as precipitation has a fuzzier relationship with elevation)
t1 <- rast('C:/a/Ecological_Sites/GIS/Climate/PRISM2010/CorrectedT/biot0112/w001001.adf')
#elevation grid matching resolution of climate grid
e1 <- rast('D:/GIS/climate/alt/alt/w001001.adf')
#elevation grid of higher resolution and cropped to an area of interest
e2 <- rast('D:/GIS/DEM/R12_DEM.tif')
e2 <- crop(e2, ext(-250000, -50000,400000, 650000))


enhanceRast <- function(t1,e1,e2){
  #create new extent to crop analysis
  expts <- data.frame(x = c(ext(e2)[1],ext(e2)[1],ext(e2)[2],ext(e2)[2]), y = c(ext(e2)[3],ext(e2)[4],ext(e2)[3],ext(e2)[4]))
  #extend a little
  expts <- expts + matrix(c(-10000,-10000,10000,10000,-10000,10000,-10000,10000),ncol = 2)
  #convert to spatVect and project to get new extent
  expts <- vect(expts, geom=c('x','y'), crs=crs(e2))
  expts <- project(expts, e1) 
  #apply new extent to crop analysis
  e1 <- crop(e1, ext(expts))
  t1 <- crop(t1, ext(expts))
  #filter out bogus elevations
  e1[e1 > 9000] <- NA; e1[e1 < -500] <- NA
  #ensure grids match extent and resolution
  e1 <- project(e1, t1)
  #ensure that grid is numeric not factor
  e1 <- e1+0
  
  names(t1) <- 't1';names(e1) <- 'e1';names(e2) <- 'e2'
  #find neighborhood means
  t5km <- aggregate(t1, fact=10, fun='mean', na.rm=T)
  e5km <- aggregate(e1, fact=10, fun='mean', na.rm=T)
  t5km <- resample(t5km, t1, method='near'); names(t5km)<- 't5km'
  e5km <- resample(e5km, t1, method='near'); names(e5km)<- 'e5km'
  #difference from neighborhood means
  tdif <- t1-t5km
  edif <- e1-e5km
  #average lapse rates per elevation in a aggregate neighborhood with higher weights where elevation differences are greatest, then smoothed and resampled to target resolution
  wts <- edif^2
  rate <- tdif/(edif+0.1)
  rate.wts <- wts*rate
  rate.sum <- aggregate(rate.wts, fact=10, fun='sum', na.rm=T)
  wts.sum <- aggregate(wts, fact=10, fun='sum', na.rm=T)
  rate.5km <- rate.sum/(wts.sum+0.001)
  rate.5km <- focal(rate.5km, na.rm=T, fun="median")
  t1.90 <- project(t1, e2) 
  e1.90 <- project(e1, e2) 
  rate.90 <- project(rate.5km, e2) 
  #apply lapse rate to high resolution DEM
  new.90 <- t1.90 + (e2 - e1.90)*rate.90
  return(new.90)}

tgs.90 <- enhanceRast(t1,e1,e2)
 
names(tgs.90) <- 'tgs.90'
tgs.90
plot(tgs.90)
