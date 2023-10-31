
# title         : examplesRSAGA.R
# purpose       : combining R and SAGA for Digital Soil Mapping applications;
# reference     : Dobos, E., Hengl, T., 2009. Soil Mapping Applications. In: Hengl, T. and Reuter, H.I. (Eds), Geomorphometry: Geomorphometry: Concepts, Software, Applications. Developments in Soil Science, vol. 33, Elsevier, 461-479 pp. [http://dx.doi.org/10.1016/S0166-2481(08)00020-2]
# producer      : Prepared by T. Hengl
# last update   : In Amsterdam, NL, May 2010.
# inputs        : baranja.txt -- 60 field observations of soil properties, elevations.txt -- 6367 point-sampled heights, SMU??.asc -- maps of soil mapping units (as indicators);
# outputs       : Simulated DEMs, corresponding land surface parameters and predicted soil variables;
# remarks 1     : data available for donwload at [http://geomorphometry.org/content/some-examples-rsaga];

library(maptools)
library(rgdal)
library(gstat)
library(RSAGA)
# rsaga.env(path="C:/Program Files/saga_vc")

***************
# 1. Error propagation and geomorphometry (both can be run via R now):

# Simulate a DEM from points, then derive a slope map 20 times and then derive average values
# First, simulate a DEM from point observations using Sequential Gaussian Simulations in gstat:

# Import the point measurements of heights to generate a DEM:

elevations <- read.delim("elevations.txt") 
coordinates(elevations)=~X+Y
# spplot(elevations)

# Import the grid definition:

gridmaps <- readGDAL("SMU1.asc")
names(gridmaps) <- "SMU1"

# Derive area in km^2:
maparea <- (gridmaps@bbox["x","max"]-gridmaps@bbox["x","min"])*(gridmaps@bbox["y","max"]-gridmaps@bbox["y","min"])/1e+06

# Fit a variogram for elevations:
elevations.or <- variogram(Z~1, elevations) 
elevations.ovgm <- fit.variogram(elevations.or, vgm(1, "Sph", 1000, 1)) 
plot(elevations.or, elevations.ovgm, plot.nu=F, pch="+")

# Make 50 realizations of a DEM using Sequential Gaussian Simulations:
DEM.sim <- krige(Z~1, elevations, gridmaps, elevations.ovgm, nmax=40, nsim=50)  # takes cca 2 mins!
# plot in a loop:
for (i in 1:length(DEM.sim@data)) {
      image(as.image.SpatialGridDataFrame(DEM.sim[i]), col=terrain.colors(16), asp=1)
}

# Write the simulated DEMs to SAGA GIS format:
for(i in 1:length(DEM.sim@data)) {
   writeGDAL(DEM.sim[i], paste("DEM", i, ".sdat", sep=""), "SAGA")
}

# Now, derive SLOPE, PLANC, PROFC and SINS maps in SAGA 50 times:
# ESRI wrapper is used to get the maps directly in ArcInfo format;

for (i in 1:length(DEM.sim@data)) {
  rsaga.slope(method="poly2zevenbergen", in.dem=paste("DEM",as.character(i),sep=""), out.slope=paste("SLOPE",as.character(i),sep=""), show.output.on.console=FALSE)
   rsaga.plan.curvature(method="poly2zevenbergen", in.dem=paste("DEM",as.character(i),sep=""), out.hcurv=paste("PLANC",as.character(i),sep=""), show.output.on.console=FALSE)
   rsaga.profile.curvature(method="poly2zevenbergen", in.dem=paste("DEM",as.character(i),sep=""), out.vcurv=paste("PROFC",as.character(i),sep=""), show.output.on.console=FALSE)
   rsaga.hillshade(method="standard", azimuth=315, declination=45, exaggeration=4, in.dem=paste("DEM",as.character(i),sep=""), out.grid=paste("SOLIN",as.character(i),sep=""), show.output.on.console=FALSE)
}
# This produces a lot of maps!!
# Create lists of the maps:

slope.list <- dir(path=getwd(), pattern=glob2rx("SLOPE*.sgrd"), full.names=TRUE)
slope.list[1:5]

# Average values per pixel:
# rsaga.get.usage("geostatistics_grid", 5)
rsaga.geoprocessor(lib="geostatistics_grid", module=5, param=list(GRIDS=paste(slope.list, collapse=";"), MEAN="SLOPEavg.sgrd", STDDEV="SLOPEstd.sgrd", STDDEVLO="SLOPElow.sgrd", STDDEVHI="SLOPEhig.sgrd"))

# Optional: generate a DEM using the Thin Plate Spline interpolation:
writeOGR(elevations, "elevations.shp", "elevations", "ESRI Shapefile")
# rsaga.get.usage("grid_spline", 1)
pixsize <- gridmaps@grid@cellsize[1]
rsaga.geoprocessor(lib="grid_spline", module=1, param=list(GRID="DEMtps.sgrd", SHAPES="elevations.shp", FIELD=1, RADIUS=sqrt(maparea)*1000/3, SELECT=1, MAXPOINTS=30, TARGET=0, USER_CELL_SIZE=pixsize, USER_X_EXTENT_MIN=gridmaps@bbox[1,1]+pixsize/2, USER_X_EXTENT_MAX=gridmaps@bbox[2,1]+pixsize/2, USER_Y_EXTENT_MIN=gridmaps@bbox[1,2]-pixsize/2, USER_Y_EXTENT_MAX=gridmaps@bbox[2,2]-pixsize/2))

# Load back results to R and visualize various correlation plots:
gridmaps$SLOPEavg <- readGDAL("SLOPEavg.sdat")$band1
gridmaps$SLOPEstd <- readGDAL("SLOPEstd.sdat")$band1
gridmaps$DEM <- readGDAL("DEMtps.sdat")$band1
# Correlation between slopes and elevations:
scatter.smooth(gridmaps$DEM, gridmaps$SLOPEavg, span=9/10, col="red", xlab="var(DEM)", ylab="AVG(SLOPE)")

# Derive all landform parameters needed:
rsaga.local.morphometry(method="poly2zevenbergen", in.dem="DEMtps", out.slope="SLOPE", out.hcurv="PLANC", out.vcurv="PROFC", show.output.on.console=FALSE)
rsaga.hillshade(method="standard", azimuth=315, declination=45, exaggeration=4, in.dem="DEMtps", out.grid="SINS", show.output.on.console=FALSE)
rsaga.wetness.index(in.dem="DEMtps", out.wetness.index="TWI", out.carea="CAREA", out.cslope="CSLOPE", show.output.on.console=FALSE)

***************

# 2. Spatial interpolation:

baranja <- read.delim("baranja.txt")
coordinates(baranja) <-~X+Y
bubble(baranja, "SOLUM")
writeOGR(baranja, "baranja.shp", "baranja", "ESRI Shapefile")

rsaga.esri.to.sgrd(in.grids=paste("SMU", 1:9, ".asc", sep=""), out.sgrds=paste("SMU", 1:9, ".sgrd", sep=""), in.path=getwd())

# read grids to R:
grid.list <- c("SLOPE", "PROFC", "PLANC", "SINS", "CAREA", "CSLOPE", "TWI", paste("SMU", 1:9, sep=""))
for(i in 1:length(grid.list)){ gridmaps@data[,grid.list[i]] <- readGDAL(set.file.extension(grid.list[i], ".sdat"), silent=TRUE)$band1 }

baranja.ov <- overlay(gridmaps, baranja)
baranja.ov@data <- cbind(baranja@data, baranja.ov@data)
# regression model:

summary(solum.fit <- lm(SOLUM~DEM+SLOPE+PLANC+PROFC+TWI+SMU1+SMU2+SMU3+SMU4+SMU5+SMU7+SMU8+SMU9, data=baranja.ov))
solum.step <- step(solum.fit)
summary(solum.step)
pred.list <- names(attr(attr(solum.step$model, "terms"), "dataClasses"))[-1]

# predict values in SAGA using only regression model;
# rsaga.get.usage("geostatistics_grid", 4)
rsaga.geoprocessor(lib="geostatistics_grid", module=4, param=list(GRIDS=paste(set.file.extension(pred.list,".sgrd"), collapse=";"), SHAPES="baranja.shp", ATTRIBUTE=0, TABLE="regout.dbf", RESIDUAL="solum_res.shp", REGRESSION="SOLUM_reg.sgrd", INTERPOL=0))

# Fit the variogram for residuals:
plot(variogram(solum.step$call$formula, baranja.ov))
SOLUM.rv <- variogram(solum.step$call$formula, baranja.ov)
SOLUM.rvgm <- fit.variogram(SOLUM.rv, vgm(0, "Exp", 500, 200))
# singluar model!
# plot(SOLUM.rv, SOLUM.rvgm, plot.nu=T)

# Regression-kriging in SAGA:
# rsaga.get.usage("geostatistics_kriging", 8)
rsaga.geoprocessor(lib="geostatistics_kriging", module=8, param=list(GRIDS=paste(set.file.extension(pred.list,".sgrd"), collapse=";"), GRID="SOLUM_rk.sgrd", SHAPES="baranja.shp", FIELD=0, BLOG=FALSE, MODEL=1, TARGET=0, NUGGET=0, SILL=200, RANGE=500, INTERPOL=0, USER_CELL_SIZE=pixsize, USER_X_EXTENT_MIN=gridmaps@bbox[1,1]+pixsize/2, USER_X_EXTENT_MAX=gridmaps@bbox[2,1]+pixsize/2, USER_Y_EXTENT_MIN=gridmaps@bbox[1,2]-pixsize/2, USER_Y_EXTENT_MAX=gridmaps@bbox[2,2]-pixsize/2))
# Reports error!

# Optional, ordinary kriging:
# rsaga.get.usage("geostatistics_kriging", 6)
rsaga.geoprocessor(lib="geostatistics_kriging", module=6, param=list(GRID="SOLUM_ok.sgrd", VARIANCE="SOLUM_okvar.sgrd", SHAPES="baranja.shp", FIELD=0, MODEL=1, NUGGET=0, SILL=200, RANGE=500, TARGET=0, NUGGET=0, SILL=200, RANGE=500, USER_CELL_SIZE=pixsize, USER_X_EXTENT_MIN=gridmaps@bbox[1,1]+pixsize/2, USER_X_EXTENT_MAX=gridmaps@bbox[2,1]+pixsize/2, USER_Y_EXTENT_MIN=gridmaps@bbox[1,2]-pixsize/2, USER_Y_EXTENT_MAX=gridmaps@bbox[2,2]-pixsize/2))
# Reports error!

# Delete temporary files:
# unlink("*.sgrd"); unlink("*.sdat"); unlink("*.mgrd")

# end of script;