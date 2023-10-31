
# title         : baranjahill.R
# purpose       : extraction of stream networks from a DEM using error propagation;
# reference     : Hengl, T. 2009. A Practical Guide to Geostatistical Mapping of Environmental Variables, 2nd Edt. EUR 22904 EN Scientific and Technical Research series, Office for Official Publications of the European Communities, Luxemburg, 264 pp. ISBN: 978-92-79-06904-8
# producer      : Prepared by T. Hengl
# last update   : In Amsterdam, NL, 1 Aug 2009.
# inputs        : point map with field sampled elevations (6367 points);
# outputs       : a probability map of finding a stream given the uncertainty in the model;
# remarks 1     : this exercises generates a large amount of maps! please make sure you have enough space on your hard disk;

# ------------------------------------------------------------
# Initial settings:
# ------------------------------------------------------------

library(maptools)
library(gstat)
library(geoR)
library(rgdal)
library(lattice)
library(RSAGA)
library(spatstat)
myenv <- rsaga.env(workspace="D:/DAY4", path="C:/Program Files/saga_vc")


# MGI Balkans zone 6:
gk_6 <- "+proj=tmerc +lat_0=0 +lon_0=18 +k=0.9999 +x_0=6500000 +y_0=0 +ellps=bessel +units=m +towgs84=550.499,164.116,475.142,5.80967,2.07902,-11.62386,0.99999445824"

# ------------------------------------------------------------
# Import of data to R:
# ------------------------------------------------------------
 
download.file("http://www.spatial-analyst.net/DATA/elevations.zip", destfile=paste(getwd(), "elevations.zip", sep="/"))
for(j in list(".shp", ".shx", ".dbf")){
fname <- zip.file.extract(file=paste("elevations", j, sep=""), zipname="elevations.zip")
file.copy(fname, paste("./elevations", j, sep=""), overwrite=TRUE)
}
unlink("elevations.zip")
list.files(getwd(), recursive=T, full=F)

# import the map to R:
elevations <- readShapePoints("elevations.shp", proj4string=CRS(gk_6))
str(elevations)
names(elevations@data) <- "Z"

# ------------------------------------------------------------
# Variogram modelling:
# ------------------------------------------------------------

range(elevations$Z)
# sub-sample - geoR can not deal with large datasets!
sel <- runif(length(elevations@data[[1]]))<0.05
Z.geo <- as.geodata(elevations[sel,"Z"])
# histogram:
plot.geodata(Z.geo, qt.col=grey(runif(4)))
# Variogram modelling (target variable):
par(mfrow=c(1,2))
# anisotropy:
plot(variog4(Z.geo, max.dist=1000, messages=FALSE), lwd=2)
# fit variogram using likfit:
Z.svar <- variog(Z.geo, max.dist=1000, messages=FALSE)
Z.vgm <- variofit(Z.svar, messages=FALSE, ini=c(var(Z.geo$data), 1000), fix.nugget=T, nugget=0)
Z.vgm
# confidence bands:
env.model <- variog.model.env(Z.geo, obj.var=Z.svar, model=Z.vgm)
plot(Z.svar, envelope=env.model); lines(Z.vgm, lwd=2);
dev.off()

# spacing between contours:
bin.Z <- (max(elevations$Z)-min(elevations$Z))*(length(elevations$Z))^(-1/3)
round(bin.Z, 0)

# ------------------------------------------------------------
# Geostatistical simulations:
# ------------------------------------------------------------

# prepare an empty grid:
demgrid <- spsample(elevations, type="regular", cellsize=c(30,30))
gridded(demgrid) <- TRUE
fullgrid(demgrid) <- TRUE
gridcell <- demgrid@grid@cellsize[1]

# locs <- pred_grid(c(demgrid@bbox[1,1]+gridcell/2, demgrid@bbox[1,2]-gridcell/2), c(demgrid@bbox[2,1]+gridcell/2, demgrid@bbox[2,2]-gridcell/2), by=gridcell)
# conditional simulations:
# Z.sims <- krige.conv(Z.geo, locations=locs, krige=krige.control(obj.m=Z.vgm), output=output.control(n.predictive=1)) 
# this is computationally very intensive; geoR is simply not fit to work with large data!

# copy the values fitted in geoR:
Z.ovgm <- vgm(psill=Z.vgm$cov.pars[1], model="Mat", range=Z.vgm$cov.pars[2], nugget=Z.vgm$nugget, kappa=1.2)
# *Kappa is artificially set at 1.2! to produce smooth surface!
Z.ovgm

# conditional simulations in gstat:
N.sim <- 50
DEM.sim <- krige(Z~1, elevations, demgrid, Z.ovgm, nmax=30, nsim=N.sim)
gridded(DEM.sim) <- TRUE
fullgrid(DEM.sim) <- TRUE
# plot 4 realizations:
spplot(DEM.sim[1:4], col.regions=grey(seq(0,1,0.025)))

# Cross-section at y=5,073,012:
cross.s <- data.frame(X=seq(demgrid@bbox[1,1]+gridcell/2, demgrid@bbox[1,2]-gridcell/2, gridcell), Y=rep(5073012, demgrid@grid@cells.dim[1]))
coordinates(cross.s) <-~X+Y
# proj4string(cross.s) <- elevations@proj4string
cross.ov <- overlay(DEM.sim, cross.s)
plot(cross.ov@coords[,1], cross.ov@data[[1]], type="l", xlab="X", ylab="Z", col="grey")
for(i in 2:N.sim-1){
lines(cross.ov@coords[,1], cross.ov@data[[i]], col="grey")
}
lines(cross.ov@coords[,1], cross.ov@data[[N.sim]], lwd=2)

# ------------------------------------------------------------
# Generation of stream network:
# ------------------------------------------------------------

# write simulated DEMs to SAGA format:
for(i in 1:N.sim){
write.asciigrid(DEM.sim[i], paste("DEM", i, ".asc", sep=""), na.value=-1)
}
# get a list of files:
dem.list <- list.files(getwd(), pattern="DEM[[:digit:]]*.asc")
rsaga.esri.to.sgrd(in.grids=dem.list, out.sgrd=set.file.extension(dem.list, ".sgrd"), in.path=getwd(), show.output.on.console=FALSE, env=myenv)
unlink(dem.list)
# learn about parameters needed:
# rsaga.get.usage("ta_channels", 0)

# empty list:
stream.list <- list(rep(NA, N.sim))
for (i in 1:N.sim) {
# First, filter the spurious sinks:
  rsaga.geoprocessor(lib="ta_preprocessor", module=2, param=list(DEM=paste("DEM", i, ".sgrd", sep=""), RESULT="DEMflt.sgrd", MINSLOPE=0.05), show.output.on.console=FALSE, env=myenv)
# Second, extract the channel network:
rsaga.geoprocessor(lib="ta_channels", module=0, param=list(ELEVATION="DEMflt.sgrd", CHNLNTWRK=paste("channels", i, ".sgrd", sep=""), CHNLROUTE="channel_route.sgrd", SHAPES="channels.shp", INIT_GRID="DEMflt.sgrd", DIV_CELLS=3, MINLEN=40), show.output.on.console=FALSE, env=myenv)
stream.list[[i]] <- readOGR("channels.shp", "channels")
proj4string(stream.list[[i]]) <- elevations@proj4string
}

# plot all derived streams at top of each other:
stream.plot <- as.list(rep(NA, N.sim))
for(i in 1:N.sim){
stream.plot[[i]] <- list("sp.lines", stream.list[[i]])
}
lines.plt <- spplot(DEM.sim[1], col.regions=grey(seq(0.5,1,0.025)), scales=list(draw=T), sp.layout=stream.plot, main="100 streams")

# read all predicted stream networks as grids:
streamgrid.list <- list.files(getwd(), pattern="channels[[:digit:]]*.sgrd")
rsaga.sgrd.to.esri(in.sgrds=streamgrid.list, out.grids=set.file.extension(streamgrid.list, ".asc"), out.path=getwd(), prec=0, show.output.on.console=FALSE, env=myenv)
# read all grids to R:
streamgrid <- readGDAL(set.file.extension(streamgrid.list[[1]], ".asc"))
streamgrid@data[[1]] <- ifelse(streamgrid$band1<0, 0, 1)  
for(i in 2:length(streamgrid.list)){
  tmp <- readGDAL(set.file.extension(streamgrid.list[[i]], ".asc"))
  streamgrid@data[[i]] <- ifelse(tmp$band1<0, 0, 1)
}  
names(streamgrid) <- set.file.extension(streamgrid.list, ".asc")
proj4string(streamgrid) <- elevations@proj4string

# sum the streams for all simulations:
streamgrid$pr <- rowSums(streamgrid@data, na.rm=T, dims=1)/length(streamgrid@data)
stream.plt <- spplot(streamgrid["pr"], col.regions=grey(rev((1:59)/60)), scales=list(draw=T), sp.layout=list("sp.lines", stream.list[[1]]), main="Stream (probability)")
print(lines.plt, split=c(1,1,2,1), more=TRUE)
print(stream.plt, split=c(2,1,2,1), more=FALSE)

write.asciigrid(streamgrid["pr"], "stream_pr.asc", na.value=-1)

# ------------------------------------------------------------
# Evaluation of the propagated uncertainty:
# ------------------------------------------------------------

# mean and sd for simulated elevations:
rsaga.geoprocessor(lib="geostatistics_grid", module=5, param=list(GRIDS=paste(set.file.extension(dem.list, ".sgrd"), collapse=";"), MEAN="DEM_avg.sgrd", STDDEV="DEM_std.sgrd"), show.output.on.console=FALSE, env=myenv)
# slope:
rsaga.esri.wrapper(rsaga.slope, method="poly2zevenbergen", in.dem="DEM_avg.sgrd", out.slope="SLOPE.sgrd", prec=3, clean.up=F, env=myenv)
# residual analysis:
rsaga.geoprocessor(lib="geostatistics_grid", 0, param=list(INPUT="DEM_avg.sgrd", MEAN="tmp.sgrd", STDDEV="tmp.sgrd", RANGE="tmp.sgrd", DEVMEAN="tmp.sgrd", PERCENTILE="tmp.sgrd", RADIUS=5, DIFF="DIFMEAN.sgrd"), env=myenv)

rsaga.sgrd.to.esri(in.sgrds=c("DEM_avg.sgrd","DEM_std.sgrd", "DIFMEAN.sgrd"), out.grids=c("DEM_avg.asc","DEM_std.asc", "DIFMEAN.asc"), out.path=getwd(), prec=2, env=myenv)
# read back to R:
gridmaps <- readGDAL("DEM_avg.asc")
names(gridmaps) <- "DEM"
gridmaps$std <- readGDAL("DEM_std.asc")$band1
gridmaps$SLOPE <- readGDAL("SLOPE.asc")$band1
gridmaps$DIFMEAN <- readGDAL("DIFMEAN.asc")$band1

# correlation plots:
par(mfrow=c(1,2))
scatter.smooth(gridmaps$SLOPE, gridmaps$std, span=18/19, col="grey", xlab="Slope", ylab="DEM error", pch=19)
scatter.smooth(streamgrid$pr[streamgrid$pr>0], gridmaps$DIFMEAN[streamgrid$pr>0], span=18/19, col="grey", xlab="Stream probability", ylab="Difference from mean value", pch=19)

# ------------------------------------------------------------
# Optional: selection of the cell size
# ------------------------------------------------------------

pixel.range <- c(20, 30, 40, 50, 60, 80, 100)
# generate DEMs using splines:
for(i in 1:length(pixel.range)){
rsaga.geoprocessor(lib="grid_spline", module=1, param=list(GRID=paste("DEMpix", pixel.range[i], ".sgrd", sep=""), SHAPES="elevations.shp", FIELD=1, RADIUS=sqrt(areaSpatialGrid(demgrid))/3, SELECT=1, MAXPOINTS=10, TARGET=0, 
USER_CELL_SIZE=pixel.range[i], USER_X_EXTENT_MIN=demgrid@bbox[1,1]+pixel.range[i]/2, USER_X_EXTENT_MAX=demgrid@bbox[1,2]-pixel.range[i]/2, USER_Y_EXTENT_MIN=demgrid@bbox[2,1]+pixel.range[i]/2, USER_Y_EXTENT_MAX=demgrid@bbox[2,2]-pixel.range[i]/2), env=myenv)
}

# estimate drainage map for each DEM:
for(i in 1:length(pixel.range)){
# filter the spurious sinks:
rsaga.geoprocessor(lib="ta_preprocessor", module=2, param=list(DEM=paste("DEMpix", pixel.range[i], ".sgrd", sep=""), RESULT="DEMflt.sgrd", MINSLOPE=0.05), show.output.on.console=FALSE, env=myenv)
# minimum length:
min.len <- round(sqrt(areaSpatialGrid(demgrid))/(pixel.range[i]*3.5), 0)
# extract the channel network:
rsaga.geoprocessor(lib="ta_channels", module=0, param=list(ELEVATION="DEMflt.sgrd", CHNLNTWRK=paste("chnlntwrk_pix", pixel.range[i], ".sgrd", sep=""), CHNLROUTE="tmp.sgrd", SHAPES=paste("channels_pix", i, ".shp", sep=""), INIT_GRID="DEMflt.sgrd", DIV_CELLS=3, MINLEN=min.len), show.output.on.console=FALSE, env=myenv)
# buffer distance to actual streams (use the finest grid cell size):
rsaga.geoprocessor(lib="grid_gridding", module=0, param=list(GRID="stream_pix.sgrd", INPUT=paste("channels_pix", i, ".shp", sep=""), FIELD=1, LINE_TYPE=0, USER_CELL_SIZE=pixel.range[1], USER_X_EXTENT_MIN=demgrid@bbox[1,1]+pixel.range[1]/2, USER_X_EXTENT_MAX=demgrid@bbox[1,2]-pixel.range[1]/2, USER_Y_EXTENT_MIN=demgrid@bbox[2,1]+pixel.range[1]/2, USER_Y_EXTENT_MAX=demgrid@bbox[2,2]-pixel.range[1]/2), show.output.on.console=FALSE, env=myenv)
# extract a buffer distance map:
rsaga.geoprocessor(lib="grid_tools", module=10, param=list(SOURCE="stream_pix.sgrd", DISTANCE="tmp.sgrd", ALLOC="tmp.sgrd", BUFFER=paste("buffer_pix", pixel.range[i], ".sgrd", sep=""), DIST=2000, IVAL=pixel.range[1]), show.output.on.console=FALSE, env=myenv)
}

# read the maps to R:
rsaga.sgrd.to.esri(in.sgrds="chnlntwrk_pix20.sgrd", out.grids="chnlntwrk_pix20.asc", out.path=getwd(), prec=1, env=myenv)
griddrain <- readGDAL("chnlntwrk_pix20.asc")
names(griddrain) <- "chnlntwrk"
for(i in 2:length(pixel.range)){
rsaga.sgrd.to.esri(in.sgrds=paste("buffer_pix", pixel.range[i], ".sgrd", sep=""), out.grids=paste("buffer_pix", pixel.range[i], ".asc", sep=""), out.path=getwd(), prec=1, env=myenv)
griddrain@data[paste("buffer_pix", pixel.range[i], sep="")] <- readGDAL(paste("buffer_pix", pixel.range[i], ".asc", sep=""))$band1
}
str(griddrain@data)

# summary statistics:
stream.dist <- as.list(rep(NA, length(pixel.range)))
mean.dist <- c(0, rep(NA, length(pixel.range)-1))
for(i in 2:length(pixel.range)){
stream.dist[[i]] <- summary(griddrain@data[!is.na(griddrain$chnlntwrk),i])
mean.dist[i] <- stream.dist[[i]][4]
}
stream.dist
# final plot:
plot(pixel.range, mean.dist, xlab="Grid cell size", ylab="Error of mapping a stream", pch=19) 
lines(pixel.range, mean.dist) 

# ------------------------------------------------------------
# Optional: Export of grids to Google Earth from R:
# ------------------------------------------------------------

# reproject the grid:
streamgrid.ll <- spTransform(streamgrid["pr"], CRS("+proj=longlat +datum=WGS84"))
streamgrid.ll@bbox
# grid cell size:
corrf <- (1 + cos((streamgrid.ll@bbox[1,"max"] + streamgrid.ll@bbox[2,"min"])/2*pi/180))/2
geogrd.cell <- corrf*(streamgrid.ll@bbox[1,"max"] - streamgrid.ll@bbox[1,"min"]) / streamgrid@grid@cells.dim[1]
geogrd.cell

# new grid:
geoarc <- spsample(streamgrid.ll, type="regular", cellsize=c(geogrd.cell,geogrd.cell))
gridded(geoarc) <- TRUE
gridparameters(geoarc)
geoarc.grid <- SpatialGridDataFrame(geoarc@grid, data=data.frame(rep(1, length(geoarc@grid.index))), proj4string=streamgrid.ll@proj4string)
# resample values:
library(akima)
streamgrid.llgrd <- interp(x=streamgrid.ll@coords[,1], y=streamgrid.ll@coords[,2], z=streamgrid.ll$pr, xo=seq(geoarc.grid@bbox[1,"min"],geoarc.grid@bbox[1,"max"], length=geoarc.grid@grid@cells.dim[[1]]), yo=seq(geoarc.grid@bbox[2,"min"],geoarc.grid@bbox[2,"max"], length=geoarc.grid@grid@cells.dim[[2]]), linear=TRUE, extrap=FALSE)
# convert to sp class:
streamgrid.llgrd <- as(as.im(streamgrid.llgrd), "SpatialGridDataFrame")
proj4string(streamgrid.llgrd) <- CRS("+proj=longlat +datum=WGS84")
# mask the "0" values:
streamgrid.llgrd$pr <- ifelse(streamgrid.llgrd$v<0.05, NA, streamgrid.llgrd$v)

# generate a KML ground overlay:
streamgrid.kml <- GE_SpatialGrid(streamgrid.llgrd)
png(file="stream.png", width=streamgrid.kml$width, height=streamgrid.kml$height, bg="transparent")
par(mar=c(0,0,0,0), xaxs="i", yaxs="i")
image(as.image.SpatialGridDataFrame(streamgrid.llgrd["pr"]), col=grey(rev(seq(0,0.95,0.05))), xlim=streamgrid.kml$xlim, ylim=streamgrid.kml$ylim)
kmlOverlay(streamgrid.kml, "stream.kml", "stream.png", name="Stream probability")
dev.off()

# ------------------------------------------------------------
# Optional: extracion of streams in GRASS GIS
# ------------------------------------------------------------

library(spgrass6) # version => 0.6-1 (http://spatial.nhh.no/R/Devel/spgrass6_0.6-1.zip)
# Location of your GRASS installation:
loc <- initGRASS("C:/GRASS", home=tempdir())
loc
# Import the ArcInfo ASCII file to GRASS:
parseGRASS("r.in.gdal") # commmand description
execGRASS("r.in.gdal", flags="o", parameters=list(input="DEM_avg.asc", output="DEM"))
# re-write the environmental parameters:
execGRASS("g.region", parameters=list(rast="DEM"))
gmeta6()
# extract the drainage network:
execGRASS("r.watershed", flags=c("m", "overwrite"), parameters=list(elevation="DEM", stream="stream", threshold=as.integer(50)))
# thin the raster map so it can be converted to vectors:
execGRASS("r.thin", parameters=list(input="stream", output="streamt"))
# convert to vectors:
execGRASS("r.to.vect", parameters=list(input="streamt", output="streamt", feature="line"))
# read to R:
streamt <- readVECT6("streamt")
plot(streamt)

# remove all temporary files:
save.image(".RData")
unlink("*.hgrd")
unlink("*.sgrd")
unlink("*.sdat")
unlink("DEM**.***")
unlink("channels**.***")


# end of script!
