
# title         : streams_error.R
# purpose       : extraction of stream networks from a DEM using error propagation technique;
# reference     : Hengl, T., Heuvelink, G.B.M., van Loon, E.E., 2010. On the uncertainty of stream networks derived from elevation data: the error propagation approach. Hydrology and Earth System Sciences, (special issue Geomorphometry 2009)
# producer      : Prepared by T. Hengl
# last update   : In Amsterdam, NL, 18 Jun 2010.
# inputs        : point map with field sampled elevations;
# outputs       : a probability map of finding a stream and associated uncertainty;
# remarks 1     : this exercises generates a large amount of maps! please make sure you have enough space on your hard disk;
# remarks 2     :  first open IE and login to the website; otherwise you will not have an access to the zip files!

# ------------------------------------------------------------
# Initial settings:
# ------------------------------------------------------------

rm(list=ls())
library(maptools)
library(gstat)
library(geoR)
library(rgdal)
library(lattice)
library(RSAGA)
setInternet2(use = TRUE)

# Project settings

# study area name:
study.name <- "zlatibor" 
# CRS:
CRS_string <- "+proj=tmerc +lat_0=0 +lon_0=21 +k=0.9999 +x_0=7500000 +y_0=0 +ellps=bessel +towgs84=574.027,170.175,401.545,4.88786,-0.66524,-13.24673,0.99999311067 +units=m"
stream.name <- "streams"

# study.name <- "BaranjaHill" 
# CRS_string <- "+proj=tmerc +lat_0=0 +lon_0=18 +k=0.9999 +x_0=6500000 +y_0=0 +ellps=bessel +units=m +towgs84=550.499,164.116,475.142,5.80967,2.07902,-11.62386,0.99999445824"
# stream.name <- "StreamsLine"

# elevations file name:
elev.name <- "elevations"
gridcell <- 30
# number of simulations:
N.sim <- 100

# ------------------------------------------------------------
# Import of data to R:
# ------------------------------------------------------------
 
download.file(paste("http://geomorphometry.org/system/files/", study.name, ".zip", sep=""), destfile=paste(getwd(), "/", study.name, ".zip", sep=""), quiet=TRUE, cacheOK=FALSE)
fname <- zip.file.extract(file=paste(getwd(), "/", elev.name, ".txt", sep=""), zipname=paste(study.name, ".zip", sep=""))
file.copy(fname, paste("./", elev.name, ".txt", sep=""), overwrite=TRUE)
for(j in list(".shp", ".shx", ".dbf")){
fname <- zip.file.extract(file=paste(getwd(), "/", stream.name, j, sep=""), zipname=paste(study.name, ".zip", sep=""))
file.copy(fname, paste("./", stream.name, j, sep=""), overwrite=TRUE)
}
unlink(paste(study.name, ".zip", sep=""))
list.files(getwd(), recursive=T, full=F)

# import the map to R:
elevations <- read.delim("elevations.txt")
true.stream <- readShapeLines(paste(stream.name, ".shp", sep=""), proj4string=CRS(CRS_string))
# rename the target variable:
coordinates(elevations) <- ~X+Y
names(elevations@data)[1] <- "Z"
proj4string(elevations) <- CRS(CRS_string)
study.size <- sqrt((elevations@bbox[1,2]-elevations@bbox[1,1])*(elevations@bbox[2,2]-elevations@bbox[2,1]))/3
# plot width/height:
p.width <- 7
p.height <- (elevations@bbox[2,2]-elevations@bbox[2,1])/(elevations@bbox[1,2]-elevations@bbox[1,1])*p.width

# ------------------------------------------------------------
# Variogram modelling:
# ------------------------------------------------------------

# sub-sample to cca. 1000 points - geoR can not deal with large datasets!
sel <- runif(length(elevations@data[[1]]))< 1000/length(elevations@data[[1]])
Z.geo <- as.geodata(elevations[sel,"Z"])
# histogram:
pdf(paste("Fig_", study.name, "_geoR_plot.pdf", sep=""), width=1.2*p.width, height=p.height)
plot(Z.geo, qt.col=grey(runif(4)))
dev.off()
# Variogram modelling (target variable):
pdf(paste("Fig_", study.name, "_vgm_plot.pdf", sep=""), width=1.7*7, height=0.9*7)
par(mfrow=c(1,2))
# anisotropy:
plot(variog4(Z.geo, max.dist=study.size, messages=FALSE), lwd=2)
title(main=study.name) 
# fit variogram using likfit:
Z.svar <- variog(Z.geo, max.dist=study.size, messages=FALSE)
Z.vgm <- variofit(Z.svar, messages=FALSE, ini=c(var(Z.geo$data), study.size), fix.nugget=T, nugget=0)
# confidence bands:
env.model <- variog.model.env(Z.geo, obj.var=Z.svar, model=Z.vgm)
plot(Z.svar, envelope=env.model); lines(Z.vgm, lwd=2);
dev.off()

# ------------------------------------------------------------
# Geostatistical simulations:
# ------------------------------------------------------------

# prepare an empty grid:
demgrid <- spsample(elevations, type="regular", cellsize=c(gridcell,gridcell))
gridded(demgrid) <- TRUE
fullgrid(demgrid) <- TRUE

# copy the values fitted in geoR:
Z.ovgm <- vgm(psill=Z.vgm$cov.pars[1], model="Mat", range=Z.vgm$cov.pars[2], nugget=Z.vgm$nugget, kappa=1.2)
# *Kappa is artificially set at 1.2! to produce smooth surface!

# conditional simulations in gstat:
vt.gt <- gstat(id=c("Z"), formula=Z~1, data=elevations, model=Z.ovgm, nmax=30)
DEM.sim <- predict.gstat(vt.gt, demgrid, nsim=N.sim, debug.level=-1)
fullgrid(DEM.sim) <- TRUE
# plot 4 realizations:
pdf(paste("Fig_", study.name, "_DEM_sims.pdf", sep=""), width=1.7*p.width, height=0.5*p.height)
spplot(DEM.sim[1:4], col.regions=grey(seq(0.2,1,0.025)))
dev.off()

# ------------------------------------------------------------
# Generation of stream network:
# ------------------------------------------------------------

# write simulated DEMs to SAGA format:
for(i in 1:N.sim){
write.asciigrid(DEM.sim[i], paste("DEM", i, ".asc", sep=""), na.value=-1)
}
# get a list of files:
dem.list <- list.files(getwd(), pattern="DEM[[:digit:]]*.asc")
rsaga.esri.to.sgrd(in.grids=dem.list, out.sgrd=set.file.extension(dem.list, ".sgrd"), in.path=getwd(), show.output.on.console=FALSE)
unlink(dem.list)
# learn about parameters needed:
# rsaga.get.usage("ta_channels", 0)

# empty list:
stream.list <- list(rep(NA, N.sim))
for (i in 1:N.sim) {
# First, filter the spurious sinks:
  rsaga.geoprocessor(lib="ta_preprocessor", module=2, param=list(DEM=paste("DEM", i, ".sgrd", sep=""), RESULT="DEMflt.sgrd", MINSLOPE=0.05), show.output.on.console=FALSE)
# Second, extract the channel network:
rsaga.geoprocessor(lib="ta_channels", module=0, param=list(ELEVATION="DEMflt.sgrd", CHNLNTWRK=paste("channels", i, ".sgrd", sep=""), CHNLROUTE="channel_route.sgrd", SHAPES="channels.shp", INIT_GRID="DEMflt.sgrd", DIV_CELLS=3, MINLEN=40), show.output.on.console=FALSE)
stream.list[[i]] <- readOGR("channels.shp", "channels", verbose=FALSE)
proj4string(stream.list[[i]]) <- elevations@proj4string
}

# plot all derived streams at top of each other:
stream.plot <- as.list(rep(NA, N.sim))
for(i in 1:N.sim){
stream.plot[[i]] <- list("sp.lines", stream.list[[i]])
}
lines.plt <- spplot(DEM.sim[1], col.regions=grey(seq(0.4,1,0.025)), scales=list(draw=T), sp.layout=stream.plot, main=study.name)
pdf(paste("Fig_", study.name, "_100streams.pdf", sep=""), width=p.width, height=p.height)
print(lines.plt)
dev.off()

# read all predicted stream networks as grids:
streamgrid.list <- list.files(getwd(), pattern="channels[[:digit:]]*.sgrd")
rsaga.sgrd.to.esri(in.sgrds=streamgrid.list, out.grids=set.file.extension(streamgrid.list, ".asc"), out.path=getwd(), prec=0, show.output.on.console=FALSE)
# read all grids to R:
streamgrid <- readGDAL(set.file.extension(streamgrid.list[[1]], ".asc"))
streamgrid@data[[1]] <- ifelse(streamgrid$band1<0, 0, 1)  
for(i in 2:length(streamgrid.list)){
  tmp <- readGDAL(set.file.extension(streamgrid.list[[i]], ".asc"), silent=TRUE)
  streamgrid@data[[i]] <- ifelse(tmp$band1<0, 0, 1)
}  
names(streamgrid) <- set.file.extension(streamgrid.list, ".asc")
proj4string(streamgrid) <- elevations@proj4string
 
# sum the streams for all simulations:
streamgrid$pr <- rowSums(streamgrid@data, na.rm=T, dims=1)/length(streamgrid@data)
streamgrid$sd <- ifelse(streamgrid$pr==0, 0, -streamgrid$pr*log(streamgrid$pr)-(1-streamgrid$pr)*log(1-streamgrid$pr))
pdf(paste("Fig_", study.name, "_streams_errormap.pdf", sep=""), width=p.width, height=p.height)
stream.plt <- spplot(streamgrid["sd"], col.regions=grey(rev((1:59)/60)), scales=list(draw=T), sp.layout=list("sp.lines", true.stream), main=paste("Stream (error map) for", study.name))
print(stream.plt)
dev.off()

write.asciigrid(streamgrid["pr"], paste(study.name, "stream_pr.asc", sep=""), na.value=-1)
write.asciigrid(streamgrid["sd"], paste(study.name, "stream_sd.asc", sep=""), na.value=-1)

# ------------------------------------------------------------
# Evaluation of the propagated uncertainty:
# ------------------------------------------------------------

# mean and sd for simulated elevations:
rsaga.geoprocessor(lib="geostatistics_grid", module=5, param=list(GRIDS=paste(set.file.extension(dem.list, ".sgrd"), collapse=";"), MEAN="DEM_avg.sgrd", STDDEV="DEM_std.sgrd"), show.output.on.console=FALSE)
# slope:
rsaga.esri.wrapper(rsaga.slope, method="poly2zevenbergen", in.dem="DEM_avg.sgrd", out.slope="SLOPE.sgrd", prec=3, clean.up=F)
# residual analysis:
rsaga.geoprocessor(lib="geostatistics_grid", 0, param=list(INPUT="DEM_avg.sgrd", MEAN="tmp.sgrd", MIN="tmp.sgrd", MAX="tmp.sgrd", STDDEV="tmp.sgrd", RANGE="tmp.sgrd", DEVMEAN="tmp.sgrd", PERCENTILE="tmp.sgrd", RADIUS=5, DIFF="DIFMEAN.sgrd"))
# read back to R:
gridmaps <- readGDAL("DEM_avg.sdat")
names(gridmaps) <- "DEM"
gridmaps$std <- readGDAL("DEM_std.sdat")$band1
gridmaps$SLOPE <- readGDAL("SLOPE.sdat")$band1
gridmaps$DIFMEAN <- readGDAL("DIFMEAN.sdat)$band1

# correlation plots:
streamgrid$DIFMEAN.c <- as.factor(round(gridmaps$DIFMEAN/5, 0)) # round by 5
levels(streamgrid$DIFMEAN.c) <- as.character(as.integer(levels(streamgrid$DIFMEAN.c))*5)
# er.c <- aggregate(streamgrid$sd[streamgrid$pr>0], by=list(streamgrid$DIFMEAN.c[streamgrid$pr>0]), FUN=sum)
pdf(paste("Fig_", study.name, "_corplots_relief.pdf", sep=""), width=5, height=7)
boxplot(streamgrid$sd[streamgrid$pr>0] ~ streamgrid$DIFMEAN.c[streamgrid$pr>0], col="grey", ylab="Stream error (H)", xlab="Difference from the mean value", main=study.name)
dev.off()

# ------------------------------------------------------------
# Optional: Export of streams / probability map to Google Earth:
# ------------------------------------------------------------

# export all streams individually:
stream.list.ll <- list(rep(NA, N.sim))
for (i in 1:N.sim) {
stream.list.ll[[i]] <- spTransform(stream.list[[i]], CRS("+proj=longlat +datum=WGS84"))
writeOGR(stream.list.ll[[i]], paste("stream_sim", i, ".kml", sep=""), paste("stream_sim", i, sep=""), "KML")
# replace color:
tmp.txt <- readLines(paste("stream_sim", i, ".kml", sep=""))
tmp.txt <- gsub(pattern=paste("ff0000ff",sep=""), replacement=paste("ff000000",sep=""), tmp.txt, fixed=T)
write(tmp.txt, paste("stream_sim", i, ".kml", sep=""))
}

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

# remove all temporary files:
save.image(paste(study.name, ".RData", sep=""))
unlink("*.hgrd")
unlink("*.sgrd")
unlink("*.sdat")
unlink("DEM**.***")
unlink("channels**.***")


# end of script!
