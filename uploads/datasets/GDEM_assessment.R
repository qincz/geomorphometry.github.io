
# title         : GDEM_assessment.R
# purpose       : Quick assessment of GDEM using LiDAR data from four regions;
# reference     : Hengl, T., Reuter, H.I., 2011. How accurate and usable is GDEM? A statistical assessment of GDEM using LiDAR data. Geomorphometry 2011 Conference proceedings, Sept 7-11, Redlands University.
# producer      : Prepared by T. Hengl
# last update   : In Wageningen, NL, Sept 2011.
# inputs        : GDEM and LiDAR DEMs for four case studies (data available at www.geomorphometry.org); 
# outputs       : regression models, tests of significance, visualizations; SAGA GIS version 2.0.7!

# ------------------------------------------------------------
# Initial settings and data loading
# ------------------------------------------------------------

library(maptools)
library(gstat)
library(rgdal)
library(RODBC)
library(RSAGA)
library(lattice)
# Get the paths for FWTools [http://fwtools.maptools.org/]
fw.path <- utils::readRegistry("SOFTWARE\\WOW6432NODE\\FWTools")$Install_Dir
fw.path
# On 32bit machine use:
# fw.path <- utils::readRegistry("SOFTWARE\\FWTools")$Install_Dir
gdalwarp <- shQuote(normalizePath(file.path(fw.path, "bin/gdalwarp.exe")))
gdal_translate <- shQuote(normalizePath(file.path(fw.path, "bin/gdal_translate.exe")))

# Download the datasets:
download.file("http://geomorphometry.org/sites/default/files/GDEM_assessment.zip", destfile=paste(getwd(), "GDEM_assessment.zip", sep="/"))
unzip("GDEM_assessment.zip")
# list all files:
LDEM.list <- dir(pattern="*LDEM.asc")
GDEM.list <- dir(pattern="*GDEM.asc")
study.list <- list(rep(NA, length(LDEM.list)))
# read maps to R:
grid.list <- list(rep(NA, length(LDEM.list)))

for(j in 1:length(LDEM.list)){
    grid.list[[j]] <- readGDAL(LDEM.list[[j]])
    study.list[[j]] <- strsplit(LDEM.list[j], "_")[[1]][1]
    names(grid.list[[j]]@data) <- paste(study.list[[j]], "_LDEM", sep="") 
    grid.list[[j]]@data[, paste(study.list[[j]], "_GDEM", sep="")] <- readGDAL(GDEM.list[[j]])$band1
}
 
# check structure:
str(grid.list[[j]]@data)
# print stdev:
for(j in 1:length(LDEM.list)){ print(sd(grid.list[[j]]@data[,1]), na.rm=TRUE) }
# size of the study areas:
for(j in 1:length(LDEM.list)){ print(areaSpatialGrid(grid.list[[j]])/1e6) }

# ------------------------------------------------------------
# Regression modelling
# ------------------------------------------------------------

# correlation plots [http://addictedtor.free.fr/graphiques/RGraphGallery.php?graph=159]
library(IDPmisc)
ipairs(grid.list[[2]]@data[,1:2], cex=0.8, ztransf = function(x){x[x<1] <- 1; log2(x)})
ipairs(grid.list[[1]]@data[,1:2], cex=0.8, ztransf = function(x){x[x<1] <- 1; log2(x)})
# simple plot:
par(mfrow=c(1,2), mar=c(2,2,1,1))
scatter.smooth(grid.list[[2]]@data[,2], grid.list[[2]]@data[,1], xlab="", ylab="", main="", cex=.5, col="grey", pch=19, asp=1)
text(locator(1), "Calabria")
scatter.smooth(grid.list[[1]]@data[,2], grid.list[[1]]@data[,1], xlab="", ylab="", main="", cex=.5, col="grey", pch=19, asp=1)
text(locator(1), "Boschoord")

# amount of variation explained by GDEM:
lm.list <- list(rep(NA, length(LDEM.list)))
formula.list <- list(rep(NA, length(LDEM.list)))
for(j in 1:length(LDEM.list)){
 formula.list[[j]] <- as.formula(paste(names(grid.list[[j]]@data)[1], "~", names(grid.list[[j]]@data)[2]))
 lm.list[[j]] <- lm(formula.list[[j]], grid.list[[j]]@data)
}
# print all R-squares:
for(j in 1:length(LDEM.list)){ print(summary(lm.list[[j]])$adj.r.squared) }

# RMSE:
for(j in 1:length(LDEM.list)){ print(sqrt(sum((grid.list[[j]]@data[,1] - grid.list[[j]]@data[,2])^2)/length(grid.list[[j]]@data[,1]))) }

# ------------------------------------------------------------
# Comparison of stream networks
# ------------------------------------------------------------

stream_LDEM.list <- list(rep(NA, length(LDEM.list)))
stream_GDEM.list <- list(rep(NA, length(LDEM.list)))
for (j in 1:length(LDEM.list)) {
# write to SAGA grid:
writeGDAL(grid.list[[j]][1], paste(names(grid.list[[j]]@data)[1], ".sdat", sep=""), "SAGA")
writeGDAL(grid.list[[j]][2], paste(names(grid.list[[j]]@data)[2], ".sdat", sep=""), "SAGA")
# filter the spurious sinks:
rsaga.geoprocessor(lib="ta_preprocessor", module=2, param=list(DEM=paste(names(grid.list[[j]]@data)[1], ".sgrd", sep=""), DEM_PREPROC=paste(names(grid.list[[j]]@data)[1], "f.sgrd", sep=""), METHOD=1, THRESHOLD=FALSE), show.output.on.console=FALSE)
rsaga.geoprocessor(lib="ta_preprocessor", module=2, param=list(DEM=paste(names(grid.list[[j]]@data)[2], ".sgrd", sep=""), DEM_PREPROC=paste(names(grid.list[[j]]@data)[2], "f.sgrd", sep=""), METHOD=1, THRESHOLD=FALSE), show.output.on.console=FALSE)
# extract the channel network:
rsaga.geoprocessor(lib="ta_channels", module=0, param=list(ELEVATION=paste(names(grid.list[[j]]@data)[1], "f.sgrd", sep=""), CHNLNTWRK="channels.sgrd", CHNLROUTE="channel_route.sgrd", SHAPES=paste("channels_", names(grid.list[[j]]@data)[1], ".shp", sep=""), INIT_GRID=paste(names(grid.list[[j]]@data)[1], "f.sgrd", sep=""), DIV_CELLS=3, MINLEN=40), show.output.on.console=FALSE)
# buffer distance from streams:
rsaga.geoprocessor(lib="grid_tools", module=10, param=list(SOURCE="channels.sgrd", DISTANCE=paste("dist_", names(grid.list[[j]]@data)[1], ".sgrd", sep=""), ALLOC="tmp.sgrd", BUFFER="tmp.sgrd", DIST=sqrt(areaSpatialGrid(grid.list[[j]]))/3, IVAL=100))
stream_LDEM.list[[j]] <- readOGR(paste("channels_", names(grid.list[[j]]@data)[1], ".shp", sep=""), paste("channels_", names(grid.list[[j]]@data)[1], sep=""), verbose=FALSE)
grid.list[[j]]@data[, paste("dist_", study.list[[j]], "_LDEM", sep="")] <- readGDAL(paste("dist_", names(grid.list[[j]]@data)[1], ".sdat", sep=""))$band1
rsaga.geoprocessor(lib="ta_channels", module=0, param=list(ELEVATION=paste(names(grid.list[[j]]@data)[2], "f.sgrd", sep=""), CHNLNTWRK="channels.sgrd", CHNLROUTE="channel_route.sgrd", SHAPES=paste("channels_", names(grid.list[[j]]@data)[2], ".shp", sep=""), INIT_GRID=paste(names(grid.list[[j]]@data)[2], "f.sgrd", sep=""), DIV_CELLS=3, MINLEN=40), show.output.on.console=FALSE)
# buffer distance from streams:
rsaga.geoprocessor(lib="grid_tools", module=10, param=list(SOURCE="channels.sgrd", DISTANCE=paste("dist_", names(grid.list[[j]]@data)[2], ".sgrd", sep=""), ALLOC="tmp.sgrd", BUFFER="tmp.sgrd", DIST=sqrt(areaSpatialGrid(grid.list[[j]]))/3, IVAL=100))
stream_GDEM.list[[j]] <- readOGR(paste("channels_", names(grid.list[[j]]@data)[2], ".shp", sep=""), paste("channels_", names(grid.list[[j]]@data)[2], sep=""), verbose=FALSE)
grid.list[[j]]@data[, paste("dist_", study.list[[j]], "_GDEM", sep="")] <- readGDAL(paste("dist_", names(grid.list[[j]]@data)[2], ".sdat", sep=""))$band1
}

# amount of variation explained by the regression model (GDEM):
dlm.list <- list(rep(NA, length(LDEM.list)))
dformula.list <- list(rep(NA, length(LDEM.list)))
for(j in 1:length(LDEM.list)){
 dformula.list[[j]] <- as.formula(paste(names(grid.list[[j]]@data)[3], "~", names(grid.list[[j]]@data)[4]))
 dlm.list[[j]] <- lm(dformula.list[[j]], grid.list[[j]]@data)
}
# print all R-squares:
for(j in 1:length(LDEM.list)){ print(summary(dlm.list[[j]])$adj.r.squared) }

# RMSE for stream networks:
for(j in 1:length(LDEM.list)){ print(sqrt(sum((grid.list[[j]]@data[,1] - grid.list[[j]]@data[,2])^2)/length(grid.list[[j]]@data[,1]))) }

# ------------------------------------------------------------
# Surface roughness
# ------------------------------------------------------------

# generate 1000 random samples and fit a variogram:
vgm_LDEM.list <- list(rep(NA, length(LDEM.list)))
var_LDEM.list <- list(rep(NA, length(LDEM.list)))
vgm_GDEM.list <- list(rep(NA, length(LDEM.list)))
var_GDEM.list <- list(rep(NA, length(LDEM.list)))
for(j in 1:length(LDEM.list)){
rnd.pnt <- spsample(grid.list[[j]][1], 1000, type="random")
zL.pnt <- overlay(grid.list[[j]][1], rnd.pnt)
zG.pnt <- overlay(grid.list[[j]][2], rnd.pnt)
vgmL.formula <- as.formula(paste(names(grid.list[[j]]@data)[1], "~ 1"))
vgmG.formula <- as.formula(paste(names(grid.list[[j]]@data)[2], "~ 1"))
var_LDEM.list[[j]] <- variogram(vgmL.formula, zL.pnt)
var_GDEM.list[[j]] <- variogram(vgmG.formula, zG.pnt)
vgm_LDEM.list[[j]] <- fit.variogram(var_LDEM.list[[j]], vgm(nugget=0, model="Mat", range=sqrt(areaSpatialGrid(grid.list[[j]]))/3, psill=var(z.pnt@data[,1]), kappa=1.2))
vgm_GDEM.list[[j]] <- fit.variogram(var_GDEM.list[[j]], vgm(nugget=0, model="Mat", range=sqrt(areaSpatialGrid(grid.list[[j]]))/3, psill=var(z.pnt@data[,1]), kappa=1.2))
}

# plot a variogram:
plot(var_LDEM.list[[2]], vgm_LDEM.list[[2]]) 
plot(var_GDEM.list[[2]], vgm_GDEM.list[[2]])

# nugget variation for LDEM and GDEM:
vgm_LDEM.list[[2]]
vgm_GDEM.list[[2]]

# ------------------------------------------------------------
# What do people think?
# ------------------------------------------------------------

# download data from the website:
geomorph.conn <- odbcConnect(dsn="geomorphometry.org", DBMSencoding="UTF-8")
# get results of poll:
poll_text.tbl <- sqlQuery(geomorph.conn, query=paste("SELECT nid, chtext, chvotes FROM poll_choices", sep=""))
str(poll_text.tbl)
# plot a pie chart:
pie.data <- poll_text.tbl[poll_text.tbl$nid==1422,"chvotes"]
names(pie.data) <- c("SRTM still better", "Inaccurate", "Difficult to obtain", "Cell size needs correction", "Better than SRTM")
# poll_text.tbl[poll_text.tbl$nid==1422,"chtext"] 
pie(pie.data, col=c(grey(.8),grey(.2),grey(.5),grey(.95),grey(.7),"white"))
title(main=paste("GDEM usability on ", format(Sys.time(), "%B %d %Y"), "  (N=", sum(pie.data),")", sep=""), cex.main=1.2, font.main=1) 

# end of script;