
# title         : fishcamp.R
# purpose       : Pre-processing of elevation data, extraction of landform units and predictive mapping of soil mapping units
# reference     : Hengl, T. 2009. A Practical Guide to Geostatistical Mapping of Environmental Variables, 2nd Edt. EUR 22904 EN Scientific and Technical Research series, Office for Official Publications of the European Communities, Luxemburg, 264 pp. ISBN: 978-92-79-06904-8
# producer      : Prepared by T. Hengl
# last update   : In Amsterdam, NL, 8 Aug 2009.
# inputs        : study area fishcamp (37.46353 N; 119.6119 W) comprises: (1) LIDAR points, (2) contour lines, (3) SRTM DEM; These can be obtained from [http://www.spatial-analyst.net/DATA/];
# outputs       : DEMs generated from points, stream networks, DEM-parameters, predicted soil mapping units and landforms;
# remarks 1     : Guides to install and set-up SAGA are available at [http://spatial-analyst.net/wiki/index.php?title=Software]
# remarks 2     : local coordinates system used - UTM zone 11N with North American Datum 83;

library(maptools)
library(rgdal)
library(gstat)
library(spatstat)
library(mda) # kappa statistics
library(vcd)
library(RSAGA); rsaga.env()

# ------------------------------------------------------------
# Download of the maps:
# ------------------------------------------------------------

# Download the datasets:
download.file("http://www.spatial-analyst.net/DATA/fishcamp.zip", destfile=paste(getwd(), "fishcamp.zip", sep="/"))
# LiDAR points:
for(j in list(".shp", ".shx", ".dbf")){
fname <- zip.file.extract(file=paste("lidar", j, sep=""), zipname="fishcamp.zip")
file.copy(fname, paste("./lidar", j, sep=""), overwrite=TRUE)
}
# contour lines:
for(j in list(".shp", ".shx", ".dbf")){
fname <- zip.file.extract(file=paste("contours", j, sep=""), zipname="fishcamp.zip")
file.copy(fname, paste("./contours", j, sep=""), overwrite=TRUE)
}
# streams:
for(j in list(".shp", ".shx", ".dbf")){
fname <- zip.file.extract(file=paste("streams", j, sep=""), zipname="fishcamp.zip")
file.copy(fname, paste("./streams", j, sep=""), overwrite=TRUE)
}
# SRTM DEM:
fname <- zip.file.extract(file="DEMSRTM1.asc", zipname="fishcamp.zip")
file.copy(fname, "./DEMSRTM1.asc", overwrite=TRUE)
# soil map:
fname <- zip.file.extract(file="soilmu.asc", zipname="fishcamp.zip")
file.copy(fname, "./soilmu.asc", overwrite=TRUE)

# delete temp file:
unlink("fishcamp.zip")
list.files(getwd(), recursive=TRUE, full=FALSE)

# import LiDAR points:
lidar <- readOGR("lidar.shp", "lidar")
str(lidar@data) # 273,028 points (a very large point dataset);
# GRIDS:
grids25m <- readGDAL("DEMSRTM1.asc")
names(grids25m) <- "DEMSRTM1"
grids5m <- readGDAL("soilmu.asc")
names(grids5m) <- "soilmu"
grids5m$soilmu.c <- as.factor(grids5m$soilmu)
spplot(grids5m, col.regions=rainbow(8))
# estimate the pixel size:
pixelsize <- round(2*sqrt(areaSpatialGrid(grids25m)/length(lidar$Z)), 0)
pixelsize
# set-up the correct CRS:
proj4string(lidar) <- CRS("+init=epsg:26911")
proj4string(lidar)
proj4string(grids25m) <- CRS("+init=epsg:26911")
proj4string(grids5m) <- CRS("+init=epsg:26911")

# coordinates of the center:
grids25m.ll <- spTransform(grids25m, CRS("+proj=longlat +ellps=WGS84"))
grids25m.ll@bbox
clon <- mean(grids25m.ll@bbox[1,])
clat <- mean(grids25m.ll@bbox[2,])

# ------------------------------------------------------------
# DEM generation and variogram modelling:
# ------------------------------------------------------------

# variogram modelling

# sub-sample:
lidar.sample <- lidar[runif(length(lidar$Z))<0.05,]
varmap.plt <- plot(variogram(Z~1, lidar.sample, map=TRUE, cutoff=50*pixelsize, width=pixelsize), col.regions=grey(rev(seq(0,1,0.025))))
# obvious anisotropy!
Z.svar <- variogram(Z~1, lidar.sample, alpha=c(45,135)) # cutoff=50*dem.pixelsize
Z.vgm <- fit.variogram(Z.svar, vgm(psill=var(lidar.sample$Z), "Gau", sqrt(areaSpatialGrid(grids25m))/4, nugget=0, anis=c(p=135, s=0.6)))
vgm.plt <- plot(Z.svar, Z.vgm, plot.nu=F, cex=2, pch="+", col="black")
Z.vgm
# plot the two variograms next to each other:
print(varmap.plt, split=c(1,1,2,1), more=T)
print(vgm.plt, split=c(2,1,2,1), more=F)

# Generate an initial DEM:
rsaga.geoprocessor(lib="grid_gridding", module=0, param=list(GRID="DEM5LIDAR.sgrd", INPUT="lidar.shp", FIELD=0, LINE_TYPE=0, USER_CELL_SIZE=pixelsize, USER_X_EXTENT_MIN=grids5m@bbox[1,1]+pixelsize/2, USER_X_EXTENT_MAX=grids5m@bbox[1,2]-pixelsize/2, USER_Y_EXTENT_MIN=grids5m@bbox[2,1]+pixelsize/2, USER_Y_EXTENT_MAX=grids5m@bbox[2,2]-pixelsize/2))
# Filter artefacts --- first run the residual analysis and mask out spikes:
rsaga.geoprocessor(lib="geostatistics_grid", 0, param=list(INPUT="DEM5LIDAR.sgrd", MEAN="tmp.sgrd", STDDEV="tmp.sgrd", RANGE="tmp.sgrd", DEVMEAN="tmp.sgrd", PERCENTILE="tmp.sgrd", RADIUS=5, DIFF="dif_lidar.sgrd")) # all features >7 pixels remain!
# read back to R and mask out all areas:
rsaga.sgrd.to.esri(in.sgrd=c("dif_lidar.sgrd", "DEM5LIDAR.sgrd"), out.grids=c("dif_lidar.asc", "DEM5LIDAR.asc"), out.path=getwd(), prec=1)
grids5m$DEM5LIDAR <- readGDAL("DEM5LIDAR.asc")$band1
grids5m$dif <- readGDAL("dif_lidar.asc")$band1
lim.dif <- quantile(grids5m$dif, c(0.025,0.975), na.rm=TRUE)
lim.dif
grids5m$DEM5LIDARf <- ifelse(grids5m$dif<=lim.dif[[1]]|grids5m$dif>=lim.dif[[2]], NA, grids5m$DEM5LIDAR)
summary(grids5m$DEM5LIDARf)[7]/length(grids5m@data[[1]])  # 15% pixels masked out!
write.asciigrid(grids5m["DEM5LIDARf"], "DEM5LIDARf.asc", na.value=-1)
rsaga.esri.to.sgrd(in.grids="DEM5LIDARf.asc", out.sgrds="DEM5LIDARf.sgrd", in.path=getwd())

# Filter the missing values (close gaps):
rsaga.geoprocessor(lib="grid_tools", module=7, param=list(INPUT="DEM5LIDARf.sgrd", RESULT="DEM5LIDARf.sgrd")) # we write to the same file!
rsaga.sgrd.to.esri(in.sgrd="DEM5LIDARf.sgrd", out.grids="DEM5LIDARf.asc", out.path=getwd(), prec=1)
grids5m$DEM5LIDARf <- readGDAL("DEM5LIDARf.asc")$band1


# Generate DEM from contours:
rsaga.geoprocessor(lib="grid_spline", module=1, param=list(GRID="DEM25TPS.sgrd", SHAPES="contours.shp", TARGET=0, SELECT=1, MAXPOINTS=10, USER_CELL_SIZE=25, USER_X_EXTENT_MIN=grids5m@bbox[1,1]+pixelsize/2, USER_X_EXTENT_MAX=grids5m@bbox[1,2]-pixelsize/2, USER_Y_EXTENT_MIN=grids5m@bbox[2,1]+pixelsize/2, USER_Y_EXTENT_MAX=grids5m@bbox[2,2]-pixelsize/2))
# Deepen drainage route:
rsaga.geoprocessor(lib="ta_preprocessor", module=1, param=list(DEM="DEM25TPS.sgrd", DEM_PREPROC="DEM25TPSf.sgrd", METHOD=0))

# create empty grid:
rsaga.geoprocessor(lib="grid_tools", module=23, param=list(GRID="DEM25LIDAR.sgrd", M_EXTENT=0, XMIN=grids5m@bbox[1,1]+pixelsize/2, YMIN=grids5m@bbox[2,1]+pixelsize/2, NX=grids25m@grid@cells.dim[1], NY=grids25m@grid@cells.dim[2], CELLSIZE=25))
# resample to 25 m:
rsaga.geoprocessor(lib="grid_tools", module=0, param=list(INPUT="DEM5LIDARf.sgrd", GRID="DEM25LIDAR.sgrd", GRID_GRID="DEM25LIDAR.sgrd", METHOD=2, KEEP_TYPE=FALSE, SCALE_UP_METHOD=5))
# read to R:
rsaga.sgrd.to.esri(in.sgrd=c("DEM25LIDAR.sgrd", "DEM25TPS.sgrd"), out.grids=c("DEM25LIDAR.asc", "DEM25TPS.asc"), out.path=getwd(), prec=1)
grids25m$DEM25LIDAR <- readGDAL("DEM25LIDAR.asc")$band1
grids25m$DEM25TPS <- readGDAL("DEM25TPS.asc")$band1
# difference between the LiDAR DEM and topo DEM:
plot(grids25m$DEM25LIDAR, grids25m$DEM25TPS, asp=1)
sqrt(sum((grids25m$DEM25LIDAR - grids25m$DEM25TPS)^2)/length(grids25m$DEM25LIDAR))

# ------------------------------------------------------------
# Extraction of Land Surface Parameters:
# ------------------------------------------------------------

# Topographic Wetness Index:
rsaga.geoprocessor(lib="ta_hydrology", module=15, param=list(DEM="DEM5LIDARf.sgrd", C="catharea.sgrd", GN="catchslope.sgrd", CS="modcatharea.sgrd", SB="TWI.sgrd", T=10))
# valley depth:
rsaga.geoprocessor(lib="ta_morphometry", module=14, param=list(DEM="DEM5LIDARf.sgrd", HO="tmp.sgrd", HU="VDEPTH.sgrd", NH="tmp.sgrd", SH="tmp.sgrd", MS="tmp.sgrd", W=12, T=120, E=4))  # what do these parameters mean?!
# incoming solar radiation:
rsaga.geoprocessor(lib="ta_lighting", module=2, param=list(ELEVATION="DEM5LIDARf.sgrd", INSOLAT="INSOLAT.sgrd", DURATION="durat.sgrd", LATITUDE=clat, HOUR_STEP=2, TIMESPAN=2, DAY_STEP=5))  # time-consuming!
# convergence index:
rsaga.geoprocessor(lib="ta_morphometry", module=2, param=list(ELEVATION="DEM5LIDARf.sgrd", RESULT="CONVI.sgrd", RADIUS=3, METHOD=0, SLOPE=FALSE))  

# read maps back to R:
LSP.list <- c("TWI.asc", "VDEPTH.asc", "INSOLAT.asc", "CONVI.asc")
rsaga.sgrd.to.esri(in.sgrds=set.file.extension(LSP.list, ".sgrd"), out.grids=LSP.list, prec=1, out.path=getwd())
for(i in 1:length(LSP.list)){
  grids5m@data[strsplit(LSP.list[i], ".asc")[[1]]] <- readGDAL(LSP.list[i])$band1
}


# ------------------------------------------------------------
# Unsupervised extraction of landforms:
# ------------------------------------------------------------

# Principal components: 
pc.dem <- prcomp(~DEM5LIDARf+TWI+VDEPTH+INSOLAT+CONVI, scale=TRUE, grids5m@data)
biplot(pc.dem, arrow.len=0.1, xlabs=rep(".", length(pc.dem$x[,1])), main="PCA biplot")
# LSPs are relatively uncorrelated;

# fuzzy k-means:

# Determine number of clusters
demdata <- as.data.frame(pc.dem$x)
wss <- (nrow(demdata)-1)*sum(apply(demdata,2,var))
for (i in 2:20) {wss[i] <- sum(kmeans(demdata, centers=i)$withinss)}
plot(1:20, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")
# does not converge :(

kmeans.dem <- kmeans(demdata, 12)
grids5m$kmeans.dem <- kmeans.dem$cluster
grids5m$landform <- as.factor(kmeans.dem$cluster)
summary(grids5m$landform)
spplot(grids5m["landform"], col.regions=rainbow(12))
write.asciigrid(grids5m["kmeans.dem"], "landform.asc", na.value=-1)


# Export to GE:
rsaga.esri.to.sgrd(in.grids="landform.asc", out.sgrd="landform.sgrd", in.path=getwd())
rsaga.geoprocessor(lib="pj_proj4", 2, param=list(SOURCE_PROJ=paste('"', proj4string(grids5m), '"', sep=""), TARGET_PROJ="\"+proj=longlat +datum=WGS84\"", SOURCE="landform.sgrd", TARGET="landform_ll.sgrd", TARGET_TYPE=0, INTERPOLATION=0))
# export to PNG:
rsaga.geoprocessor(lib="io_grid_image", 0, param=list(GRID="landform_ll.sgrd", FILE="landform.png"))
# read back to R:
rsaga.sgrd.to.esri(in.sgrds="landform_ll.sgrd", out.grids="landform_ll.asc", prec=1, out.path=getwd())
landform.ll <- readGDAL("landform_ll.asc")
proj4string(landform.ll) <- CRS("+proj=longlat +datum=WGS84")
landform.kml <- GE_SpatialGrid(landform.ll)
kmlOverlay(landform.kml, kmlfile="landform.kml", imagefile="landform.png", name="Landform classes (12)")

# ------------------------------------------------------------
# Fitting variograms for different landform classes:
# ------------------------------------------------------------

lidar.sample.ov <- overlay(grids5m["landform"], lidar.sample)
lidar.sample.ov$Z <- lidar.sample$Z

landform.no <- length(levels(lidar.sample.ov$landform))
landform.vgm <- as.list(rep(NA, landform.no))
landform.par <- data.frame(landform=as.factor(levels(lidar.sample.ov$landform)), Nug=rep(NA, landform.no), Sill=rep(NA, landform.no), range=rep(NA, landform.no))

for(i in 1:length(levels(lidar.sample.ov$landform))) {
   tmp <- subset(lidar.sample.ov, lidar.sample.ov$landform==levels(lidar.sample.ov$landform)[i])
   landform.vgm[[i]] <- fit.variogram(variogram(Z~1, tmp, cutoff=50*pixelsize), vgm(psill=var(tmp$Z), "Gau", sqrt(areaSpatialGrid(grids25m))/4, nugget=0))
   landform.par$Nug[i] <- round(landform.vgm[[i]]$psill[1], 1)
   landform.par$Sill[i] <- round(landform.vgm[[i]]$psill[2], 1)
   landform.par$range[i] <- round(landform.vgm[[i]]$range[2], 1)
}
landform.par
   
# compare the fitted variogram parameters:


for(i in 1:length(landforms.vgm[[1]])) {
    landforms.vgm$Nug[noquote(i)] <- get(paste("Z.vgm",i,sep=""))$psill[1]
    landforms.vgm$Sill[noquote(i)] <- get(paste("Z.vgm",i,sep=""))$psill[2]
    landforms.vgm$range[noquote(i)] <- get(paste("Z.vgm",i,sep=""))$range[2]
}

landforms.vgm
# this is interesting! the landforms seems to differ largely

# ------------------------------------------------------------
# Spatial prediction of soil mapping units:
# ------------------------------------------------------------

# convert the raster map to polygon map:
rsaga.esri.to.sgrd(in.grids="soilmu.asc", out.sgrd="soilmu.sgrd", in.path=getwd())
rsaga.geoprocessor(lib="shapes_grid", module=6, param=list(GRID="soilmu.sgrd", SHAPES="soilmu.shp", CLASS_ALL=1))
# convert the polygon to line map:
rsaga.geoprocessor(lib="shapes_lines", module=0, param=list(POLYGONS="soilmu.shp", LINES="soilmu_l.shp"))
# derive the buffer map using the shape file:
rsaga.geoprocessor(lib="grid_gridding", module=0, param=list(GRID="soilmu_r.sgrd", INPUT="soilmu_l.shp", FIELD=0, LINE_TYPE=0, TARGET_TYPE=0, USER_CELL_SIZE=pixelsize, USER_X_EXTENT_MIN=grids5m@bbox[1,1]+pixelsize/2, USER_X_EXTENT_MAX=grids5m@bbox[1,2]-pixelsize/2, USER_Y_EXTENT_MIN=grids5m@bbox[2,1]+pixelsize/2, USER_Y_EXTENT_MAX=grids5m@bbox[2,2]-pixelsize/2))
# buffer distance:
rsaga.geoprocessor(lib="grid_tools", module=10, param=list(SOURCE="soilmu_r.sgrd", DISTANCE="soilmu_dist.sgrd", ALLOC="tmp.sgrd", BUFFER="tmp.sgrd", DIST=sqrt(areaSpatialGrid(grids25m))/3, IVAL=pixelsize))
# surface specific points (medial axes!):
rsaga.geoprocessor(lib="ta_morphometry", module=3, param=list(ELEVATION="soilmu_dist.sgrd", RESULT="soilmu_medial.sgrd", METHOD=1))
# read to R:
rsaga.sgrd.to.esri(in.sgrds="soilmu_medial.sgrd", out.grids="soilmu_medial.asc", prec=0, out.path=getwd())
grids5m$soilmu_medial <- readGDAL("soilmu_medial.asc")$band1
# generate the training pixels:
grids5m$weight <- abs(ifelse(grids5m$soilmu_medial>=0, 0, grids5m$soilmu_medial))
dens.weight <- as.im(as.image.SpatialGridDataFrame(grids5m["weight"]))
# image(dens.weight)
training.pix <- rpoint(length(grids5m$weight)/10, f=dens.weight)
# plot(training.pix)
training.pix <- data.frame(x=training.pix$x, y=training.pix$y , no=1:length(training.pix$x))
coordinates(training.pix) <- ~x+y
# spplot(grids5m["weight"], col.regions=grey(rev(seq(0,0.95,0.05))), sp.layout=list("sp.points", pch="+", training.pix, col="yellow"))
writeOGR(training.pix, "training_pix.shp", "training.pix", "ESRI Shapefile")

# Run MLR:
training.pix.ov <- overlay(grids5m, training.pix)
library(nnet)
# fit the model:
mlr.soilmu <- multinom(soilmu.c~DEM5LIDARf+TWI+VDEPTH+INSOLAT+CONVI, training.pix.ov)
# summary(mlr.soilmu)
grids5m$soilmu.mlr <- predict(mlr.soilmu, newdata=grids5m)
spplot(grids5m["soilmu.mlr"], col.regions=rainbow(length(levels(grids5m$soilmu.c))))
grids5m$soilmumlr <- as.integer(grids5m$soilmu.mlr)
write.asciigrid(grids5m["soilmumlr"], "soilmu_mlr.asc", na.value=-1)

# the kappa statistics (the whole map!)
sel <- !is.na(grids5m$soilmu.c)
Kappa(confusion(grids5m$soilmu.c[sel], grids5m$soilmu.mlr[sel]))
agreementplot(confusion(grids5m$soilmu.c, grids5m$soilmu.mlr))

# ------------------------------------------------------------
# Extraction of memberships:
# ------------------------------------------------------------

# mask-out classes with <5 points:
mask.c <- as.integer(attr(summary(training.pix.ov$soilmu.c[summary(training.pix.ov$soilmu.c)<5]), "names"))
mask.c
# fuzzy exponent:
fuzzy.e <- 1.2
# extract the class centroids:
class.c <- aggregate(training.pix.ov@data[c("DEM5LIDARf", "TWI", "VDEPTH", "INSOLAT", "CONVI")], by=list(training.pix.ov$soilmu.c), FUN="mean")
class.sd <- aggregate(training.pix.ov@data[c("DEM5LIDARf", "TWI", "VDEPTH", "INSOLAT", "CONVI")], by=list(training.pix.ov$soilmu.c), FUN="sd")
# derive distances in feature space:
distmaps <- as.list(levels(grids5m$soilmu.c)[mask.c])
tmp <- rep(NA, length(grids5m@data[[1]]))
for(c in (1:length(levels(grids5m$soilmu.c)))[mask.c]){
distmaps[[c]] <- data.frame(DEM5LIDARf=tmp, TWI=tmp, VDEPTH=tmp, INSOLAT=tmp, CONVI=tmp)
for(j in list("DEM5LIDARf", "TWI", "VDEPTH", "INSOLAT", "CONVI")){
distmaps[[c]][j] <- ((grids5m@data[j]-class.c[c,j])/class.sd[c,j])^2
}
}
# sum up distances per class:
distsum <- data.frame(tmp)
for(c in (1:length(levels(grids5m$soilmu.c)))[mask.c]){
distsum[paste(c)] <- sqrt(rowSums(distmaps[[c]], na.rm=T, dims=1))
}
distsum[1] <- NULL
str(distsum)
totsum <- rowSums(distsum^(-2/(fuzzy.e-1)), na.rm=T, dims=1)
# derive the fuzzy membership:
for(c in (1:length(levels(grids5m$soilmu.c)))[mask.c]){
grids5m@data[paste("mu_", c, sep="")] <- (distsum[paste(c)]^(-2/(fuzzy.e-1))/totsum)[,1]
}
spplot(grids5m[c("mu_1","mu_2","mu_3","mu_4","mu_5","mu_6")], at=seq(0,1,0.05), col.regions=grey(rev(seq(0,0.95,0.05))))
# this is different of course from the exercises 3!
write.asciigrid(grids5m["mu_4"], "mu_4.asc", na.value=-1)
write.asciigrid(grids5m["mu_5"], "mu_5.asc", na.value=-1)


# end of script!