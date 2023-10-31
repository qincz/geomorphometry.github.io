
# title         : Boschoord.R
# purpose       : Automated supervised and usupervised extraction of geomorphological features using digital elevation data: a case study in Drenthe;
# reference     : Seijmonsbergen, A.C., Hengl, T., Anders, N.S., 2010. Automated extraction of geomorphological features using DEMs. In: Smith, M.J., Paron, P. and Griffiths, J. eds. "Geomorphological Mapping: a professional handbook of techniques and applications." Developments in Earth Surface Processes, Elsevier.
# producer      : Prepared by T. Hengl
# last update   : In Amsterdam, NL, 15 Jun 2009.
# inputs        : ahn5m.img - LiDAR-based DEM; BOD50.img - soil mapping units 1:50k; geomorph.img - Geomorphological map; TOP10vlak.shp - land use areas; TOP10vlak.shp - location buildings; hoogte_16ef.shp - field measurements of elevation;
# outputs       : extracted classes and summary statistics; various plots and images;


##############################################################
#   TABLE OF CONTENT:
#   - Initial settings
#    1. Import and preparation of maps:
#    2. Derivation of DEM-parameters
#    3. Multinomial Logistic Regression
#    4. Unsupervised classification
#
##############################################################

download.file("http://www.geomorphometry.org/sites/default/files/Boschoord.zip", destfile=paste(getwd(), "Boschoord.zip", sep="/"), mode='wb')
unzip("Boschoord.zip")

##############################################################
#   Initial settings:
##############################################################

library(maptools)
library(rgdal)
library(RSAGA)
library(gstat)
library(spatstat)
library(mda) # kappa statistics
library(vcd)
ILWIS <- "C:\\Progra~1\\N52\\Ilwis35\\IlwisClient.exe -C"
NL_RD <- "+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.999908 +x_0=155000 +y_0=463000 +ellps=bessel +towgs84=565.237,50.0087,465.658,-0.406857,0.350733,-1.87035,4.0812 +units=m +no_defs"

##############################################################
#    1. Import and preparation of maps:
#
##############################################################

##############################################################
#    1.1 Import of maps
##############################################################

grids5m <- readGDAL("ahn5m.img", p4s=NL_RD)
names(grids5m) <- "ahn5m"
# ~1 milion pixels! we will only work with 25 maps
str(grids5m)
# import to SAGA GIS format:
rsaga.geoprocessor(lib="io_gdal", 0, param=list(GRIDS="ahn5m.sgrd", FILE="ahn5m.img"))

grids25m.list <- c("geomorph.img", "BOD50.img")
grids25m <- readGDAL(grids25m.list[1], p4s=NL_RD)
names(grids25m) <- strsplit(grids25m.list[1], ".img")[[1]]
for(i in 2:length(grids25m.list)){
  grids25m@data[strsplit(grids25m.list[i], ".img")[[1]]] <- readGDAL(grids25m.list[i], p4s=NL_RD)$band1
}
# mask "Bebouwing" and "Water", merge some classes:
grids25m$geomorph.c <- as.factor(ifelse(grids25m$geomorph==15|grids25m$geomorph==9, NA, ifelse(grids25m$geomorph==6, 4,ifelse(grids25m$geomorph==8, 0, grids25m$geomorph)))) 
summary(grids25m$geomorph.c)
grids25m$geomorph <- as.integer(grids25m$geomorph.c)
grids25m$BOD50 <- as.factor(grids25m$BOD50)
proj4string(grids25m) <- CRS(NL_RD) 
gridcell <- grids25m@grid@cellsize[[1]]
# write to ILWIS:
# writeGDAL(grids25m["geomorph"], "ilwis/geomorph.mpr", "ILWIS") 
# write.asciigrid(grids25m["geomorph"], "geomorph_eng.asc", na.value=-1)

# Geomorphological map:
geomorph <- readOGR("geomorph.shp", "geomorph")
proj4string(geomorph) <- CRS(NL_RD)
# determine the coordinates of the centre:
geomorph.ll <- spTransform(geomorph, CRS("+proj=longlat +ellps=WGS84"))
geomorph.ll@bbox
clon <- mean(geomorph.ll@bbox[1,])
clat <- mean(geomorph.ll@bbox[2,])

geomorph.utm <- spTransform(geomorph, CRS("+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
geomorph.utm@bbox
# writeOGR(geomorph.utm, "geomorph_utm.shp", "geomorph_utm", "ESRI Shapefile")

##############################################################
#    1.2 Derivation of the TOPO-DEM
##############################################################

# Read the field measured evalations:
hoogte_16ef <- readOGR("hoogte_16ef.shp", "hoogte_16ef")
# set the correct proj4 string:
proj4string(hoogte_16ef) <- CRS(NL_RD)
# determine the search radius:
search.r <- sqrt(diff(hoogte_16ef@bbox[1,])^2 + diff(hoogte_16ef@bbox[2,])^2)/4
summary(hoogte_16ef$JAAR[hoogte_16ef$JAAR<2009]) # years from 1959 to 1966;
hist(hoogte_16ef$HOOGTE, col="grey", breaks=25)
sd(hoogte_16ef$HOOGTE)
# plot the variogram map:
plot(variogram(HOOGTE~1, hoogte_16ef, map=TRUE, cutoff=sqrt(areaSpatialGrid(grids25m))/2, width=16*gridcell), col.regions=bpy.colors())
# estimate the semivariances:
HOOGTE.svar <- variogram(HOOGTE~1, hoogte_16ef, alpha=c(22.5,112.5))  
plot(HOOGTE.svar, pch="+", cex=2) # obvious anisotropy!
# initial variogram model:
HOOGTE.ivgm <- vgm(nugget=0, model="Lin", anis=c(p=35, s=0.4))
# fit the variogram:
HOOGTE.vgm <- fit.variogram(HOOGTE.svar, model=HOOGTE.ivgm)
plot(HOOGTE.svar, HOOGTE.vgm, pch="+", cex=2)  # linear variogram model!
HOOGTE.vgm

# Generate a DEM:
grids25m$DEM25TOPO <- krige(HOOGTE~1, hoogte_16ef, HOOGTE.vgm, newdata=grids25m, nmax=20)$var1.pred
spplot(grids25m["DEM25TOPO"], col.regions=terrain.colors(35)) 
write.asciigrid(grids25m["DEM25TOPO"], "DEM25TOPO.asc", na.value=-999)
rsaga.esri.to.sgrd(in.grids="DEM25TOPO.asc", out.sgrds="DEM25TOPO.sgrd", in.path=getwd())
# Compare with the Thin plate splines (SAGA):
# rsaga.geoprocessor(lib="grid_spline", module=1, param=list(GRID="DEM25TOPO.sgrd", SHAPES="hoogte_16ef.shp", FIELD=4, TARGET=0, SELECT=1, RADIUS=search.r/3, REGUL=1.5, MAXPOINTS=40, USER_CELL_SIZE=gridcell, USER_X_EXTENT_MIN=grids25m@bbox[1,1]+gridcell/2, USER_X_EXTENT_MAX=grids25m@bbox[1,2]-gridcell/2, USER_Y_EXTENT_MIN=grids25m@bbox[2,1]+gridcell/2, USER_Y_EXTENT_MAX=grids25m@bbox[2,2]-gridcell/2))
# writeGDAL(grids25m["DEM25TOPO"], "ilwis/DEM25TOPO.mpr", "ILWIS")

##############################################################
#    1.3 Filtering of the LIDAR DEM
##############################################################

# first run the residual analysis and mask out spikes:
rsaga.geoprocessor(lib="geostatistics_grid", 0, param=list(INPUT="ahn5m.sgrd", MEAN="tmp.sgrd", STDDEV="tmp.sgrd", RANGE="tmp.sgrd", DEVMEAN="tmp.sgrd", PERCENTILE="tmp.sgrd", RADIUS=14, DIFF="dif_ahn5m.sgrd")) # all features >14 pixels remain!
# then maks out all artificial borders/roads (derive representativness): 
rsaga.geoprocessor(lib="geostatistics_grid", 1, param=list(INPUT="ahn5m.sgrd", RESULT="rps_ahn5m.sgrd", RADIUS=7, EXPONENT=1.4))   # arbitrary parameters!
# read back to R and mask out all areas:
rsaga.sgrd.to.esri(in.sgrd=c("dif_ahn5m.sgrd", "rps_ahn5m.sgrd"), out.grids=c("dif_ahn5m.asc", "rps_ahn5m.asc"), out.path=getwd(), prec=1)
grids5m$dif <- readGDAL("dif_ahn5m.asc")$band1
grids5m$rps <- readGDAL("rps_ahn5m.asc")$band1
grids5m$ahn5mf <- ifelse(grids5m$dif<=-110|grids5m$dif>=110|grids5m$rps<=6.5, NA, grids5m$ahn5m)
summary(grids5m$ahn5mf)  # 10% pixels masked out!
write.asciigrid(grids5m["ahn5mf"], "ahn5mf.asc", na.value=-1)

# Up-scale the filtered 5 m LiDAR DEM:
rsaga.esri.to.sgrd(in.grids="ahn5mf.asc", out.sgrds="ahn5mf.sgrd", in.path=getwd())
rsaga.geoprocessor(lib="grid_tools", module=23, param=list(GRID="DEM25LIDAR.sgrd", M_EXTENT=1, CELLSIZE=gridcell, XMIN=grids25m@bbox[1,1]+gridcell/2, XMAX=grids25m@bbox[1,2]-gridcell/2, YMIN=grids25m@bbox[2,1]+gridcell/2, YMAX=grids25m@bbox[2,2]-gridcell/2))  # empty new grid;
rsaga.geoprocessor(lib="grid_tools", module=0, param=list(INPUT="ahn5mf.sgrd", GRID="DEM25LIDAR.sgrd", GRID_GRID="DEM25LIDAR.sgrd", METHOD=2, KEEP_TYPE=F, SCALE_UP_METHOD=3))
# close the gaps (fill in the missing pixels!)
rsaga.geoprocessor(lib="grid_tools", module=7, param=list(INPUT="DEM25LIDAR.sgrd", RESULT="DEM25LIDAR.sgrd"))
rsaga.geoprocessor("grid_calculus", module=1, param=list(INPUT="DEM25LIDAR.sgrd", RESULT="DEM25LIDAR.sgrd", FORMUL="a/100")) # convert to meters!

rsaga.sgrd.to.esri(in.sgrd="DEM25LIDAR.sgrd", out.grids="DEM25LIDAR.asc", out.path=getwd(), prec=1)
grids25m$DEM25LIDAR <- readGDAL("DEM25LIDAR.asc")$band1
writeGDAL(grids25m["lidarDEM"], "ilwis/DEM25LIDAR.mpr", "ILWIS")
# difference between the two DEMs:
grids25m$DEM.dif <- grids25m2$DEM25TOPO-grids25m$DEM25LIDAR
spplot(grids25m["DEM.dif"], col.regions=bpy.colors(25)) 
writeGDAL(grids25m["DEM.dif"], "ilwis/DEM_dif.mpr", "ILWIS")


##############################################################
#    1.5 Optional: Download of SRTM DEM
##############################################################

# http://srtm.csi.cgiar.org/SELECTION/inputCoord.asp
download.file("ftp://anon:password@137.73.24.202/array1/portal/srtm41/srtm_data_geotiff/srtm_38_02.zip", destfile=paste(getwd(), "srtm_38_02.zip", sep="/"))
fname <- zip.file.extract(file="srtm_38_02.tif", zipname="srtm_38_02.zip")
file.copy(fname, "./srtm_38_02.tif", overwrite=TRUE)
unlink("srtm_38_02.zip")
DEMsrtm.info <- GDALinfo("srtm_38_02.tif")
# import map to SAGA format:
rsaga.geoprocessor(lib="io_gdal", 0, param=list(GRIDS="srtm_38_02.sgrd", FILE="srtm_38_02.tif"))
unlink("srtm_38_02.tif")
# create an empty grid:
rsaga.geoprocessor(lib="grid_tools", module=23, param=list(GRID="DEMsrtm.sgrd", M_EXTENT=1, XMIN=grids25m@bbox[1,1]+50, YMIN=grids25m@bbox[2,1]+50, XMAX=grids25m@bbox[1,2]-50, YMAX=grids25m@bbox[2,2]-50, CELLSIZE=100)) # empty grid!
# resample to the same study area:
rsaga.geoprocessor(lib="pj_proj4", 6, param=list(SOURCE_PROJ=as.character(attr(DEMsrtm.info, "projection")), TARGET_PROJ=NL_RD, SOURCE="srtm_38_02.sgrd", TARGET="DEMsrtm.sgrd", TARGET_TYPE=2, GET_GRID_GRID="DEMsrtm.sgrd", INTERPOLATION=2))

##############################################################
#    2. Derivation of DEM-parameters
#
##############################################################

for(DEM.i in c("DEM25LIDAR.sgrd", "DEM25TOPO.sgrd")){
# filter the spurious sinks:
rsaga.geoprocessor(lib="ta_preprocessor", module=2, param=list(DEM=DEM.i, RESULT="DEM25LIDAR_f.sgrd"))
# SAGA wetness index:
rsaga.geoprocessor(lib="ta_hydrology", module=15, param=list(DEM="DEM25LIDAR_f.sgrd", C="tmp.sgrd", GN="tmp.sgrd", CS="tmp.sgrd", SB="twi.sgrd", T=120))  # floating point
# valley depth:
rsaga.geoprocessor(lib="ta_morphometry", module=14, param=list(DEM=DEM.i, HO="tmp.sgrd", HU="vdepth.sgrd", NH="tmp.sgrd", SH="tmp.sgrd", MS="tmp.sgrd", W=12, T=120, E=4))  # what do these parameters mean?!
# MRVBF:
rsaga.geoprocessor(lib="ta_morphometry", module=8, param=list(DEM=DEM.i, MRVBF="mrvbf.sgrd", MRRTF="tmp.sgrd", T_SLOPE=4, T_PCTL_V=0.6, T_PCTL_R=0.55, P_SLOPE=1.2, P_PCTL=0.8))  
# Difference From Mean value / Percentile:
rsaga.geoprocessor(lib="geostatistics_grid", 0, param=list(INPUT=DEM.i, MEAN="tmp.sgrd", STDDEV="tmp.sgrd", RANGE="tmp.sgrd", DEVMEAN="tmp.sgrd", DIFF="dfm.sgrd", PERCENTILE="perc.sgrd", RADIUS=80)) # wide search radius!
# Convergence index:
rsaga.geoprocessor(lib="ta_morphometry", module=2, param=list(ELEVATION=DEM.i, RESULT="convi.sgrd", RADIUS=3, METHOD=0, SLOPE=TRUE))  
rsaga.geoprocessor(lib="grid_tools", module=7, param=list(INPUT="convi.sgrd", RESULT="convi.sgrd"))
DEM.list <- c("twi.asc", "vdepth.asc", "mrvbf.asc", "dfm.asc", "perc.asc", "convi.asc")
rsaga.sgrd.to.esri(in.sgrds=set.file.extension(DEM.list, ".sgrd"), out.grids=DEM.list, prec=1, out.path=getwd())
for(i in 1:length(DEM.list)){
  DEM.name <- strsplit(strsplit(DEM.i, "25")[[1]][2], ".sgrd")[[1]]
  grids25m@data[paste(strsplit(DEM.list[i], ".asc")[[1]], "_", DEM.name, sep="")] <- readGDAL(DEM.list[i])$band1
}
}
# this takes few minutes...

str(grids25m@data)
hist(grids25m$twi_LIDAR, col="grey")
hist(grids25m$vdepth_LIDAR, col="grey")
hist(grids25m$dfm_LIDAR, col="grey")

##############################################################
#    3. Multinomial Logistic Regression
# 
##############################################################

##############################################################
#    3.1 Selection of training pixels
##############################################################

# convert the polygon to line map:
rsaga.geoprocessor(lib="shapes_lines", module=0, param=list(POLYGONS="geomorph.shp", LINES="geomorph_l.shp"))
# Derive the buffer map using the shape file:
rsaga.geoprocessor(lib="grid_gridding", module=0, param=list(GRID="geomorph_r.sgrd", INPUT="geomorph_l.shp", FIELD=0, LINE_TYPE=0, TARGET_TYPE=0, USER_CELL_SIZE=gridcell, USER_X_EXTENT_MIN=grids25m@bbox[1,1]+gridcell/2, USER_X_EXTENT_MAX=grids25m@bbox[1,2]-gridcell/2, USER_Y_EXTENT_MIN=grids25m@bbox[2,1]+gridcell/2, USER_Y_EXTENT_MAX=grids25m@bbox[2,2]-gridcell/2))
# buffer distance:
rsaga.geoprocessor(lib="grid_tools", module=10, param=list(SOURCE="geomorph_r.sgrd", DISTANCE="geomorph_dist.sgrd", ALLOC="tmp.sgrd", BUFFER="tmp.sgrd", DIST=sqrt(areaSpatialGrid(grids25m))/3, IVAL=gridcell))
# surface specific points (medial axes!):
rsaga.geoprocessor(lib="ta_morphometry", module=3, param=list(ELEVATION="geomorph_dist.sgrd", RESULT="geomorph_medial.sgrd", METHOD=1))
# read to R:
rsaga.sgrd.to.esri(in.sgrds="geomorph_medial.sgrd", out.grids="geomorph_medial.asc", prec=0, out.path=getwd())
grids25m$geomorph_medial <- readGDAL("geomorph_medial.asc")$band1
# Generate the training pixels:
grids25m$weight <- abs(ifelse(grids25m$geomorph_medial>=0, 0, grids25m$geomorph_medial))
dens.weight <- as.im(as.image.SpatialGridDataFrame(grids25m["weight"]))
training.pix <- rpoint(length(grids25m$weight)/20, f=dens.weight)
training.pix <- data.frame(x=training.pix$x, y=training.pix$y , no=1:length(training.pix$x))
coordinates(training.pix) <- ~x+y
spplot(grids25m["weight"], col.regions=grey(rev(seq(0,0.95,0.05))), sp.layout=list("sp.points", pch="+", training.pix, col="yellow"))
writeOGR(training.pix, "training_pix.shp", "training.pix", "ESRI Shapefile")

##############################################################
#    3.2 Model fitting and prediction
##############################################################

# DEM25LIDAR
# Run MLR:
training.pix.ov <- overlay(grids25m, training.pix)
library(nnet)
# fit the model:
mlr.geomorph <- multinom(geomorph.c~DEM25LIDAR+twi_LIDAR+vdepth_LIDAR+mrvbf_LIDAR+dfm_LIDAR+perc_LIDAR+convi_LIDAR, training.pix.ov)
summary(mlr.geomorph)
grids25m$geomorph.mlr <- predict(mlr.geomorph, newdata=grids25m)
spplot(grids25m["geomorph.mlr"], col.regions=rainbow(length(levels(grids25m$geomorph.c))))
grids25m$geomorphmlr <- as.integer(grids25m$geomorph.mlr)
write.asciigrid(grids25m["geomorphmlr"], "geomorph_mlr.asc", na.value=-1)
# writeGDAL(grids25m["geomorphmlr"], "ilwis/geomorph_mlr.mpr", "ILWIS")

# DEM25TOPO
mlr2.geomorph <- multinom(geomorph.c~DEM25TOPO+twi_TOPO+vdepth_TOPO+mrvbf_TOPO+dfm_TOPO+perc_TOPO+convi_TOPO, training.pix.ov)
grids25m$geomorph.mlr2 <- predict(mlr2.geomorph, newdata=grids25m)
grids25m$geomorphmlr2 <- as.integer(grids25m$geomorph.mlr2)
write.asciigrid(grids25m["geomorphmlr2"], "geomorph_mlr2.asc", na.value=-1)

# the kappa statistics (the whole map!)
sel <- !is.na(grids25m$geomorph.c)
Kappa(confusion(grids25m$geomorph.c[sel], grids25m$geomorph.mlr[sel]))
agreementplot(confusion(grids25m$geomorph.c, grids25m$geomorph.mlr))
Kappa(confusion(grids25m$geomorph.c[sel], grids25m$geomorph.mlr2[sel]))


##############################################################
#    4. Unsupervised classification
# 
##############################################################

##############################################################
#    4.1 Optimal class selection
##############################################################

# fuzzy k-means:
pc.dem <- prcomp(~DEM25LIDAR+twi+vdepth+mrvbf+dfm+perc+convi, scale=TRUE, grids25m@data)
biplot(pc.dem, arrow.len=0.1, xlabs=rep(".", length(pc.dem$x[,1])), main="PCA biplot")
# twi and mrvbf; perc and dfm are highly correlated; use PCs instead;

# Determine number of clusters
demdata <- as.data.frame(pc.dem$x)
wss <- (nrow(demdata)-1)*sum(apply(demdata,2,var))
for (i in 2:20) {wss[i] <- sum(kmeans(demdata, centers=i)$withinss)}
plot(1:20, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")
# does not converge :(

# fuzzy k-means:
kmeans.dem <- kmeans(demdata, length(levels(grids25m$geomorph.c)))
grids25m$kmeans.dem <- kmeans.dem$cluster
grids25m$landform <- as.factor(kmeans.dem$cluster)
summary(grids25m$landform)
spplot(grids25m["landform"], col.regions=rainbow(length(levels(grids25m$geomorph.c))))
write.asciigrid(grids25m["kmeans.dem"], "landform.asc", na.value=-1)

# end of script!