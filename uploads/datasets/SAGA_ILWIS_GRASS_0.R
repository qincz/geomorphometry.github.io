
# R + SAGA/ILWIS/GRASS integration
# http://geomorphometry.org/R.asp

######## R script ###############

# DEM processing and extraction of channel networks;

library(maptools)
library(rgdal)
library(RSAGA) 
# ! first download and install SAGA GIS [http://www.saga-gis.org], ILWIS GIS [https://52north.org/download/Ilwis/52n-Ilwis-v3-05-02.zip] and GRASS GIS [http://grass.itc.it] to your machine first!
rsaga.env(path="C:/Progra~1/saga_vc")
ILWIS <- "C:\\Progra~1\\N52\\Ilwis35\\IlwisClient.exe -C"
MGI_Z6 <- "+proj=tmerc +lat_0=0 +lon_0=18 +k=0.9999 +x_0=6500000 +y_0=0 +ellps=bessel +towgs84=550.499,164.116,475.142,5.80967,2.07902,-11.62386,0.99999445824 +units=m"
# The coordinate system is "MGI Zone 6"; see [http://spatial-analyst.net/wiki/index.php?title=MGI_/_Balkans_coordinate_systems]
setwd("C:/tmp")

# obtain the data:
download.file("http://geomorphometry.org/system/files/BaranjaHill.zip", destfile=paste(getwd(), "BaranjaHill.zip", sep="/"))
fname <- zip.file.extract(file="DEM25m.asc", zipname="BaranjaHill.zip")
file.copy(fname, "./DEM25m.asc", overwrite=TRUE)
list.files(pattern="asc")

# view the data in ILWIS:
dem25m <- readGDAL("DEM25m.asc")
proj4string(dem25m) <- CRS(MGI_Z6)
writeGDAL(dem25m[1], "dem25m.mpr", "ILWIS")
# set the right coordinate system!
download.file("http://spatial-analyst.net/CRS/gk_6.csy", destfile=paste(getwd(), "gk_6.csy", sep="/"))
shell(cmd=paste(ILWIS, "setcsy dem25m.grf gk_6.csy -force"), wait=F)
shell(cmd=paste(ILWIS, "open dem25m.mpr -noask"), wait=F)
# make a relief view in ILWIS:
shell(cmd=paste(ILWIS, "run C:\\Progra~1\\Ilwis3\\Scripts\\Hydro-DEM\\dem_visualization dem25m.mpr dem25m_c.mpr"), wait=F)

# load to SAGA and derive drainage network:
rsaga.esri.to.sgrd(in.grids="dem25m.asc", out.sgrds="dem25m.sgrd", in.path=getwd())
# First, fill the spurious sinks:
rsaga.get.modules("ta_preprocessor")
rsaga.get.usage("ta_preprocessor", 2)
rsaga.geoprocessor(lib="ta_preprocessor", module=2, param=list(DEM="dem25m.sgrd", RESULT="dem25m_f.sgrd", MINSLOPE=0.05))
# Second, extract the channel network:
rsaga.get.modules("ta_channels")
rsaga.get.usage("ta_channels", 0)
rsaga.geoprocessor(lib="ta_channels", module=0, param=list(ELEVATION="dem25m.sgrd", CHNLNTWRK="channel_ntwrk.sgrd", CHNLROUTE="channel_route.sgrd", SHAPES="channels2.shp", INIT_GRID="dem25m.sgrd", DIV_CELLS=3, MINLEN=40))
channels <- readOGR("channels.shp", "channels")
spplot(channels["LENGTH"], col.regions=bpy.colors())
# derive Topographic Wetness Index:
rsaga.geoprocessor(lib="ta_hydrology", module=15, param=list(DEM="dem25m.sgrd", C="catharea.sgrd", GN="catchslope.sgrd", CS="modcatharea.sgrd", SB="TWI.sgrd", T=10))

# Export of grids to Google Earth using SAGA GIS:
rsaga.geoprocessor(lib="pj_proj4", 2, param=list(SOURCE_PROJ=paste('"', proj4string(dem25m), '"', sep=""), TARGET_PROJ="\"+proj=longlat +datum=WGS84\"", SOURCE="TWI.sgrd", TARGET="TWI_ll.sgrd", TARGET_TYPE=0, INTERPOLATION=0))
# export to PNG:
rsaga.geoprocessor(lib="io_grid_image", 0, param=list(GRID="TWI_ll.sgrd", FILE="TWI.png"))
# read back to R:
rsaga.sgrd.to.esri(in.sgrds="TWI_ll.sgrd", out.grids="TWI_ll.asc", prec=1, out.path=getwd())
TWI.ll <- readGDAL("TWI_ll.asc")
proj4string(TWI.ll) <- CRS("+proj=longlat +datum=WGS84")
TWI.kml <- GE_SpatialGrid(TWI.ll)
kmlOverlay(TWI.kml, kmlfile="TWI.kml", imagefile="TWI.png", name="Topographic Wetness Index")

# Optional: export to Google Earth using ILWIS:
dem25m.ll <- spTransform(dem25m[1], CRS("+proj=longlat +datum=WGS84"))
corrf <- (1+cos((dem25m.ll@bbox[2,2]+dem25m.ll@bbox[2,1])/2*pi/180))/2
geogrd.cell <- corrf*(dem25m.ll@bbox[1,2]-dem25m.ll@bbox[1,1])/dem25m@grid@cells.dim[1]
geoarc <- spsample(dem25m.ll, type="regular", cellsize=c(geogrd.cell,geogrd.cell))
gridded(geoarc) <- TRUE
gridparameters(geoarc)
gridparameters(dem25m)

# resample the map (Bilinear) to the new geographic grid: 
shell(cmd=paste(ILWIS, " crgrf geoarc.grf ",geoarc@grid@cells.dim[[2]]," ",geoarc@grid@cells.dim[[1]]," -crdsys=LatlonWGS84 -lowleft=(",geoarc@grid@cellcentre.offset[[1]],",",geoarc@grid@cellcentre.offset[[2]],") -pixsize=",geoarc@grid@cellsize[[1]],sep=""), wait=F)
shell(cmd=paste(ILWIS, "dem25m_ll_c.mpr = MapResample(dem25m_c.mpr, geoarc, BiLinear)"), wait=F)
shell(cmd=paste(ILWIS, "open dem25m_ll_c.mpr -noask"))
shell(cmd=paste(ILWIS, "export tiff(dem25m_ll_c.mpr, dem25m_c.tif)"), wait=F)

# generate a KML (ground overlay):
dem25m.kml <- GE_SpatialGrid(geoarc)
kmlOverlay(dem25m.kml, "dem25m_c.kml", "dem25m_c.tif", name="Shaded relief in ILWIS")

# export of extracted channels to Google Earth:
proj4string(channels) <- CRS(MGI_Z6)
channels.ll <- spTransform(channels[3], CRS("+proj=longlat +datum=WGS84"))
writeOGR(channels.ll, "channels.kml", "channels", "KML")

# extract the drainage network in GRASS GIS:
library(spgrass6) # version => 0.6-1 (http://spatial.nhh.no/R/Devel/spgrass6_0.6-1.zip)
# Location of your GRASS installation:
loc <- initGRASS("C:/GRASS", home=tempdir())
loc
# Import the ArcInfo ASCII file to GRASS:
parseGRASS("r.in.gdal")  # commmand description
execGRASS("r.in.gdal", flags="o", parameters=list(input="DEM25m.asc", output="DEM"))
execGRASS("g.region", parameters=list(rast="DEM"))
gmeta6()
# extract the drainage network:
execGRASS("r.watershed", flags=c("m", "overwrite"), parameters=list(elevation="DEM", stream="stream", threshold=as.integer(50)))
# thin the raster map so it can be converted to vectors:
execGRASS("r.thin", parameters=list(input="stream", output="streamt"))
# convert to vectors:
execGRASS("r.to.vect", parameters=list(input="streamt", output="streamt", feature="line"))
streamt <- readVECT6("streamt")
plot(streamt)

######## R script ############### 
