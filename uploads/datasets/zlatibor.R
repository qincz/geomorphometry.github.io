# Script to simulate the control to generate DEM and assess the error of the height measurements
# prepared for the paper "Geostatistical simulation of DEMs using auxiliary data"
# Developed by T. Hengl (spatial-analyst.net)
# In Ispra, Italy, March 2007.
# Inputs: control.txt - 1020 precise measurements (photogrametric + spot heights); elevations.txt - 2051 points (contours + spot heights); dem10m_tin.asc - 100x150 pixels 30m DEM derived from contour lines in ArcGIS
# Coordinate system Gauss Kruger zone 7 (gk_7);

# Import the point control (table) in R:

control <- read.delim("control.txt")
elevations <- read.delim("elevations.txt")
SRTMDEM <- read.delim("SRTMDEM.txt")
library(sp)
# convert to spatial
coordinates(control)=~X+Y
coordinates(elevations)=~X+Y
coordinates(SRTMDEM)=~X+Y
spplot(elevations["Z"], col.regions=bpy.colors(), at = seq(850,1170,10))
spplot(control["h"], col.regions=bpy.colors(), at = seq(850,1170,10))
spplot(SRTMDEM["SRTM"])

# Import the raster maps control - demz and errors - deme

library(rgdal)
dem = readGDAL("dem30m.asc")
dem$z = readGDAL("dem30m.asc")$band1
dem$dstream = readGDAL("dstream.asc")$band1
dem$stdev = readGDAL("stdev.asc")$band1
dem$shade = readGDAL("shade.asc")$band1
dem$NDVI = readGDAL("NDVI.asc")$band1
spplot(dem["z"], col.regions=bpy.colors(), at = seq(850,1150,5))
spplot(dem["dstream"], col.regions=bpy.colors())

srtmdem = readGDAL("SRTMDEM.asc")
srtmdem$srtm = readGDAL("SRTMDEM.asc")$band1

# Overlay points over demz and deme

control.ov = overlay(dem, control)
control$z = control.ov$z
control$dstream = control.ov$dstream
control$stdev = control.ov$stdev
control$shade = control.ov$shade
control$NDVI = control.ov$NDVI
control$delta = control$z-control$h
sel = !is.na(control$delta)
spplot(control["delta"], col.regions=bpy.colors())
summary(control["delta"])
hist(control$delta, col="grey")

elevations.ov = overlay(dem, elevations)
elevations$z = elevations.ov$z
elevations$dstream = elevations.ov$dstream
elevations$stdev = elevations.ov$stdev
elevations$shade = elevations.ov$shade
elevations$NDVI = elevations.ov$NDVI
sel = !is.na(elevations$z)

# 1. Geostatistical simulations to directly generate a DEM from control

# 1.1 No auxiliary maps available

library(gstat)
plot(variogram(Z~1, elevations))
elevations.or = variogram(Z~1, elevations)
elevations.ovgm = fit.variogram(elevations.or, vgm(1, "Sph", 1000, 1))
plot(elevations.or, elevations.ovgm, plot.nu=F, pch="+")
elevations.ovgm
hist(elevations$Z, col="grey")

elevations.ok = krige(Z~1, elevations, dem, elevations.ovgm, nmax=80)
spplot(elevations.ok[1], col.regions=bpy.colors(), at = seq(850,1200,5))

elevations.sim = krige(Z~1, elevations, dem, elevations.ovgm, nmax=80, nsim = 2)
spplot(elevations.sim[1], col.regions=bpy.colors(), at = seq(850,1200,5))
summary(elevations.sim[1])

write.asciigrid(elevations.ok[1], "elevations_ok.asc")
write.asciigrid(elevations.sim[1], "elevations_sim002.asc")

# 1.2 Auxiliary maps available
# First check if the predictors are significant


summary(lm(Z~dstream+NDVI, elevations))
scatter.smooth(elevations$dstream, elevations$Z, span=2/3)
scatter.smooth(elevations$NDVI, elevations$Z, span=2/3)

plot(variogram(Z~dstream, elevations[sel,]))
elevations.rev = variogram(Z~dstream, elevations[sel,])
elevations.rvgm = fit.variogram(elevations.rev, vgm(1, "Sph", 1000, 1))
plot(elevations.rev, elevations.rvgm, plot.nu=F, pch="+")
elevations.rvgm

# predictions
elevations.rk = krige(Z~dstream, elevations[sel,], dem, elevations.rvgm, nmin=50, nmax=200)
spplot(elevations.rk[1], col.regions=bpy.colors(), at = seq(850,1200,5), sp.layout=list("sp.points", pch="+", col="black", elevations))

#simulations
elevations.rksim = krige(Z~dstream, elevations[sel,], dem, elevations.rvgm, nmin=20, nmax=80, nsim = 2)
spplot(elevations.rksim[1], col.regions=bpy.colors(), at = seq(850,1200,5))
summary(elevations.rksim[1])

write.asciigrid(elevations.rk[2], "elevations_rkvar.asc")
write.asciigrid(elevations.rksim[2], "elevations_rksim002.asc")

# 2. Geostatistical simulations of the errors

# 2.1 No auxiliary maps available

plot(variogram(delta~1, control))
delta.or = variogram(delta~1, control)
delta.ovgm = fit.variogram(delta.or, vgm(1, "Exp", 1000, 1))
plot(delta.or, delta.ovgm, plot.nu=F, pch="+")
delta.ovgm

# predictions
delta.ok = krige(delta~1, control, dem, delta.ovgm, nmax = 80)
spplot(delta.ok[1], col.regions=bpy.colors(), at = seq(-5,5,0.2), sp.layout=list("sp.points", pch="+", col="black", elevations))
# simulations
delta.sim = krige(delta~1, control, dem, delta.ovgm, nmax = 50, nsim = 4)
spplot(delta.sim[2], col.regions=bpy.colors(), at = seq(-5,5,0.2))

write.asciigrid(delta.ok[1], "delta_ok.asc")
write.asciigrid(delta.sim[2], "delta_sim002.asc")

# 2.2 Auxiliary maps available
# First check the nature of correlation and then fit the model

scatter.smooth(log(control$stdev), abs(control$delta), span=2/3)
scatter.smooth(sqrt(control$shade), abs(control$delta), span=2/3)
scatter.smooth(control$NDVI, abs(control$delta), span=2/3)
summary(lm(abs(delta)~NDVI+log(stdev)+sqrt(shade), control))

plot(variogram(abs(delta)~NDVI+log(stdev)+sqrt(shade), control))
delta.rev = variogram(abs(delta)~NDVI+log(stdev)+sqrt(shade), control)
delta.rvgm = fit.variogram(delta.rev, vgm(1, "Exp", 1000, 1))
plot(delta.rev, delta.rvgm, plot.nu=F, pch="+")
delta.rvgm

delta.rk = krige(abs(delta)~NDVI+log(stdev)+sqrt(shade), control, dem, delta.rvgm, nmax = 80)
spplot(delta.rk[1], col.regions=bpy.colors(), at = seq(0,3,0.2), sp.layout=list("sp.points", pch="+", col="black", elevations))
summary(delta.rk[1])

write.asciigrid(delta.rk[1], "delta_rk.asc")

# 3.1 Modelling SRTM DEM - downscaling to 30 m resolution

# fit the regression and variogram models

summary(lm(SRTM~dstream+NDVI, SRTMDEM))

plot(variogram(SRTM~1, SRTMDEM))
SRTM.or = variogram(SRTM~1, SRTMDEM)
SRTM.ovgm = fit.variogram(SRTM.or, vgm(1, "Bes", 1000, 1))
plot(SRTM.or, SRTM.ovgm, plot.nu=T, pch="+")
SRTM.ovgm 

# We use the original variogram fitted for punctual measurements

sel = !is.na(SRTMDEM$z)

SRTM.rk = krige(SRTM~dstream+NDVI, SRTMDEM[sel,], dem, elevations.rvgm, nmin=50, nmax=200)
spplot(SRTM.rk[1], col.regions=bpy.colors(), at = seq(850,1200,5))

SRTM.rksim = krige(SRTM~dstream+NDVI, SRTMDEM[sel,], dem, elevations.rvgm, nmin=50, nmax=200, nsim=2)
spplot(SRTM.rksim[1], col.regions=bpy.colors(), at = seq(850,1200,5))


write.asciigrid(SRTM.rk[1], "SRTM_rk.asc")

# 3.2 Error assessment for the SRTM DEM

# overlay the maps

SRTM.ov = overlay(dem, SRTMDEM)
SRTMDEM$z = SRTM.ov$z
SRTMDEM$dstream = SRTM.ov$dstream
SRTMDEM$NDVI = SRTM.ov$NDVI

SRTMDEM.ov = overlay(srtmdem, control)
control$srtm = SRTMDEM.ov$srtm
control$delta2 = control$srtm-control$h
spplot(control["delta2"], col.regions=bpy.colors())
summary(control["delta2"])
hist(control$delta2, col="grey")

plot(variogram(delta2~1, control))
delta2.or = variogram(delta2~1, control)
delta2.ovgm = fit.variogram(delta2.or, vgm(1, "Exp", 1000, 1))
plot(delta2.or, delta2.ovgm, plot.nu=F, pch="+")
delta2.ovgm

# simulations
delta2.sim = krige(delta2~1, control, dem, delta2.ovgm, nmax = 50, nsim = 4)
spplot(delta2.sim[2], col.regions=bpy.colors(), at = seq(-25,25,0.5))

write.asciigrid(delta2.sim[2], "delta2_sim002.asc")

scatter.smooth(log(control$stdev), control$delta2, span=2/3)
scatter.smooth(sqrt(control$shade), control$delta2, span=2/3)
scatter.smooth(control$NDVI, control$delta, span=2/3)
summary(lm(delta2~NDVI+stdev+shade, control))

plot(variogram(delta2~NDVI+stdev+shade, control))
delta2.rev = variogram(delta2~NDVI+stdev+shade, control)
delta2.rvgm = fit.variogram(delta2.rev, vgm(1, "Exp", 1000, 1))
plot(delta2.rev, delta2.rvgm, plot.nu=F, pch="+")
delta2.rvgm

# simulations with auxiliary maps
delta2.rksim = krige(delta2~NDVI+stdev+shade, control, dem, delta2.rvgm, nmax = 200, nsim = 4)
spplot(delta2.rksim[1], col.regions=bpy.colors(), at = seq(-25,25,0.5))

write.asciigrid(delta2.rksim[2], "delta2_rksim001.asc")


# Put back the simulated values to the control table

delta.rksim1 = as(delta.rksim[1], "SpatialGridDataFrame")
control.ovd = overlay(delta.rksim1, control)
control$sim1 = control.ovd$sim1
summary(control["sim1"])

write.table(control, file="control_delta.txt")