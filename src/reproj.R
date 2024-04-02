library(raster)

ele <- raster("data/rasters/terrain/ele.tif")

dists <- stack(list.files("data/rasters/sites/distances", full.names = T, pattern = ".tif$"))
dists_wgs <- projectRaster(dists, ele)

for (i in 1:nlayers(dists_wgs)) {
    writeRaster(dists_wgs[[i]], paste("data/rasters/site_dist/", names(dists_wgs)[i], ".tif", sep=""), overwrite=TRUE)
}