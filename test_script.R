library("lidR")
library("future")
library("raster")
library("sf")
library("dplyr")
library("rgeos")
remove(list = ls())

col <- height.colors(50)
ori_las <- readLAS("D:\\PROJECT_TREES\\PIRONGIA\\laser\\CL2_BE33_2021_1000_1019.las", filter="-drop_classification 6 7 18")
las <- classify_noise(ori_las, sor(15,7))
las_denoise <- filter_poi(las, Classification != LASNOISE)

n_las <- normalize_height(las_denoise, knnidw())
filtered_las = filter_poi(n_las, Z>=1, Z<=1.5)
pd = rasterize_density(filtered_las,5)
plot(pd)
writeRaster(pd, "D:\\PROJECT_TREES\\isotree_trail\\pt_density_100CM_150CM.tif", overwrite=TRUE)

summary(las)

getWS <- function(x){
  y <- 1.6125 + x^0.665
  y[x<2] <- 2
  return(y)
}

chm <- rasterize_canopy(n_las, 1, pitfree(subcircle = 0.2))

# to smooth the chm
fill.na <- function(x, i=5) { if (is.na(x)[i]) { return(mean(x, na.rm = TRUE)) } else { return(x[i]) }}
w <- matrix(1, 3, 3)

smooth_chm <- terra::focal(chm, w, fun = mean, na.rm = TRUE)
tree_tops <- locate_trees(n_las, lmf(ws = getWS, hmin = 5))
tree_tops2d <- sf::st_zm(tree_tops)
plot(smooth_chm, col = col)
plot(sf::st_geometry(tree_tops), add = TRUE, pch = 3)
tree <- segment_trees(n_las, silva2016(smooth_chm, tree_tops, max_cr_factor = 0.3, exclusion = 0.25, ID="treeID"))
#tree <- segment_trees(n_las, dalponte2016(smooth_chm, tree_tops, th_tree = 20 , max_cr = 20, ID="treeID"))

#tree <- segment_trees(n_las, li2012(dt1 = 4, dt2 = 8, R = 0, Zu = 40, hmin = 20, speed_up = 8))
crowns <- crown_metrics(  tree,  type = "convex",  func = .stdtreemetrics, attribute = "treeID", geom = "convex")
ttop <- gCentroid( spgeom = methods::as( object = crowns, Class = "Spatial" ) , byid = TRUE )



st_write(tree_tops2d, "D:\\PROJECT_TREES\\isotree_trail\\tree_tops_1019_silva_5plus.shp", append=FALSE)
st_write(crowns, "D:\\PROJECT_TREES\\isotree_trail\\crowns_1019_silva_5plus.shp", append=FALSE)
st_write(st_as_sf(ttop), "D:\\PROJECT_TREES\\isotree_trail\\ttop_1019_silva_5plus.shp", append=FALSE)
writeLAS(tree, "D:\\PROJECT_TREES\\isotree_trail\\trees_1019_silva_5plus.las")
plot(crowns['Z'], add=TRUE)


