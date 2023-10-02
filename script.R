library("lidR")
library("future")
library("raster")
library("sf")
library("dplyr")
remove(list = ls())

col <- height.colors(50)
ori_las <- readLAS("D:\\PROJECT_TREES\\REQ201107_Lars\\laser\\CL2_BE33_2021_1000_0504.las", filter="-drop_classification 6 7 18")

las <- classify_noise(ori_las, sor(15,7))
las_denoise <- filter_poi(las, Classification != LASNOISE)

n_las <- normalize_height(las_denoise, knnidw())

getWS <- function(x){
  y <- 2.6125 + x^0.665
  y[x<2] <- 2
  y[x>40] <- 15
  return(y)
}

chm <- rasterize_canopy(n_las, 0.3, pitfree(subcircle = 0.2))

# to smooth the chm
fill.na <- function(x, i=5) { if (is.na(x)[i]) { return(mean(x, na.rm = TRUE)) } else { return(x[i]) }}
w <- matrix(1, 3, 3)

smooth_chm <- terra::focal(chm, w, fun = mean, na.rm = TRUE)
tree_tops <- locate_trees(n_las, lmf(ws = getWS, hmin = 18))
tree_tops2d <- sf::st_zm(tree_tops)
plot(smooth_chm, col = col)
plot(sf::st_geometry(tree_tops), add = TRUE, pch = 3)
tree <- segment_trees(n_las, silva2016(smooth_chm, tree_tops, max_cr_factor = 0.3, exclusion = 0.25, ID="treeID"))

#tree <- segment_trees(n_las, li2012(dt1 = 1.5, dt2 = 2, R = 2, Zu = 30, hmin = 20))
crowns <- crown_metrics(  tree,  type = "convex",  func = .stdtreemetrics, attribute = "treeID", geom = "convex")

st_write(crowns, "D:\\PROJECT_TREES\\isotree_trail\\crowns_0504.shp", append=FALSE)


