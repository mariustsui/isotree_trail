library("lidR")
library("future")
library("raster")
library("sf")
library("dplyr")
remove(list = ls())
gc()

col <- height.colors(50)
ori_las <- readLAS("D:\\PROJECT_TREES\\CL2_BE33_2021_1000_0504.las", filter="-drop_classification 6 7 18")

las <- classify_noise(ori_las, sor(15,7))
las_denoise <- filter_poi(las, Classification != LASNOISE)

n_las <- normalize_height(las_denoise, knnidw())


tree_li2012 <- segment_trees(n_las, li2012(dt1 = 4, dt2 = 8, R = 0, Zu = 40, hmin = 20, speed_up = 8))

crowns_li2012 <- crown_metrics(  tree_li2012,  type = "convex",  func = .stdtreemetrics, attribute = "treeID", geom = "convex")
plot(crowns_li2012["Z"], main = "Crown area (convex hull)", col = height.colors(20))

st_write(crowns_li2012, "D:\\PROJECT_TREES\\isotree_trail\\crowns_0504_li2012.shp", append=FALSE)

plot(tree_li2012, color="treeID", colorPalette = pastel.colors(500))