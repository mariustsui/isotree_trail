library("lidR")
library("future")
library("raster")
library("sf")
library("dplyr")
remove(list = ls())

col <- height.colors(50)
ori_las <- readLAS("D:\\PROJECT_TREES\\CL2_BE33_2021_1000_1019.las", filter="-drop_classification 6 7 18")
n_las <- normalize_height(ori_las, knnidw())

tree_li2012 <- segment_trees(n_las, li2012(dt1 = 3, dt2 = 1, R = 10, Zu = 15, hmin = 18))

crowns_li2012 <- delineate_crowns(  tree_li2012,  type = "convex",  func = .stdmetrics)
plot(crowns_li2012["area"], main = "Crown area (convex hull)")
crowns_li2012 <- st_as_sf(crowns_li2012)
st_write(crowns_li2012, "D:\\PROJECT_TREES\\isotree_trail\\crowns_li2012.shp", append=FALSE)