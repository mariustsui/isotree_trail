library("lidR")
library("future")
library("raster")
library("sf")
library("sp")
library("dplyr")
remove(list = ls())
gc()

# Set the minimum height of a tree to be detected ( in meters )
hmin = 20

# read all the point cloud (.las) into a lasCatalog
las_catg <- readLAScatalog("D:\\PROJECT_TREES\\REQ201107_Lars\\laser\\", filter="-drop_classification 6 7 18")
#las_catg <- readLAScatalog("D:\\PROJECT_TREES\\LAS\\", filter="-drop_classification 6 7 18")
las_catg
las_check(las_catg)

# Locate Trees from LAS

start <- Sys.time()

f_treecrown <- function(las){
  
  las <- classify_noise(las, sor(15,7))
  las_denoise <- filter_poi(las, Classification != LASNOISE)
  n_las <- normalize_height(las_denoise, knnidw())
  tree <- segment_trees(n_las, li2012(dt1 = 4, dt2 = 8, R = 0, Zu = 40, hmin = 20, speed_up = 8))
  crowns <- crown_metrics(  tree,  type = "convex",  func = .stdtreemetrics, attribute = "treeID", geom = "convex")
  diff <- Sys.time() - start
  print(paste0("Lapsed for ", format(round(diff, 2), nsmall = 2)))
  return(crowns)
  
}

opt    <- list(automerge = TRUE)   # catalog_apply will merge the outputs into a single object
output = NULL
output <- catalog_map(las_catg, f_treecrown, .options = opt)
op <- dplyr::bind_rows(output)
types <- st_geometry_type(op)
types_df <- data.frame(types)
my_labelled_sf_object <- merge(op, types_df,by.x=0, by.y=0, all.x=TRUE)
my_filtered_sf_object <- my_labelled_sf_object[my_labelled_sf_object$types == "POLYGON",]
st_write(my_filtered_sf_object , paste0("D:\\PROJECT_TREES\\R_OUTPUT\\CROWNS",'_',hmin,'PLUS_li2012.shp' ), append=FALSE)


