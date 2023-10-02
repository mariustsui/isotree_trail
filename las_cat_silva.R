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
  # function to define the window size for locate_trees
  getWS <- function(x){
    y <- 2.6125 + x^0.665
    y[x<2] <- 2
    y[x>40] <- 15
    return(y)
  }
  
  # function to smooth the derived chm by focal mean
  fill.na <- function(x, i=5) { if (is.na(x)[i]) { return(mean(x, na.rm = TRUE)) } else { return(x[i]) }}
  w <- matrix(1, 3, 3)
  
  f_treecrown <- function(las){
    
    las <- classify_noise(las, sor(15,7))
    las_denoise <- filter_poi(las, Classification != LASNOISE)
    nlas <- normalize_height(las_denoise, knnidw())
    tree_tops <- locate_trees(nlas, lmf(ws = getWS, hmin = hmin))
    chm <- rasterize_canopy(nlas, 0.3, pitfree(subcircle = 0.2))
    smooth_chm <- terra::focal(chm, w, fun = mean, na.rm = TRUE)
    tree <- segment_trees(nlas, silva2016(smooth_chm, tree_tops, max_cr_factor = 0.3, exclusion = 0.25, ID="treeID"))
    crowns <- crown_metrics(  tree,  type = "convex",  func = .stdtreemetrics, attribute = "treeID", geom = "convex")
    return(crowns)
  }
  
  opt    <- list(automerge = TRUE)   # catalog_apply will merge the outputs into a single object
  output = NULL
  output <- catalog_map(las_catg, f_treecrown, .options = opt)
  names(output) <- basename(las_catg$filename)
  
  lapply( seq_along(output), function(i){
    single_sf <- dplyr::bind_rows(output[i])
    types <- st_geometry_type(single_sf)
    types_df <- data.frame(types)
    my_labelled_sf_object <- merge(single_sf, types_df,by.x=0, by.y=0, all.x=TRUE)
    my_filtered_sf_object <- my_labelled_sf_object[my_labelled_sf_object$types == "POLYGON",]
    st_write(my_filtered_sf_object, paste0("D:\\PROJECT_TREES\\R_OUTPUT\\CROWNS\\",substring(names(output)[[i]], 1, 23),'_',hmin,'PLUS.shp' ), append=FALSE)
  })
  
  single_sf <- dplyr::bind_rows(output)
  types <- st_geometry_type(single_sf)
  types_df <- data.frame(types)
  my_labelled_sf_object <- merge(single_sf, types_df,by.x=0, by.y=0, all.x=TRUE)
  my_filtered_sf_object <- my_labelled_sf_object[my_labelled_sf_object$types == "POLYGON",]
  st_write(my_filtered_sf_object, paste0("D:\\PROJECT_TREES\\R_OUTPUT\\CROWNS",'_',hmin,'PLUS_SILVA.shp' ), append=FALSE)
  
  
  