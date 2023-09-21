library("lidR")
library("future")
library("raster")
library("sf")
library("sp")
library("rgdal")
remove(list = ls())

# read all the point cloud (.las) into a lasCatalog
  las_catg <- readLAScatalog("D:\\PROJECT_TREES\\REQ201107_Lars\\laser\\", filter="-drop_classification 6 7 18")

# Locate Trees from LAS
  # function to define the window size for locate_trees
  getWS <- function(x){
    y <- 2.6125 + x^0.665
    y[x<2] <- 2
    y[x>40] <- 15
    return(y)
  }
  
  fill.na <- function(x, i=5) { if (is.na(x)[i]) { return(mean(x, na.rm = TRUE)) } else { return(x[i]) }}
  w <- matrix(1, 3, 3)
  
  f_treecrown <- function(las){
    nlas <- normalize_height(las, knnidw())
    tree_tops <- locate_trees(nlas, lmf(ws = getWS, hmin = 18))
    chm <- rasterize_canopy(nlas, 0.3, pitfree(subcircle = 0.2))
    smooth_chm <- terra::focal(chm, w, fun = mean, na.rm = TRUE)
    tree <- segment_trees(nlas, silva2016(smooth_chm, tree_tops, max_cr_factor = 0.3, exclusion = 0.25, ID="treeID"))
    crowns <- delineate_crowns(  tree,  type = "convex",  func = .stdtreemetrics)
    return(crowns)
  }
  
  opt    <- list(automerge = TRUE)   # catalog_apply will merge the outputs into a single object
  output <- catalog_map(las_catg, f_treecrown, .options = opt)
  output <- st_as_sf(output)
  st_write(output, "D:\\PROJECT_TREES\\R_OUTPUT\\CROWNS.shp", layer = "Tree_Crowns",append=FALSE)
  
  