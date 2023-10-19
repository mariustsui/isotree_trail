library("lidR")
library("future")
library("raster")
library("sf")
library("sp")
library("dplyr")
library("rgeos")
remove(list = ls())
gc()

# Set the minimum height of a tree to be detected ( in meters )

hmin <- 1
loop_index <- 0
forest_metrics <- matrix(nrow = 11, ncol=4)
colnames(forest_metrics) <- c("Number of Trees", "Trees per hectare", "Crown area(m2)", "Crown Cover")
rownames(forest_metrics) <- c("1M", "1.5M", "2M", "5M", "10M", "15M", "20M", "25M", "30M", "35M", "40M")

# read all the point cloud (.las) into a lasCatalog
  
  las_catg <- readLAScatalog("D:\\PROJECT_TREES\\MAUNGATAUTARI\\laser\\", filter="-drop_classification 6 7 18")
  
  #las_catg <- readLAScatalog("D:\\PROJECT_TREES\\PIRONGIA\\laser\\", filter="-drop_classification 6 7 18")
  #las_catg <- readLAScatalog("D:\\PROJECT_TREES\\LAS\\", filter="-drop_classification 6 7 18")
  
  opt_chunk_buffer(las_catg) <- 30
  opt_chunk_size(las_catg)   <- 1000      
  
  las_catg
  las_check(las_catg)

# Locate Trees from LAS
  # function to define the window size for locate_trees
  getWS <- function(x){
    y <- 1.6125 + x^0.665
    y[x<2] <- 2
    return(y)
  }
  
  # function to smooth the derived chm by focal mean
  fill.na <- function(x, i=5) { if (is.na(x)[i]) { return(mean(x, na.rm = TRUE)) } else { return(x[i]) }}
  w <- matrix(1, 3, 3)
  
  f_treecrown <- function(chunk, hmin){
    
    # The chunk argument is a LAScluster object. The users do not need to know how it works.
    # readLAS will load the region of interest (chunk) with a buffer around it, taking advantage of
    # point cloud indexation if possible. The filter and select options are propagated automatically
    las <- readLAS(chunk)
    if (is.empty(las)) return(NULL)
    
    las <- classify_noise(las, sor(15,7))
    las_denoise <- filter_poi(las, Classification != LASNOISE)
    nlas <- normalize_height(las_denoise, knnidw())
    tree_tops <- locate_trees(nlas, lmf(ws = getWS, hmin = hmin))
    
    bbox <- st_bbox(chunk)
    tree_tops <- sf::st_crop(tree_tops, bbox)
    
    
    chm <- rasterize_canopy(nlas, 1, pitfree(subcircle = 0.2))
    smooth_chm <- terra::focal(chm, w, fun = mean, na.rm = TRUE)
    tree <- segment_trees(nlas, silva2016(smooth_chm, tree_tops, max_cr_factor = 0.3, exclusion = 0.25, ID="treeID"))
    crowns <- crown_metrics(  tree,  type = "convex",  func = .stdtreemetrics, attribute = "treeID", geom = "convex")
    diff <- Sys.time() - start
    print(paste0("Lapsed for ", format(round(diff, 2), nsmall = 2)))
    return(crowns)
  }
  
  start <- Sys.time()
  print(paste0("Project started at ", start))
  
  loop_index <- loop_index + 1
  opt    <- list(need_buffer = TRUE, automerge = TRUE)   # catalog_apply will merge the outputs into a single object
  
  output = NULL
  output <- catalog_apply(las_catg, f_treecrown, hmin = hmin, .options = opt)
  
  single_sf <- dplyr::bind_rows(output)
  types <- st_geometry_type(single_sf)
  types_df <- data.frame(types)
  my_labelled_sf_object <- merge(single_sf, types_df,by.x=0, by.y=0, all.x=TRUE)
  my_filtered_sf_object <- my_labelled_sf_object[my_labelled_sf_object$types == "POLYGON",]
  
  
  st_write(my_filtered_sf_object, paste0("D:\\PROJECT_TREES\\MAUNGATAUTARI\\R_OUTPUT\\CROWNS",'_',hmin,'PLUS_SILVA.shp' ), append=FALSE)
  #st_write(my_filtered_sf_object, paste0("D:\\PROJECT_TREES\\isotree_trail\\ALLCROWNS",'_',hmin,'PLUS_SILVA.shp' ), append=FALSE)
  
  
  a_ha <- units::set_units(st_area(las_catg), ha)
  a_m2 <- units::set_units(st_area(las_catg), m^2)


  
  df <- as.data.frame((forest_metrics))
  
  h_list <- c(1, 1.5, 2, 5, 10, 15, 20, 25, 30, 35, 40)
  i <- 1
  for (h in h_list){

    tree <- my_filtered_sf_object[my_filtered_sf_object$Z>=h,]
    df[i, ] <- c(nrow(tree), round(nrow(tree)/ as.numeric(a_ha), 2), round(sum(tree$convhull_area), 2), round( (sum(tree$convhull_area)/as.numeric(a_m2)*100), 2))  
    ttop <- gCentroid( spgeom = methods::as( object = tree, Class = "Spatial" ) , byid = TRUE )
    
    st_write(tree, paste0("D:\\PROJECT_TREES\\MAUNGATAUTARI\\R_OUTPUT\\CROWNS\\CROWNS",'_',h*100,'CM_PLUS_SILVA.shp' ), append=FALSE)
    #st_write(tree, paste0("D:\\PROJECT_TREES\\isotree_trail\\CROWNS",'_',h,'PLUS_SILVA.shp' ), append=FALSE)    
    st_write(st_as_sf(ttop), paste0("D:\\PROJECT_TREES\\MAUNGATAUTARI\\R_OUTPUT\\TREETOPS\\TTOPS",'_',h*100,'CM_PLUS_SILVA.shp' ), append=FALSE)
    #st_write(st_as_sf(ttop), paste0("D:\\PROJECT_TREES\\isotree_trail\\TTOPS",'_',h,'PLUS_SILVA.shp' ), append=FALSE)

    
    i <- i+1
  }
  
  write.csv(df, "D:\\PROJECT_TREES\\MAUNGATAUTARI\\R_OUTPUT\\FOREST_METRICS.csv")
  #write.csv(df, "D:\\PROJECT_TREES\\isotree_trail\\FOREST_METRICS.csv")
  