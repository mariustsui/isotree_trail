library("lidR")
library("future")
library("raster")
library("sf")
library("sp")
library("dplyr")
remove(list = ls())
gc()

# Set the minimum height of a tree to be detected ( in meters )

hmin <- c(5)
loop_index <- 0
forest_metrics <- matrix(nrow = 6, ncol=4)
colnames(forest_metrics) <- c("Number of Trees", "Trees per hectare", "Crown area(m2)", "Crown Cover")
rownames(forest_metrics) <- c("5M", "10M", "15M", "20M", "30M", "40M")

# read all the point cloud (.las) into a lasCatalog
  

  las_catg <- readLAScatalog("D:\\PROJECT_TREES\\REQ201107_Lars\\laser\\", filter="-drop_classification 6 7 18")
  #las_catg <- readLAScatalog("D:\\PROJECT_TREES\\LAS\\", filter="-drop_classification 6 7 18")
  
  
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
  
  f_treecrown <- function(las, hmin){
    
    las <- classify_noise(las, sor(15,7))
    las_denoise <- filter_poi(las, Classification != LASNOISE)
    nlas <- normalize_height(las_denoise, knnidw())
    tree_tops <- locate_trees(nlas, lmf(ws = getWS, hmin = hmin))
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
  for (h in hmin){
    loop_index <- loop_index + 1
    opt    <- list(automerge = TRUE)   # catalog_apply will merge the outputs into a single object
    output = NULL
    output <- catalog_map(las_catg, f_treecrown, hmin = h, .options = opt)
    
    single_sf <- dplyr::bind_rows(output)
    types <- st_geometry_type(single_sf)
    types_df <- data.frame(types)
    my_labelled_sf_object <- merge(single_sf, types_df,by.x=0, by.y=0, all.x=TRUE)
    my_filtered_sf_object <- my_labelled_sf_object[my_labelled_sf_object$types == "POLYGON",]
    
    
    st_write(my_filtered_sf_object, paste0("D:\\PROJECT_TREES\\R_OUTPUT\\CROWNS",'_',h,'PLUS_SILVA.shp' ), append=FALSE)
    #st_write(my_filtered_sf_object, paste0("D:\\PROJECT_TREES\\isotree_trail\\ALLCROWNS",'_',h,'PLUS_SILVA.shp' ), append=FALSE)
    
    
    a_ha <- units::set_units(st_area(las_catg), ha)
    a_m2 <- units::set_units(st_area(las_catg), m^2)
    tph <- nrow(my_filtered_sf_object) / (as.numeric(a_ha))
    cc <- sum(my_filtered_sf_object$convhull_area) / (as.numeric(a_m2))
    forest_metrics[loop_index, 1] <- nrow(my_filtered_sf_object)
    forest_metrics[loop_index, 2] <- round(tph, 2)
    forest_metrics[loop_index, 3] <- sum(my_filtered_sf_object$convhull_area)
    forest_metrics[loop_index, 4] <- round(cc, 2)
  }
  
  df <- as.data.frame((forest_metrics))
  
  
  tree_5 <- my_filtered_sf_object
  tree_10 <- tree_5[tree_5$Z >= 10,]
  st_write(tree_10, paste0("D:\\PROJECT_TREES\\R_OUTPUT\\CROWNS",'_10','PLUS_SILVA.shp' ), append=FALSE)
  df[2, ] <- c(nrow(tree_10), round(nrow(tree_10)/ as.numeric(a_ha), 2), round(sum(tree_10$convhull_area), 2), round( (sum(tree_10$convhull_area)/as.numeric(a_m2)*100), 2))
  tree_15 <- tree_5[tree_5$Z >= 15,]
  st_write(tree_15, paste0("D:\\PROJECT_TREES\\R_OUTPUT\\CROWNS",'_15','PLUS_SILVA.shp' ), append=FALSE)
  df[3, ] <- c(nrow(tree_15), round(nrow(tree_15)/ as.numeric(a_ha), 2), round(sum(tree_15$convhull_area), 2), round( (sum(tree_15$convhull_area)/as.numeric(a_m2)*100), 2))
  tree_20 <- tree_5[tree_5$Z >= 20,]
  st_write(tree_20, paste0("D:\\PROJECT_TREES\\R_OUTPUT\\CROWNS",'_20','PLUS_SILVA.shp' ), append=FALSE)
  df[4, ] <- c(nrow(tree_20), round(nrow(tree_20)/ as.numeric(a_ha), 2), round(sum(tree_20$convhull_area), 2), round( (sum(tree_20$convhull_area)/as.numeric(a_m2)*100), 2))
  tree_30 <- tree_5[tree_5$Z >= 30,]
  st_write(tree_30, paste0("D:\\PROJECT_TREES\\R_OUTPUT\\CROWNS",'_30','PLUS_SILVA.shp' ), append=FALSE)
  df[5, ] <- c(nrow(tree_30), round(nrow(tree_30)/ as.numeric(a_ha), 2), round(sum(tree_30$convhull_area), 2), round( (sum(tree_30$convhull_area)/as.numeric(a_m2)*100), 2))
  tree_40 <- tree_5[tree_5$Z >= 40,]
  st_write(tree_40, paste0("D:\\PROJECT_TREES\\R_OUTPUT\\CROWNS",'_40','PLUS_SILVA.shp' ), append=FALSE)
  df[6, ] <- c(nrow(tree_40), round(nrow(tree_40)/ as.numeric(a_ha), 2), round(sum(tree_40$convhull_area), 2), round( (sum(tree_40$convhull_area)/as.numeric(a_m2)*100), 2))
  
  write.csv(df, "D:\\PROJECT_TREES\\R_OUTPUT\\FOREST_METRICS.csv")
  #write.csv(df, "D:\\PROJECT_TREES\\isotree_trail\\FOREST_METRICS.csv")
  