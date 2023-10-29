library("lidR")
library("sf")
library("raster")
remove(list = ls())
gc()

las_catg <- readLAScatalog("D:\\PROJECT_TREES\\MAUNGATAUTARI\\laser\\", filter="-drop_classification 6 7 18")
#las_catg <- readLAScatalog("D:\\PROJECT_TREES\\LAS\\", filter="-drop_classification 6 7 18")


hMin = 30
hMax = 30

opt_chunk_buffer(las_catg) <- 30
opt_chunk_size(las_catg)   <- 1000  

f_pointdensity <- function(chunk, hmin, hmax){
  las <- readLAS(chunk)
  if (is.empty(las)) return(NULL)
  
  las <- classify_noise(las, sor(15,7))
  las_denoise <- filter_poi(las, Classification != LASNOISE)
  nlas <- normalize_height(las_denoise, knnidw())
  filtered_las = filter_poi(nlas, Z>=hmin)
  pt_density = rasterize_density(filtered_las, 5)
  return(pt_density)
}

opt    <- list(need_buffer = TRUE, automerge = TRUE)   # catalog_apply will merge the outputs into a single object
pt_den <- catalog_apply(las_catg, f_pointdensity, hmin = hMin, hmax = hMax, .options = opt)  
writeRaster(pt_den, paste0("D:\\PROJECT_TREES\\MAUNGATAUTARI\\R_OUTPUT\\DENSITYMAP\\DENSITY",'_',hMin*100,'CM_PLUS.tif' ), overwrite=TRUE)

