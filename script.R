library("lidR")
library("future")

col <- height.colors(50)
ori_las <- readLAS("D:\\PROJECT_TREES\\CL2_BE33_2021_1000_1019.las", filter="-drop_classification 6 7 18")

n_las <- normalize_height(ori_las, knnidw())

getWS <- function(x){
  y <- 2.6125 + x^0.665
  y[x<2] <- 2
  y[x>40] <- 15
  return(y)
}

chm <- rasterize_canopy(n_las, 0.3, pitfree(subcircle = 0.3))

# to smooth the chm
fill.na <- function(x, i=5) { if (is.na(x)[i]) { return(mean(x, na.rm = TRUE)) } else { return(x[i]) }}
w <- matrix(1, 3, 3)

smooth_chm <- terra::focal(chm, w, fun = mean, na.rm = TRUE)

tree_tops <- locate_trees(n_las, lmf(getWS))

tree_dal <- segment_trees(n_las, dalponte2016(smooth_chm, tree_tops))
tree_li2012 <- segment_trees(n_las, li2012())

plot(tree_dal, bg = "white", size = 4, color = "treeID") # visualize trees
plot(tree_li2012, bg = "white", size = 4, color = "treeID") # visualize trees

#plot(smooth_chm, col = col)
#plot(sf::st_geometry(tree_tops), add = TRUE, pch = 3)