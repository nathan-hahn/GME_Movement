#### Landscape Metrics ####

library(raster)
library(tidyverse)
#remotes::install_github("r-spatialecology/landscapemetrics")
library(landscapemetrics)
library(sf)
library(doParallel)

#' Generate Landscape Metrics 
#' 
#' Step 1: Import list of buffered steps for processing
#' Step 2: Parallel processing of steps in landscapemetrics package
#' Step 3: Export results as a dataframe - to be bound to tracking df later

# set wd - for ssh
setwd('~/Dropbox/CSU/GME_Movement')

##### Re-import data #####

# add split
split <- readRDS('./Staging Areas/movdata/ddm_buffer.RDS')

# add ev data
# Estes landcover
gse <- raster('./spatial data/change03_181_reclassMara_2019-11-22.tif')

# Reclassify to Natural-Ag cover
# 0 = nodata, 1 = ag, 2 = natural cover
m <- c(0, 0, 1, 1, 2, 2, 3, 2, 4, 0)
rclmat <- matrix(m, ncol=2, byrow=TRUE)
gse <- reclassify(gse, rcl = rclmat)


##### Extract in Parallel #####

# spin up cluster
cores <- 6
cl <- makeCluster(cores, output="") #output should make it spit errors
registerDoParallel(cl)

split <- split[91:119]


# ddm-buffer metrics - 5 hours for list of 30
system.time({
l <- length(split)
ddm.metrics <- foreach::foreach(i=1:l) %dopar% {
  metrics <- landscapemetrics::sample_lsm(gse, y = split[[i]], size = 5201, shape = 'square',
                             what = c('lsm_c_ed', 'lsm_c_pland', 'lsm_p_para', 'lsm_l_contag', 'lsm_l_lsi'),
                             plot_id = split[[i]]$uid)
}
})

stopCluster(cl)

t <- data.table::rbindlist(ddm.metrics)
write.csv(t, './Staging Areas/movdata/ddm_metrics_r4.csv')

