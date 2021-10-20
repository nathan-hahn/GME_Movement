#### Landscape Metrics ####

library(raster)
library(tidyverse)
library(landscapemetrics)
library(sf)
library(doParallel)

#' Generate Landscape Metrics 
#' 
#' Step 1: Import staging area relocations and spatial data
#' Step 2: Calculate buffer distances for mean step and daily distance moved
#' Step 3: Filter spatial data rasters by staging relocations (buffer)
#' Step 4: Apply landscape metrics to filtered spatial data


# set wd - for ssh
setwd('~Dropbox/CSU/GME_Movement')

# Import spatial data
# Forest cover - hansen
forest <- raster('./spatial data/Forest_hansen/hansen_cover60_reclass_2019_2.tif')
m <- c(NA, 0, 1, 1)
rclmat <- matrix(m, ncol=2, byrow=TRUE)
forest <- reclassify(forest, rcl = rclmat)

# Estes landcover
gse <- raster('./spatial data/change03_181_reclassMara_2019-11-22.tif')

# Reclassify to Natural-Ag cover
# 0 = nodata, 1 = ag, 2 = natural cover
m <- c(0, 0, 1, 1, 2, 2, 3, 2, 4, 0)
rclmat <- matrix(m, ncol=2, byrow=TRUE)
gse <- reclassify(gse, rcl = rclmat)

# Import staging area relocations
df <- as.data.frame(data.table::fread('./Staging Areas/movdata/GMEcollars_004_stageclassified_20211019'))
df$V1 <- NULL

##### 2: Buffer Distances #####

# mean step length (by ag phase) - 312m
step <- df %>%
  group_by(ag.window) %>%
  summarise(mean.step = mean(dist))
step

# mean daily distance moved (by ag phase) - 5201m
ddm <- df %>%
  group_by(ag.window, dayBurst) %>%
  summarise(ddm = sum(dist)) %>%
  group_by(ag.window) %>%
  summarise(m.ddm = mean(ddm))
ddm

##### 3: Calculate Metrics in Buffers #####

# filter to relocs of interest
#df <- df[df$vote == 1,] 
#sample.df <- df[1:20000,] # for testing

# check points
#plot(gse)
#points(sample.df$x, sample.df$y)

# create spatial df
sf <- st_as_sf(sample.df, coords = c('x','y'), crs = 32736)

##### Extract in Parallel #####
# spin up cluster
cl <- makeCluster(7)
registerDoParallel(cl)

# split dataframe into equal chunks
d <- (dim(df)[1])/7
split <- split(sf, (seq(nrow(sf))-1) %/% d)

# # step-buffer metrics
# step.gse <- foreach(i=1:length(split), .combine=rbind) %do% {
#   step.metrics <- sample_lsm(gse, y = split.200[[i]], size = 312, shape = 'square',
#                              what = c('lsm_c_ed','lsm_p_para','lsm_l_contag'))
# }
# 
# step.gse

# ddm-buffer metrics
system.time({
ddm.metrics <- foreach(i=1:length(split), .combine=rbind) %do% {
  metrics <- sample_lsm(gse, y = split[[i]], size = 5201, shape = 'square',
                             what = c('lsm_c_ed','lsm_p_para','lsm_l_contag', 'lsm_l_lsi'))
}
})

write.csv(ddm.metrics, './Staging Areas/dd_metrics_20211020')

