#### Landscape Metrics ####

library(raster)
library(tidyverse)
#remotes::install_github("r-spatialecology/landscapemetrics")
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
setwd('~/Dropbox/CSU/GME_Movement')

# Import staging area relocations
df <- as.data.frame(data.table::fread('./Staging Areas/movdata/GMEcollars_004_stageclassified_20211019'))
df$V1 <- NULL

##### 2: Buffer Distances #####

# mean step length (by ag phase) - 312m
step <- df %>%
  group_by(ag.window) %>%
  summarise(mean.step = mean(dist))
step

step <- 312

# mean daily distance moved (by ag phase) - 5201m
ddm <- df %>%
  group_by(ag.window, dayBurst) %>%
  summarise(ddm = sum(dist)) %>%
  group_by(ag.window) %>%
  summarise(m.ddm = mean(ddm))
ddm

ddm <- 5201

##### 3: Generate buffers #####

# filter to relocs of interest
#df <- df[df$vote == 1,] 
#sample.df <- df[1:20000,] # for testing

# check points
#plot(gse)
#points(sample.df$x, sample.df$y)

# create spatial df
sf <- st_as_sf(df, coords = c('x','y'), crs = 32736)

# create buffers
sf.ddm <- st_buffer(sf, dist = ddm)

# split dataframe into equal chunks
#d <- (dim(sf.ddm)[1])/10
split <- split(sf.ddm, (seq(nrow(sf.ddm))-1) %/% 10000)

# save as rds file
saveRDS(split, './Staging Areas/movdata/ddm_buffer.RDS')

# clear environment to free memory - if using Rstudio, also need to restart r session
rm(list = ls())

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
cores <- 2
cl <- makeCluster(cores, output="") #output should make it spit errors
registerDoParallel(cl)

# # step-buffer metrics
# step.metrics <- foreach(i=1:length(split), .combine=rbind) %do% {
#   step.metrics <- sample_lsm(gse, y = split.200[[i]], size = 312, shape = 'square',
#                              what = c('lsm_c_ed','lsm_p_para','lsm_l_contag'))
# }

# write.csv(step.metrics, './Staging Areas/step_metrics_20211020')

# ddm-buffer metrics
system.time({
ddm.metrics <- foreach(i=1:length(split)) %dopar% {
  metrics <- landscapemetrics::sample_lsm(gse, y = split[[i]], size = 5201, shape = 'square',
                             what = c('lsm_c_ed', 'lsm_c_pland', 'lsm_p_para', 'lsm_l_contag', 'lsm_l_lsi'))
}
})

stopCluster(cl)

write.csv(ddm.metrics, './Staging Areas/dd_metrics_20211020')

