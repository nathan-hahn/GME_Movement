#'----
#'title: Raster extraction
#'author: Nathan Hahn
#'date: 
#'----

# Set environment to EAT
Sys.setenv(TZ="Africa/Nairobi") 

library(dplyr)
library(sp)
library(raster)

##########################################################################

####Prep data####

# load RData
df <- readRDS("./movdata/GMEcollars_002_clean_2020-06-25.rds")

# summarize relocs by individual  
df %>%
  group_by(id) %>%
  summarize(n.reloc = n()) 

# remove rows with NA
df <- df[complete.cases(df[,1:2]),]

library(raster)
# system CRS - for some reason doesn't work if 32736???
study.area <- '+proj=utm +init=epsg:32736'

# get rasters from spatial data folder and set CRS
dist2ag <- raster("./spatial data/dist2ag_estes_32736_2019-11-21.tif")
dist2permwater <- raster("./spatial data/dist2permanent_water_20200624.tif")
dist2seasonalwater <- raster("./spatial data/dist2seasonal_water_20200624.tif")
dist2water <- raster("./spatial data/dist2merged_water_20200624.tif")
slope <- raster("./spatial data/slope_estes_32736_2020-05-12.tif")
lc <- raster("./spatial data/change03_181_reclassMara_2019-11-22.tif")
gHM <- raster("./spatial data/gHM_estes_32736_2020-05-12.tif")
pa <- raster("./spatial data/GSEr_estes_32736_2019-11-21.tif") # protection status
pa <- reclassify(pa, cbind(NA, 1)) # set no data to 'not protected' (1)

# downsampled sentinel classification to 30 meters to keep consistent with other proximity layers. extent different
dist2forest <- raster("./spatial data/dist2forest_tiedmen_sentinelpts_32736_30.tif") 

# plot and check raster distributions
rasters <- list(slope=slope, dist2ag=dist2ag, dist2forest = dist2forest,
               dist2water=dist2water, gHM = gHM, protection = pa, landcover = lc)


s <- stack(dist2ag, dist2water, slope, pa, lc)

# # rasters
# par(mfrow = c(3,2))
# for (i in 1:length(rasters)){
#  plot(rasters[[i]], main = names(rasters[[i]]))
# }
# 
# # histograms
# par(mfrow = c(3,2))
# for (i in 1:length(rasters)) {
#  hist(rasters[[i]], breaks = 100, main = names(rasters[[i]]))
# }

## Extract raster covariates

# create matrix for used points
used <- matrix(1, nrow = nrow(df), ncol = 10)
# create spatial points 
locs <- SpatialPointsDataFrame(as.matrix(df[c("x","y")]), data = df, 
                               proj4string = crs(study.area))

# ~20 minutes
system.time({
used[,2] <- extract(dist2ag, locs)
used[,3] <- extract(dist2water, locs)
used[,4] <- extract(dist2permwater, locs)
used[,5] <- extract(dist2seasonalwater, locs)
used[,6] <- extract(slope, locs)
used[,7] <- extract(gHM, locs)
used[,8] <- extract(pa, locs)
used[,9] <- extract(lc, locs)
used[,10] <- extract(dist2forest, locs)
})

# check
head(used)
summary(used)


## Extract raster covariates (in parallel) - testing
# library(snow)
# 
# locs2 <- SpatialPointsDataFrame(as.matrix(df[c("x","y")]), data = df,
#                                 proj4string = crs(study.area))
# # Extract
# nCores <- detectCores() - 1
# beginCluster(n = nCores)
# system.time(
#   used2 <- extract(s, locs2)
# )
# endCluster()


####Add to Tracking DF####

# create data frame
mode(used) = "numeric"
used2 <- as.data.frame(used)
used2$ID <- as.character(locs@data$id)
colnames(used2) <- c("used", "dist2ag", "dist2water", "slope", "gHM", "pa", "lc.estes", "dist2forest", "merge_id")
head(used2)

# unstandardized data frame
used.df <- cbind(df, used2)
head(used.df)


test <- subset(used.df, id != merge_id)
nrow(test) # should be zero
used.df$merge_id <- NULL

##########################################################################################

####Add season variable####
#' Must extend season clustering to cover full extent of tracking data
#' Assign season (wet or dry) based on season window dates indentified by seasonal NDVI mixture model. Season is added as a column to the main ele mov dataframe with covariates
#' Wet = 1
#' Dry = 2 

# use csv with all season windows combined
windows <- read.csv("~/Dropbox (Personal)/CSU/MEP/Data/spatial data/season_windows.csv")
windows$start <- as.POSIXct(strptime(windows$start, format = "%Y-%m-%d %H:%M:%S", tz="Africa/Nairobi" ))
windows$end <- as.POSIXct(strptime(windows$end, format = "%Y-%m-%d %H:%M:%S", tz="Africa/Nairobi" ))

# cut, setting breaks as the start dates and the max end date for season
rng <- cut(used.df$date,
           breaks = c(windows$start, max(windows$end)),
           include.lowest = T)

# add rng to ele dates
test <- cbind(used.df$Fixtime, rng)
head(test)

# convert to data.frame. rng window values show up as the factor level 1 - 47
# Assign season value based on factor level. 1-23 are wet windows (1), 24-47 are dry windows (2)
test <- as.data.frame(rng) %>%
  mutate(rng = as.numeric(rng)) %>%
  mutate(season = ifelse(rng < 24, 1, 2)) %>%
  dplyr::select(season)

# attach to ele mov dataframe
used.df <- bind_cols(used.df, test) 


####Export####

# write to rds file
outfile <- paste0("./movdata/GMEcollars_002_used_", Sys.Date(), ".rds")
saveRDS(used.df, outfile)

outfile <- paste0("./movdata/GMEcollars_002_used_", Sys.Date(), ".csv")
write.csv(used.df, outfile)

# write raster stack
saveRDS(s, "./spatial data/GME_rasterStack_20200602.rds")

########################################################################################################
########################################################################################################