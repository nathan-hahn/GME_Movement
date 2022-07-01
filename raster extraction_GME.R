#'----
#'title: Raster extraction
#'author: Nathan Hahn
#'date: 
#'----

# Set environment to EAT
Sys.setenv(TZ="Africa/Nairobi") 

library(dplyr)
library(lubridate)
library(sp)
library(raster)
library(terra)
library(parallel)

##########################################################################

####Prep data####

# load RData
df <- readRDS("./movdata/GMEcollars_004_clean_2021-10-04.rds")
df$date <- ymd_hms(df$date, tz = "Africa/Nairobi")
# summarize relocs by individual  
df %>%
  group_by(id) %>%
  summarize(n.reloc = n()) 

# remove rows with NA
df <- df[complete.cases(df[,1:2]),]

library(raster)
# system CRS 
study.area <- '+proj=utm +init=epsg:32736'

# get rasters from spatial data folder and set CRS
dist2ag <- rast("./spatial data/dist2ag_estes_32736_2019-11-21.tif")
dist2agedge <- rast("./spatial data/dist2agedge_estes_20200629.tif")
dist2permwater <- rast("./spatial data/dist2permanent_water_20200624.tif")
dist2seasonalwater <- rast("./spatial data/dist2seasonal_water_20200624.tif")
dist2water <- rast("./spatial data/dist2merged_water_20200624.tif")
dist2forest <- rast("./spatial data/dist2forest_hansen_cover60_32736_30.tif") 
slope <- rast("./spatial data/slope_estes_32736_2020-05-12.tif")
#lc <- raster("./spatial data/change03_181_reclassMara_2019-11-22.tif")
lc <- rast("./spatial data/Tiedman/sentinel2018-3yr-GSE-03-2022_32736_agmask.tif") #SME - march2022
gHM <- rast("./spatial data/gHM_estes_32736_2020-05-12.tif")
pa <- rast("./spatial data/GSEr_estes_32736_2020-10-20.tif") # protection status
pa <- classify(pa, cbind(NA, 1)) # set no data to 'not protected' (1)
pa <- classify(pa, cbind(0, 1))

# raster list
r.list <- list(dist2ag, dist2agedge, dist2permwater, dist2seasonalwater, dist2water, dist2forest, 
               gHM, slope, pa, lc)

# plot and check raster distributions
# 
# par(mfrow = c(4,3))
# for (i in 1:10){
#  plot(s[[i]], main = names(s[[i]]))
# }

# # histograms
# par(mfrow = c(3,2))
# for (i in 1:length(rasters)) {
#  hist(rasters[[i]], breaks = 100, main = names(rasters[[i]]))
# }

## Extract raster covariates

# create matrix for used points
# used <- matrix(1, nrow = nrow(df), ncol = 11)
# # create spatial points 
# locs <- SpatialPointsDataFrame(as.matrix(df[c("x","y")]), data = df, 
#                                proj4string = crs(study.area))

# # ~35 minutes
# system.time({
# used[,2] <- extract(dist2ag, locs)
# used[,3] <- extract(dist2agedge, locs) 
# used[,4] <- extract(dist2water, locs)
# used[,5] <- extract(dist2permwater, locs)
# used[,6] <- extract(dist2seasonalwater, locs)
# used[,7] <- extract(slope, locs)
# used[,8] <- extract(gHM, locs)
# used[,9] <- extract(pa, locs)
# used[,10] <- extract(lc, locs)
# used[,11] <- extract(dist2forest, locs)
# })
# 
# # Extract raster covariates (in parallel) - testing
# library(doParallel)
# cores<- 7
# cl <- makeCluster(cores, output="") #output should make it spit errors
# registerDoParallel(cl)
# 
# 
# # create spatial dataframe
# locs <- SpatialPointsDataFrame(as.matrix(df[c("x","y")]), data = df,
#                                 proj4string = crs(study.area))
# 
# # Extract - 45 minutes for 1.5 million relocs
# system.time({
# extract <- foreach::foreach(i = 1:10) %dopar% 
#   raster::extract(r.list[[i]], locs)
# })
# 
# used <- do.call(cbind, extract)
# 
# stopCluster(cl)

cores = 7
vect <- st_as_sf(df, coords = c('x','y'), crs = study.area) %>%
  terra::vect()
r.extract <- mclapply(r.list, function(x)
  terra::extract(x, vect)[,2], mc.cores = cores)
used <- do.call(cbind, r.extract) # bind results into an extracted dataframe for each individual

# check
head(used)
summary(used)

####Add to Tracking DF####

# create data frame
mode(used) = "numeric"
used2 <- as.data.frame(used)
used2$ID <- as.character(df$id)
colnames(used2) <- c("dist2ag", "dist2agedge", "dist2permwater", "dist2seasonalwater", "dist2water", 
                     "dist2forest", "gHM", "slope", 'pa', 'lc.estes', "merge_id")
head(used2)

# unstandardized data frame
used.df <- cbind(df, used2)
head(used.df)

test <- subset(used.df, id != merge_id)
nrow(test) # should be zero
used.df$merge_id <- NULL
used.df$used <- 1

##### Adjust ag edge #####
used.df$dist2agedge <- ifelse(used.df$lc.estes == 1, -(used.df$dist2agedge), used.df$dist2agedge)


##########################################################################################


####Export####

# write to rds file
outfile <- paste0("./movdata/GMEcollars_004_used_", Sys.Date(), ".rds")
saveRDS(used.df, outfile)

# write to csv file
outfile <- paste0("./movdata/GMEcollars_004_used_", Sys.Date(), ".csv")
write.csv(used.df, outfile)

########################################################################################################
########################################################################################################