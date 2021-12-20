#### Logistic Regression Test - Staging Areas ####
library(tidyverse)
library(lme4)
library(terra)
library(sf)
library(MuMIn)
library(mapview)
library(tmap)

Sys.setenv(TZ='Africa/Nairobi')

###########################################################################################################
## ///////////////////////////////////////////// Data Prep ///////////////////////////////////////////// ##
###########################################################################################################

#### Load Data ####
df <- as.data.frame(data.table::fread('./movdata/GMEcollars_004_stageclassified_20211019.csv', tz = ''))
df$V1 <- NULL
df$V1 <- NULL
df$...1 <- NULL
df$used <- NULL

df$vote <- as.factor(df$vote)
table(df$vote)

# filter to relocs of interest -- fixes during ag use periods
df <- df %>%
  filter(ag.window.ext == 1)

#### Extract Covariates ####
# covs of interest are prop.forest, prop.ag, drainage (1/0), dist2water, gHM, slope

## Remove existing covariates
df <- df %>%
  select(-c('dist2ag', 'dist2agedge', 'dist2permwater', 'dist2seasonalwater', 'dist2water', 'dist2forest', 'gHM', 'slope', 'pa'))

##### Import covariate layers #####

## Raster layer prep
prop.ag.250 <- rast("./spatial data/estes_ag_pct_250.tif")
prop.ag.1500 <- rast("./spatial data/estes_ag_pct_1500.tif")
prop.forest.250 <- rast("./spatial data/hansen_forest_pct_250.tif")
prop.forest.1500 <- rast("./spatial data/hansen_forest_pct_1500.tif")
dist2ag <- rast("./spatial data/dist2ag_estes_32736_2019-11-21.tif")
dist2agedge <- rast("./spatial data/dist2agedge_estes_20200629.tif")
dist2water <- rast("./spatial data/dist2merged_water_20200624.tif")
dist2forest <- rast("./spatial data/dist2forest_hansen_cover60_32736_30.tif") 
slope <- rast("./spatial data/slope_estes_32736_2020-05-12.tif")
lc <- rast("./spatial data/change03_181_reclassMara_2019-11-22.tif")
gHM <- rast("./spatial data/gHM_estes_32736_2020-05-12.tif")
prop.settlement.250 <- rast("spatial data/estes_settlement_pct_250.tif")
prop.settlement.1500 <- rast("spatial data/estes_settlement_pct_1500.tif")
dist2paedge <- rast("./spatial data/dist2paedge_estes_32736_20211118.tif")

pa <- rast("./spatial data/GSEr_dissolved_estes_32736_20211118.tif") # binary protection status dissolved (0/1)
pa <- terra::classify(pa, cbind(NA, 0)) # set no data to 'not protected' (0)
plot(pa)

## Polygon layer prep
# get shapefile for drains buffer
drains <- st_read("./spatial data/drains/drains_estes_20211117/drains_estes_-2021-11-17.shp", 
                  layer="drains_estes_-2021-11-17", crs = 4326) %>%
  st_transform(crs = 32736) %>%
  filter(RIV_ORD <= 7)
  
# buffer by 1000m
drains.1000 <- st_buffer(drains, dist = 1000) %>%
  st_union() 
plot(drains.1000)

##### Extract #####

# create spatial dataframe
study.area <- 32736
locs.sf <- st_as_sf(df, coords = c('x','y'), crs = study.area)
locs <- locs.sf %>% terra::vect()

# create matrix for used points
used <- matrix(1, nrow = nrow(df), ncol = 17)

# ~60 seconds
system.time({
used[,2] <- terra::extract(dist2ag, locs)[,2]
used[,3] <- terra::extract(dist2agedge, locs)[,2]
used[,4] <- terra::extract(dist2water, locs)[,2]
used[,5] <- as.numeric(st_intersects(locs.sf, drains.1000)) # point to polygon
used[,5][is.na(used[,5])] <- 0 # convert NAs to 0's
used[,6] <- terra::extract(slope, locs)[,2]
used[,7] <- terra::extract(gHM, locs)[,2]
used[,8] <- terra::extract(pa, locs)[,2]
used[,9] <- terra::extract(lc, locs)[,2]
used[,10] <- terra::extract(dist2forest, locs)[,2]
used[,11] <- terra::extract(prop.ag.250, locs)[,2]
used[,12] <- terra::extract(prop.ag.1500, locs)[,2]
used[,13] <- terra::extract(prop.forest.250, locs)[,2]
used[,14] <- terra::extract(prop.forest.1500, locs)[,2]
used[,15] <- terra::extract(dist2paedge, locs)[,2]
used[,16] <- terra::extract(prop.settlement.250, locs)[,2]
used[,17] <- terra::extract(prop.settlement.1500, locs)[,2]
})


# check
head(used)
summary(used)

# create data frame
mode(used) = "numeric"
used2 <- as.data.frame(used)
used2$uid <- as.numeric(locs.sf$uid)
colnames(used2) <- c("used","dist2ag", "dist2agedge", "dist2water", 'drains1000', 'slope', 'gHM', 'pa', 'lc', 'dist2forest', 
                     'prop.ag.250', 'prop.ag.1500', 'prop.forest.250', 'prop.forest.1500', 'dist2paedge', 'prop.settlement.250', 
                     'prop.settlement.1500', 'merge_id')
head(used2)

# check for duplicates
length(unique(used2$merge_id)) == nrow(used2)

# unstandardized data frame
used.df <- cbind(df, used2)
head(used.df)

test <- subset(used.df, uid != merge_id)
nrow(test) # should be zero
used.df$merge_id <- NULL

# check for duplicates
length(unique(used.df$uid)) == nrow(used.df)

##### Adjust ag edge #####
used.df$dist2agedge <- ifelse(used.df$lc.estes == 1, -(used.df$dist2agedge), used.df$dist2agedge)
used.df$dist2paedge <- ifelse(used.df$pa == 1, -(used.df$dist2paedge), used.df$dist2paedge)

summary(used.df$dist2agedge)
summary(used.df$dist2paedge)

# check for duplicates
length(unique(used.df$uid)) == nrow(used.df)

##### Add Forest Proximity Index #####
## Pre-computed 

# read in results - tagged by uid
# prox.index.250 <- read.csv('spatial data/Forest_hansen/distanceloop/output/proximity_index_table_250.csv') %>%
#   select(prox.index.250 = prox.index, prox.index.mean, uid)
# 
# used.df <- merge(used.df, prox.index.250, 'uid', all.x = TRUE)
# summary(used.df$prox.index.250)


##### Save output #####
write.csv(used.df, './Staging Areas/movdata/movdat_004_lsdv.csv')


###########################################################################################################
## ///////////////////////////////////////// Model Fitting ///////////////////////////////////////////// ##
###########################################################################################################

#### Model Fitting ####

##### Read in data #####
used.df <- as.data.frame(data.table::fread('./Staging Areas/movdata/movdat_004_lsdv.csv'))
used.df$V1 <- NULL

## **remove relocations within ag**
ds.st.sub <- used.df %>%
  filter(lubridate::hour(date) %in% c(6:18)) #%>%
  #filter(dist2ag > 0)

# check for duplicates
length(unique(ds.st.sub$uid)) == nrow(ds.st.sub)

## Load proximity index 
prox.index <- read.csv('spatial data/Forest_hansen/distanceloop/output/proximity_index_table_250.csv')
ds.st.sub <- merge(ds.st.sub, prox.index, c('uid', 'x', 'y'), all.x = T)

# update with search radius
names(ds.st.sub)[names(ds.st.sub) == 'prox.index'] <- 'prox.index.250'

## adjust proximity index values
ds.st.sub$prox.index.250 <- ifelse(is.na(ds.st.sub$prox.index.250), 0.00001, ds.st.sub$prox.index.250)

# ds.st.sub$prox.norm = (ds.st.sub$prox.index.250 - min(ds.st.sub$prox.index.250)) / (max(ds.st.sub$prox.index.250) - min(ds.st.sub$prox.index.250))
# 
# ds.st.sub %>% group_by(vote) %>% summarise(mean = mean(prox.index.250),
#                                            mean.norm = mean(prox.norm),
#                                            mean.scale = mean(scale(prox.index.250)),
#                                            mean.log = mean(log(prox.index.250)))

## standardize covs
covariates <- c("dist2ag", "dist2agedge", "dist2water", 'slope', 'dist2forest', 'dist2paedge')
ds.st.sub <- ds.st.sub %>%
  dplyr::select(uid, subject_name, burst, x, y, date, vote, all_of(covariates), drains1000, gHM, prop.ag.250, prop.ag.1500,
                prop.forest.250, prop.forest.1500, prop.settlement.250, prop.settlement.1500, prox.index.250, pa) %>%
  mutate(forest = if_else(dist2forest == 0, 1, 0)) %>%
  mutate_at(covariates, .funs = scale) %>%
  mutate(drains1000 = as.factor(drains1000)) %>%
  mutate(pa = as.factor(pa)) %>%
  #mutate(prox.index.250 = log(prox.index.250)) %>%
  droplevels()


##### test autocorrelation variogram #####
library(gstat)
t <- st_as_sf(ds.st.sub, coords = c('x','y'), crs = 32736)
t <- as_Spatial(t)

t.var<-gstat::variogram(vote~1, t)
plot(t.var)

##### test spatial scales for moving window metrics #####
mod.ag.step <- glmer(vote ~ prop.ag.250 + (1|subject_name),
                data = ds.st.sub, family = binomial)
mod.ag.daily <- glmer(vote ~ prop.ag.1500 + (1|subject_name),
                      data = ds.st.sub, family = binomial)
AICc(mod.ag.step, mod.ag.daily) # ag daily scale is better by AICc, ag step model is rank deficient and ag drops from model

mod.for.step <- glmer(vote ~ prop.forest.250 + (1|subject_name),
                 data = ds.st.sub, family = binomial)
mod.for.daily <- glmer(vote ~ prop.forest.1500 + (1|subject_name),
                       data = ds.st.sub, family = binomial)
mod.for.dist <- glmer(vote ~ dist2forest + (1|subject_name),
                      data = ds.st.sub, family = binomial)

AICc(mod.for.step, mod.for.daily, mod.for.dist) # forest step scale is better by AICc

# check proximity index
mod.prox.250 <- glmer(vote ~ log(prox.index.250) + (1|subject_name),
                      data = ds.st.sub, family = binomial)


## correlation between forest prop and proximity index
t <- lm(log(prox.index.250) ~ prop.forest.250, data = ds.st.sub)

##### Fit candidate models #####

## Global model
mod.global.sub <- glmer(vote ~ prop.ag.1500 + prop.forest.250 + drains1000 + slope + gHM + pa + (1|subject_name), 
                    data = ds.st.sub, family = binomial)

## Human footprint model
mod.hm.sub <- glmer(vote ~ prop.ag.1500 + gHM + pa + (1|subject_name),
                data = ds.st.sub, family = binomial)

## Natural features model
mod.nat.sub <- glmer(vote ~ prop.forest.250 + drains1000 + slope + (1|subject_name),
                 data = ds.st.sub, family = binomial)

## Strongest predictors
mod.str.sub <- glmer(vote ~ prop.ag.1500 + prop.forest.250 + gHM + (1|subject_name),
                     data = ds.st.sub, family = binomial)

mod.str2.sub <- glmer(vote ~ prop.ag.1500 + prop.forest.250 + (1|subject_name),
                     data = ds.st.sub, family = binomial)

mod.str3.sub <- glmer(vote ~ gHM + (1|subject_name),
                                    data = ds.st.sub, family = binomial)

mod.str4.sub <- glmer(vote ~ gHM*prop.forest.250 + (1|subject_name),
                      data = ds.st.sub, family = binomial)

mod.global.interaction <- glmer(vote ~ prop.ag.1500 + gHM*prop.forest.250 + log(prox.index.250) + drains1000 + slope + pa + (1|subject_name), 
                          data = ds.st.sub, family = binomial)

AICc(mod.global.sub, mod.hm.sub, mod.nat.sub, mod.str.sub, mod.str2.sub, mod.str3.sub, mod.str4.sub, mod.global.interaction)

summary(mod.global.interaction)
sjPlot::plot_model(mod.global.interaction)

#### Grid-Based Model Fitting ####

##### Create Grid #####
## create grid -- 250m with extent based on gHM Estes (1000m res)
# use model dataframe to create terra vect object
locs <- used.df %>%
  filter(lc.estes != 1) %>%
  st_as_sf(coords = c('x','y'), crs = 32736) %>%
  terra::vect()
# split into stage/no.stage dataframes
locs.stage <- locs[locs$vote == 1,]
locs.nostage <- locs[locs$vote == 0,]

# reference raster - 250m or 1000m
r <- terra::rast(locs, extent = ext(gHM), resolution = 250)

r.stage <- terra::rasterize(locs.stage, r, fun = 'sum')
r.nostage <- terra::rasterize(locs.nostage, r, fun = 'sum')

stack <- raster::stack(raster::raster(r.stage), raster::raster(r.nostage))

prop.dens <- function(x,y){
  x = ifelse(is.na(x) & !is.na(y), 0, x) # replace NAs for staging raster with 0 if y has a value
  sum = x + y
  prop.dens = x/sum
  return(prop.dens)
}

rel.staging.occ <- raster::overlay(stack, fun = prop.dens)
rel.staging.occ <- rast(rel.staging.occ) # keep everything in terra
plot(rel.staging.occ)


weights <- function(x,y){
  x = ifelse(is.na(x) & !is.na(y), 0, x) # replace NAs for staging raster with 0 if y has a value
  sum = x + y
  return(sum)
}

occ.weights <- raster::overlay(stack, fun = weights)
occ.weights <- rast(occ.weights)
plot(occ.weights)

# check results
gme <- sf::st_read('~/Dropbox (Personal)/CSU/GME_Movement/spatial data/GSE/GSE_2020.shp')
YlOrRd <- RColorBrewer::brewer.pal(9, 'YlOrRd')
library(mapview)
mapview(gme, zcol = 'pa_status', col.regions = RColorBrewer::brewer.pal(3, 'Greens')) + 
  mapview(raster::raster(r.nostage), col.regions = YlOrRd, layer.name = 'non-staging relocations count') + 
  mapview(raster::raster(r.stage), col.regions = YlOrRd, layer.name = 'staging relocations count') + 
  mapview(raster::raster(rel.staging.occ), col.regions = YlOrRd, layer.name = 'relative density of staging')

# check with staging relocs
mapview(rel.staging.occ, col.regions = YlOrRd, layer.name = 'relative density of staging', alpha = 0.6) 
  #mapview(st_as_sf(used.df[used.df$vote == 1,], coords = c('x','y'), crs = 32736), col.region = 'red', cex = 1.5, layer.name = 'staging relocs')

## save the output -- check resolution size
terra::writeRaster(rel.staging.occ, './Staging Areas/staging_reldensity_250m.tif', overwrite = T) 

## relative density of staging to non-staging
library(tmap)
gme <- sf::st_read('~/Dropbox (Personal)/CSU/GME_Movement/spatial data/GSE/GSE_2020.shp')
dens.r <- raster::raster('./Staging Areas/staging_reldensity_250m.tif')
rel.staging.occ <- tm_shape(gme) + tm_polygons(col = 'pa_status', palette = RColorBrewer::brewer.pal(3, 'Greens'), title = 'serengeti-mara') +
  tm_shape(dens.r, raster.downsample = FALSE) + tm_raster(palette = RColorBrewer::brewer.pal(5, 'YlOrRd'), title = 'Relative Staging Occurance')
rel.staging.occ

tmap_save(rel.staging.occ, "Rel Staging Map.png", dpi = 300)



##### Extract Relative Staging Occurrence #####
#' Extract staging occurrence cell values for each relocation
#' Filter to relocations of interest (staging relocs and non-staging relocs within ag days and not in ag)

# create spatial dataframe
study.area <- 32736
locs.sf <- st_as_sf(ds.st.sub, coords = c('x','y'), crs = study.area)
locs <- locs.sf %>% terra::vect()

# extract to new column
#rel.staging.occ <- terra::rast('./Staging Areas/staging_reldensity_1000m.tif')
#occ.weights <- terra::rast('...')
ds.st.sub$rel.staging.occ <- terra::extract(rel.staging.occ, locs)[,2]
ds.st.sub$occ.weights <- terra::extract(occ.weights, locs)[,2]
ds.st.sub <- ds.st.sub[!is.na(ds.st.sub$rel.staging.occ),]
summary(ds.st.sub$rel.staging.occ)
summary(ds.st.sub$occ.weights)

#ds.st.sub$occ.weights <- (ds.st.sub$occ.weights-min(ds.st.sub$occ.weights))/(max(ds.st.sub$occ.weights)-min(ds.st.sub$occ.weights))

hist(ds.st.sub$rel.staging.occ)

##### test autocorrelation variogram #####
library(gstat)
t <- st_as_sf(ds.st.sub, coords = c('x','y'), crs = 32736)
t <- as_Spatial(t)

hscat(rel.staging.occ~1, t, 1:9*10000)


t.var<-gstat::variogram(rel.staging.occ~1, t, cloud = T)
plot(t.var)

TheVariogramModel <- vgm(psill=50, model="Exp", nugget=0.02, range=100000)
FittedModel <- fit.variogram(t.var, model=TheVariogramModel)
plot(t.var, model=FittedModel)


##### Fit Regression Model #####

mod.global.sub.1000 <- glmer(rel.staging.occ ~ prop.ag.1500 + prop.forest.250 + drains1000 + slope + gHM + pa + (1|subject_name), 
                            weights = log(occ.weights), # weights supplied as total counts that proportions arise from
                            data = ds.st.sub,
                            family = binomial)
summary(mod.global.sub.1000)

mod.global.sub.250.int <- glmer(rel.staging.occ ~ prop.ag.1500 + gHM*prop.forest.250 + drains1000 + slope + pa + (1|subject_name), 
                             weights = log(occ.weights), # weights supplied as total counts that proportions arise from
                             data = ds.st.sub,
                             family = binomial)
summary(mod.global.sub.250.int)


sjPlot::plot_models(mod.global.sub.1000, mod.global.sub.250)

plot(mod.global.sub.1000)
plot(mod.global.sub.1000.xy)


