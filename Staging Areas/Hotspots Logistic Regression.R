##### Import covariate layers #####

library(terra)
library(tidyverse)
library(sf)
library(lme4)
library(MuMIn)
library(spdep)
library(data.table)

## Raster layer prep
prop.ag.250 <- rast("./spatial data/estes_ag_pct_250.tif")
prop.ag.1500 <- rast("./spatial data/estes_ag_pct_1500.tif")
prop.forest.250 <- rast("./spatial data/hansen_forest_pct_250.tif")
prop.forest.1500 <- rast("./spatial data/hansen_forest_pct_1500.tif")
dist2ag <- rast("./spatial data/dist2ag_estes_32736_2022-01-03.tif")
dist2agedge <- rast("./spatial data/dist2agedge_estes_20200629.tif")
dist2water <- rast("./spatial data/dist2merged_water_20200624.tif")
dist2forest <- rast("./spatial data/dist2forest_hansen_cover60_32736_30.tif") 
slope <- rast("./spatial data/slope_estes_32736_2020-05-12.tif")
lc <- rast("./spatial data/change03_181_reclassMara_2019-11-22.tif")
gHM <- rast("./spatial data/gHM_estes_32736_2020-05-12.tif")
prop.settlement.250 <- rast("spatial data/estes_settlement_pct_250.tif")
prop.settlement.1500 <- rast("spatial data/estes_settlement_pct_1500.tif")
dist2paedge <- rast("./spatial data/dist2paedge_estes_32736_20211118.tif")
ndviCOV <- rast("./spatial data/NDVICoV_estes_32736-2022-02-03.tif")

pa <- rast("./spatial data/GSEr_dissolved_estes_32736_20220120.tif") # 3 - PA, 2 - cc
pa <- terra::classify(pa, cbind(NA, 0)) # set no data to 'not protected' (0)
plot(pa)

rel.staging.occ <- rast("./Staging Areas/rel.staging.occ_20220128.tif")
weights.occ <- rast("./Staging Areas/occ.weights_20220128.tif")

## Polygon layer prep
# get shapefile for drains buffer
drains <- st_read("./spatial data/drains/drains_estes_20211117/drains_estes_-2021-11-17.shp", 
                  layer="drains_estes_-2021-11-17", crs = 4326) %>%
  st_transform(crs = 32736) %>%
  filter(RIV_ORD <= 7)

# buffer by 250m
drains <- st_buffer(drains, dist = 250) %>%
  st_union() %>%
  terra::vect()


## Resample
stack.30 <- list(prop.ag.250, prop.ag.1500, prop.forest.250, prop.forest.1500, slope, dist2paedge)
stack.250 <- lapply(stack.30, terra::aggregate, fact = 8)
stack.250[[7]] <- terra::aggregate(pa, fact = 8, fun = 'min')

# NDVI predictability already at 250m res
stack.250[[8]] <- ndviCOV

# gHM at 1000m res
stack.250[[9]] <- terra::disagg(gHM, fact = 4)

# drains rasterized at 250m res
drains <- terra::rasterize(drains, stack.250[[1]])
stack.250[[10]] <- drains

names(stack.250) <- c('prop.ag.250','prop.ag.1500','prop.forest.250','prop.forest.1500','slope','dist2paedge','pa', 'ndviCOV','gHM','drains')


##### Make Staging Hotspots Grid #####

## Import Movement Data 
used.df <- as.data.frame(data.table::fread('./Staging Areas/movdata/movdat_004_lsdv_20220109.csv'))
used.df$V1 <- NULL


##### Create Grid #####
## create grid -- 250m with extent based on gHM Estes (1000m res)
# use model dataframe to create terra vect object
locs <- used.df %>%
  filter(ag.window.ext == 1) %>%
  filter(lc.estes != 1) %>% 
  st_as_sf(coords = c('x','y'), crs = 32736) %>%
  terra::vect()
# split into stage/no.stage dataframes
locs.stage <- locs[locs$vote.ag == 1,]
locs.nostage <- locs[locs$vote.ag == 0,]

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

# writeRaster(rel.staging.occ, './Staging Areas/rel.staging.occ_hmm_20220128.tif')
# writeRaster(occ.weights, './Staging Areas/occ.weights_hmm_20220128.tif')

##### Make Individual Movement Grids #####

## **remove relocations within ag** 
stagelocs <- used.df %>%
  filter(ag.window.ext == 1) %>%
  filter(lc.estes != 1)

# create spatial points in terra vec format
study.area <- 32736
locs.sf <- st_as_sf(stagelocs, coords = c('x','y'), crs = study.area)
locs <- locs.sf %>% terra::vect()
# create locs list - split by individual
locs <- split(locs, stagelocs$subject_name)

# create base raster - based on the extent and res of the staging hotspot raster 
r.base <- rast(xmin = 480885, xmax = 864635, ymin = 9583755, ymax = 9938255, resolution = 250)

# create raster with counts - each cell value is number of relocations within it for the individual 
r.count <- lapply(locs, rasterize, y = r.base, fun = length)
names(r.count) <- sort(unique(stagelocs$subject_name))

# reclassify to 0/1 used/unused pixels
r.used <- lapply(r.count, terra::classify, cbind(0, Inf, 1))

# convert individual movement grids - generates an xy coordinate for every used cell for each individual
r.points <- lapply(r.used, FUN = function (x) {
  r <- raster::raster(x)
  points <- as.data.frame(raster::rasterToPoints(r))
  colnames(points) <- c('x','y','used')
  return(points)
})
# add subject names as list names
names(r.points) <- sort(unique(stagelocs$subject_name))

##### Extract covariates #####
## Covariates are extracted using cell centers of the r.base raster

# get xy coordinates for r.base
r.base.points <- raster::rasterToPoints(raster::raster(r.base))
r.base.locs <- as.data.frame(r.base.points) %>%
  st_as_sf(coords = c('x','y'), crs = study.area) %>%
  terra::vect()

# create matrix for used points
used <- matrix(1, nrow = nrow(r.base.locs), ncol = 12)

# ~60 seconds
system.time({
  used[,1] <- terra::extract(rel.staging.occ, r.base.locs)[,2]
  used[,2] <- terra::extract(occ.weights, r.base.locs)[,2]
  used[,3] <- terra::extract(stack.250$prop.ag.250, r.base.locs) [,2]
  used[,4] <- terra::extract(stack.250$prop.ag.1500, r.base.locs)[,2]
  used[,5] <- terra::extract(stack.250$prop.forest.250, r.base.locs)[,2]
  used[,6] <- terra::extract(stack.250$prop.forest.1500, r.base.locs)[,2]
  used[,7] <- terra::extract(stack.250$slope, r.base.locs)[,2]
  used[,8] <- terra::extract(stack.250$dist2paedge, r.base.locs)[,2]
  used[,9] <- terra::extract(stack.250$pa, r.base.locs)[,2]
  used[,10] <- terra::extract(stack.250$ndviCOV, r.base.locs)[,2]
  used[,11] <- terra::extract(stack.250$gHM, r.base.locs)[,2]
  used[,12] <- terra::extract(stack.250$drains, r.base.locs)[,2]

})


# check
head(used)
summary(used)

# create data frame
mode(used) = "numeric"
used2 <- as.data.frame(used)
colnames(used2) <- c('rel.staging.occ','occ.weights','prop.ag.250','prop.ag.1500','prop.forest.250','prop.forest.1500','slope','dist2paedge','pa','ndviCOV','gHM','drains250')
head(used2)

# unstandardized data frame
used.df <- cbind(r.base.points, used2)
head(used.df)

# match covariate values with individual grids - adds covariate values to each individual grid table
match.df <- lapply(r.points, function (x) merge(x, used.df, by = c('x','y'), all.x = T))

# add subject_name covariate to each cell
names(match.df) <- sort(unique(stagelocs$subject_name))
match.df <- Map(cbind, match.df, subject_name = names(match.df))

# create final modeling dataframe
hotspots <- do.call(rbind, match.df)
dim(hotspots)
head(hotspots)
rownames(hotspots) <- NULL

summary(hotspots)

#write.csv(hotspots, './Staging Areas/movdata/staging_hotspots_datatable.csv')
hotspots <- read.csv('./Staging Areas/movdata/staging_hotspots_datatable.csv')

## Clean modeling dataframe

# remove cell positions that have no staging (cells in ag)
hotspots <- filter(hotspots, !is.na(occ.weights))

# fix drains
hotspots$drains250 <- ifelse(is.na(hotspots$drains250), 0, hotspots$drains250)

# add sex
sex <- stagelocs %>% group_by(subject_name, subject_sex) %>% tally() %>% select(-n)
hotspots <- merge(hotspots, sex, by = 'subject_name', all.x = T)

## standardize covs
covariates <- c('slope', 'dist2paedge', 'ndviCOV')
hotspots <- hotspots %>%
  mutate_at(covariates, .funs = scale) %>%
  mutate(drains250 = as.factor(drains250)) %>% #update to drains250 later!!
  mutate(pa = as.factor(pa)) %>%
  droplevels()


##### Calculate Autocovariate - by individual #####

library(sp)
library(spdep)

## Create sp dataframe
# dividing utm's by 1000 speeds up AC calculations
hotspots$x.1 <- hotspots$x/1000
hotspots$y.1 <- hotspots$y/1000
dataSp <- hotspots
coordinates(dataSp) <- ~ x.1+y.1

# list of spatial dataframes by individual 
split <- split(dataSp, hotspots$subject_name)

## Define a function to find the search radius and create an autocovariate for each individual in the list
autocov_apply <- function(x, type){
  # find search radius for the individual
  k1 <- knn2nb(knearneigh(coordinates(x)))
  d1 <- nbdists(k1,coordinates(x))
  crit.thresh <- max(unlist(nbdists(k1,coordinates(x))))
  
  # define autocovariate parameters
  type=type;                  # weights independent of distance 
  nbsize=crit.thresh;  # distance radius based on advice to select the lowest value at which all points have neighbors, this is the number of cells or pixels
  wstyle="B";           # valid weighting scheme with symmetrical weighting
  z = x$rel.staging.occ # response variable 
  
  # apply autocov function from spdep
  AC <- autocov_dist(z, coordinates(x), nbs = nbsize, zero.policy = T, style=wstyle, type = type)
  return(AC)
}

# apply function wit equal weighting
AC.list.equal <- lapply(split, autocov_apply, type = 'one')
AC.equal <- unlist(AC.list.equal)
# apply function with inverse weighting
AC.list.inverse <- lapply(split, autocov_apply, type = 'inverse')
AC.inverse <- unlist(AC.list.inverse)


hotspots$AC.equal <- scale(AC.equal, center = F, scale = T)
hotspots$AC.inverse <- scale(AC.inverse, center = F, scale = T)

##### Fit AC models #####

# define model statements
f1.ac <- rel.staging.occ ~ prop.forest.250 + AC.inverse + (1|subject_name)
f2.ac <- rel.staging.occ ~ prop.forest.250 + prop.ag.1500 + AC.inverse + (1|subject_name)
f3.ac <- rel.staging.occ ~ prop.forest.250 + prop.ag.1500 + gHM + AC.inverse + (1|subject_name)
f4.ac <- rel.staging.occ ~ prop.forest.250 + prop.ag.1500 + gHM + slope + AC.inverse + (1|subject_name)
f5.ac <- rel.staging.occ ~ prop.forest.250 + prop.ag.1500 + gHM + slope + dist2paedge + AC.inverse + (1|subject_name)
f6.ac <- rel.staging.occ ~ prop.forest.250 + prop.ag.1500 + gHM + slope + dist2paedge + drains250 + AC.equal + (1|subject_name)
f7.ac <- rel.staging.occ ~ prop.forest.250 + prop.ag.1500 + gHM + dist2paedge + AC.inverse + (1|subject_name)
f8.ac <- rel.staging.occ ~ prop.forest.250 + slope + drains250 + AC.inverse + (1|subject_name)
f9.ac <- rel.staging.occ ~ prop.ag.1500 + gHM + dist2paedge + AC.inverse + (1|subject_name) # very poor performance without forest
f10.safety1.ac <- rel.staging.occ ~ gHM*prop.forest.250 + prop.ag.1500 + gHM + dist2paedge + AC.inverse + (1|subject_name)
f11.safety2.ac <- rel.staging.occ ~ prop.forest.250 + gHM*prop.ag.1500 + gHM + dist2paedge + AC.inverse + (1|subject_name)

m1.ac <- glmer(f1.ac, 
                            weights = log(occ.weights), # weights supplied as total counts that proportions arise from
                            data = hotspots,
                            family = binomial)
summary(m1.ac)

m2.ac <- glmer(f2.ac, 
               weights = log(occ.weights), # weights supplied as total counts that proportions arise from
               data = hotspots,
               family = binomial)
summary(m2.ac)


m3.ac <- glmer(f3.ac, 
               weights = log(occ.weights), # weights supplied as total counts that proportions arise from
               data = hotspots,
               family = binomial)
summary(m3.ac)

m4.ac <- glmer(f4.ac, 
               weights = log(occ.weights), # weights supplied as total counts that proportions arise from
               data = hotspots,
               family = binomial)
summary(m4.ac)

m5.ac <- glmer(f5.ac, 
               weights = log(occ.weights), # weights supplied as total counts that proportions arise from
               data = hotspots,
               family = binomial)
summary(m5.ac)

m6.ac <- glmer(f6.ac, 
               weights = log(occ.weights), # weights supplied as total counts that proportions arise from
               data = hotspots,
               family = binomial)
summary(m6.ac)

m7.ac <- glmer(f7.ac, 
               weights = log(occ.weights), # weights supplied as total counts that proportions arise from
               data = hotspots,
               family = binomial)
summary(m7.ac)

m8.ac <- glmer(f8.ac, 
               weights = log(occ.weights), # weights supplied as total counts that proportions arise from
               data = hotspots,
               family = binomial)
summary(m8.ac)

m9.ac <- glmer(f9.ac, 
               weights = log(occ.weights), # weights supplied as total counts that proportions arise from
               data = hotspots,
               family = binomial)
summary(m9.ac)

m10.ac <- glmer(f10.safety1.ac, 
               weights = log(occ.weights), # weights supplied as total counts that proportions arise from
               data = hotspots,
               family = binomial)
summary(m10.ac)

m11.ac <- glmer(f11.safety2.ac, 
                weights = log(occ.weights), # weights supplied as total counts that proportions arise from
                data = hotspots,
                family = binomial)
summary(m11.ac)sum

# model table with AICc
selection <- model.sel(list(m1.ac, m2.ac, m3.ac, m4.ac, m5.ac, m6.ac, m7.ac, m8.ac, m9.ac, m10.ac), rank = 'AICc')
selection

# VIF test of top model
vif <- car::vif(m6.ac)
vif # nothing above 1.5

## Plot out results (Odds Ratios)
library(sjPlot)
tab_model(m6.ac, rm.terms = 'AC.inverse') 
plot_model(m6.ac, rm.terms = 'AC.inverse')
plot_model(m6.ac, type = "pred", terms = 'prop.forest.250')
plot_model(m6.ac, type = "pred", terms = 'gHM')
plot_model(m6.ac, type = "pred", terms = "dist2paedge")


f6 <- rel.staging.occ ~ prop.forest.250 + prop.ag.1500 + gHM + slope + dist2paedge + drains250 + (1|subject_name)
m6 <- glmer(f6.ac, 
               weights = log(occ.weights), # weights supplied as total counts that proportions arise from
               data = hotspots,
               family = binomial)
summary(m6)


# save top model
#saveRDS(m6.ac, './Staging Areas/hotspots_hmm_m6.RDS')
#m6.ac <- readRDS('./Staging Areas/hotspots_hmm_m6.RDS')

##### Compare HMM and Tortuosity Models #####
m6.ac.tor <- readRDS('./Staging Areas/hotspots_tor_m6.RDS')
m6.ac.hmm <- readRDS('./Staging Areas/hotspots_hmm_m6.RDS')

plot_models(m6.ac.hmm, m6.ac.tor, m.labels = c('Occurance of HMM-defined staging', 'Occurance of Tortuosity-defined staging'))


##### Maps #####

## map of staging and ag-use relocs
# define staging and ag-day reloc datasets
stage.relocs <- filter(used.df, vote.ag == 1) %>%
  st_as_sf(coords = c('x','y'), crs = 32736)
ag.relocs <- used.df %>% 
  filter(vote.ag == 0 & ag.window.ext == 1) %>%
  st_as_sf(coords = c('x','y'), crs = 32736) 

# tmap
library(tmap)
stage.points.hmm <- tm_shape(gme) + tm_polygons(col = 'pa_status', palette = RColorBrewer::brewer.pal(3, 'Greens'), title = 'Protected Areas') +
  tm_shape(ag.relocs) + tm_dots(size = 0.01, alpha = 0.03, col = '#FEE0D2') + # ag-day locs
  tm_shape(stage.relocs) + tm_dots(size = 0.01, alpha = 0.02, col = '#DE2D26') + # staging locs 
  # add manual legend for tm_dot layers
  tm_add_legend(title = 'Agricultural Use Locations', type = 'symbol', labels = c('ag-day', 'staging'), col = c('#FEE0D2', '#DE2D26')) 
# # add map title
# tm_layout(main.title = 'Staging Relocations (HMM)',
#           main.title.position = 'center')
stage.points.hmm

# tmap_save(stage.points.hmm, "Staging Relocs Map HMM 20220216.png", dpi = 300)
