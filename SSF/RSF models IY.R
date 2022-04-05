##### Resource Selection Model (Third Order Selection) #####

library(amt)
library(tidyverse)
library(sjPlot)
library(terra)
library(sf)
library(parallel)
library(tictoc)
library(lubridate)
library(mclust)

theme_set(  theme_bw() + # set theme with no legend of strip text
              theme(#panel.grid.major = element_blank(),
                    #panel.grid.minor = element_blank(), 
                    strip.background = element_blank(),
                    panel.border = element_rect(colour = "black"),
                    strip.text = element_text(size = 12),
                    legend.text = element_text(size = 10),
                    axis.line = element_line(colour = "black"), axis.title=element_text(size=14), axis.text = element_text(size=12))
)


set.seed(1)

##### 1a. Prepare Data #####

# read data 
#movdata <- readRDS('./SSF/eledata_expanded.rds')
#movdata <- readRDS('./SSF/eledata_allmara.RDS')
#movdata <- readRDS("./movdata/GMEcollars_004_usedClust_2021-10-05.rds") 
movdata <- readRDS("./movdata/GMEcollars_004_usedClust_2022-04-04.rds")
movdata <- movdata[!movdata$subject_name %in% c('Shamba','Courtney','David','Pepper','Harakati'),]
movdata <- movdata[movdata$fixType != 'irregular',]
movdata$uid <- 1:nrow(movdata)

## Remove relocations beyond landcover prediction layer
# keeps >80% of data for Lowana. 50% for Harakati, 95% Tai, 92% for Jer
ymin = 9691532.88
movdata <- filter(movdata, y >= (ymin+1000)) # add buffer

# create date object
movdata$date <- as.POSIXct(movdata$date) # check still in EAT
movdata$tod <- ifelse(hour(movdata$date) %in% c(6:17), 1,0)
# add new variable for individual-year (name-year) that can be used to group and fit models by. Years use april cutoff
movdata$subject_year <- paste(movdata$subject_name, movdata$year.cuts, sep = '-')

# check - 
t <- movdata %>% group_by(subject_year)
tt <- t %>% summarise(mean = round(mean(ag.used),2))

##### Define Seasons #####

# load season windows (Ecoscope colab export)
library(lubridate)
seasons <- read.csv('./spatial data/season_time_windows_wetdrytrans_20112021.csv')
seasons$start <- ymd_hms(seasons$start)
seasons$end <- ymd_hms(seasons$end)
# transition to dry -- based on staging paper, transition and dry season are best for detecting crop raiding season
seasons$season <- ifelse(seasons$season != 'wet', 'dry', 'wet')

# cut the dataframe by the start date of each season and label it with season name. 
# NOTE: switch the label to `unique_season` to check that season dates and labels are matching correctly for your dataset
rng.name <- cut(movdata$date, breaks = seasons$start, include.lowest = T, labels = head(seasons$season, -1))
# check
levels(rng.name)
# create variables
movdata$season <- rng.name


# make a move object for the amt package - we will use the id, which is a
#         combination of the animal's name and collar number
track <- make_track(movdata, .x = x, .y = y, .t = date, subject_name = subject_name, subject_year = subject_year, 
                    subject_sex = subject_sex, tactic.season = tactic.season, season = season, uid = uid, tod = tod)

##### 1b. Generate random samples #####
#' We will use amt functions to make random samples for our models. First, we will generate random points for the rsf
#' model using the `random_points` function. We will use an mcp home range estimator to generate the random points.
#' We are only using one individual, but the code is set up to scale with multiple individuals using the nest and unest
#' functions, as long as you identify the individuals using a column called `id`.

## make random points for RSF
# set cores
cores = 7

# create random points
{tic()
rand.factor = 2 # 30:1 unused:used
avail <- split(track, track$subject_year) # split variable (subject_name or subject_year)
avail <- parallel::mclapply(avail, function(x) {
  set.seed(1)
  rp <- x %>% amt::random_points(n = nrow(.)*rand.factor, hr = 'mcp') #
  rp$subject_name = unique(x$subject_name)
  rp$subject_year = unique(x$subject_year) # update this by the split variable
  rp$subject_sex = unique(x$subject_sex)
  rp$tactic.season = unique(x$tactic.season)
  return(rp) }, mc.cores = cores)
toc()}

avail <- do.call(rbind, avail)
avail <- filter(avail, case_ == FALSE)

# create day/night for available points
set.seed(5)
tod <- sample(c(0,1), nrow(avail), replace = TRUE, prob = c(0.5, 0.5))
avail$tod <- tod

# combine used and available
used <- track %>% mutate(case_ = "TRUE") %>%
  select(c(case_,x_,y_,subject_name,subject_year,subject_sex,tactic.season,tod))
rsf.df <- as.data.frame(rbind(used, avail))

##### Covariate Extraction #####

#' We now have used (case_ == TRUE) and available (case_ == FALSE) points for each individual. Before we can begin 
#' modeling, we have to extract covariates. We will use the amt package function to do this. 

# import layers - note that to use the amt functions they need to be raster objects, not terra objects
library(terra)
# import layers
prop.ag <- rast("./spatial data/Estes_ag_pct_1500.tif")	
prop.forest <- rast("./spatial data/hansen_forest_pct_250.tif")		
slope <- rast('./spatial data/slope_estes_32736_2020-05-12.tif')
ndviCoV <- rast('./spatial data/NDVICoV_estes_32736-2022-02-03.tif')
prop.settlement.250 <- rast("spatial data/google settlements/estes_settlement_pct_250.tif")
prop.settlement.1500 <- rast("spatial data/google settlements/estes_settlement_pct_1500.tif")
gHM <- rast("./spatial data/gHM_estes_32736_2020-05-12.tif")
roadsPrimary <- rast("./spatial data/Roads landDx/road_primary_landDx.tif")
roadsSecondary <- rast("./spatial data/Roads landDx/road_secondary_landDx.tif")
dist2roads.primary <- rast("./spatial data/Roads landDx/dist2road_primary_landDx.tif")
dist2roads.secondary <- rast("./spatial data/Roads landDx/dist2road_secondary_landDx.tif")
dist2ag <- rast("./spatial data/Tiedman/dist2ag_tiedman.tif")
dist2paedge <- rast("./spatial data/dist2paedge_estes_32736_20211118.tif")
pa.1000 <- rast("./spatial data/GSE/GSE Focal Rasters/protected_focal_1000.tif")
lu.1000 <- rast("./spatial data/GSE/GSE Focal Rasters/lu_focal_1000.tif")
up.1000 <- rast("./spatial data/GSE/GSE Focal Rasters/up_focal_1000.tif")
dist2settlement <- rast("./spatial data/google settlements/estes_dist2settlement_32736.tif")


# tiedman layer
# 1 = ag
# 2 = cover20
# 3 = cover2070
# 4 = cover70
# 5 = degraded
#lc <- rast("./spatial data/Tiedman/agmask_reclass.tif") # maraonly
lc <- rast("./spatial data/Tiedman/sentinel2018-3yr-GSE-03-2022_32736_agmask.tif") #SME - march2022
ag <- terra::clamp(lc, lower = 1, upper = 1, values = FALSE)
cover20 <- terra::clamp(lc, lower = 2, upper = 2, values = FALSE)
cover2070 <- terra::clamp(lc, lower = 3, upper = 3, values = FALSE)
cover70 <- terra::clamp(lc, lower = 4, upper = 4, values = FALSE)

# drains
# get shapefile for drains buffer
drains <- st_read("./spatial data/drains/drains_estes_20211117/drains_estes_-2021-11-17.shp", 
                  layer="drains_estes_-2021-11-17", crs = 4326) %>%
  st_transform(crs = 32736) %>%
  filter(RIV_ORD <= 7) %>%
# buffer by 250m
  st_buffer(dist = 250) %>%
  st_union() %>%
  terra::vect() %>%
  terra::rasterize(.,lc) # rasterize at 10m sentinel 

# pa (factor)
# Needs to potentially be updated a bit based on KS email responses 
pa <- st_read('./spatial data/GSE/GSEv_undissolved_temp_32736.shp') %>%
  terra::vect()
# rasterize mask layer
landuse <- terra::rasterize(pa, slope, field = 'pa_level') 
# NA to 1 (unprotected land)
landuse[is.na(landuse)] <- 1

## Create covariate list
# rasters must be in a list for mcapply
r.list <- list(slope, ndviCoV, ag, cover20, cover2070, cover70, drains, 
               prop.settlement.250, prop.settlement.1500, gHM, 
               dist2roads.primary, dist2roads.secondary, dist2ag, 
               landuse, dist2paedge, pa.1000,lu.1000,up.1000,dist2settlement)

## Extract
study.area <- 32736

## Extraction in parallel - 40 minutes for mara only
# need to split up to avoid ram issues creating large terra::vect() objects. can't handle 10mil+ datasets
split <- split(rsf.df, rsf.df$subject_name)
cores = 7
used <- NULL
{tic()
  for(i in 1:length(split)) {
    # create vect
    vect <- split[[i]] %>% 
      st_as_sf(coords = c('x_','y_'), crs = study.area) %>%
      terra::vect()
    r.extract <- mclapply(r.list, function(x)
      terra::extract(x, vect)[,2], mc.cores = cores)
    used[[i]] <- do.call(cbind, r.extract) # bind results into an extracted dataframe for each individual
  }
toc()}

used2 <- do.call(rbind, used)

# check
head(used2)
summary(used2)

# create data frame
mode(used2) = "numeric"
used2 <- as.data.frame(used2)
colnames(used2) <- c('slope','ndviCoV','ag','cover20','cover2070','cover70','drains',
                     'prop.settlement.250','prop.settlement.1500','gHM',
                     'dist2roads.primary','dist2roads.secondary','dist2ag',
                     'landuse','dist2paedge','pa.1000','lu.1000','up.1000','dist2settlement')
head(used2)

# unstandardized data frame
rsf.ext <- cbind(rsf.df, used2)

## Adjust dist2edge metrics
rsf.ext$dist2paedge <- ifelse(rsf.ext$landuse > 1, -rsf.ext$dist2paedge, rsf.ext$dist2paedge)


# ##### Add Ag Availability for Homerange
# ## Ag availability by home range
# # 1. Calculate mcp homeranges
# # 2. Calculate proportion of ag within each homerange polygon
# # 3. Attach individual-year ag availability values to the dataset
# 
# df <- rsf.ext %>% ungroup() %>% filter(case_ == TRUE) %>% select(subject_year, x_, y_)
# sp::coordinates(df) <- ~x_+y_
# split <- split(df, df$subject_year)
# mcp <- lapply(split, function(x) {
#   x <- adehabitatHR::mcp(x)
#   x <- st_as_sf(x) %>%
#     terra::vect()
#   return(x)})
# 
# # get mean ag in homerange (proportion 1/0)
# ag2 <- terra::classify(ag, cbind(NA, 0), right=FALSE)
# {tic()
#   cores = 7
#   ag.ext <- mclapply(mcp, function(x) terra::extract(ag2, x, fun = mean)[,2])
#   toc()}
# 
# t <- as.data.frame(do.call(rbind, ag.ext))
# t$subject_year <- rownames(t)
# colnames(t) <- c('ag.avail', 'subject_year')
# 
# write.csv(t, './SSF/ag.homerange.IY.datatable_20220312.csv', row.names = FALSE)

# read in ag available table
t <- read.csv('./SSF/ag.homerange.IY.datatable_20220312.csv')

rsf.ext <- merge(rsf.ext, t, by = 'subject_year') 

# save extracted df and compress
write.csv(rsf.ext, './SSF/RSF_dataset_extracted_20220405.csv')
system('gzip ./SSF/RSF_dataset_extracted_20220405.csv')

# drop old dataframes
remove(r.extract)
remove(split)
remove(used)
remove(rsf.df)
remove(r.list)
remove(used2)

##### Save Point #####
# Can make the following code a seperate file for model fitting
rsf.ext <- as.data.frame(data.table::fread('./SSF/RSF_dataset_extracted_20220403.csv.gz'))
rsf.ext$subject_year <- as.character(rsf.ext$subject_year)
rsf.ext$tactic.season <- as.factor(rsf.ext$tactic.season)
rsf.ext$V1 <- NULL


## Prep Dataframe

## scale some covariates before we fit it
rsf.ext <- filter(rsf.ext, ndviCoV > -13 & ndviCoV < 54)

t <- rsf.ext %>%
  mutate_at(.vars = c('slope', 'ndviCoV'), .funs = scale) %>%
  #mutate_at(.vars = c('gHM','prop.settlement.250','prop.settlement.1500'), .funs = scale)
  mutate_at(.vars = c('ag','cover20','cover2070','cover70','drains'), .funs = function(x) if_else(!is.na(x), 1, 0)) %>%
  mutate_at(.vars = c('ag','cover20','cover2070','cover70', 'drains', 'landuse'), .funs = as.factor)

rsf.ext$case_ <- as.factor(rsf.ext$case_)
rsf.ext$subject_name <- as.factor(rsf.ext$subject_name)

# filter out NDVI cov outliers
rsf.ext <- filter(rsf.ext, ndviCoV < 20) # removes three major outliers from Marima dataset

# add day/night
rsf.ext$tod 


# check dataframe
summary(rsf.ext)

## filter individuals (>1000 points)
t <- rsf.ext %>% group_by(subject_year) %>% 
  filter(case_ == 'TRUE') %>%
  tally() %>%
  filter(n > 1000) # set min number of points for subject_name or subject_year

rsf.ext<- filter(rsf.ext, subject_year %in% t$subject_year)
rsf.ext <- filter(rsf.ext, !(subject_year %in% c('Caroline-2014','Tressa-2015','Chelsea-2016','Jinomoja-2018','Jinomoja-2019','Jinomoja-2020')))


## add extra metadata
t <- movdata %>% group_by(subject_year, region, tactic.season) %>% rename(tactic.season.new = tactic.season) %>% 
  tally() %>% select(-n)

tt <- merge(rsf.ext, t, by = 'subject_year')
dim(rsf.ext)

# check ag
t <- rsf.ext %>% group_by(subject_year) %>% 
  summarise(mean = mean(as.numeric(as.character(ag))))

# # fit pop model
# m <- glm(case_ ~ ag + cover20 + cover2070 + cover70 + slope + ndviCoV + drains + scale(gHM), data = rsf.ext, family = binomial)
# summary(m)
# m1 <- rsf.ext %>% fit_logit(case_ ~ ag + cover2070 + cover70 + slope + ndviCoV + drains + gHM)
# m2 <- rsf.ext %>% fit_logit(case_ ~ ag + cover70 + slope + ndviCoV + gHM)
# m3 <- rsf.ext %>% fit_logit(case_ ~ ag + prop.settlement.1500)
# m4 <- rsf.ext %>% fit_logit(case_ ~ cover20 + cover2070 + cover70 + slope + ndviCoV + drains)
# m5 <- rsf.ext %>% fit_logit(case_ ~ ag*cover70 + slope + ndviCoV + gHM)
# m.250 <- glm(case_ ~ ag + cover20 + cover2070 + cover70 + slope + ndviCoV + drains + scale(prop.settlement.250), data = rsf.ext, family = binomial)
# m.1500 <- glm(case_ ~ ag + cover20 + cover2070 + cover70 + slope + ndviCoV + drains + scale(prop.settlement.1500), data = rsf.ext, family = binomial)
# 
# tab <- model.sel(m, m.250, m.1500)
# tab
# 
# # m.rsf <- glmer(case_ ~ ag + cover20 + cover2070 + cover70 + slope + ndviCoV + (1|subject_name),
# #               data = rsf.ext, family = binomial)
# 
# library(MuMIn)
# tab <- model.sel(m, m1$model, m2$model, m3$model, m4$model, m5$model, m.250, m.1500)
# View(tab)
# 
# # check out the summary of top model
# summary(m)
# 
# # viz the covariates
# #sjPlot::plot_model(m, transform = NULL, df_method='wald')


##### Fit Individual-level RSF #####
## Run in parallel 
# fit model
# 
# rsf.ext$dist2lu <- ifelse(rsf.ext$dist2lu > 5000, 5000, rsf.ext$dist2lu)
# rsf.ext$dist2lu <- ifelse(rsf.ext$dist2lu < -5000, -5000, rsf.ext$dist2lu)
# rsf.ext$dist2up <- ifelse(rsf.ext$dist2up > 5000, 5000, rsf.ext$dist2up)
# rsf.ext$dist2up <- ifelse(rsf.ext$dist2up < -5000, -5000, rsf.ext$dist2up)

{tic()
  split <- split(rsf.ext, rsf.ext$subject_year) # subject_name or subject_year
  m.ind.rsf <- lapply(split, function(x)
    glm(case_ ~ ag + cover20 + cover2070 + cover70 + slope + ndviCoV + drains + scale(prop.settlement.250) +
          scale(dist2roads.primary) + scale(dist2roads.secondary) + scale(dist2paedge), data = x, family = 'binomial'))
toc()}

remove(split)
#saveRDS(m.ind.rsf, paste0('./SSF/RSF_modfit_20220404.RDS'))

# get table of estimates for each individual
id.est.rsf <- lapply(m.ind.rsf, function(x) 
  tab <- x %>% broom::tidy() %>%
    mutate(
      ymin = estimate - 1.96 * std.error,
      ymax = estimate + 1.96 * std.error
    ))
id.est.rsf <- dplyr::bind_rows(id.est.rsf, .id = "subject_year")

# test variance inflation factors
id.vif.rsf <- lapply(m.ind.rsf, function(x) {
  t <- car::vif(x)
  vif <- as.data.frame(cbind(names(t), unname(t))) })
t <- dplyr::bind_rows(id.vif.rsf, .id = "subject_year")
t <- as.data.frame(t)

# add other ele metadata (sex, ag tactic, etc.)
t <- rsf.ext %>% group_by(subject_year, subject_name, subject_sex, tactic.season, tactic.season.new, region, ag.avail) %>% tally() %>% dplyr::select(-n)
id.est.rsf <- merge(id.est.rsf, t, by = 'subject_year')

#id.est.rsf <- filter(id.est.rsf, subject_name != 'Mtembezi')

#saveRDS(id.est.rsf, paste0('./SSF/RSF_estimates_20220404.RDS'))

# remove intercept
id.est.rsf <- filter(id.est.rsf, term != '(Intercept)')
#id.est.rsf <- filter(id.est.rsf, term != 'ag1')
##### plot results #####

# coeff estimates for all individuals
tt <- filter(id.est.rsf, subject_name %in% c('Ivy','Fitz', 'Chelsea', 'Lowana'))
ggplot(tt, aes(x = term, y = estimate, color = tactic.season.new)) + 
  geom_pointrange(aes(ymin = ymin, ymax = ymax),
                  position = position_dodge(width = 0.6)) +
  xlab('covariate') + ggtitle('rsf individual-level estimates') +
  labs(color = 'ag tactic') + 
  coord_flip() + facet_wrap(~subject_name)

# coeff estimates by sex
ggplot(id.est.rsf, aes(x = as.factor(term), y = estimate, color = subject_sex)) +
  geom_boxplot() +
  xlab('covariate') + ggtitle('rsf individual-level estimates by sex') +
  labs(color = 'sex') + 
  coord_flip()

# coeff estimates by tactic (IY)
ggplot(id.est.rsf, aes(x = as.factor(term), y = estimate, color = factor(tactic.season.new))) +
  geom_boxplot() +
  xlab('covariate') + ggtitle('rsf individual-level estimates by ag tactic') +
  labs(color = 'ag tactic (year)') + 
  coord_flip() 

##### Model Function Response with Ag #####
ag.est <- filter(id.est.rsf, term == 'ag1')
plot(ag.est$ag.avail, ag.est$estimate)
ag.est$ag.avail <- log(ag.est$ag.avail+0.0001)
ag.est$subject_name <- as.factor(ag.est$subject_name)
m.funct <- lme4::lmer(estimate ~ ag.avail + (1|subject_name), data = ag.est)
sjPlot::plot_model(m.funct, terms = "ag.avail [all]", type = 'pred', show.data = TRUE, title = '', axis.labels = '', dot.size = 1.5)


##### Random Forest Clustering #####
# Use tactic season classifications as testing for clustering

# create RF table - pivot wider
wide <- id.est.rsf %>% select(-c("std.error","statistic","p.value","ymin","ymax","subject_name")) %>%
  pivot_wider(names_from = term, values_from = estimate) #%>%
  #filter(ag1 != 0) 
#wide$tactic.season <- ifelse(wide$tactic.season == 3, 2, wide$tactic.season)

colnames(wide) <- c("subject_year","subject_sex","tactic.season","tactic.season.new","region","ag.avail","ag","cover20","cover2070","cover70","slope",
                    "ndviCoV","drains","prop.settlement.250","dist2roads.primary","dist2roads.secondary","dist2paedge")

## Data partition
set.seed(5)
ind <- sample(2, nrow(wide), replace = TRUE, prob = c(0.8, 0.2))
train <- wide[ind==1,]
test <- wide[ind==2,]
dim(train)
dim(test)

## Train
rf <- randomForest::randomForest(as.factor(tactic.season.new) ~ ag + cover20 + cover2070 + cover70 + slope + ndviCoV + drains + 
                                   prop.settlement.250 + dist2roads.primary + dist2roads.secondary + dist2paedge + region, 
                                 data = train, ntree = 1000, proximity=TRUE)
rf

randomForest::varImpPlot(rf)
randomForest::importance(rf)


# Predict - Confusion matrix
p1 <- predict(rf, test)
caret::confusionMatrix(p1, as.factor(test$tactic.season))


##### Random Forest Regression #####
# Use tactic season classifications as response var in a RF regression to 
# identify most important covariates for distinguishing clusters

wide$subject_name <- stringr::str_extract(wide$subject_year, "[^-]+")

rf <- randomForest::randomForest(as.numeric(as.factor(subject_year)) ~ ag + cover20 + cover2070 + cover70 + slope + ndviCoV + drains + 
                                   prop.settlement.250 + dist2roads.primary + dist2roads.secondary + dist2paedge + region, 
                                 data = wide, ntree = 1000, importance=TRUE)

rf

randomForest::varImpPlot(rf)
randomForest::importance(rf)


##### Gaussian Mixture Clustering #####
library(mclust)

# prepare data - remove non-informative variables for clustering
wide.clust <- wide %>% select(region, slope, dist2paedge, prop.settlement.250, ndviCoV, drains, dist2roads.primary, ag, cover20, cover70)

clust <- mclust::Mclust(wide.clust)

plot(clust, what = "classification",
     xlab = "max",
     ylab = "cluster")


# dimension reduction
drclust <- mclust::MclustDR(clust)

plot(drclust)

# classify
wide$classification <- clust$classification













