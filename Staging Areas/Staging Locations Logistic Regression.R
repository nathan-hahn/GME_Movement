#### Logistic Regression - Staging Locations ####
library(tidyverse)
library(lme4)
library(terra)
library(sf)
library(MuMIn)
library(mapview)
library(tmap)
library(lubridate)

theme_set(  theme_bw() + # set theme with no legend of strip text
              theme(panel.grid.major = element_blank(),
                    strip.background = element_blank(),
                    panel.border = element_rect(colour = "black"),
                    strip.text = element_text(size = 12),
                    legend.text = element_text(size = 10),
                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title=element_text(size=14), axis.text = element_text(size=12))
)

Sys.setenv(TZ='Africa/Nairobi')

###########################################################################################################
## ///////////////////////////////////////////// Data Prep ///////////////////////////////////////////// ##
###########################################################################################################

#### Load Data ####
#df <- as.data.frame(data.table::fread('./movdata/GMEcollars_004_stageclassified_20211019.csv', tz = ''))
df <- as.data.frame(data.table::fread('./Staging Areas/movdata/GMEcollars_004_stageclassified_20220109.csv', tz=''))
df$V1 <- NULL
df$V1 <- NULL
df$...1 <- NULL
df$used <- NULL

# # for nsd**
nsd <- as.data.frame(data.table::fread('./Staging Areas/movdata/GMEcollars_004_stageclassified_NSD_20211019.csv', tz='')) %>%
  select(uid, vote.ag, r2n)
df$vote.ag <- NULL # drop HMM vote
df <- merge(nsd, df, by = 'uid', all.x = T)

# for tor**
# tor <- as.data.frame(data.table::fread('./Staging Areas/movdata/GMEcollars_004_stageclassified_TOR_20211019.csv', tz='')) %>%
#   select(uid, vote.ag, tor)
# df.vote.ag <- NULL #  drop HMM vote
# df <- merge(tor, df, by = 'uid', all.x = T)

# correct pardamat ag to non-ag
df$ag.used <- ifelse(df$site == 'mep' & df$pa == 2, 0, df$ag.used)

df$vote.ag <- as.factor(df$vote.ag)
table(df$vote.ag)

# filter to relocs of interest -- fixes during ag use periods
# df <- df %>%
#   filter(ag.window.ext == 1)

#### Extract Covariates ####
# covs of interest are prop.forest, prop.ag, drainage (1/0), dist2water, gHM, slope

## Remove existing covariates
df <- df %>%
  dplyr::select(-c('dist2ag', 'dist2agedge', 'dist2permwater', 'dist2seasonalwater', 'dist2water', 'dist2forest', 'gHM', 'slope', 'pa'))

##### Import covariate layers #####

## Raster layer prep
prop.ag.250 <- rast("./spatial data/estes_ag_pct_250.tif")
prop.ag.1500 <- rast("./spatial data/estes_ag_pct_1500.tif")
prop.forest.250 <- rast("./spatial data/hansen_forest_pct_250.tif")
prop.forest.1500 <- rast("./spatial data/hansen_forest_pct_1500.tif")
dist2ag <- rast("./spatial data/dist2ag_estes_32736_2022-01-03.tif")
#dist2ag <- rast("./spatial data/dist2ag_estes_32736_2019-11-21.tif")
dist2agedge <- rast("./spatial data/dist2agedge_estes_20200629.tif")
dist2water <- rast("./spatial data/dist2merged_water_20200624.tif")
dist2forest <- rast("./spatial data/dist2forest_hansen_cover60_32736_30.tif") 
slope <- rast("./spatial data/slope_estes_32736_2020-05-12.tif")
lc <- rast("./spatial data/change03_181_reclassMara_2019-11-22.tif")
gHM <- rast("./spatial data/gHM_estes_32736_2020-05-12.tif")
prop.settlement.250 <- rast("spatial data/estes_settlement_pct_250.tif")
prop.settlement.1500 <- rast("spatial data/estes_settlement_pct_1500.tif")
dist2paedge <- rast("./spatial data/dist2paedge_estes_32736_20211118.tif")

pa <- rast("./spatial data/GSEr_dissolved_estes_32736_20220120.tif") # 3 - PA, 2 - cc
pa <- terra::classify(pa, cbind(NA, 0)) # set no data to 'not protected' (0)
plot(pa)

## Polygon layer prep
# get shapefile for drains buffer
drains <- st_read("./spatial data/drains/drains_estes_20211117/drains_estes_-2021-11-17.shp", 
                  layer="drains_estes_-2021-11-17", crs = 4326) %>%
  st_transform(crs = 32736) %>%
  filter(RIV_ORD <= 7)

# buffer by 1000m
drains <- st_buffer(drains, dist = 250) %>%
  st_union() 
#plot(drains.1000)

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
  used[,5] <- as.numeric(st_intersects(locs.sf, drains)) # point to polygon
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
colnames(used2) <- c("used","dist2ag", "dist2agedge", "dist2water", 'drains250', 'slope', 'gHM', 'pa', 'lc', 'dist2forest', 
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

##### Add Forest Proximity Index #####
## Pre-computed 

# read in results - tagged by uid
setwd('./spatial data/Forest_hansen/distanceloop/output')
prox.index.1000 <- list.files(pattern = "1000.RDS") %>%
  map_dfr(readRDS) %>%
  dplyr::select(prox.index.1000 = prox.index, uid = point.id) %>%
  group_by(uid) %>%
  summarise(prox.index.1000 = mean(prox.index.1000))

head(prox.index.1000)

used.df <- merge(used.df, prox.index.1000, 'uid', all.x = TRUE)

# check for duplicates
length(unique(used.df$uid)) == nrow(used.df)

setwd("/Users/nhahn/Dropbox/CSU/GME_Movement")


##### Define Seasons #####

# load season windows (Ecoscope colab export)
library(lubridate)
seasons <- read.csv('./spatial data/season_time_windows_wetdrytrans_20112021.csv')
seasons$start <- ymd_hms(seasons$start)
seasons$end <- ymd_hms(seasons$end)
used.df$date <- ymd_hms(used.df$date)

# cut the dataframe by the start date of each season and label it with season name. 
# NOTE: switch the label to `unique_season` to check that season dates and labels are matching correctly for your dataset
rng.name <- cut(used.df$date, breaks = seasons$start, include.lowest = T, labels = head(seasons$season, -1))

# check
levels(rng.name)

# create variables
used.df$season <- rng.name

# check for duplicates
length(unique(used.df$uid)) == nrow(used.df)

##### Save output #####
#write.csv(used.df, './Staging Areas/movdata/movdat_004_lsdv_20220109.csv')
#write.csv(used.df, './Staging Areas/movdata/movdat_004_lsdv_nsd_20220109.csv')
#write.csv(used.df, './Staging Areas/movdata/movdat_004_lsdv_tor_20220109.csv')



###########################################################################################################
## ///////////////////////////////////////// Stage Timing ////////////////////////////////////////////// ##
###########################################################################################################

used.df <- as.data.frame(data.table::fread('./Staging Areas/movdata/movdat_004_lsdv_20220109.csv', tz=''))
used.df$V1 <- NULL

stage <- used.df %>%
  filter(vote.ag == 1)

# staging by season
prop.table(table(stage$season))

# staging by season -- with individual-years
stage.season <- stage %>% group_by(subject_name, year.cuts, season) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n),
         sd = sd(n)) %>%
  filter(!is.na(season))

boxplot(stage.season$freq ~ stage.season$season, xlab = 'Season', ylab = "Percent of Stages in Year")


ggplot(stage.season, aes(x = season, y = freq)) + geom_boxplot() + 
  geom_point(aes(size = n/100), position = position_dodge2(width = .6), color = 'darkgrey') +
  ylab('Percent of Stages in Year') + xlab('Season')


## Stage start/end times
stage.time <- stage %>% 
  group_by(subject_name, year.cuts, stage.event.index) %>%
  mutate(hour = hour(date)) %>%
  mutate(stage.start = min(hour),
         stage.end = max(hour)) %>%
  group_by(stage.event.index, season) %>%
  summarise(stage.start = min(stage.start),
            stage.end = min(stage.end)) %>%
  filter(!is.na(season))

library(reshape2)
melt <- melt(stage.time, id.vars = c('stage.event.index', 'season'))

boxplot(melt$value ~ melt$variable)  

ggplot(melt, aes(x = variable, y = value)) + geom_boxplot() + geom_point(shape = 15, position = position_dodge2(width = 0.5))

ggplot(melt, aes(x = value, fill = variable)) + geom_histogram(position = "dodge", bins = 10) +
  scale_x_continuous(breaks=seq(8,18,1)) +
  #facet_wrap(. ~ season) + 
  xlab('Hour of Day') + ylab('Number of Stages')

ggsave('./Staging Areas/plots/stage_time_hist_seasonsfacet.tiff', dpi = 300)



t <- stage %>% 
  group_by(stage.event.index) %>%
  tally()
summary(t$n)


## Activity Budgets of timing
stage.ac <- used.df %>%
  mutate(vote.ag = as.numeric(levels(vote.ag))[vote.ag]) %>%
  # tag full days for ag stage
  group_by(dayBurst) %>%
  mutate(vote.ag = ifelse(sum(vote.ag) > 1, 1, 0)) %>%
  # identify non-ag, ag non-stage, and ag stage days
  mutate(stage2 = as.numeric(vote.ag) + ag.window.ext) %>% 
  mutate(stage2 = as.factor(stage2)) %>% 
  mutate(stage2 = recode_factor(stage2, '0'='non-ag use day','1'='ag use day/non-stage','2'='ag use day/stage'))

# calculate pct encamped 
stage.ac <- stage.ac %>% 
  filter(!is.na(ag.window.ext)) %>%
  mutate(hour = hour(date)) %>%
  group_by(subject_name, hour, stage2) %>%
  mutate(enc.day = sum(viterbi==1), meander.day = sum(viterbi==2), dw.day = sum(viterbi==3)) %>%
  mutate(pct = round((enc.day/n()), 3))

# calculate 95% CIs
stage.ac <- stage.ac %>%
  group_by(hour, stage2) %>%
  summarise(mean.pct = mean(pct, na.rm = TRUE),
            sd = sd(pct, na.rm = TRUE),
            n = n(),
            se = sd / sqrt(n),
            lwr.ci = mean.pct - qt(1 - (0.05 / 2), n - 1) * se,
            upr.ci = mean.pct + qt(1 - (0.05 / 2), n - 1) * se) 

# plot it!!
ggplot(stage.ac, aes(x = hour, y = mean.pct)) + 
  geom_pointrange(aes(ymin = lwr.ci, ymax = upr.ci, shape = stage2)) +
  xlab("Hour of Day") + ylab("Percent of Encamped Fixes") + labs(shape = '') + 
  # add sunrise/sunset
  geom_vline(xintercept=6
             ,color="dark grey", linetype="dashed", size=0.3) + 
  geom_vline(xintercept=18,
             color="dark grey", linetype="dashed", size=0.3)
#ggsave('./Staging Areas/plots/activitybudgets_HMM.tiff', dpi = 300)





###########################################################################################################
## ///////////////////////////////////////// Model Fitting ///////////////////////////////////////////// ##
###########################################################################################################

#### Model Fitting ####

##### Create Model Dataframe #####
# used.df <- as.data.frame(data.table::fread('./Staging Areas/movdata/movdat_004_lsdv_20220109.csv'))
# used.df$V1 <- NULL

## **remove relocations within ag**
stagelocs <- used.df %>%
  #filter(lubridate::hour(date) %in% c(6:18)) #%>%
  filter(ag.window.ext == 1) %>%
  filter(dist2ag > 0)

# check for duplicates
length(unique(stagelocs$uid)) == nrow(stagelocs)

## adjust proximity index values
stagelocs$prox.index.1000 <- ifelse(is.na(stagelocs$prox.index.1000), 0.00001, stagelocs$prox.index.1000)

## standardize covs
covariates <- c("dist2ag", "dist2agedge", "dist2water", 'slope', 'dist2forest', 'dist2paedge')
stagelocs <- stagelocs %>%
  dplyr::select(uid, subject_name, subject_sex, burst, x, y, date, vote.ag, all_of(covariates), drains250, gHM, prop.ag.250, prop.ag.1500,
                prop.forest.250, prop.forest.1500, prop.settlement.250, prop.settlement.1500, prox.index.1000, pa, season) %>%
  mutate(forest = if_else(dist2forest == 0, 1, 0)) %>%
  mutate_at(covariates, .funs = scale) %>%
  mutate(drains250 = as.factor(drains250)) %>%
  mutate(pa = as.factor(pa), season = as.factor(season)) %>%
  droplevels()

# set dummy reference levels
stagelocs$season <- relevel(stagelocs$season, ref = 'wet') # season dummy
stagelocs$pa <- relevel(stagelocs$pa, ref = '3') # pa dummy -- fully protected

##### test autocorrelation variogram #####
# library(gstat)
# t <- st_as_sf(stagelocs, coords = c('x','y'), crs = 32736)
# t <- as_Spatial(t)
# 
# t.var<-gstat::variogram(vote~1, t)
# plot(t.var)

##### test spatial scales for moving window metrics #####
mod.ag.step <- glmer(vote ~ prop.ag.250 + (1|subject_name),
                     data = stagelocs, family = binomial)
mod.ag.daily <- glmer(vote ~ prop.ag.1500 + (1|subject_name),
                      data = stagelocs, family = binomial)
AICc(mod.ag.step, mod.ag.daily) # ag daily scale is better by AICc, ag step model is rank deficient and ag drops from model

mod.for.step <- glmer(vote ~ prop.forest.250 + (1|subject_name),
                      data = stagelocs, family = binomial)
mod.for.daily <- glmer(vote ~ prop.forest.1500 + (1|subject_name),
                       data = stagelocs, family = binomial)
mod.for.dist <- glmer(vote ~ dist2forest + (1|subject_name),
                      data = stagelocs, family = binomial)

AICc(mod.for.step, mod.for.daily, mod.for.dist) # forest step scale is better by AICc

# check proximity index
mod.prox.1000 <- glmer(vote.ag ~ log(prox.index.1000) + (1|subject_name),
                       data = stagelocs, family = binomial)


## correlation between forest prop and proximity index
t <- lm(log(prox.index.250) ~ prop.forest.250, data = stagelocs)

##### Fit candidate models #####

## Global model
mod.global.sub <- glmer(vote.ag ~ prop.ag.1500 + prop.forest.250 + drains250 + slope + gHM + dist2paedge + season + (1|subject_name), 
                        data = stagelocs, family = binomial)

## Human footprint model - danger
mod.hm.sub <- glmer(vote.ag ~ prop.ag.1500 + gHM + dist2paedge + (1|subject_name),
                    data = stagelocs, family = binomial)

## Natural features model - safety
mod.nat.sub <- glmer(vote.ag ~ prop.forest.250 + drains250 + slope + (1|subject_name),
                     data = stagelocs, family = binomial)

## Strongest predictors
mod.str.sub <- glmer(vote.ag ~ prop.ag.1500 + prop.forest.250 + gHM + (1|subject_name),
                     data = stagelocs, family = binomial)

mod.str2.sub <- glmer(vote.ag ~ prop.ag.1500 + prop.forest.250 + (1|subject_name),
                      data = stagelocs, family = binomial)

mod.str3.sub <- glmer(vote.ag ~ gHM + (1|subject_name),
                      data = stagelocs, family = binomial)

# staging relative to safety (-gHM) and cover (+forest)
mod.global.safetyint <- glmer(vote.ag ~ prop.ag.1500 + gHM*prop.forest.250 + drains250 + slope + dist2paedge + (1|subject_name), 
                                data = stagelocs, family = binomial)

# staging relative to food availability (+ag) and safety (-gHM)/(+forest)
mod.global.foodint <- glmer(vote.ag ~ gHM*prop.ag.1500 + prop.forest.250 + drains250 + slope + dist2paedge + (1|subject_name), 
                                data = stagelocs, family = binomial)
mod.global.foodint2 <- glmer(vote.ag ~ prop.ag.1500*prop.forest.250 + gHM + drains250 + slope + dist2paedge + (1|subject_name), 
                            data = stagelocs, family = binomial)

t <- AICc(mod.global.sub, mod.hm.sub, mod.nat.sub, mod.str.sub, mod.str2.sub, mod.str3.sub, mod.global.safetyint, mod.global.foodint, mod.global.foodint2)

t <- AICc(mod.global.sub, mod.global.safetyint, mod.global.foodint, mod.global.foodint2)
t

summary(mod.global.sub)
sjPlot::tab_model(mod.global.sub, transform = 'plogis') # NULL/plogis

sjPlot::plot_model(mod.global.sub, type = "pred", terms = c("prop.ag.1500 [all]")) + xlab('Proportion of Agriculture (1500m)') + ylab('Probability of Staging')

sjPlot::plot_model(mod.global.sub, type = "pred", terms = c("prop.forest.250 [all]")) + xlab('Proportion of Forest (250m)') + ylab('Probability of Staging')

sjPlot::plot_model(mod.global.sub, type = "pred", terms = c("gHM [all]")) + xlab('Human Modification Index') + ylab('Probability of Staging')


##### Male/Female Model #####


male <- filter(stagelocs, subject_sex == "male")
# update scaling for the subsetted dataset
male$dist2paedge <- male$dist2paedge * attr(male$dist2paedge, 'scaled:scale') + attr(male$dist2paedge, 'scaled:center')
male$dist2paedge <- scale(male$dist2paedge)

female <- filter(stagelocs, subject_sex == "female")
# update scaling f or the subsetted dataset
female$dist2paedge <- female$dist2paedge * attr(female$dist2paedge, 'scaled:scale') + attr(female$dist2paedge, 'scaled:center')
female$dist2paedge <- scale(female$dist2paedge)


## Global model
mod.global.male <- glmer(vote.ag ~ prop.ag.1500 + prop.forest.250 + slope + gHM + dist2paedge + season + (1|subject_name), 
                        data = male, family = binomial)

summary(mod.global.male)

mod.global.female <- glmer(vote.ag ~ prop.ag.1500 + prop.forest.250 + slope + gHM + dist2paedge + season + (1|subject_name), 
                         data = female, family = binomial)


summary(mod.global.female)

plot_models(mod.global.male, mod.global.female)

