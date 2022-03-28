##### Resource Selection Model (Third Order Selection) #####

library(amt)
library(tidyverse)
library(sjPlot)
library(terra)
library(sf)
library(parallel)
library(tictoc)

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
movdata <- readRDS("./movdata/GMEcollars_004_usedClust_2021-10-05.rds") 
movdata <- movdata[!movdata$subject_name %in% c('Shamba','Courtney','David','Pepper','Harakati'),]
movdata <- movdata[movdata$fixType != 'irregular',]
movdata$uid <- 1:nrow(movdata)

## Remove relocations beyond landcover prediction layer
# keeps >80% of data for Lowana. 50% for Harakati, 95% Tai, 92% for Jer
ymin = 9691532.88
movdata <- filter(movdata, y >= (ymin+1000)) # add buffer

# create date object
movdata$date <- as.POSIXct(movdata$date) # check still in EAT

# add new variable for individual-year (name-year) that can be used to group and fit models by. Years use april cutoff
movdata$subject_year <- paste(movdata$subject_name, movdata$year.cuts, sep = '-')

# make a move object for the amt package - we will use the id, which is a
#         combination of the animal's name and collar number
track <- make_track(movdata, .x = x, .y = y, .t = date, subject_name = subject_name, subject_year = subject_year, subject_sex = subject_sex, tactic.season = tactic.season, uid = uid)
t1 <- split(track, track$subject_year) # subject_name or subject_year

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
rand.factor = 30 # 30:1 unused:used
rsf.df <- do.call(rbind, t1)
rsf.df <- split(rsf.df, rsf.df$subject_year) # split variable (subject_name or subject_year)
rsf.df <- parallel::mclapply(rsf.df, function(x) {
  set.seed(1)
  rp <- x %>% amt::random_points(n = nrow(.)*rand.factor, hr = 'mcp') # 
  rp$subject_name = unique(x$subject_name)
  rp$subject_year = unique(x$subject_year) # update this by the split variable
  rp$subject_sex = unique(x$subject_sex)
  rp$tactic.season = unique(x$tactic.season)
  return(rp) }, mc.cores = cores)
toc()}

t <- rsf.df

rsf.df <- do.call(rbind, rsf.df)

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
prop.settlement.250 <- rast("spatial data/estes_settlement_pct_250.tif")
prop.settlement.1500 <- rast("spatial data/estes_settlement_pct_1500.tif")
gHM <- rast("./spatial data/gHM_estes_32736_2020-05-12.tif")
roadsPrimary <- rast("./spatial data/Roads landDx/road_primary_landDx.tif")
roadsSecondary <- rast("./spatial data/Roads landDx/road_secondary_landDx.tif")
dist2roads.primary <- rast("./spatial data/Roads landDx/dist2road_primary_landDx.tif")
dist2roads.secondary <- rast("./spatial data/Roads landDx/dist2road_secondary_landDx.tif")
dist2ag <- rast("./spatial data/Tiedman/dist2ag_tiedman.tif")

# tiedman layer
# 1 = ag
# 2 = cover20
# 3 = cover2070
# 4 = cover70
# 5 = degraded
#lc <- rast("./spatial data/Tiedman/agmask_reclass.tif") # maraonly
lc <- rast("./spatial data/Tiedman/sentinel2018-3yr-GSE-02-2022_32736_agmask.tif") #SME - temp
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


## Create covariate list
# rasters must be in a list for mcapply
r.list <- list(slope, ndviCoV, ag, cover20, cover2070, cover70, drains, 
               prop.settlement.250, prop.settlement.1500, gHM, 
               roadsPrimary, roadsSecondary, dist2roads.primary, dist2roads.secondary, dist2ag)
## Extract
study.area <- 32736

## Extraction in parallel - 40 minutes
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
    extract <- mclapply(r.list, function(x)
      terra::extract(x, vect)[,2], mc.cores = cores)
    used[[i]] <- do.call(cbind, extract) # bind results into an extracted dataframe for each individual
  }
toc()}

used2 <- do.call(rbind, used)

# check
head(used2)
summary(used2)

# create data frame
mode(used2) = "numeric"
used2 <- as.data.frame(used2)
colnames(used2) <- c('slope','ndviCoV','ag','cover20','cover2070','cover70','drains','prop.settlement.250','prop.settlement.1500','gHM','roadsPrimary','roadsSecondary','dist2roads.primary','dist2roads.secondary','dist2ag')
head(used2)

# unstandardized data frame
rsf.ext <- cbind(rsf.df, used2)

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

write.csv(t, './SSF/ag.homerange.IY.datatable_20220312.csv', row.names = FALSE)

# read in ag available table
t <- read.csv('./SSF/ag.homerange.IY.datatable_20220312.csv')

rsf.ext <- merge(rsf.ext, t, by = 'subject_year') 

# save extracted df
write.csv(rsf.ext, './SSF/RSF_dataset_extracted_20220325.csv')

# drop old dataframes
remove(extract)
remove(used)
remove(rsf.df)
remove(rsf.sf)
remove(r.list)
remove(used2)

##### Save Point #####
# Can make the following code a seperate file for model fitting
rsf.ext <- as.data.frame(data.table::fread('./SSF/RSF_dataset_extracted_20220325.csv'))
rsf.ext$subject_year <- as.character(rsf.ext$subject_year)
rsf.ext$tactic.season <- as.factor(rsf.ext$tactic.season)
rsf.ext$V1 <- NULL


## Prep Dataframe

## scale some covariates before we fit it
rsf.ext <- filter(rsf.ext, ndviCoV > -13 & ndviCoV < 54)

rsf.ext <- rsf.ext %>%
  mutate_at(.vars = c('slope', 'ndviCoV'), .funs = scale) %>%
  #mutate_at(.vars = c('gHM','prop.settlement.250','prop.settlement.1500'), .funs = scale)
  mutate_at(.vars = c('ag','cover20','cover2070','cover70','drains','roadsPrimary','roadsSecondary'), .funs = function(x) if_else(!is.na(x), 1, 0)) %>%
  mutate_at(.vars = c('ag','cover20','cover2070','cover70', 'drains','roadsPrimary','roadsSecondary'), .funs = as.factor)

rsf.ext$case_ <- as.factor(rsf.ext$case_)
rsf.ext$subject_name <- as.factor(rsf.ext$subject_name)

# filter out NDVI cov outliers
rsf.ext <- filter(rsf.ext, ndviCoV < 20) # removes three major outliers from Marima dataset
t <- filter(rsf.ext, ndviCoV > -6) # removes major outliers from Serengeti eles
# check dataframe
summary(rsf.ext)

## filter individuals (>1000 points)
t <- rsf.ext %>% group_by(subject_year) %>% 
  filter(case_ == 'TRUE') %>%
  tally() %>%
  filter(n > 1000) # set min number of points for subject_name or subject_year

rsf.ext<- filter(rsf.ext, subject_year %in% t$subject_year)
rsf.ext <- filter(rsf.ext, !(subject_year %in% c('Caroline-2014','Tressa-2015','Chelsea-2016','Jinomoja-2018','Jinomoja-2019','Jinomoja-2020')))


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

{tic()
  split <- split(rsf.ext, rsf.ext$subject_year) # subject_name or subject_year
  m.ind.rsf <- lapply(split, function(x)
    glm(case_ ~ ag + cover20 + cover2070 + cover70 + slope + ndviCoV + drains + scale(prop.settlement.250) + 
          scale(dist2roads.primary) + scale(dist2roads.secondary), data = x, family = 'binomial'))
toc()}

remove(split)
#saveRDS(m.ind.rsf, paste0('./SSF/RSF_modfit_', region,'_20220325.RDS'))

# get table of estimates for each individual
id.est.rsf <- lapply(m.ind.rsf, function(x) 
  tab <- x %>% broom::tidy() %>%
    mutate(
      ymin = estimate - 1.96 * std.error,
      ymax = estimate + 1.96 * std.error
    ))
id.est.rsf <- dplyr::bind_rows(id.est.rsf, .id = "subject_year")


# add other ele metadata (sex, ag tactic, etc.)
t <- rsf.ext %>% group_by(subject_year, subject_name, subject_sex, tactic.season, ag.avail) %>% tally() %>% dplyr::select(-n)
id.est.rsf <- merge(id.est.rsf, t, by = 'subject_year')

#saveRDS(id.est.rsf, paste0('./SSF/RSF_estimates_', region, '_20220325.RDS'))

# remove intercept
id.est.rsf <- filter(id.est.rsf, term != '(Intercept)')
id.est.rsf <- filter(id.est.rsf, term != 'ag1')
##### plot results #####

# coeff estimates for all individuals
#tt <- filter(id.est.rsf, subject_name %in% c('Ivy','Lucy', 'Chelsea', 'Alina'))
ggplot(id.est.rsf, aes(x = term, y = estimate, color = tactic.season)) + 
  geom_pointrange(aes(ymin = ymin, ymax = ymax),
                  position = position_dodge(width = 0.6)) +
  xlab('covariate') + ggtitle('rsf individual-level estimates') +
  labs(color = 'ag tactic') + 
  coord_flip() #+ facet_wrap(~subject_name)

# coeff estimates by sex
ggplot(id.est.rsf, aes(x = as.factor(term), y = estimate, color = subject_sex)) +
  geom_boxplot() +
  xlab('covariate') + ggtitle('rsf individual-level estimates by sex') +
  labs(color = 'sex') + 
  coord_flip()

# coeff estimates by tactic (IY)
ggplot(id.est.rsf, aes(x = as.factor(term), y = estimate, color = factor(tactic.season))) +
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
sjPlot::plot_model(m.funct, type = 'pred', show.data = TRUE, title = '', axis.labels = '', dot.size = 1.5)

