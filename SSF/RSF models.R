##### Resource Selection Model (Third Order Selection) #####

library(amt)
library(tidyverse)
library(sjPlot)
library(terra)
library(sf)
library(tictoc)

set.seed(69)

##### 1a. Prepare Data #####

# read data
#movdata <- readRDS('./SSF/eledata_expanded.rds')
movdata <- readRDS('./SSF/eledata_allmara.RDS')
movdata <- movdata[!movdata$subject_name %in% c('Shamba','Courtney','David','Pepper'),]

# create date object
movdata$date <- as.POSIXct(movdata$date) # check still in EAT

# add new variable for individual-year (name-year) that can be used to group and fit models by. Years use april cutoff
movdata$subject_year <- paste(movdata$subject_name, movdata$year.cuts, sep = '-')

# make a move object for the amt package - we will use the id, which is a
#         combination of the animal's name and collar number
track <- make_track(movdata, .x = x, .y = y, .t = date, id = id, subject_name = subject_name)
t1 <- split(track, track$id)

##### 1b. Generate random samples #####
#' We will use amt functions to make random samples for our models. First, we will generate random points for the rsf
#' model using the `random_points` function. We will use an mcp home range estimator to generate the random points.
#' We are only using one individual, but the code is set up to scale with multiple individuals using the nest and unest
#' functions, as long as you identify the individuals using a column called `id`.

## make random points for RSF
# set cores
cores = 6

# create random points - 1.5hrs for maraonly dataset
{tic()
rand.factor = 5 # 5:1 unused:used
rsf.randpoints <- do.call(rbind, t1)
rsf.randpoints <- split(rsf.randpoints, rsf.randpoints$subject_name) # split variable (subject_name or subject_year)
rsf.randpoints <- parallel::mclapply(rsf.randpoints, function(x) {
  set.seed(69)
  rp <- x %>% amt::random_points(n = nrow(.)*rand.factor, hr = 'mcp') # 
  rp$subject_name = unique(x$subject_name) # update this by the split variable
  return(rp) }, mc.cores = cores)
toc()}

# create dataframe for used points
rsf.usedpoints <- lapply(t1, function(x)
  x %>% mutate(case_ = 'TRUE') %>% select(c(case_, x_, y_, subject_name))) # signify our used points using the amt 'case' nomenclature

# create a single dataframe
rsf.randpoints <- do.call(rbind, rsf.randpoints)
rsf.usedpoints <- do.call(rbind, rsf.usedpoints)
rsf.df <- rbind(rsf.randpoints, rsf.usedpoints)

t <- movdata %>% group_by(subject_name, subject_sex, tactic.agg) %>% tally() %>% dplyr::select(-c(n))
rsf.df <- inner_join(x = t, y = rsf.df, by = 'subject_name')

# remove old dataframes
remove(rsf.randpoints)
remove(rsf.usedpoints)

#saveRDS(rsf.df, './SSF/rsf.df.allmara.RDS')

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

# tiedman layer
# 1 = ag
# 2 = cover20
# 3 = cover2070
# 4 = cover70
# 5 = degraded
lc <- rast("./spatial data/Tiedman/agmask_reclass.tif")
ag <- terra::clamp(lc, lower = 1, upper = 1, values = FALSE)
cover20 <- terra::clamp(lc, lower = 2, upper = 2, values = FALSE)
cover2070 <- terra::clamp(lc, lower = 3, upper = 3, values = FALSE)
cover70 <- terra::clamp(lc, lower = 4, upper = 4, values = FALSE)

#drains
# get shapefile for drains buffer
drains <- st_read("./spatial data/drains/drains_estes_20211117/drains_estes_-2021-11-17.shp", 
                  layer="drains_estes_-2021-11-17", crs = 4326) %>%
  st_transform(crs = 32736) %>%
  filter(RIV_ORD <= 7)
# buffer by 250m
drains <- st_buffer(drains, dist = 250) %>%
  st_union() %>%
  terra::vect() %>%
  terra::rasterize(.,ag) # rasterize at 10m sentinel 


## Create covariate list
# rasters must be in a list for mcapply
r.list <- list(slope, ndviCoV, ag, cover20, cover2070, cover70, drains, prop.settlement.250, prop.settlement.1500, gHM)

## Extract

# create terra vect object
study.area <- 32736
rsf.sf <- rsf.df %>%
  st_as_sf(coords = c('x_','y_'), crs = study.area) %>%
  terra::vect() 

# extract in parallel - 212 seconds with allmara dataset
{tic()
  cores = 6
  extract <- mclapply(r.list, function(x)
    terra::extract(x, rsf.sf)[,2], mc.cores = cores)
  used <- do.call(cbind, extract) # bind results into an extracted dataframe
toc()}

# check
head(used)
summary(used)

# create data frame
mode(used) = "numeric"
used <- as.data.frame(used)
colnames(used) <- c('slope','ndviCoV','ag','cover20','cover2070','cover70','drains','prop.settlement.250','prop.settlement.1500','gHM')
head(used)

# unstandardized data frame
rsf.ext <- cbind(rsf.df, used)

# scale some covariates before we fit it
rsf.ext <- rsf.ext %>%
  mutate_at(.vars = c('slope', 'ndviCoV'), .funs = scale) %>%
  mutate_at(.vars = c('ag','cover20','cover2070','cover70','drains'), .funs = function(x) if_else(!is.na(x), 1, 0)) %>%
  mutate_at(.vars = c('ag','cover20','cover2070','cover70', 'drains'), .funs = as.factor)

# filter out NDVI cov outliers
rsf.ext <- rsf.ext[rsf.ext$ndviCoV < 20,] # removes three major outliers from Marima dataset

#' Our covariates are now extracted. Before model fitting, check the data frame to make sure we did everything 
#' correctly!

summary(rsf.ext)

# drop old dataframes
#remove(extract)
#remove(used)

##### Fit Population-Level RSF Model #####
#' To make things simple, we are going to forego model selection and fit the same model for each of the RSF, SSF, and 
#' iSSF so we can see how covariate estimates may differ at different orders of selection

# fit model
rsf.ext$case_ <- as.factor(rsf.ext$case_)
rsf.ext$subject_name <- as.factor(rsf.ext$subject_name)

m <- glm(case_ ~ ag + cover20 + cover2070 + cover70 + slope + ndviCoV + drains + gHM, data = rsf.ext, family = binomial)
summary(m)

m1 <- rsf.ext %>% fit_logit(case_ ~ ag + cover2070 + cover70 + slope + ndviCoV + drains + gHM)
m2 <- rsf.ext %>% fit_logit(case_ ~ ag + cover70 + slope + ndviCoV + gHM)
m3 <- rsf.ext %>% fit_logit(case_ ~ ag + prop.settlement.1500)
m4 <- rsf.ext %>% fit_logit(case_ ~ cover20 + cover2070 + cover70 + slope + ndviCoV + drains)
m5 <- rsf.ext %>% fit_logit(case_ ~ ag*cover70 + slope + ndviCoV + gHM)

# m.rsf <- glmer(case_ ~ ag + cover20 + cover2070 + cover70 + slope + ndviCoV + (1|subject_name),
#               data = rsf.ext, family = binomial)

library(MuMIn)
tab <- model.sel(m$model, m1$model, m2$model, m3$model, m4$model, m5$model)
View(tab)

# check out the summary of top model
summary(m)

# viz the covariates
#sjPlot::plot_model(m$model, transform = NULL, df_method='wald')


##### Fit Individual-level RSF #####
## Run in parallel 

# fit model
{tic()
  cores = 6
  m.ind.rsf <- split(rsf.ext, rsf.ext$subject_name)
  m.ind.rsf <- mclapply(m.ind.rsf, function(x) 
    glm(case_ ~ ag + cover20 + cover2070 + cover70 + slope + ndviCoV + drains + gHM, data = x, family = 'binomial'), 
    mc.cores = cores)
  names(m.ind.rsf) <- unique(rsf.ext$subject_name)
toc()}


# get table of estimates for each individual
id.est.rsf <- lapply(m.ind.rsf, function(x) 
  tab <- x %>% broom::tidy() %>%
    mutate(
      ymin = estimate - 1.96 * std.error,
      ymax = estimate + 1.96 * std.error
    ))
id.est.rsf <- dplyr::bind_rows(id.est.rsf, .id = "subject_name")

# add other ele metadata (sex, ag tactic, etc.)
t <- movdata %>% group_by(subject_name, subject_sex, tactic.agg) %>% tally() %>% dplyr::select(-n)
id.est.rsf <- merge(id.est.rsf, t, by = 'subject_name')

id.est.rsf <- filter(id.est.rsf, !subject_name %in% c('Courtney','David','Pepper','Chelsea'))

##### plot results #####

# coeff estimates for all individuals
ggplot(id.est.rsf, aes(x = as.factor(term), y = estimate, color = subject_name)) + 
  geom_pointrange(aes(ymin = ymin, ymax = ymax),
                  position = position_dodge(width = 0.3)) +
  xlab('covariate') + ggtitle('rsf individual-level estimates') +
  labs(color = 'subject_name') + 
  coord_flip()

# coeff estimates by sex
ggplot(id.est.rsf, aes(x = as.factor(term), y = estimate, color = subject_sex)) +
  geom_boxplot() +
  xlab('covariate') + ggtitle('rsf individual-level estimates by sex') +
  labs(color = 'sex') + 
  coord_flip()

# coeff estimates by tactic (aggregate)
ggplot(id.est.rsf, aes(x = as.factor(term), y = estimate, color = factor(tactic.agg))) +
  geom_boxplot() +
  xlab('covariate') + ggtitle('rsf individual-level estimates by ag tactic') +
  labs(color = 'ag tactic (aggregate)') + 
  coord_flip() 
