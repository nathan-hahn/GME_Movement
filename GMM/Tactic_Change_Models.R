##### Tactic Change Model #####
# Fit and evaluate tactic change models
# Calculates tactic change and logistic model candidates for eval

library(tidyverse)
library(lubridate)
library(adehabitatHR)
library(raster)

source('GME_functions.R')

##### Prep Data #####
# Import cluster-classified data and cluster model object
movdata <- readRDS('./GMM/movdata/GMEcollars_003_usedClust_2020-10-30.rds')
movdata$individual.year <- paste(movdata$subject_name, movdata$year.cuts, sep = '.')
mod <- readRDS('./GMM/results/mSelect_2020-10-30.rds')
dist2ag <- raster("./spatial data/dist2ag_estes_32736_2019-11-21.tif")

##### Home Range Calcs #####
## Build homeranges for each individual-year and extract the areas

# prep spdf
t <-  filter(movdata, !is.na(x)) %>%
  dplyr::select(x, y, subject_name, year.cuts) 
coordinates(t) <-~x+y

# build MCP's by individual + year
split <- split(t, list(t$subject_name, t$year.cuts), drop = TRUE)
mcp <- lapply(split, mcp)

# extract areas and attach to movdata
split.area <- list()
for(i in 1:length(split)){
  split.area[[i]] <- cbind(split[[i]], newcol = mcp[[i]]$area)
  colnames(split.area[[i]]@data) <- c("subject_name", "year.cuts", "year.mcp.area")
}
mcp.areas <- do.call(rbind, split.area)
movdata$year.mcp.area <- mcp.areas@data$year.mcp.area


## Calculate centroid to ag distance
# get centroids to check distance to ag
centroids <- lapply(mcp, getSpPPolygonsLabptSlots)
centroids <- as.data.frame(do.call(rbind, centroids))
centroids$individual.year <- names(mcp)
colnames(centroids) <- c('x', 'y', 'individual.year')

# # create spatial points 
locs <- SpatialPointsDataFrame(as.matrix(centroids[c("x","y")]), data = centroids, 
                               proj4string = crs('+proj=utm +init=epsg:32736'))


centroids$centroid.dist <- extract(dist2ag, locs)
centroids$x <- NULL
centroids$y <- NULL
movdata <- inner_join(movdata, centroids, by = c("individual.year"))


## Calculate mean daily displacement
movdata$ymd <- as.Date(movdata$date, "%Y/%m/%d")
stat <- function(x) c(sum = sum(x))
ag.daily <- as.data.frame(aggregate(dist ~ ymd + subject_name, movdata, stat)) # aggregated stats by day

# bind daily displacement to full dataset
displacement <- inner_join(movdata, ag.daily, by = c("ymd", "subject_name")) %>%
  rename(daily.disp = dist.y) %>%
  mutate(tactic.season = as.integer(tactic.season), year.cuts = as.character(year.cuts)) %>%
  group_by(subject_name, year.cuts, tactic.season) %>%
  summarise(mu.daily.disp = mean(daily.disp, na.rm = TRUE)) %>%
  group_by(subject_name) %>%
  mutate(m.lag = dplyr::lag(mu.daily.disp, n = 1)) 

movdata <- inner_join(movdata, displacement, by = c("subject_name", "year.cuts", "tactic.season"))


##### Cut Points #####
# get cutpoints of the GMM to calculate duration over thresholds
cut.points <- vector()
cut.points[1] <- cluster_cutpoint(mod, comp = c(1,2), grph=T)
cut.points[2] <- cluster_cutpoint(mod, comp = c(2,3), grph = T)
cut.points[3] <- cluster_cutpoint(mod, comp = c(3,4), grph = T)
cut.points


##### Ag Use Variation ######
# get the max, min, amplitude values for each individual-year -- from moving window
amplitude <- movdata %>%
  filter(!is.na(mean)) %>%
  group_by(subject_name, year.cuts) %>%
  mutate(year.mean = mean(ag.used, na.rm = TRUE),
         year.max = max(mean, na.rm = TRUE),
         year.min = min(mean, na.rm = TRUE),
         year.amp = year.max - year.min,
         duration = sum(mean >= cut.points[3]), # hours spent above habitual threshold
         year.begin = min(date),
         year.end = max(date)) %>%
  group_by(subject_name, tactic.agg, site, year.cuts, tactic.season, year.begin, year.end, year.mean, year.max, year.min, year.amp, duration,
           year.mcp.area, mu.daily.disp, subject_sex, subject_ageClass, m.lag, centroid.dist) %>%
  # filter to individual years with at least 1 month's worth of fixes (does not account for NAs)
  tally() %>% filter(n > 500) %>% droplevels()
#View(amplitude)



##### Temporal Tactic Changes ######

change <- amplitude %>%
  #filter(subject_name == "Olchoda") %>%
  group_by(subject_name, tactic.agg, subject_sex, subject_ageClass) %>%
  mutate(tactic.prev = lag(tactic.season)) %>%
  mutate(tactic.change = if_else(tactic.season != dplyr::lag(tactic.season), 1, 0)) %>%
  mutate(tactic.direction = tactic.season - dplyr::lag(tactic.season)) %>%
  dplyr::select(subject_name, year.cuts, tactic.season, tactic.prev, tactic.change, tactic.direction) %>%
  drop_na() %>%
  droplevels()


## Tactic Change Rates

change.sum <- change %>% 
  filter(subject_name %in% change$subject_name) %>%
  group_by(tactic.season) %>%
  summarise(year.n = length(year.cuts),
            change = sum(tactic.change)/year.n,
            sd = sd(tactic.change),
            se = sd/sqrt(year.n),
            lwr = change - se,
            upr = change + se) %>%
  mutate(tactic.season = recode_factor(as.factor(tactic.season),
                                       "1" = "Rare",
                                       "2" = "Sporadic",
                                       "3" = "Seasonal",
                                       "4" = "Habitual"))

                                                                              
head(change.sum)

# plot
ggplot(change.sum, aes(factor(tactic.season), change)) + geom_pointrange(
  aes(ymin = lwr, ymax = upr)) + 
  xlab("Tactic Class") + ylab("Tactic Change Rate")

# plot by sex
change.sex <- change %>% 
  filter(subject_name %in% change$subject_name) %>%
  group_by(subject_sex, tactic.season) %>%
  summarise(year.n = length(year.cuts),
            change = sum(tactic.change)/year.n,
            sd = sd(tactic.change),
            se = sd/sqrt(year.n),
            lwr = change - se,
            upr = change + se) %>%
  mutate(tactic.season = recode_factor(as.factor(tactic.season),
                                      "1" = "Rare",
                                      "2" = "Sporadic",
                                      "3" = "Seasonal",
                                      "4" = "Habitual"))

# plot
ggplot(change.sex, aes(factor(tactic.season), change, color = subject_sex)) + geom_pointrange(
  aes(ymin = lwr, ymax = upr), position = position_dodge(0.4)) + 
  xlab("Tactic Class") + ylab("Tactic Change Rate")

# store tactic change info 
change.df <- amplitude %>%
  group_by(subject_name, tactic.agg, subject_sex, subject_ageClass) %>%
  mutate(tactic.prev = lag(tactic.season)) %>%
  mutate(tactic.change = if_else(tactic.season != dplyr::lag(tactic.season), 1, 0)) %>%
  mutate(tactic.direction = tactic.season - dplyr::lag(tactic.season)) %>%
  #dplyr::select(subject_name, year.cuts, tactic.season, tactic.prev, tactic.change, tactic.direction) %>%
  #drop_na() %>%
  droplevels()

#write.csv(change.df, './GMM/results/individual_year_classification 20210216.csv')

## Are ele characteristics driving tactic choice? 
## Model the current years tactic as a function of prev. year tactic, peak NDVI, home range, ageClass, sex
library(lme4)
library(nnet)
library(MuMIn)

mod.df <- change.df %>%
  ungroup() %>%
  
  # format vectors
  mutate_at(c("tactic.prev", "tactic.change", "tactic.season", "subject_sex", "subject_ageClass"), as.factor) %>%
  mutate(subject_ageClass = recode_factor(subject_ageClass, 
                                          "adult" = "young adult",
                                          "mature" = "old adult")) %>%
  dplyr::select(subject_name, tactic.change, subject_sex, subject_ageClass, year.mcp.area, mu.daily.disp, m.lag, tactic.season, tactic.prev) %>%
  
  # drop years with no previous tactic (first years)
  drop_na() %>% droplevels() %>%
  
  # log homerange data
  mutate(year.mcp.area = log(year.mcp.area))


# fit models
m1<- glmer(tactic.change ~ tactic.prev + subject_sex*subject_ageClass + (1|subject_name), data = mod.df, family = 'binomial')
m2 <- glmer(tactic.change ~ tactic.prev + subject_sex + subject_ageClass + (1|subject_name), data = mod.df, family = 'binomial')
m3 <- glmer(tactic.change ~ tactic.prev + (1|subject_name), data = mod.df, family = 'binomial')
m4 <- glmer(tactic.change ~ subject_sex + subject_ageClass + (1|subject_name), data = mod.df, family = 'binomial')
m5 <- glmer(tactic.change ~ subject_sex*subject_ageClass + (1|subject_name), data = mod.df, family = 'binomial')
m6 <- glmer(tactic.change ~ subject_sex*subject_ageClass + year.mcp.area + (1|subject_name), data = mod.df, family = 'binomial')
m7 <- glmer(tactic.change ~ subject_sex*subject_ageClass + year.mcp.area + tactic.prev + (1|subject_name), data = mod.df, family = 'binomial')
m8 <- glmer(tactic.change ~ year.mcp.area + tactic.prev + (1|subject_name), data = mod.df, family = 'binomial')
m9 <- glmer(tactic.change ~ subject_sex*subject_ageClass + log(mu.daily.disp) + (1|subject_name), data = mod.df, family = 'binomial')
m10 <- glmer(tactic.change ~ subject_sex*subject_ageClass + log(mu.daily.disp) + year.mcp.area + (1|subject_name), 
             data = mod.df, family = 'binomial')
# AIC table w/ likelihoods
knitr::kable(AICc(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10))

aicc <- AICc(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10)
mods <- list(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10)
t <- lapply(mods, logLik)
t <- do.call(rbind, t)

AIC.table <- cbind(aicc, t)
AIC.table %>% arrange(AICc)


# coeff + CIs for top model
summary(m6)
confint(m6)

# diagnostics
hist(resid(m6))

## Get odd ratios from top model
# odds <- exp(summary(m6)$coefficients)
# odds.CI <- exp(confint(m6))
# odds <- cbind(OR = odds[,1], odds.CI[2:6,])

## Plot sex/age interaction using prediction from top model - probabilities or odds
library(DescTools)
plot.dat <- m6@frame %>% 
  mutate(prob = predict(m6, newdata = ., type = "response")) %>%
  group_by(subject_sex, subject_ageClass) %>%
  summarise(n = n(),
            mean = MeanCI(prob, method = "boot", type = "norm", R=1000)[1],
            lwr.ci = MeanCI(prob, method = "boot", type = "norm", R=1000)[2],
            upr.ci = MeanCI(prob, method = "boot", type = "norm", R=1000)[3]) 

# plot.dat.odds <- plot.dat %>%
#   mutate(mean = exp(mean), lwr.ci = exp(lwr.ci), upr.ci = exp(upr.ci))
# plot.dat.odds

# plot
levels(plot.dat$subject_sex) <- c('Female', 'Male')
levels(plot.dat$subject_ageClass) <- c('Young adult', 'Mature adult')

ggplot(plot.dat, aes(x = subject_ageClass, y = mean, shape = subject_sex)) + 
  geom_pointrange(aes(ymin = lwr.ci, ymax = upr.ci)) + 
  geom_line(aes(group = subject_sex), linetype = "dashed") + 
  xlab("Age Class") + ylab("Probability of Switching Tactics") + labs(shape = "Sex")



