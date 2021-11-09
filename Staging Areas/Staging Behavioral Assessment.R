#### Staging Use Stats ####
library(tidyverse)
library(lubridate)

## Load the data
gme.stage <- as.data.frame(data.table::fread('./movdata/GMEcollars_004_stageclassified_20211019.csv'))

gme.stage <- gme.stage %>%
  dplyr::select(uid, subject_name, year.cuts, tactic.season, tactic.agg, x, y, date, dist, id, subject_sex, region, 
                subject_ageClass, ag.used, viterbi, ag.window, dayBurst, ag.window.ext, vote)

head(gme.stage)

#' ##### Staging Area Data Summaries

#' Check the hourly distribution of staging relocations

# check distribution of stage relocs over the course of the day
t <- filter(gme.stage, vote == 1) %>%
  group_by(hour(date)) %>% tally()
ggplot(t, aes(as.factor(`hour(date)`), n)) + geom_bar(stat = 'identity') +
  xlab('hour of day') + ylab('n stage relocations')

#' The total number of stage events in the dataset is 5692.

## index ag days
eventIndex <- ifelse(gme.stage$ag.window.ext == 1, gme.stage$dayBurst, 0)
# assign event index
gme.stage$ag.day.index <- eventIndex

## index staging events
eventFlag <- ifelse(gme.stage$vote == 1, TRUE, FALSE)
eventIndex <- inverse.rle(within.list(rle(eventFlag),
                                      values[values] <- seq_along(values[values])))
gme.stage$stage.event.index <- eventIndex

length(unique(gme.stage$stage.event.index))

#' Calculate distribution of staging frequency prior to crop raiding among individuals. 

# summarise overall % staging distribution
stage.summary <- gme.stage %>%
  group_by(subject_name) %>%
  summarise(n.stage = length(unique(stage.event.index)),
            n.raid = length(unique(ag.day.index)),
            pct.stage = n.stage/n.raid) 

summary(stage.summary$pct.stage)

#' What is the distribution of staging lengths? Is capped at 12 hrs by parameter space

length <- gme.stage %>%
  filter(vote == 1) %>%
  group_by(stage.event.index) %>%
  tally()

boxplot(length$n)

hist(length$n)

#' Summarise staging frequency distribution by individual and tactic. Use yearly tactics (n=101) 
#' Note the Rare group with 100% stage percentages have very small numbers of 
#' stages. See the end of the script for data summary.

# summarise % staging distribution by individual-year tactic
stage.summary <- gme.stage %>%
  group_by(subject_name, tactic.season, year.cuts) %>%
  summarise(n.stage = length(unique(stage.event.index)),
            n.raid = length(unique(ag.day.index)),
            pct.stage = n.stage/n.raid) %>%
  filter(n.raid > 5)

boxplot(stage.summary$pct.stage ~ stage.summary$tactic.season,
        main = 'distribution of staging frequency by individuals in the 4 tactics',
        xlab = 'Tactics (Rare to Habitual)', ylab = 'Pct. of Raids with a Stage')

#' Summarise inter-annual variation in staging by individuals - coefficient of variation

cv.summary <- stage.summary %>%
  group_by(subject_name) %>%
  summarise(cv = raster::cv(pct.stage)/100)
hist(cv.summary$cv, main = 'Individual Inter-Annual Variation in Staging Pct', xlab = 'Coefficient of Variation')  

#' How does the length of the ag bout affect staging frequency? 
#' 1. Calculate unique id for ag.windows
#' 2. Calculate the number of stage events within the ag.window events
#' 3. Plot relationship between ag.window length in days (x) and number of stage events (y) to see if it's a linear relationship. Plot as boxplot

# 1. index ag days

eventIndex <- gme.stage %>%
  # group by id to avoid merges between individuals
  group_by(id) %>%
  # add flag and unique number for each ag day - note the number is only unique to the individual at this point
  mutate(eventFlag = if_else(ag.window.ext == 1, TRUE, FALSE)) %>%
  mutate(eventIndex = inverse.rle(within.list(rle(eventFlag),
                                                      values[values] <- seq_along(values[values])))) %>%
  # create unique index by pasting the event number and individual id
  mutate(eventIndex = paste(id, eventIndex, sep = '-')) %>%
  select(id, ag.window.ext, eventFlag, eventIndex)

# assign event index
gme.stage$ag.window.ext.index <- eventIndex$eventIndex


# 2. # stage events within ag.window events
window <- gme.stage %>%
  filter(ag.window.ext == 1) %>%
  mutate(stage.event.index = na_if(stage.event.index,0)) %>%
  group_by(ag.window.ext.index) %>%
  summarise(window.length = round(difftime(max(date), min(date), units = 'days'),0),
    n.stage.window = length(unique(na.omit(stage.event.index))),
    stage.pct = n.stage.window/as.numeric(window.length),
    id = unique(id)) # exclude NAs, which are non-staging relocs
window
#window <- window[window$ag.window.index != 0,] # all non-ag days tagged as 0, so remove them
window$window.length <- as.numeric(window$window.length)

# check
#View(window)

# plot and correlation test
boxplot(window$stage.pct ~ window$window.length)

cor.test(window$stage.pct, window$window.length, method = 'spearman')


#' ##### Plot the data on a map

##### Plot Staging Events #####
library(sf)
gme <- sf::st_read('~/Dropbox (Personal)/CSU/GME_Movement/spatial data/GSE/GSE_2020.shp') 
stage.relocs <- filter(gme.stage, vote == 1)

ggplot(data = gme) + geom_sf() + coord_sf(datum=st_crs(32736)) + 
  geom_point(data = stage.relocs, aes(x, y, color = 'red'), size = 0.2, alpha = 0.2) +
  labs(color = 'staging relocations')

#' Table of staging events and ag bouts for each individual-year
stage.summary <- gme.stage %>%
  group_by(subject_name, tactic.season, year.cuts) %>%
  summarise(n.stage = length(unique(stage.event.index)),
            n.raid = length(unique(ag.day.index)),
            pct.stage = n.stage/n.raid,
            )
(stage.summary)

#' ##### Model Staging Frequency #####

#' % ag use for the year, age, sex, hr size, daily displacement, avg. length ag bout for year, (1|ID)

## Homerange - MCP
## Build homeranges for each individual-year and extract the areas
library(adehabitatHR)
# prep spdf
t <-  filter(gme.stage, !is.na(x)) %>%
  dplyr::select(x, y, subject_name, year.cuts) 
coordinates(t) <-~x+y

# build MCP's by individual + year
split <- split(t, list(t$subject_name, t$year.cuts), drop = TRUE)
mcp <- lapply(split, adehabitatHR::mcp, percent=95)

# extract areas and attach to movdata
split.area <- list()
for(i in 1:length(split)){
  split.area[[i]] <- cbind(split[[i]], newcol = mcp[[i]]$area)
  colnames(split.area[[i]]@data) <- c("subject_name", "year.cuts", "year.mcp.area")
}
mcp.areas <- do.call(rbind, split.area) 
mcp.areas <- mcp.areas@data %>% 
  group_by(subject_name, year.cuts) %>%
  summarise(year.mcp.area = unique(year.mcp.area)) 

gme.stage <- inner_join(gme.stage, mcp.areas, by = c('subject_name', 'year.cuts'))


## Get daily displacement
gme.stage <- gme.stage %>%
  group_by(dayBurst) %>%
  mutate(dailyDist = sum(dist)) %>% ungroup()

## Get daily displacement only during ag use
ag.mu.dailyDist <- gme.stage %>%
  group_by(subject_name, year.cuts, ag.window.ext) %>%
  summarise(ag.mu.dailyDist = mean(dailyDist)) %>%
  filter(ag.window.ext == 1) %>%
  dplyr::select(-ag.window.ext)

gme.stage <- inner_join(gme.stage, ag.mu.dailyDist, by = c('subject_name', 'year.cuts'))

## Mean bout length by year
bout.length <- gme.stage %>%
  filter(ag.window.ext == 1) %>%
  mutate(stage.event.index = na_if(stage.event.index,0)) %>%
  group_by(ag.window.ext.index) %>%
  summarise(window.length = round(difftime(max(date), min(date), units = 'days'),0),
            id = unique(id), year.cuts = unique(year.cuts)) %>%
  group_by(id, year.cuts) %>%
  summarise(mean.bout = median(as.numeric(window.length)))

gme.stage <- inner_join(gme.stage, bout.length, by = c('id', 'year.cuts'))

## Summarise all data
mod.df <- gme.stage %>%
  group_by(subject_name, year.cuts) %>%
  summarise(
    # add year-constant variables
    tactic.season = unique(tactic.season),
    subject_sex = unique(subject_sex),
    subject_ageClass = unique(subject_ageClass),
    mean.bout = unique(mean.bout),
    year.mcp.area = unique(year.mcp.area),
    # calc year variables
    mean.ag = mean(ag.used),
    mu.dailyDist = mean(dailyDist),
    sd.dailyDist = sd(dailyDist),
    ag.mu.dailyDist = mean(ag.mu.dailyDist),
    n.stage = length(unique(stage.event.index)),
    n.raid = length(unique(ag.day.index)),
    #n.day = length(unique(dayBurst)),
    pct.stage = n.stage/n.raid,
  ) %>% 
  ungroup()
mod.df <- as.data.frame(mod.df)

## Format data
mod.df$subject_name <- as.factor(mod.df$subject_name)
mod.df$subject_sex <- as.factor(mod.df$subject_sex)
mod.df$subject_ageClass <- as.factor(mod.df$subject_ageClass)
mod.df$mean.bout <- as.numeric(mod.df$mean.bout)

mod.df <- mod.df[mod.df$n.raid >= 5,] # 5 stages is the first quartile

## Fit Model ##
library(lme4)
mod.global <- lmer(pct.stage ~ subject_sex + subject_ageClass + mean.bout + log(year.mcp.area) + log(mu.dailyDist) + mean.ag + (1|subject_name),
           data = mod.df, REML = FALSE)

summary(mod.global)
plot(mod.global)

# sex/age model
mod.1 <- lmer(pct.stage ~ subject_sex + subject_ageClass + (1|subject_name), data = mod.df, REML = FALSE)
summary(mod.1)
plot(mod.1)

# movement behaviors model
mod.2 <- lmer(pct.stage ~ log(mu.dailyDist)+log(ag.mu.dailyDist)+log(year.mcp.area)+(1|subject_name), data = mod.df, REML = FALSE)
summary(mod.2)
plot(mod.2)

# ag use behaviors model
mod.3 <- lmer(pct.stage ~ mean.ag + mean.bout + (1|subject_name), data = mod.df, REML = FALSE)
summary(mod.3)
plot(mod.3)

# ag use + movement model
mod.4 <- lmer(pct.stage ~ mean.ag + mean.bout + log(year.mcp.area) + log(mu.dailyDist) + log(ag.mu.dailyDist) + (1|subject_name), data = mod.df, REML = FALSE)
summary(mod.4)

# ag use + sex/age model
#mod.5 <- lmer(pct.stage ~ mean.ag + mean.bout + subject_sex + subject_ageClass + (1|subject_name), data = mod.df, REML = FALSE)

# combo best performers
mod.5 <- lmer(pct.stage ~ mean.ag + subject_sex + subject_ageClass + (1|subject_name), data = mod.df, REML = FALSE)
mod.6 <- lmer(pct.stage ~ mean.ag + log(year.mcp.area) + log(mu.dailyDist) + log(ag.mu.dailyDist) + (1|subject_name), data = mod.df, REML = FALSE)
mod.7 <- lmer(pct.stage ~ log(mu.dailyDist) + (1|subject_name), data = mod.df, REML = FALSE)


library(MuMIn)
# AICc
t <- as.data.frame(AICc(mod.global, mod.1, mod.2, mod.3, mod.4, mod.5, mod.6, mod.7))
t$deltaAICc <- qpcR::akaike.weights(t$AICc)$deltaAIC
t$AICc_weight <- round(qpcR::akaike.weights(t$AICc)$weights, 2)
t$model_likelihood <- round(qpcR::akaike.weights(t$AICc)$rel.LL, 2)
t$loglik <- c(logLik(mod.global), logLik(mod.1), logLik(mod.2), logLik(mod.3), logLik(mod.4), logLik(mod.5), logLik(mod.6), logLik(mod.7))
t$mod.formula <- c(mod.global@call$formula, mod.1@call$formula, mod.2@call$formula, mod.3@call$formula, mod.4@call$formula, 
 mod.5@call$formula, mod.6@call$formula, mod.7@call$formula)
t <- t[order(t$AICc),]
t

View(t)

# Final table
t$mod.formula <- as.character(t$mod.formula)
write.csv(t, './Staging Areas/stage_behavior_modtable.csv')

## Deeper eval of mod.6

#'
plot(mod.2) 
summary(mod.2)

confint(mod.2)

sjPlot::tab_model(mod.2)
sjPlot::plot_model(mod.2)




