#### Staging Use Stats ####
library(tidyverse)
library(lubridate)

Sys.setenv(TZ='Africa/Nairobi')

## Load the data
#gme.stage <- as.data.frame(data.table::fread('./movdata/GMEcollars_004_stageclassified_20211019.csv', tz=''))
gme.stage <- as.data.frame(data.table::fread('./Staging Areas/movdata/GMEcollars_004_stageclassified_20220109.csv', tz=''))


gme.stage <- gme.stage %>%
  dplyr::select(uid, subject_name, year.cuts, tactic.season, tactic.agg, x, y, date, dist, id, subject_sex, region, 
                subject_ageClass, ag.used, viterbi, ag.window, dayBurst, ag.window.ext, vote.ag, site, pa)

head(gme.stage)

# correct pardamat ag to non-ag
gme.stage$ag.used <- ifelse(gme.stage$site == 'mep' & gme.stage$pa == 2, 0, gme.stage$ag.used)

#' ##### Staging Area Data Summaries

#' Check the hourly distribution of staging relocations

# check distribution of stage relocs over the course of the day
t <- filter(gme.stage, vote.ag == 1) %>%
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
  group_by(subject_name, subject_sex) %>%
  summarise(n.stage = length(unique(stage.event.index)),
            n.raid = length(unique(ag.day.index)),
            pct.stage = n.stage/n.raid) 

summary(stage.summary$pct.stage)
tapply(stage.summary$pct.stage,stage.summary$subject_sex, summary)

hist(stage.summary$pct.stage, xlab = 'staging percentage (individual years)', main = 'Distribution of staging percentage')

#' What is the distribution of staging lengths? Is capped at 12 hrs by parameter space

length <- gme.stage %>%
  filter(vote.ag == 1) %>%
  group_by(stage.event.index) %>%
  tally()

boxplot(length$n)

hist(length$n, xlab = 'staging length (hours)', main = 'Histogram of staging length', breaks = 12)

summary(length$n)

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
  dplyr::select(id, ag.window.ext, eventFlag, eventIndex)

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

# plot
boxplot(window$stage.pct ~ window$window.length)

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

remove(t)

## Get daily displacement
gme.stage <- gme.stage %>%
  group_by(dayBurst) %>%
  mutate(dailyDist = sum(dist)) %>% ungroup()

## Get daily displacement only during ag use
ag.mu.dailyDist <- gme.stage %>%
  #filter(vote != 1) %>% # filter out staging events to prevent conflation of staging with movement
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
    weights = n.raid
  ) %>% 
  ungroup()
mod.df <- as.data.frame(mod.df)

## Format data
mod.df$subject_name <- as.factor(mod.df$subject_name)
mod.df$subject_sex <- as.factor(mod.df$subject_sex)
mod.df$subject_ageClass <- as.factor(mod.df$subject_ageClass)
mod.df$mean.bout <- as.numeric(mod.df$mean.bout)

#mod.df <- mod.df[mod.df$n.raid >= 5,] # 5 stages is the first quartile

## Fit Model ##
library(lme4)
mod.global <- lmer(pct.stage ~ subject_sex + subject_ageClass + mean.bout + log(year.mcp.area) + log(mu.dailyDist) + mean.ag + (1|subject_name),
           data = mod.df, REML = FALSE)

summary(mod.global)
plot(mod.global)

mod.global.bi <- glmer(pct.stage ~ subject_sex + subject_ageClass + mean.bout + log(year.mcp.area) + log(mu.dailyDist) + mean.ag + (1|subject_name), 
                  data = mod.df,
                  weights = weights,
                  family = binomial)

# sex/age model
mod.1 <- lmer(pct.stage ~ subject_sex + subject_ageClass + (1|subject_name), data = mod.df, REML = FALSE)
summary(mod.1)
plot(mod.1)

mod.1.bi <- glmer(pct.stage ~ subject_sex + subject_ageClass + (1|subject_name), 
                  data = mod.df,
                  weights = weights,
                  family = binomial)

# movement behaviors model
mod.2 <- lmer(pct.stage ~ log(mu.dailyDist)+log(ag.mu.dailyDist)+log(year.mcp.area)+(1|subject_name), data = mod.df, REML = FALSE)
summary(mod.2)
plot(mod.2)

mod.2.bi <- glmer(pct.stage ~ log(mu.dailyDist)+log(ag.mu.dailyDist)+log(year.mcp.area)+(1|subject_name), 
                  data = mod.df,
                  weights = weights,
                  family = binomial)

# ag use behaviors model
mod.3 <- lmer(pct.stage ~ mean.ag + mean.bout + (1|subject_name), data = mod.df, REML = FALSE)
summary(mod.3)
plot(mod.3)

mod.3.bi <- glmer(pct.stage ~ mean.ag + mean.bout + (1|subject_name),
                  data = mod.df,
                  weights = weights,
                  family = binomial)

# ag use + movement model
mod.4 <- lmer(pct.stage ~ mean.ag + mean.bout + log(year.mcp.area) + log(mu.dailyDist) + log(ag.mu.dailyDist) + (1|subject_name), data = mod.df, REML = FALSE)
summary(mod.4)

mod.4.bi <- glmer(pct.stage ~ mean.ag + mean.bout + log(year.mcp.area) + log(mu.dailyDist) + log(ag.mu.dailyDist) + (1|subject_name),
                  data = mod.df,
                  weights = weights,
                  family = binomial)

# ag use + sex/age model
#mod.5 <- lmer(pct.stage ~ mean.ag + mean.bout + subject_sex + subject_ageClass + (1|subject_name), data = mod.df, REML = FALSE)

# combo best performers
mod.5 <- lmer(pct.stage ~ mean.ag + subject_sex + subject_ageClass + (1|subject_name), data = mod.df, REML = FALSE)
mod.5.bi <- glmer(pct.stage ~ mean.ag + subject_sex + subject_ageClass + (1|subject_name),
                  data = mod.df,
                  weights = weights,
                  family = binomial)


mod.6 <- lmer(pct.stage ~ mean.ag + log(year.mcp.area) + log(mu.dailyDist) + log(ag.mu.dailyDist) + (1|subject_name), data = mod.df, REML = FALSE)
mod.6.bi <- glmer(pct.stage ~ mean.ag + log(year.mcp.area) + log(mu.dailyDist) + log(ag.mu.dailyDist) + (1|subject_name),
                  data = mod.df,
                  weights = weights,
                  family = binomial)

mod.7 <- lmer(pct.stage ~ log(mu.dailyDist) + (1|subject_name), data = mod.df, REML = FALSE)
mod.7.bi <- glmer(pct.stage ~ log(mu.dailyDist) + (1|subject_name),
                  data = mod.df,
                  weights = weights,
                  family = binomial)

library(MuMIn)
## NOTE: global model and mod.4 failed to converge for the binomial weighted models. Removed from model table calculations

# AICc
t <- as.data.frame(AICc(mod.1.bi, mod.2.bi, mod.3.bi, mod.5.bi, mod.6.bi, mod.7.bi))
t$deltaAICc <- qpcR::akaike.weights(t$AICc)$deltaAIC
t$AICc_weight <- round(qpcR::akaike.weights(t$AICc)$weights, 2)
t$model_likelihood <- round(qpcR::akaike.weights(t$AICc)$rel.LL, 2)
t$loglik <- c(logLik(mod.1.bi), logLik(mod.2.bi), logLik(mod.3.bi), logLik(mod.5.bi), logLik(mod.6.bi), logLik(mod.7.bi))
t$mod.formula <- c(mod.1.bi@call$formula, mod.2.bi@call$formula, mod.3.bi@call$formula, 
                   mod.5.bi@call$formula, mod.6.bi@call$formula, mod.7.bi@call$formula)
t <- t[order(t$AICc),]
t

View(t)

# summary 
summary(mod.6.bi)

sjPlot::plot_model(mod.6.bi)
sjPlot::tab_model(mod.6.bi)

# Final table
t$mod.formula <- as.character(t$mod.formula)
write.csv(t, './Staging Areas/stage_behavior_modtable_binomialadjust_20220109.csv')


#### TODO Add Model Selection Based on Many Competing Models ####


#### TODO Try Beta Regression ####
# https://hansjoerg.me/2019/05/10/regression-modeling-with-proportion-data-part-1/



