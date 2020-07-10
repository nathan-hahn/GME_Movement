##### GMM Model Visualization #####
# For visualizing and compiling GMM results. 
# Attach classifications to the dataframe
library(tidyverse)
library(lubridate)
library(mclust)


# Load cluster results table

movdat <- readRDS('./GMM/GMEcollars_002_res90_2020-07-09.rds') %>%
  mutate(subject_sex = if_else(subject_sex %in% c("male", "M"), "male", "female"))
  
mod <- readRDS('./GMM/results/modSelect.rds')

clust.result <- read.csv("./GMM/results/clustResult_GME_2020-07-09.csv")
clust.summary <- clust.result %>%
  group_by(ag.class.mean) %>%
  summarise(n = length(ag.class.mean),
            mean = mean(mean.occupancy),
            sd = sd(mean.occupancy),
            se = sd/sqrt(n),
            lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
            upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se)
(clust.summary) # overall


##### Predict tactics for each season #####
split <- split(movdat, as.character(movdat$year.cuts))
mod.predict <- list()
result <- list()
class <- list()

for(i in 1:length(split)) {
  ag.mean <- split[[i]] %>%
    mutate(year.cuts = as.character(year.cuts)) %>%
    group_by(subject_name, subject_sex, subject_age, site, region, year.cuts, cut.date) %>%
    summarise(n = n(),
              mean.occupancy = mean(ag.used))
  mod.predict[[i]] <- predict.Mclust(mod, ag.mean$mean.occupancy)
  
  # create result table for each season year
  
  result[[i]] <- ag.mean
  #result[[i]][,paste0("ag.class.",ag.mean$year.cuts[1])] <- mod.predict[[i]]$classification
  result[[i]]$ag.class.season <- mod.predict[[i]]$classification
}

names(result) <- c('2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018', '2019')

# plot predictions by season year
par(mfrow = c(3,3))
for(i in 1:length(mod.predict)){
  ag.mean <- split[[i]] %>%
    group_by(subject_name, year.cuts) %>%
    summarise(n = n(),
              mean.occupancy = mean(ag.used))
  mclust1Dplot(ag.mean$mean.occupancy, z = mod.predict[[i]]$z, classification = mod.predict[[i]]$classification,
               what = "classification")
}


##### Frequency Tables #####

## clust summary - by individual season
tactics.summary <- do.call(rbind, result) %>%
  group_by(ag.class.season) %>%
  summarise(n = n(),
            mean = mean(mean.occupancy),
            sd = sd(mean.occupancy),
            se = sd/sqrt(n),
            c.05 = mean - qt(1 - (0.05 / 2), n - 1) * se,
            c.95 = mean + qt(1 - (0.05 / 2), n - 1) * se)
tactics.summary 

## individual season classifications
tactics.seasons <- do.call(rbind, result) %>%
  group_by(ag.class.season) %>%
  select(subject_name, subject_sex, site, region, year.cuts, ag.class.season) %>%
  spread(key = year.cuts, value = ag.class.season) %>%
  mutate(mean.occupancy = clust.result$mean.occupancy, ag.class = clust.result$ag.class.mean) %>%
  relocate(subject_name, subject_sex, region, site, mean.occupancy, ag.class)
tactics.seasons

## tactic summary over time with rel freqs
season <- do.call(rbind, result) %>%
  group_by(year.cuts, cut.date, ag.class.season) %>%
  summarise(n = n()) %>%
  mutate(prop=round(n/sum(n), 2)) %>%
  ungroup() %>% mutate(cut.date = ymd(cut.date))

# plot
ggplot(data = season, aes(x = cut.date, y = prop)) + 
  geom_line(aes(color = as.factor(ag.class.season), linetype = as.factor(ag.class.season))) + 
  geom_point(aes(color = as.factor(ag.class.season)))

## tactics by sex
tactics.sex <- do.call(rbind, result) %>%
  group_by(subject_sex, ag.class.season) %>%
  summarise(n = n()) %>%
  mutate(prop=round(n/sum(n), 2))
tactics.sex


## tactics by site - frequency calculated relative to the site
tactics.region <- do.call(rbind, result) %>%
  group_by(region, ag.class.season) %>%
  summarise(n = n()) %>%
  mutate(prop=round(n/sum(n), 2))
tactics.region

tactics.site <- do.call(rbind, result) %>%
  group_by(site, ag.class.season) %>%
  summarise(n = n()) %>%
  mutate(prop=round(n/sum(n), 2))
tactics.site

##### Add seasonal and aggregated results to the dataframe

# bind lumped results to movdata
agg.clust <- clust.result %>%
  select(subject_name, tactic.agg = ag.class.mean) %>%
  inner_join(., movdat, by = c("subject_name")) 

# bind summarized results
clust.df <- do.call(rbind, result) %>%
  ungroup() %>%
  select(subject_name, year.cuts, tactic.season = ag.class.season) %>%
  inner_join(., agg.clust, by = c("subject_name", "year.cuts")) %>%
  as.data.frame()


##### Save Tables #####
# write.csv(tactics.summary, './GMM/results/cluster tables/tactics_summary.csv')
# write.csv(tactics.seasons, './GMM/results/cluster tables/tactics_season_classes.csv')
# write.csv(tactics.sex, './GMM/results/cluster tables/tactics_sex.csv')
#write.csv(tactics.site, './GMM/results/cluster tables/tactics_site.csv')
# write.csv(tactics.region, './GMM/results/cluster tables/tactics_region.csv')
saveRDS(clust.df, './GMM/GMEcollars_002_usedClust_2020-07-09.rds')


