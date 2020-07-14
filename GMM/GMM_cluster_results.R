##### GMM Model Visualization #####
# For visualizing and compiling GMM results. 
# Attach classifications to the dataframe
library(tidyverse)
library(lubridate)
library(mclust)


# Load cluster results table

movdat <- readRDS('./GMM/GMEcollars_002_res90_2020-07-14.rds') 
  
mod <- readRDS('./GMM/results/mSelect_2020-07-14.rds')

clust.result <- read.csv("./GMM/results/clustResult_GME_2020-07-14.csv")
clust.summary <- clust.result %>%
  mutate(ag.class.mean = as.factor(ag.class.mean)) %>%
  group_by(ag.class.mean) %>%
  summarise(n = length(ag.class.mean),
            mean = MeanCI(mean.occupancy, method = "boot", type = "norm", R=200)[1],
            lwr.ci = MeanCI(mean.occupancy, method = "boot", type = "norm", R=200)[2],
            upr.ci = MeanCI(mean.occupancy, method = "boot", type = "norm", R=200)[3]) 
(clust.summary) # overall

ggplot(clust.summary, aes(ag.class.mean, mean)) + geom_pointrange(
  aes(ymin = lwr.ci, ymax = upr.ci, color = ag.class.mean), 
  position = position_dodge(0.3)) 

##### Predict tactics for each season #####
split <- split(movdat, as.character(movdat$year.cuts))
mod.predict <- list()
result <- list()
class <- list()

for(i in 1:length(split)) {
  ag.mean <- split[[i]] %>%
    mutate(year.cuts = as.character(year.cuts)) %>%
    group_by(subject_name, subject_sex, subject_ageClass, site, region, year.cuts, cut.date) %>%
    summarise(n = n(),
              mean.occupancy = mean(ag.used))
  mod.predict[[i]] <- predict.Mclust(mod, ag.mean$mean.occupancy)
  
  # create result table for each season year
  
  result[[i]] <- ag.mean
  #result[[i]][,paste0("ag.class.",ag.mean$year.cuts[1])] <- mod.predict[[i]]$classification
  result[[i]]$ag.class.year <- mod.predict[[i]]$classification
}


#names(result) <- c("2011", '2012', "2013", "2014", "2015", "2016", "2017", "2018", "2019")

# plot predictions by season year
par(mfrow = c(3,3))
for(i in 1:length(mod.predict)){
  mclust1Dplot(result[[i]]$mean.occupancy, z = mod.predict[[i]]$z, classification = mod.predict[[i]]$classification,
               what = "classification")
}


##### Frequency Tables #####
library(DescTools)

tactics <- do.call(rbind, result)

# check group counts before generating statistics. Sex is evenly distributed but other variables cannot be compared across groups! 
tactics %>%  group_by(subject_ageClass) %>%
  summarise(n = n()) 

t <- filter(tactics, n > 24*30)

## clust summary - by individual season
tactics.summary <- do.call(rbind, result) %>%
  mutate(ag.class.year = as.factor(ag.class.year)) %>%
  group_by(ag.class.year) %>%
  summarise(n = n(),
            mean = MeanCI(mean.occupancy, method = "boot", type = "norm", R=200)[1],
            lwr.ci = MeanCI(mean.occupancy, method = "boot", type = "norm", R=200)[2],
            upr.ci = MeanCI(mean.occupancy, method = "boot", type = "norm", R=200)[3]) 
tactics.summary 

ggplot(tactics.summary, aes(ag.class.year, mean)) + geom_pointrange(
  aes(ymin = lwr.ci, ymax = upr.ci, color = ag.class.year), 
  position = position_dodge(0.3)) 

## individual season classifications
tactics.seasons <- do.call(rbind, result) %>%
  group_by(ag.class.year) %>%
  dplyr::select(year.cuts, subject_name, subject_sex, subject_ageClass, site, region, ag.class.year) %>%
  droplevels() %>%
  pivot_wider(names_from = year.cuts, values_from = ag.class.year) %>%
  arrange(subject_name) %>% # make sure the order of the pivot table is the same as clust.result!!
  mutate(mean.occupancy = clust.result$mean.occupancy, ag.class = clust.result$ag.class.mean) %>%
  relocate(subject_name, subject_sex, region, site, mean.occupancy, ag.class)
tactics.seasons

## tactic summary over time with rel freqs
season <- do.call(rbind, result) %>%
  filter(n > 24*30*2, year.cuts %in% c('2016','2017', '2018', '2019')) %>%
  group_by(year.cuts, cut.date, ag.class.year) %>%
  summarise(n = n()) %>%
  mutate(prop=round(n/sum(n), 2)) %>%
  ungroup() %>% mutate(cut.date = ymd(cut.date))

# plot
ggplot(data = season, aes(x = cut.date, y = prop)) + 
  geom_line(aes(color = as.factor(ag.class.year), linetype = as.factor(ag.class.year))) + 
  geom_point(aes(color = as.factor(ag.class.year)))

## tactics by sex - compare prop of sexes by ag class
tactics.sex <- do.call(rbind, result) %>%
  group_by(ag.class.year, subject_sex) %>%
  summarise(n = n()) %>%
  mutate(prop = MultinomCI(n, conf.level = 0.95, sides = 'two.sided')[,1],
         lwr.ci = MultinomCI(n, conf.level = 0.95, sides = 'two.sided')[,2],
         upr.ci = MultinomCI(n, conf.level = 0.95, sides = 'two.sided')[,3]
  )
tactics.sex

ggplot(tactics.sex, aes(ag.class.year, prop)) + geom_pointrange(
  aes(ymin = lwr.ci, ymax = upr.ci, color = subject_sex), 
    position = position_dodge(0.3)) 


## tactics by age class - compare prop of ag classes by age class. Distribution of age classes is not even so hard to compare props between age classes
tactics.age <- do.call(rbind, result) %>%
  mutate(subject_ageClass = factor(subject_ageClass, levels = c('young adult', 'adult', 'mature'))) %>%
  group_by(subject_ageClass, ag.class.year) %>%
  summarise(n = n()) %>%
  mutate(prop = MultinomCI(n, conf.level = 0.95, sides = 'two.sided')[,1],
         lwr.ci = MultinomCI(n, conf.level = 0.95, sides = 'two.sided')[,2],
         upr.ci = MultinomCI(n, conf.level = 0.95, sides = 'two.sided')[,3]
  )
tactics.age

ggplot(tactics.age, aes(ag.class.year, prop)) + geom_pointrange(
  aes(ymin = lwr.ci, ymax = upr.ci, color = subject_ageClass), 
  position = position_dodge(0.3)) 


## tactics by region - frequency calculated relative to the site
tactics.region <- do.call(rbind, result) %>%
  group_by(region, ag.class.year) %>%
  summarise(n = n()) %>%
  mutate(prop = MultinomCI(n, conf.level = 0.95, sides = 'two.sided')[,1],
         lwr.ci = MultinomCI(n, conf.level = 0.95, sides = 'two.sided')[,2],
         upr.ci = MultinomCI(n, conf.level = 0.95, sides = 'two.sided')[,3]
  )
tactics.region

ggplot(tactics.region, aes(ag.class.year, prop)) + geom_pointrange(
  aes(ymin = lwr.ci, ymax = upr.ci, color = region), 
  position = position_dodge(0.3)) 

# site-based results
tactics.site <- do.call(rbind, result) %>%
  group_by(site, ag.class.year) %>%
  summarise(n = n()) %>%
  mutate(prop = MultinomCI(n, conf.level = 0.95, sides = 'two.sided')[,1],
         lwr.ci = MultinomCI(n, conf.level = 0.95, sides = 'two.sided')[,2],
         upr.ci = MultinomCI(n, conf.level = 0.95, sides = 'two.sided')[,3]
  )
tactics.site



##### Add seasonal and aggregated results to the dataframe

# bind lumped results to movdata
agg.clust <- clust.result %>%
  dplyr::select(subject_name, tactic.agg = ag.class.mean) %>%
  inner_join(., movdat, by = c("subject_name")) 

# bind summarized results
clust.df <- do.call(rbind, result) %>%
  ungroup() %>%
  dplyr::select(subject_name, year.cuts, tactic.season = ag.class.year) %>%
  inner_join(., agg.clust, by = c("subject_name", "year.cuts")) %>%
  as.data.frame()


##### Save Tables #####
# write.csv(tactics.summary, './GMM/results/cluster tables/tactics_summary.csv')
# write.csv(tactics.seasons, './GMM/results/cluster tables/tactics_season_classes.csv')
# write.csv(tactics.sex, './GMM/results/cluster tables/tactics_sex.csv')
#write.csv(tactics.site, './GMM/results/cluster tables/tactics_site.csv')
# write.csv(tactics.region, './GMM/results/cluster tables/tactics_region.csv')
# saveRDS(clust.df, './GMM/GMEcollars_002_usedClust_2020-07-14.rds')


