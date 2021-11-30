##### Classify Tactics in New Data #####

# Script to classify tactics in new data


library(tidyverse)
library(lubridate)
library(mclust)
library(DescTools)
theme_set(  theme_bw() + # set theme with no legend of strip text
              theme(panel.grid.major = element_blank(),
                    strip.background = element_blank(),
                    panel.border = element_rect(colour = "black"),
                    strip.text = element_text(size = 12),
                    legend.text = element_text(size = 10),
                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title=element_text(size=14), axis.text = element_text(size=12))
)

source('GME_functions.R')

# Load cluster results table

movdat <- readRDS('./movdata/GMEcollars_004_usedFiltercuts_2021-10-05.rds')

mod <- readRDS('./GMM/results/mSelect_2020-10-30.rds')


##### Get cut points for the ag classes #####
cut.points <- vector()

cut.points[1] <- cluster_cutpoint(mod, comp = c(1,2), grph=T)
cut.points[2] <- cluster_cutpoint(mod, comp = c(2,3), grph = T)
cut.points[3] <- cluster_cutpoint(mod, comp = c(3,4), grph = T)

cut.points*100


##### Add aggregated tactic #####

ag.mean <- movdat %>%
  mutate(year.cuts = as.character(year.cuts)) %>%
  group_by(subject_name) %>%
  summarise(n = n(),
            mean.occupancy = mean(ag.used))
mod.predict <- predict.Mclust(mod, ag.mean$mean.occupancy)

# create result table for each season year
aggregate.result <- ag.mean
aggregate.result$ag.class.agg <- mod.predict$classification


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
  result[[i]]$ag.class.year <- mod.predict[[i]]$classification
}


# plot predictions by season year
par(mfrow = c(3,4))
for(i in 1:length(mod.predict)){
  mclust1Dplot(result[[i]]$mean.occupancy, z = mod.predict[[i]]$z, classification = mod.predict[[i]]$classification,
               what = "classification")
}


##### Frequency Tables #####
library(DescTools)

tactics.df <- do.call(rbind, result) %>%
  mutate(ag.class.year = recode_factor(as.factor(ag.class.year),
                                       "1" = "Rare",
                                       "2" = "Sporadic",
                                       "3" = "Seasonal",
                                       "4" = "Habitual")) %>%
  filter(n > 500)

# re-order factor levels
tactics.df$ag.class.year <- factor(tactics.df$ag.class.year, levels = c("Rare", "Sporadic", "Seasonal", "Habitual"))

# plot
ggplot(tactics.df, aes(ag.class.year, mean.occupancy)) + geom_boxplot() + 
  geom_point(shape = 15, color = 'dark grey', position = position_jitter(width = 0.21)) +
  xlab("Tactic Class") + ylab("Mean Agricultural Use") 


## clust summary - by individual season
tactics.summary <- tactics.df %>%
  mutate(ag.class.year = as.factor(ag.class.year)) %>%
  group_by(ag.class.year) %>%
  summarise(n = n(),
            mean = MeanCI(mean.occupancy, method = "boot", type = "norm", R=200)[1],
            lwr.ci = MeanCI(mean.occupancy, method = "boot", type = "norm", R=200)[2],
            upr.ci = MeanCI(mean.occupancy, method = "boot", type = "norm", R=200)[3]) 
tactics.summary 

ggplot(tactics.summary, aes(ag.class.year, mean)) + geom_pointrange(
  aes(ymin = lwr.ci, ymax = upr.ci)) + 
  xlab("Tactic class") + ylab("Mean Ag Occupancy") 

t <- tactics.summary %>%
  mutate(avg.season = n/9,
         prop.season = avg.season/sum(avg.season))

t <- tactics.df %>%
  group_by(year.cuts, ag.class.year) %>%
  summarise(n = n()) %>%
  group_by(year.cuts) %>%
  mutate(prop = n/sum(n)) 

ggplot(t, aes(ag.class.year, prop)) + geom_boxplot(aes(color = ag.class.year)) +
  scale_color_manual(values = c("blue", "#FF0000FF", "#80FF00FF", "#8000FFFF")) +
  xlab("Tactic Class") + ylab("Yearly Population Proportion") + labs(color = "Tactic Class")


##### Add classifications to dataframe #####

# bind aggregate tactic classification results to movdata
agg.clust <- aggregate.result %>%
  dplyr::select(subject_name, tactic.agg = ag.class.agg) %>%
  inner_join(., movdat, by = c("subject_name")) 

# bind tactic-year classification results to movdata
clust.df <- do.call(rbind, result) %>%
  ungroup() %>%
  dplyr::select(subject_name, year.cuts, tactic.season = ag.class.year) %>%
  inner_join(., agg.clust, by = c("subject_name", "year.cuts")) %>%
  mutate(year.cuts = as.factor(year.cuts), tactic.season = as.factor(tactic.season)) %>%
  as.data.frame() 

##### Export #####
saveRDS(clust.df, './movdata/GMEcollars_004_usedClust_2021-10-05.rds')




