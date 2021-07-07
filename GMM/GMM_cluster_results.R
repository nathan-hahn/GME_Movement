##### GMM Model Visualization #####
# For visualizing and compiling GMM results. 
# Attach classifications to the dataframe
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

movdat <- readRDS('./GMM/movdata/GMEcollars_003_res90_2020-10-30.rds') 
  
mod <- readRDS('./GMM/results/mSelect_2020-10-30.rds')

clust.result <- read.csv("./GMM/results/clustResult_GME_2020-10-30.csv") 
clust.summary <- clust.result %>%
  mutate(ag.class.mean = recode_factor(as.factor(ag.class.mean),
                                       "1" = "Rare",
                                       "2" = "Sporadic",
                                       "3" = "Seasonal",
                                       "4" = "Habitual")) %>%
  group_by(ag.class.mean) %>%
  summarise(n = length(ag.class.mean),
            mean = MeanCI(mean.occupancy, method = "boot", type = "norm", R=200)[1],
            lwr.ci = MeanCI(mean.occupancy, method = "boot", type = "norm", R=200)[2],
            upr.ci = MeanCI(mean.occupancy, method = "boot", type = "norm", R=200)[3]) 
(clust.summary) # overall

## Plot cluster summaries by mean ag use trends
ggplot(clust.summary, aes(ag.class.mean, mean)) + geom_pointrange(
  aes(ymin = lwr.ci, ymax = upr.ci, color = ag.class.mean), 
  position = position_dodge(0.3)) + 
  scale_color_manual(values = c("blue", "#FF0000FF", "#80FF00FF", "#8000FFFF")) +
  xlab("Tactic Class") + ylab("Mean Occupancy") + labs(color = "Tactic Class")
  
ggplot(clust.summary, aes(ag.class.mean, mean)) + geom_boxplot(
  aes(ymin = lwr.ci, ymax = upr.ci, color = ag.class.mean), 
  position = position_dodge(0.3)) + 
  scale_color_manual(values = c("blue", "#FF0000FF", "#80FF00FF", "#8000FFFF")) +
  xlab("Tactic Class") + ylab("Mean Occupancy") + labs(color = "Tactic Class")


##### Get cut points for the ag classes #####
cut.points <- vector()

cut.points[1] <- cluster_cutpoint(mod, comp = c(1,2), grph=T)
cut.points[2] <- cluster_cutpoint(mod, comp = c(2,3), grph = T)
cut.points[3] <- cluster_cutpoint(mod, comp = c(3,4), grph = T)

cut.points*100

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
par(mfrow = c(3,3))
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


## individual season classifications
tactics.seasons <- tactics.df %>%
  group_by(ag.class.year) %>%
  dplyr::select(year.cuts, subject_name, subject_sex, subject_ageClass, site, region, ag.class.year) %>%
  droplevels() %>%
  pivot_wider(names_from = year.cuts, values_from = ag.class.year) %>%
  arrange(subject_name) %>% # make sure the order of the pivot table is the same as clust.result!!
  mutate(mean.occupancy = clust.result$mean.occupancy, ag.class.aggregate = clust.result$ag.class.mean) %>%
  relocate(subject_name, subject_sex, region, site, mean.occupancy, ag.class.aggregate)
tactics.seasons

# total proportions of tactics (individual-year)
t <- table(factor(tactics.df$ag.class.year))
prop.table(t)


## tactic summary over time with rel freqs
season <- tactics.df %>%
  filter(n > 24*30*2, year.cuts %in% c('2016','2017', '2018', '2019')) %>%
  group_by(year.cuts, cut.date, ag.class.year) %>%
  summarise(n = n()) %>%
  mutate(prop=round(n/sum(n), 2)) %>%
  ungroup() %>% mutate(cut.date = ymd(cut.date))

ggplot(data = season, aes(x = cut.date, y = prop)) + 
  geom_line(aes(color = as.factor(ag.class.year), linetype = as.factor(ag.class.year))) + 
  geom_point(aes(color = as.factor(ag.class.year))) 


## Tactics by sex - compare prop of sexes by ag class
tactics.sex <- amplitude %>%
  group_by(tactic.season, subject_sex) %>%
  summarise(n = n()) %>%
  mutate(prop = MultinomCI(n, conf.level = 0.95, sides = 'two.sided')[,1],
         lwr.ci = MultinomCI(n, conf.level = 0.95, sides = 'two.sided')[,2],
         upr.ci = MultinomCI(n, conf.level = 0.95, sides = 'two.sided')[,3],
  )
tactics.sex

tactics.sex$subject_sex <- as.factor(tactics.sex$subject_sex)
levels(tactics.sex$subject_sex) <- c('Female', 'Male')

ggplot(tactics.sex, aes(tactic.season, prop)) + geom_pointrange(
  aes(ymin = lwr.ci, ymax = upr.ci, shape = subject_sex), 
    position = position_dodge(0.3)) + xlab("Tactic class") + ylab("Proportion") +
  labs(shape = 'Sex')


## Tactics by age class - compare prop of ag classes by age class. Distribution of age classes is not even so hard to compare props between age classes
tactics.age <- tactics.df %>%
  mutate(subject_ageClass = recode_factor(as.factor(subject_ageClass), 
                                          "mature adult" = "old adult")) %>%
  group_by(ag.class.year, subject_ageClass) %>%
  summarise(n = n()) %>%
  mutate(prop = MultinomCI(n, conf.level = 0.95, sides = 'two.sided')[,1],
         lwr.ci = MultinomCI(n, conf.level = 0.95, sides = 'two.sided')[,2],
         upr.ci = MultinomCI(n, conf.level = 0.95, sides = 'two.sided')[,3]
  )
tactics.age

ggplot(tactics.age, aes(ag.class.year, prop)) + geom_pointrange(
  aes(ymin = lwr.ci, ymax = upr.ci, color = subject_ageClass), 
  position = position_dodge(0.3)) + xlab("Tactic class") + ylab("Proportion") +
  labs(color = 'Age Class')


## Tactics by age class and sex
tactics.age.sex <- tactics.df %>%
  mutate(subject_ageClass = recode_factor(as.factor(subject_ageClass), 
                                          "mature adult" = "old adult")) %>%
  group_by(ag.class.year, subject_ageClass, subject_sex) %>%
  summarise(n = n()) %>%
  mutate(prop = MultinomCI(n, conf.level = 0.95, sides = 'two.sided')[,1],
         lwr.ci = MultinomCI(n, conf.level = 0.95, sides = 'two.sided')[,2],
         upr.ci = MultinomCI(n, conf.level = 0.95, sides = 'two.sided')[,3]
  )
tactics.age.sex

  
## Tactics by region - frequency calculated relative to the region
tactics.region <- tactics.df %>%
  group_by(region, ag.class.year) %>%
  summarise(n = n()) %>%
  mutate(prop = MultinomCI(n, conf.level = 0.95, sides = 'two.sided')[,1],
         lwr.ci = MultinomCI(n, conf.level = 0.95, sides = 'two.sided')[,2],
         upr.ci = MultinomCI(n, conf.level = 0.95, sides = 'two.sided')[,3]
  )
tactics.region

ggplot(tactics.region, aes(ag.class.year, prop)) + geom_pointrange(
  aes(ymin = lwr.ci, ymax = upr.ci, color = region), 
  position = position_dodge(0.3)) + xlab("Tactic class") + ylab("Proportion") +
  labs(color = 'Region')


## Tactics by site - frequency calculated relative to the site
tactics.site <- tactics.df %>%
  group_by(site, ag.class.year) %>%
  summarise(n = n()) %>%
  mutate(prop = MultinomCI(n, conf.level = 0.95, sides = 'two.sided')[,1],
         lwr.ci = MultinomCI(n, conf.level = 0.95, sides = 'two.sided')[,2],
         upr.ci = MultinomCI(n, conf.level = 0.95, sides = 'two.sided')[,3]
  )
tactics.site


##### Tactics Over Time #####
# Plot tactics over time for certain individuals

tactics <- do.call(rbind, result) %>%
  mutate(cut.date = ymd(cut.date)) %>% rename(Year = cut.date) %>%
  filter(subject_name %in% c("Fred", "Kegol", "Kiambi", "Harakati"))


ggplot(data = tactics, aes(x = Year, y = ag.class.year)) + 
  geom_line(aes(color = as.factor(subject_name), linetype = as.factor(subject_name)), position=position_dodge2(width = 0.4, padding = 2)) + 
  geom_point(aes(color = as.factor(subject_name)), position=position_dodge(width = 0.4))


##### Add seasonal and aggregated results to the dataframe #####

# bind lumped results to movdata
agg.clust <- clust.result %>%
  dplyr::select(subject_name, tactic.agg = ag.class.mean) %>%
  inner_join(., movdat, by = c("subject_name")) 

# bind summarized results
clust.df <- do.call(rbind, result) %>%
  ungroup() %>%
  dplyr::select(subject_name, year.cuts, tactic.season = ag.class.year) %>%
  inner_join(., agg.clust, by = c("subject_name", "year.cuts")) %>%
  mutate(year.cuts = as.factor(year.cuts), tactic.season = as.factor(tactic.season)) %>%
  as.data.frame() 


##### Save Tables #####
write.csv(tactics.summary, './GMM/results/cluster tables/tactics_summary.csv')
write.csv(tactics.seasons, './GMM/results/cluster tables/tactics_season_classes_20201027.csv')
write.csv(tactics.sex, './GMM/results/cluster tables/tactics_sex.csv')
write.csv(tactics.site, './GMM/results/cluster tables/tactics_site.csv')
write.csv(tactics.region, './GMM/results/cluster tables/tactics_region.csv')
saveRDS(clust.df, './GMM/movdata/GMEcollars_003_usedClust_2020-10-30.rds')

##### Mean Ag Use Regression #####
## Build homeranges for each individual-year and extract the areas
library(adehabitatLT)
# prep spdf
t <-  filter(clust.df, !is.na(x)) %>%
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
clust.df$year.mcp.area <- mcp.areas@data$year.mcp.area

## Calculate mean daily displacement
clust.df$ymd <- as.Date(clust.df$date, "%Y/%m/%d")
stat <- function(x) c(sum = sum(x))
ag.daily <- as.data.frame(aggregate(dist ~ ymd + subject_name, clust.df, stat)) # aggregated stats by day

# bind daily displacment to full dataset
displacement <- inner_join(clust.df, ag.daily, by = c("ymd", "subject_name")) %>%
  rename(daily.disp = dist.y) %>%
  mutate(tactic.season = as.factor(tactic.season), year.cuts = as.factor(year.cuts)) %>%
  group_by(subject_name, year.cuts, tactic.season) %>%
  summarise(mu.daily.disp = mean(daily.disp, na.rm = TRUE))

output <- inner_join(clust.df, displacement, by = c("subject_name", "year.cuts", "tactic.season"))

## Model Fit
# prep dataframe
mod.df <- output %>%
  group_by(subject_name, year.cuts, tactic.season, subject_sex, subject_ageClass, year.mcp.area, mu.daily.disp) %>%
  summarise(year.mean = mean(ag.used))
  
mod.df$subject_sex <- as.factor(mod.df$subject_sex) # females as reference level
mod.df$subject_ageClass <- as.factor(mod.df$subject_ageClass) 
mod.df$subject_ageClass <- relevel(mod.df$subject_ageClass, ref = "young adult") # young adults as reference level

m1 <- lmer(year.mean ~ subject_sex + subject_ageClass + log(year.mcp.area) + (1|subject_name), data = mod.df, REML = FALSE)
#summary(m1)

m2 <- lmer(year.mean ~ subject_sex + subject_ageClass + (1|subject_name), data = mod.df, REML = FALSE)
#summary(m2)

m3 <- lmer(year.mean ~ subject_sex + log(year.mcp.area) + (1|subject_name), data = mod.df, REML = FALSE)
#summary(m3)

m3.disp <- lmer(year.mean ~ subject_sex + log(year.mcp.area) + log(mu.daily.disp) + (1|subject_name), data = mod.df, REML = FALSE)
#summary(m3.disp)

m4 <- lmer(year.mean ~ subject_sex + (1|subject_name), data = mod.df, REML = FALSE)
#summary(m4)

m5 <- lmer(year.mean ~ subject_ageClass + (1|subject_name), data = mod.df, REML = FALSE)
#summary(m5)

m6 <- lmer(year.mean ~ (1|subject_name), data = mod.df, REML = FALSE)
#summary(m6)

# AICc table
knitr::kable(AICc(m1, m2, m3, m3.disp, m4, m5, m6))

# regression table for top model
sjPlot::tab_model(m3.disp)

# plot outputs
dat <- ggpredict(m3.disp, terms = c("year.mcp.area"))
plot(dat)

