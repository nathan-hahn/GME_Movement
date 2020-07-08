##### GMM Model Visualization #####
# For visualizing and compiling GMM results. 
# Attach classifications to the dataframe

# Load cluster results table

movdat <- readRDS('./GMM/res90_agmovingwindow.rds') 
  
mod <- readRDS('./GMM/results/modSelect.rds')

clust.result <- read.csv("./GMM/results/clustResult_GME_2020-06-25.csv")
clust.summary <- clust.result %>%
  group_by(ag.class.roll) %>%
  summarise(n = length(ag.class.roll),
            mean = mean(mean.occupancy),
            sd = sd(mean.occupancy),
            se = sd/sqrt(n),
            lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
            upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se)
(clust.summary)


##### Predict tactics for each season #####
split <- split(movdat, as.character(movdat$year.cuts))
mod.predict <- list()
result <- list()
class <- list()

for(i in 1:length(split)) {
  ag.mean <- split[[i]] %>%
    mutate(year.cuts = as.character(year.cuts)) %>%
    group_by(subject_name, subject_sex, year.cuts, cut.date) %>%
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


# TODO: bind seasonal and overall results to movdata
t <- do.call(rbind, result) %>%
  select(subject_name, year.cuts, ag.class.season) %>%
  



## clust summary - by individual season
tactics.season <- do.call(rbind, result) %>%
  group_by(ag.class.season) %>%
  summarise(n = n(),
            mean = mean(mean.occupancy),
            sd = sd(mean.occupancy),
            se = sd/sqrt(n),
            c.05 = mean - qt(1 - (0.05 / 2), n - 1) * se,
            c.95 = mean + qt(1 - (0.05 / 2), n - 1) * se)
tactics.season 


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
  group_by(ag.class.season, subject_sex) %>%
  summarise(n = n())


## tactics by site


