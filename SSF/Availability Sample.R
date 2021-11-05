##### Generate Availability Sample - Step Selection Function #####
## Nathan Hahn
## 2021-10-11

library(amt)
library(tidyverse)
library(ggplot2)

# read data
movdata <- readRDS("./movdata/GMEcollars_004_usedClust_2021-10-05.rds")

# subset for testing
movdata <- movdata[movdata$fixType != 'irregular',]
#movdata <- movdata[1:100000,]

# look at step length distribution
hist(movdata$dist)

# look at turn angle distribution
hist(movdata$rel.angle)

# make a move object 
track <- make_track(movdata, .x = x, .y = y, .t = date, id = id)

# create a nested dataframe - needed for multiple id
t1 <- track %>% nest(data = -"id")

# create bursts and traj - map the two functions to a nested dataset t1
t2 <- t1 %>% 
  mutate(steps = map(data, function(x) 
    x %>% track_resample(rate = minutes(60), tolerance = minutes(30)) %>% steps_by_burst() ))

t2 %>% select(id, steps) %>% unnest(cols = steps) %>% 
  ggplot(aes(ta_, fill = factor(id))) + geom_density(alpha = 0.4)

# make random steps - NOT WORKING
t3.randsteps <- t2 %>%
  unnest(cols = steps) %>%
  random_steps(n_control = 10)

t3.randsteps$data <- NULL

head(t3.randsteps)

write.csv(t3.randsteps, './SSF/movdata_004_randsteps_2021-10-12.csv')
