#### Stage Test - Net Squared Displacement ####

library(tidyverse)
library(lubridate)


# import state-classified data
gme <- as.data.frame(data.table::fread("HMM/movdata/GMEcollars_003_HMMclassified_20201206.csv"))
colnames(gme)[1] <- 'uid'
unique(gme$subject_name)

# filter out outlier steps 
split <- split(gme, gme$viterbi)
split[[1]] <- filter(split[[1]], dist <= median((split[[2]]$dist + runif(1, -200, 200)), na.rm = T))
split[[2]] <- filter(split[[2]], dist <= median((split[[3]]$dist + runif(1, -200, 200)), na.rm = T))
split[[3]] <- filter(split[[3]], dist >= median((split[[2]]$dist + runif(1, -200, 200)), na.rm = T))
gme <- do.call(rbind, split)
gme <- gme[with(gme, order(uid)), ]

## define the period for staging assessment (10am - 2pm, 5 hrs)
gme$stage.period <- ifelse(hour(gme$date) >= 10 & hour(gme$date) <= 14, 1, 0)
# unique id for each individual-day
gme$dayBurst <- paste(gme$id, as.Date(gme$date))

# 24-hour ag window identifies days where they used ag - dayburst helps group by day and individual 
gme <- gme %>%
  group_by(dayBurst) %>%
  mutate(ag.window.ext = ifelse(sum(ag.window) >= 1, 1, 0)) # if ag.window identifies ag use phase, classify that day as ag use. Helps to pick up stages at the beginning of an ag use phase that might otherwise be split up.


## assign unique id to each ag window event
# define unique id for window/non-window phases
gme <- gme[!is.na(gme$ag.window),]
eventFlag <- ifelse(gme$ag.window == 1, TRUE, FALSE)
eventIndex <- inverse.rle(within.list(rle(eventFlag), 
                                      values[values] <- seq_along(values[values])))
# assign event index
gme$ag.window.index <- eventIndex


## Viz net square displacement
test <- filter(gme, ag.window == 1, subject_name == 'Ivy')
split <- split(test, test$ag.window.index)

plot(split$`1839`$R2n)

## Calculate NSD by stage period
library(adehabitatLT)
test <- filter(gme, stage.period == 1 & subject_name == 'Ivy')
traj <- as.ltraj(cbind(test$x, test$y), date = ymd_hms(test$date), id = test$id, burst = test$dayBurst)
traj.df <- ld(traj)
test$nsd <- traj.df$R2n
day.nsd <- test %>%
  group_by(id, dayBurst, ag.window) %>%
  summarise(maxNSD = max(R2n)) %>%
  group_by(id, ag.window) %>%
  summarise_at('maxNSD', .funs = list(min=min, Q1=~quantile(., probs = 0.25),
                                    median=median, Q3=~quantile(., probs = 0.75),
                                    max=max))



##### Figure - Ratios by hour #####

# overall ratio
t <- gme %>% 
  filter(!is.na(ag.window.ext)) %>%
  mutate(hour = hour(date)) %>%
  group_by(hour, ag.window.ext) %>%
  summarise(encamped = sum(viterbi==1),
            meandering = sum(viterbi==2),
            dw = sum(viterbi==3),
            n = n()) %>%
  mutate(ratio = encamped/meandering) %>%
  filter(!is.na(ratio)) %>%
  mutate(ratio = ifelse(ratio == Inf, encamped, ratio))

t$ag.window <- as.factor(t$ag.window.ext)
levels(t$ag.window) <- c('non-window', 'ag-window')

ggplot(t, aes(x = as.factor(hour), y = ratio, color=ag.window)) + geom_point() +
  ggtitle("Encamped/Meandering Ratio by Hour and Ag Window")

# ratio with individual variation -- to get confidence intervals
t <- gme %>%
  filter(!is.na(ag.window.ext)) %>%
  mutate(hour = hour(date)) %>%
  group_by(subject_name, hour, ag.window.ext) %>%
  mutate(enc.day = sum(viterbi==1), meander.day = sum(viterbi==2), dw.day = sum(viterbi==3)) %>%
  mutate(ratio = enc.day/meander.day) %>%
  mutate(ratio = ifelse(ratio == Inf, enc.day, ratio))

# get confidence intervals
t <- t %>%
  group_by(hour, ag.window.ext) %>%
  summarise(mean.ratio = mean(ratio, na.rm = TRUE),
            sd = sd(ratio, na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd / sqrt(n),
         lwr.ci = mean.ratio - qt(1 - (0.05 / 2), n - 1) * se,
         upr.ci = mean.ratio + qt(1 - (0.05 / 2), n - 1) * se)

t$ag.window <- as.factor(t$ag.window.ext)
levels(t$ag.window) <- c('non-window', 'ag-window')

ggplot(t, aes(as.factor(hour), mean.ratio)) + geom_pointrange(
  aes(ymin = lwr.ci, ymax = upr.ci, color = ag.window), 
  position = position_dodge(0.5), size = .2) + 
  ggtitle("Encamped/Meandering Ratio by Hour and Ag Window")

##### Look at Ratios #####
## Ratio of encamped to meandering in ag and non-ag windows during the day
gme %>% 
  group_by(stage.period, ag.window) %>%
  summarise(encamped = sum(viterbi==1),
            meandering = sum(viterbi==2),
            dw = sum(viterbi==3),
            n = n()) %>%
  mutate(ratio = encamped/meandering)

## What is the quartiles for encamped/meandering ratio OUTSIDE of ag windows
t <- gme %>% 
  filter(!is.na(ag.window.ext)) %>%
  group_by(dayBurst, stage.period, ag.window.ext) %>%
  summarise(encamped = sum(viterbi==1),
            meandering = sum(viterbi==2),
            dw = sum(viterbi==3),
            n = n()) %>%
  mutate(ratio = encamped/meandering) %>%
  filter(!is.na(ratio)) %>%
  mutate(ratio = ifelse(ratio == Inf, encamped, ratio))

t %>%  group_by(stage.period, ag.window.ext) %>%
  summarise_at('ratio', .funs = list(min=min, Q1=~quantile(., probs = 0.25),
                                     median=median, Q3=~quantile(., probs = 0.75),
                                     max=max))

f <- filter(t, stage.period == 1)
boxplot(f$ratio ~ f$ag.window.ext)

stage <- t %>%
  filter(ratio >= 0.5 & stage.period == 1) %>%
  filter(ag.window.ext == 1) %>%
  filter(dw == 0)

dim(stage)


##### Tag Stagging Events by Ratio #####
## Tag staging events in the full GME dataset based on dayBurst and stage.period 
# Check if we need to filter out 

## define the period for staging assessment (10am - 2pm, 5 hrs)
gme$stage.period <- ifelse(hour(gme$date) >= 10 & hour(gme$date) <= 15, 1, 0)

t <- gme %>%
  filter(!is.na(ag.window.ext)) %>%
  group_by(dayBurst, stage.period, ag.window.ext) %>%
  mutate(enc.day = sum(viterbi==1), meander.day = sum(viterbi==2), dw.day = sum(viterbi==3), n.day = n()) %>%
  mutate(ratio = enc.day/meander.day) %>%
  mutate(ratio = ifelse(ratio == Inf, enc.day, ratio))


t %>% filter(!is.na(ratio)) %>%
  group_by(stage.period, ag.window.ext) %>%
  summarise_at('ratio', .funs = list(min=min, Q1=~quantile(., probs = 0.25),
                                     median=median, Q3=~quantile(., probs = 0.75),
                                     max=max))

# tag all staging events
t$stage.event <- ifelse(t$ratio >0.5 & t$stage.period == 1 & t$dw.day == 0, 1, 0)
# tag ag staging events
t$ag.stage.event <- ifelse(t$ratio >0.5 & t$stage.period == 1 & t$dw.day == 0 & t$ag.window.ext == 1, 1, 0)

## index staging events
eventFlag <- ifelse(t$stage.event == 1, TRUE, FALSE)
eventIndex <- inverse.rle(within.list(rle(eventFlag), 
                                      values[values] <- seq_along(values[values])))
# assign event index
t$stage.event.index <- eventIndex

stage <- filter(t, stage.event == 1)



                                                                 

plot(stage$x, stage$y, pch = 19)


ag.stage <- filter(t, ag.stage.event == 1)
plot(ag.stage$x, ag.stage$y, pch = 19)

##### Plot staging event relocs #####

library(sf)
library(mapview)
stage.sf <- st_as_sf(stage, coords = c('x', 'y'), crs = 32736)

#mapview(stage.sf, cex = 2, zcol = 'ag.stage.event', burst = TRUE)

# check dataframe
temp <- stage %>%
  dplyr::select(uid, subject_name, date, dist, viterbi, stage.period, stage.event.index)


temp2 <- stage %>%
  group_by(stage.event.index) %>%
  filter(!any(dist > 214))

temp2 %>% group_by(viterbi) %>% summarise_at('dist', .funs = list(min=min, Q1=~quantile(., probs = 0.25),
                                                                 median=median, Q3=~quantile(., probs = 0.75),
                                                                 max=max))

##### Make Polygons #####
library(sp)
library(adehabitatHR)

## All staging areas
stage <- stage %>%
  group_by(stage.event.index) %>%
  mutate(n = n()) %>%
  filter(n >= 5) %>% # remove staging events with missing data
  filter(ag.window.ext == 1) %>%
  dplyr::select(x, y, stage.event.index)
stage.sp <- stage
coordinates(stage.sp) <-~x+y

# build MCPs and attach data
mcp.stage <- mcp(stage.sp, percent = 100)
mcp.stage@data <- as.data.frame(stage.sp@data)



## Buffer and Union of staging hulls
# convert data frame to sf object
sf.object <- stage %>%
  st_as_sf(coords = c("x", "y"))

# convert MCPs to sf polygon objects - tied to the stage id
hulls <- sf.object %>%
  group_by(stage.event.index) %>%
  summarise( geometry = st_combine( geometry ) ) %>%
  st_convex_hull() 

# Buffer and union of hulls
hull.union <- hulls %>%
  st_cast("POLYGON") %>%
  ungroup() %>%
  # buffer
  st_buffer(dist = median(gme$dist)) %>% # buffer by median step length of the dataset
  # union polygons
  st_union() %>%
  # set crs
  st_set_crs(st_crs(32736))

mapview(hull.union)





##### Staging Event Stats #####
# % of ag use days that had a staging event
# % of overall time spent staging during ag use windows
# overall % of encamped/meandering during staging
# summary of displacement during staging (NSD and Displacement)


# % of ag use days that had a staging event -- how often are elephants staging before a raid?






##### Recursion #####
library(recurse)

df.recurse <- gme %>% ungroup() %>%
  filter(ag.window == 1) %>%
  dplyr::select(c('x', 'y', 't'='date', 'id'='subject_name')) %>%
  mutate(t = as.POSIXct(t), id = as.factor(id))
df.recurse <- as.data.frame(df.recurse)

# test with ag windows - best
test <- getRecursions(df.recurse, radius = 350, threshold = 100, timeunits = 'hours')

plot(test, df.recurse, legendPos = c(13, -10), pch = 19, cex = 0.5)

# test with staging polygon hulls (can only take 1 at a time)
poly <- hull.union %>% st_cast("POLYGON") %>%
  as_Spatial() %>%
  split(1:length(.))
test <- getRecursionsInPolygon(as.data.frame(df.recurse), out[[550]])
test

plot(test, df.recurse, legendPos = c(13, -10), pch = 19, cex = 0.5)



recurse.test <- lapply(out, getRecursionsInPolygon, trajectory = df.recurse)
plot(recurse.test)
