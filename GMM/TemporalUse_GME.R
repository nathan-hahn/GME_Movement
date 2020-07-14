##### Ag Temporal Use #####
library(tidyverse)
library(lubridate)
source("GME_functions.R")

##### Data Prep #####
movdat <- readRDS('./movdata/GMEcollars_002_usedFilter_2020-07-14.rds')
movdat$ag.used <- ifelse(movdat$lc.estes == 1, 1, 0) 

##### Rolling Stats #####
roll.filter <- movdat %>%
  group_by(subject_name) %>%
  tally() %>% 
  filter(n >= 24*30*3) # months

roll.df <- filter(movdat, subject_name %in% roll.filter$subject_name) %>% droplevels()
split <- split(roll.df, roll.df$subject_name)

window <- c(90*24) # in hours. Split to before and after when align = center

output.plot <- NULL
for(i in 1:length(window)){
  # calculate moving window stats - see rollstats function
  output.plot[[i]] <- window_stats(df.list = split, window = window[[i]], align = 'center')
  # remove NAs from ag windows - created by window cutoff
  #output.plot[[i]] <- filter(output.plot[[i]], is.na(ag.window) == FALSE) %>% droplevels()
}


##### Add Crop Season and Year Cuts #####
# create a df with the best rolling window (90 days, requiring 4 months of data for one elephant)
# create unique seasons by combining season.cut and cropseason
# Bi-seasonal variable - long/short rains crop seasons
roll.90 <- output.plot[[1]] %>%
  mutate(
    month = month(date),
    cropseason = case_when(
      month %in% 4:10 ~ "long",
      TRUE ~ "short"))

# Yearly cuts - set september break points
year.cuts <- ymd_hms(c("2010-04-01 00:00:00", "2011-04-01 00:00:00", "2012-04-01 00:00:00", "2013-04-01 00:00:00", "2014-04-01 00:00:00", 
                                   "2015-04-01 00:00:00", "2016-04-01 00:00:00", "2017-04-01 00:00:00", 
                                   "2018-04-01 00:00:00", "2019-04-01 00:00:00", "2020-04-01 00:00:00"), tz = 'Africa/Nairobi')
year.names <- c("2010", "2011", '2012', "2013", "2014", "2015", "2016", "2017", "2018", "2019")
# get yearly cuts using september break points
rng <- cut(roll.90$date, breaks = c(year.cuts), include.lowest = T)
rng.name <- cut(roll.90$date, breaks = c(year.cuts), include.lowest = T, labels = year.names)
roll.90$cut.date <- rng
roll.90$year.cuts <- rng.name


roll.90$roll.season <- as.factor(paste(roll.90$year.cuts, roll.90$cropseason, sep = '-'))

# export movdata with season and year cuts
saveRDS(roll.90, paste0('./GMM/GMEcollars_002_res90_', Sys.Date(),'.rds'))


##### Amplitude #####
# get the max and min values for each season
amplitude <- roll.90 %>%
  group_by(subject_name, year.cuts) %>%
  mutate(season.mean = mean(ag.used),
         season.max = max(mean),
         season.min = min(mean),
         season.diff = season.max - season.min,
         season.begin = min(date),
         season.end = max(date)) %>%
  group_by(subject_name, year.cuts, season.begin, season.end, season.mean, season.max, season.min, season.diff) %>%
  #filter to individual seasons with at least 1 month's worth of fixes (does not account for NAs)
  tally() %>% filter(n > 500) %>% droplevels()
View(amplitude)


# test model fitting on individual-year data
split <- split(amplitude, amplitude$year.cuts)


t <- Mclust(amplitude$season.mean)
plot(t, what = "classification")

t <- Mclust(split[[1]]$season.max, G = 4)
plot(t, what = "classification")

dat <- split[[13]]$season.mean
p <- predict.Mclust(t, dat)
mclust1Dplot(dat, z = p$z, classification = p$classification)


##### Visualize #####
year.cuts <- ymd_hms(c("2010-04-01 00:00:00", "2011-04-01 00:00:00", "2012-04-01 00:00:00", "2013-04-01 00:00:00", "2014-04-01 00:00:00", 
                       "2015-04-01 00:00:00", "2016-04-01 00:00:00", "2017-04-01 00:00:00", 
                       "2018-04-01 00:00:00", "2019-04-01 00:00:00", "2020-04-01 00:00:00"), tz = 'Africa/Nairobi')

season.cuts <- ymd_hms(c("2010-04-01 00:00:00", "2011-04-01 00:00:00", "2012-04-01 00:00:00", "2013-04-01 00:00:00", "2014-04-01 00:00:00", 
                       "2015-04-01 00:00:00", "2016-04-01 00:00:00", "2017-04-01 00:00:00", 
                       "2018-04-01 00:00:00", "2019-04-01 00:00:00", "2020-04-01 00:00:00", 
                       "2010-10-30 00:00:00", "2011-10-30 00:00:00", "2012-10-30 00:00:00", "2013-10-30 00:00:00", "2014-09-30 00:00:00", 
                       "2015-10-30 00:00:00", "2016-10-30 00:00:00", "2017-10-30 00:00:00", 
                       "2018-10-30 00:00:00", "2019-10-30 00:00:00", "2020-10-30 00:00:00"), 
                       tz = 'Africa/Nairobi')

res.90 <- output.plot[[1]] %>%
  filter(subject_name %in% c("Omondi")) 
res.90$ag.used <- ifelse(res.90$ag.used == 1, 0.75, NA) # adjust ag.used values for plotting for visibility
p90 <- ggplot(res.90, aes(x = date, color = subject_name)) + #, color = name
  geom_point(aes(y = ag.used), color = "grey40", size = .2, alpha = 0.5) +
  #geom_ribbon(aes(ymin = lo.95, ymax = hi.95), size = .5, alpha = 0.4) +
  geom_line(aes(y = mean), size = 1, alpha = 0.7) +
  facet_wrap(~ subject_name) +
  labs(title = "Ag Usage: Volatility and Trend",
       subtitle = "90-Day Moving Average with 95% CI Bands") +
  theme(legend.position="none")

p90 + geom_vline(xintercept=as.numeric(year.cuts, linetype=4)) 

  
  





