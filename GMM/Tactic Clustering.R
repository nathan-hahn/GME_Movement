###### Gaussian Mixture Model clustering #####
# Tactic cluster fitting and evaluation using GMMs in the mclust package
# Final output stored as a result file and tactic cluster IDs added to the tracking data

library(ggplot2)
theme_set(theme_bw() + theme(panel.border = element_rect(), panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title=element_text(size=14), axis.text = element_text(size=12))) 
library(mclust)
library(lubridate)
library(dplyr)
source("GME_functions.R")

##### Data Prep #####
movdat <- readRDS('./movdata/GMEcollars_004_usedFilter_2022-04-01.rds')
movdat$ag.used <- ifelse(movdat$lc.estes == 1, 1, 0) 

##### Add Crop-Based Year Cuts #####
# create year cuts based on crop seasons across mara and serengeti. 
# calendar years split Dec-Feb crop seasons

# Yearly cuts - set april cut points
year.cuts <- ymd_hms(c("2010-04-01 00:00:00", "2011-04-01 00:00:00", "2012-04-01 00:00:00", "2013-04-01 00:00:00", "2014-04-01 00:00:00", 
                       "2015-04-01 00:00:00", "2016-04-01 00:00:00", "2017-04-01 00:00:00", 
                       "2018-04-01 00:00:00", "2019-04-01 00:00:00", "2020-04-01 00:00:00", "2021-04-01 00:00:00"), tz = 'Africa/Nairobi')
year.names <- c("2010", "2011", '2012', "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020")
# get yearly cuts using april cut points
rng <- cut( movdat$date, breaks = c(year.cuts), include.lowest = T)
rng.name <- cut(movdat$date, breaks = c(year.cuts), include.lowest = T, labels = year.names)

# create variables
movdat$cut.date <- rng
movdat$year.cuts <- rng.name
movdat$month <- month(movdat$date)

# check
summary(movdat$year.cuts)

# export movdata with season and year cuts
saveRDS(movdat, paste0('./movdata/GMEcollars_004_usedFiltercuts_', Sys.Date(),'.rds'))


# dataframe for fitting - add custom filters here as needed (e.g. region-specific)
df <- movdat
  
##### Mean Ag use #####
ag.mean <- df %>%
  filter(!is.na(ag.used)) %>%
  group_by(subject_name) %>%
  summarise(n = n(),
            mean.occupancy = mean(ag.used))

(ag.mean)

##### Ag Use by Month #####
# calculate occupancy by month
ag.month <- df %>%
  group_by(site, subject_name, month) %>%
  summarise(n.month = n(),
            month.occupancy = mean(ag.used)) 

# max monthly occupancy
max <- ag.month %>%
  group_by(subject_name) %>%
  filter(n.month >= 5*24) %>%
  slice(which.max(month.occupancy)) %>%
  ungroup() %>%
  mutate(max.occupancy = round(month.occupancy, 5)) %>% dplyr::select(-month.occupancy) 

# bar plot of ag occupancy by month
p <- ggplot(ag.month, aes(month, month.occupancy)) + geom_bar(stat = "identity") + facet_wrap(.~subject_name) + 
  scale_x_discrete(name ="month", breaks=c(3,6,9))
p



##### 90-day Mean Ag Use - Moving window stats #####
roll.filter <- df %>%
  group_by(subject_name) %>%
  tally() %>% 
  filter(n >= 24*30*3) # months

roll.df <- filter(df, subject_name %in% roll.filter$subject_name) %>% droplevels()
split <- split(roll.df, roll.df$subject_name)

window <- c(90*24) # in hours. Split to before and after when align = center

output <- NULL
for(i in 1:length(window)){
  # calculate moving window stats - see rollstats function
  output[[i]] <- window_stats(df.list = split, window = window[[i]], align = 'center')
  
}

##### Rolling Occupancy #####
# max.max - the max 90-day ag use for each individual
max.max <- output[[1]] %>%
  group_by(subject_name) %>%
  filter(is.na(mean) == FALSE) %>%
  summarise(max.max = max(mean)) 

# Rolling stats by individual-year - used to calculate the max.mean. 
# additional stats are useful for context and to eval data before model fitting
roll.stats <- output[[1]] %>%
  filter(!is.na(mean)) %>%
  group_by(subject_name, year.cuts) %>%
  mutate(year.mean = mean(ag.used),
         year.max = max(mean), # value of interest for max.mean
         year.min = min(mean),
         year.diff = year.max - year.min,
         year.begin = min(date),
         year.end = max(date)) %>%
  group_by(subject_name, year.cuts, year.begin, year.end, year.mean, year.max, year.min, year.diff) %>%
  #filter to individual seasons with at least 1 month's worth of fixes and check that they overlap a crop season
  tally() %>% filter(n > 500) %>% droplevels()

# max.mean value calculation for each individual - mean of max values across years
max.mean <- roll.stats %>% ungroup() %>% group_by(subject_name) %>%
  summarise(max.mean = mean(year.max))

##### Tabulate Ag Use Results #####
# Join the results from each occupancy method into a single table of individuals for model fitting
result <- full_join(max, ag.mean) %>% full_join(., max.max) %>% full_join(., max.mean)

# diff value calculation - tried different equations for diff
result$diff <- result$max.mean - result$mean.occupancy

(result)

####Model Fitting####

# univariate mean
m1.mean <- Mclust(result$mean.occupancy)
plot(m1.mean, what = "classification",
     xlab = "mean ag occupancy",
     ylab = "cluster")

# univariate rolling max occupancy
t <- result$max.mean
m1.max <- Mclust(t)
plot(m1.max, what = "classification",
     xlab = "max",
     ylab = "cluster")

# bivariate mean and roll.max
# fit model with mean and roll.max occupancy data
roll.df <- cbind(result$mean.occupancy, result$max.mean)
roll.df <- roll.df[complete.cases(roll.df), ]
m2 <- Mclust(roll.df, G = 4)
plot(m2, what = "classification",
     xlab = "mean ag occupancy",
     ylab = "max ag occupancy")


##### Results #####

# select top based on BIC values and biological realism (m1 and m2)
par(mfrow = c(2,2))
plot(m1.mean, what = "BIC", main = "m1 - Mean Ag")
plot(m1.max, what = "BIC", main = "m1 - Max Rolling Ag")
plot(m2, what = "BIC", main = "m2 - Mean - Max Rolling Ag")

# add classification to results
ag.class.mean <- as.factor(m1.mean$classification)
ag.class.max <- as.factor(m1.max$classification)
ag.class.both <- as.factor(m2$classification)
ag.class.both <- factor(ag.class.both)

result$ag.class.mean <- ag.class.mean
result$ag.class.max <- ag.class.max
result$ag.class.both <- ag.class.both

# Summarize by individual
clust.result <- result %>%
  dplyr::select(-c(month)) 
(clust.result)

##### Model validation #####
### CV and consistency evaluation

set.seed(401)

# m1.mean
train.err.1 <- NULL
test.err.1 <- NULL
cv.1 <- NULL
for(i in 1:20){
  train <- sample(1:nrow(clust.result), size = round(nrow(clust.result)*(2/3)), replace = FALSE)
  
  X.train <- clust.result[train,]
  Class.train <- clust.result$ag.class.mean[train]
  
  #table(Class.train)
  
  X.test <- clust.result[-train,]
  Class.test <- clust.result$ag.class.mean[-train]
  #table(Class.test)
  
  mod1 <- MclustDA(X.train$mean.occupancy, Class.train, modelType = "EDDA")
  sum1 <- summary(mod1, newdata = X.test$mean.occupancy, newclass = Class.test)
  
  train.err.1[i] <- sum1$err
  test.err.1[i] <- sum1$err.newdata
  cv.1[i] <- cvMclustDA(mod1, nfold = 10, metric = "error")$error
  
}

# m1.max
train.err.2 <- NULL
test.err.2 <- NULL
cv.2 <- NULL

for(i in 1:20){
  train <- sample(1:nrow(clust.result), size = round(nrow(clust.result)*(2/3)), replace = FALSE)
  
  X.train <- clust.result[train,]
  Class.train <- clust.result$ag.class.max[train]
  
  X.test <- clust.result[-train,]
  Class.test <- clust.result$ag.class.max[-train]
  
  mod2 <- MclustDA(X.train$max.mean, Class.train, modelType = "EDDA")
  sum2 <- summary(mod2, newdata = X.test$max.mean, newclass = Class.test)
  
  train.err.2[i] <- sum2$err
  test.err.2[i] <- sum2$err.newdata
  cv.2[i] <- cvMclustDA(mod2, nfold = 10, metric = "error")$error
  
}


# Both
train.err.3 <- NULL
test.err.3 <- NULL
cv.3 <- NULL

for(i in 1:20) {
  train <- sample(1:nrow(clust.result), size = round(nrow(clust.result)*(3/4)), replace = FALSE)
  
  X.train <- clust.result[train,]
  Class.train <- clust.result$ag.class.both[train]
  
  X.test <- clust.result[-train,]
  Class.test <- clust.result$ag.class.both[-train]
  
  mod3 <- MclustDA(cbind(X.train$roll.max, X.train$mean.occupancy), Class.train, modelType = "EDDA")
  sum3 <- summary(mod3, newdata = cbind(X.test$roll.max, X.test$mean.occupancy), newclass = Class.test)
  
  train.err.3[i] <- sum3$err
  test.err.3[i] <- sum3$err.newdata

}

# Training Error
mean(train.err.1)
mean(train.err.2)
mean(train.err.3)

# CROSS VALIDATION RESULTS
mean(cv.1)
sd(cv.1)
mean(cv.2)
sd(cv.2)
mean(test.err.3) 
sd(test.err.3)


##### Export #####
# new used df with cluster classification. 
# has ag.class.mean and ag.class.both classifications included

# cluster results - lumped over collar lifetime
{outfile <- paste0("./GMM/results/clustResult_GME_", Sys.Date(), ".csv")
  write.csv(clust.result, outfile)}

{outfile <- paste0("./GMM/results/clustResult_GME_", Sys.Date(), ".rds")
  saveRDS(clust.result, outfile)}

# save top model for further analysis
saveRDS(m1.mean, paste0('./GMM/results/mSelect_', Sys.Date(), '.rds'))


