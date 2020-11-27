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

#### Prep Data ####
# Load filtered data for individuals that meet temporal and spatial criteria for clustering
movdat.filter <- readRDS('./GMM/movdata/GMEcollars_003_res90_2020-10-30.rds')

# dataframe for fitting - add custom filters here as needed
df <- movdat.filter 
  
##### Mean occupancy #####
ag.mean <- df %>%
  group_by(subject_name) %>%
  summarise(n = n(),
            mean.occupancy = mean(ag.used))

(ag.mean)

##### Occupancy by Month #####
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


##### Rolling Occupancy #####
max.max <- df %>%
  group_by(subject_name) %>%
  filter(is.na(mean) == FALSE) %>%
  summarise(max.max = max(mean)) 

# Rolling stats by individual-year - used to calculate the max.mean. 
# additional stats are useful for context and to eval data before model fitting
amplitude <- df %>%
  filter(!is.na(mean)) %>%
  group_by(subject_name, year.cuts) %>%
  mutate(season.mean = mean(ag.used),
         season.max = max(mean), # value of interest for max.mean
         season.min = min(mean),
         season.diff = season.max - season.min,
         season.begin = min(date),
         season.end = max(date)) %>%
  group_by(subject_name, year.cuts, season.begin, season.end, season.mean, season.max, season.min, season.diff) %>%
  #filter to individual seasons with at least 1 month's worth of fixes
  tally() %>% filter(n > 500) %>% droplevels()

# max.mean value calculation for each individual
max.mean <- amplitude %>% ungroup() %>% group_by(subject_name) %>%
  summarise(max.mean = mean(season.max))

##### Tabulate Occupancy Results #####
# Join the results from each occupancy method into a single table of individuals for model fitting
result <- full_join(max, ag.mean) %>% full_join(., max.max) %>% full_join(., max.mean)
(result)

# diff value calculation - tried different equations for diff
result$diff <- result$max.mean - result$mean.occupancy

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
saveRDS(m1, paste0('./GMM/results/mSelect_', Sys.Date(), '.rds'))


