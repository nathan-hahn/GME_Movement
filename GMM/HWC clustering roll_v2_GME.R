## Gaussian Mixture Model clustering by roll.max

library(ggplot2)
theme_set(theme_bw() + theme(panel.border = element_rect(), panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title=element_text(size=14), axis.text = element_text(size=12))) 
library(mclust)
library(lubridate)
library(dplyr)
source("GME_functions.R")

#### Prep Data ####
# Load filtered data for individuals that meet temporal and spatial criteria for clustering
movdat.filter <- readRDS('./GMM/GMEcollars_002_res90_2020-09-06.rds')

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
roll.max <- df %>%
  group_by(subject_name) %>%
  filter(is.na(mean) == FALSE) %>%
  summarise(roll.max = max(mean)) 

# Using amplitude table
amplitude <- df %>%
  filter(!is.na(mean)) %>%
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

max.mean <- amplitude %>% ungroup() %>% group_by(subject_name) %>%
  summarise(max.mean = mean(season.max))

##### Tabulate Occupancy Results #####
# Join the results from each occupancy method into a single table of individuals for model fitting
result <- full_join(max, ag.mean) %>% full_join(., roll.max) %>% full_join(., max.mean)
(result)

result$diff <- result$roll.max - result$mean.occupancy

####Model Fitting####

# fit model with mean ag occupancy

# univariate mean
m1 <- Mclust(result$mean.occupancy)
plot(m1, what = "classification",
     xlab = "mean ag occupancy",
     ylab = "cluster")


# fit model with max monthly ag occupancy
m2 <- Mclust(result$max.occupancy)
summary(m2, classification = TRUE)
plot(m2, what = "classification",
     xlab = "max monthly ag occupancy",
     ylab = "cluster")


# univariate rolling max occupancy
t <- result$max.mean
m2.roll <- Mclust(t)
plot(m2.roll, what = "classification",
     xlab = "roll max",
     ylab = "cluster")

# t <- result$roll.max
# t <- t[complete.cases(t)]
# m2.roll <- Mclust(t)
# plot(m2.roll, what = "classification",
#      xlab = "roll max",
#      ylab = "cluster")

# fit model with mean and max occupancy data
both.df <- cbind(result$mean.occupancy, result$max.occupancy)
m3 <- Mclust(both.df, G = 6)
plot(m3, what = "classification",
     xlab = "mean ag occupancy",
     ylab = "max monthly ag occupancy")

# bivariate mean and roll.max
# fit model with mean and roll.max occupancy data
roll.df <- cbind(result$mean.occupancy, result$max.mean)
roll.df <- roll.df[complete.cases(roll.df), ]
m4 <- Mclust(roll.df, G = 4)
plot(m4, what = "classification",
     xlab = "mean ag occupancy",
     ylab = "roll.max ag occupancy")





##### Results #####

# select top based on BIC values and biological realism (m1 and m2)
par(mfrow = c(2,2))
plot(m1, what = "BIC", main = "m1 - Mean Ag")
plot(m2.roll, what = "BIC", main = "m2 - Max Rolling Ag")
plot(m3, what = "BIC", main = "m3 - Mean - Max Monthly Ag")
plot(m4, what = "BIC", main = "m4 - Mean - Max Rolling Ag")

# add classification to results
ag.class.mean <- as.factor(m1$classification)
ag.class.roll <- as.factor(m2.roll$classification)
ag.class.both <- as.factor(m4$classification)
ag.class.both <- factor(ag.class.both)

result$ag.class.mean <- ag.class.mean
result$ag.class.roll <- ag.class.roll
result$ag.class.both <- ag.class.both

# Summarize by individual
clust.result <- result %>%
  dplyr::select(-c(month)) 
(clust.result)

##### Model validation #####
### CV and consistency evaluation

# Mean
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

# Roll.max
train.err.2 <- NULL
test.err.2 <- NULL
cv.2 <- NULL

for(i in 1:20){
  train <- sample(1:nrow(clust.result), size = round(nrow(clust.result)*(2/3)), replace = FALSE)
  
  X.train <- clust.result[train,]
  Class.train <- clust.result$ag.class.roll[train]
  
  X.test <- clust.result[-train,]
  Class.test <- clust.result$ag.class.roll[-train]
  
  mod2 <- MclustDA(X.train$roll.max, Class.train, modelType = "EDDA")
  sum2 <- summary(mod2, newdata = X.test$roll.max, newclass = Class.test)
  
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
mean(cv.2)
mean(test.err.3) # cannot get CV error directly for bivariate


##### Export #####
# new used df with cluster classification. 
# has ag.class.mean and ag.class.both classifications included

# cluster results - lumped over collar lifetime
{outfile <- paste0("./GMM/results/clustResult_GME_", Sys.Date(), ".csv")
  write.csv(clust.result, outfile)}


saveRDS(m1, './GMM/results/mSelect_2020-07-14.rds')


