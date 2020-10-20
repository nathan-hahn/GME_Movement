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
movdat.filter <- readRDS('./GMM/GMEcollars_002_res90_2020-09-03.rds')

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


##### Tabulate Occupancy Results #####
# Join the results from each occupancy method into a single table of individuals for model fitting
result <- full_join(max, ag.mean) %>% full_join(., roll.max) 
(result)

result$diff <- result$roll.max - result$mean.occupancy

####Model Fitting####

# fit model with mean ag occupancy

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

# fit model with rolling max occupancy
t <- result$roll.max
t <- t[complete.cases(t)]
m2.roll <- Mclust(t)
plot(m2.roll, what = "classification",
     xlab = "roll max",
     ylab = "cluster")

# fit model with mean and max occupancy data
both.df <- cbind(result$mean.occupancy, result$max.occupancy)
m3 <- Mclust(both.df, G = 6)
plot(m3, what = "classification",
     xlab = "mean ag occupancy",
     ylab = "max monthly ag occupancy")

# fit model with mean and roll.max occupancy data
roll.df <- cbind(result$mean.occupancy, result$roll.max)
roll.df <- roll.df[complete.cases(roll.df), ]
m4 <- Mclust(roll.df)
plot(m4, what = "classification",
     xlab = "mean ag occupancy",
     ylab = "roll.max ag occupancy")



diff.df <- cbind(result$roll.max, result$diff)
diff.df <- diff.df[complete.cases(diff.df), ]
m5 <- Mclust(diff.df)
plot(m5, what = "classification",
     xlab = "mean ag occupancy",
     ylab = "diff")

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

##### Model validation/inspecition #####

### bootstrap - m3 vs m1

boot.m4 <- MclustBootstrap(m4, type = "bs")
summary(boot.m4, what = "ci")
par(mfrow=c(2,4))
plot(boot.m4, what = "mean")

boot.m1 <- MclustBootstrap(m1, type = "bs")
summary(boot.m1, what = "ci")
plot(boot.m1, what = "mean")

### crossvalidation - Mean Ag, Mean=Max, Mean=Roll.Max

# Mean
set.seed(1992)
train <- sample(1:nrow(clust.result), size = round(nrow(clust.result)*(2/3)), replace = FALSE)
X.train <- clust.result[train,]
Class.train <- clust.result$ag.class.mean[train]
#table(Class.train)

X.test <- clust.result[-train,]
Class.test <- clust.result$ag.class.mean[-train]
#table(Class.test)

mod1 <- MclustDA(X.train$mean.occupancy, Class.train, modelType = "EDDA")
sum1 <- summary(mod1, newdata = X.test$mean.occupancy, newclass = Class.test)



# Roll.max
set.seed(54)
train <- sample(1:nrow(clust.result), size = round(nrow(clust.result)*(2/3)), replace = FALSE)
X.train <- clust.result[train,]
Class.train <- clust.result$ag.class.roll[train]
#table(Class.train)

X.test <- clust.result[-train,]
Class.test <- clust.result$ag.class.roll[-train]
#table(Class.test)

mod2 <- MclustDA(X.train$roll.max, Class.train, modelType = "EDDA")
sum2 <- summary(mod2, newdata = X.test$roll.max, newclass = Class.test)


# Both
set.seed(112)
train <- sample(1:nrow(clust.result), size = round(nrow(clust.result)*(2/3)), replace = FALSE)
X.train <- clust.result[train,]
Class.train <- clust.result$ag.class.both[train]
#table(Class.train)

X.test <- clust.result[-train,]
Class.test <- clust.result$ag.class.both[-train]
#table(Class.test)

mod3 <- MclustDA(cbind(X.train$roll.max, X.train$mean.occupancy), Class.train, modelType = "EDDA")
sum3 <- summary(mod3, newdata = cbind(X.test$roll.max, X.test$mean.occupancy), newclass = Class.test)

sum1
sum2
sum3

cvMclustDA(mod1)
cvMclustDA(mod2)
cvMclustDA(mod3)







# Mean Ag, G = 4
trainData <- cbind(clust.result$mean.occupancy)
trainClass <- clust.result$ag.class.mean

models <- c("E", "V")
tab <- matrix(NA, nrow = length(models), ncol = 6)
rownames(tab) <- models
colnames(tab) <- c("model","BIC", "CV_err", "se", "G", "mod")
for(i in seq(models)) {
  tab[i,1] <- models[i]
  mod <- MclustDA(trainData, trainClass,
                  modelType = "EDDA", modelNames = models[i])
  tab[i,2] <- mod$bic
  t <- cvMclustDA(mod, nfold = 10, verbose = FALSE)
  tab[i,3] <- t$error
  tab[i,4] <- t$se
  tab[i,5] <- "4"
  tab[i,6] <- "m1"
}
tab.m1 <- as.data.frame(tab)
tab.m1

# G = 4
trainData <- cbind(clust.result$mean.occupancy, clust.result$max.occupancy)
trainClass <- clust.result$ag.class.max

models <- mclust.options()$emModelNames
tab <- matrix(NA, nrow = length(models), ncol = 6)
rownames(tab) <- models
colnames(tab) <- c("model","BIC", "CV_err", "se", "G", "mod")
for(i in seq(models)) {
  tab[i,1] <- models[i]
  mod <- MclustDA(trainData, trainClass,
                  modelType = "EDDA", modelNames = models[i])
  tab[i,2] <- mod$bic
  t <- cvMclustDA(mod, nfold = 10, verbose = FALSE)
  tab[i,3] <- t$error
  tab[i,4] <- t$se
  tab[i,5] <- "3"
  tab[i,6] <- "m3"
}
tab.m3 <- as.data.frame(tab)

# roll, G = 4
trainData <- cbind(clust.result$mean.occupancy, clust.result$roll.max)
trainData <- trainData[complete.cases(trainData), ]
trainClass <- clust.result$ag.class.roll
trainClass <- trainClass[complete.cases(trainClass)]

models <- mclust.options()$emModelNames
tab <- matrix(NA, nrow = length(models), ncol = 6)
rownames(tab) <- models
colnames(tab) <- c("model","BIC", "CV_err", "se", "G", "mod")
for(i in seq(models)) {
  tab[i,1] <- models[i]
  mod <- MclustDA(trainData, trainClass,
                  modelType = "EDDA", modelNames = models[i])
  tab[i,2] <- mod$bic
  t <- cvMclustDA(mod, nfold = 10, verbose = FALSE)
  tab[i,3] <- t$error
  tab[i,4] <- t$se
  tab[i,5] <- "6"
  tab[i,6] <- "m4"
}

tab.m4 <- as.data.frame(tab)
tab.m4

# create single table and filter to 5 largest BIC values
tab.cv <- rbind(tab.m1, tab.m3, tab.m4) %>%
  mutate(BIC = as.numeric(as.character(BIC)), CV_err = as.numeric(as.character(CV_err))) %>%
  top_n(n = 5, BIC) %>% arrange(desc(BIC))
tab.cv

 
##### Fit Top Model to Data #####
cov.df <- cbind(clust.result$mean.occupancy, clust.result$roll.max)
cov.df <- cov.df[complete.cases(cov.df), ]
m.max <- Mclust(data = cov.df, modelNames = 'VVE', G = 4)

par(mfrow = c(1,2))
plot(m.max,  what = "classification",
     xlab = "Mean Agricultural Occupancy",
     ylab = "Max Rolling Agricultural Occupancy")

plot(m.max,  what = "uncertainty",
     xlab = "Mean Agricultural Occupancy",
     ylab = "Max Rolling Agricultural Occupancy")


cov.df <- clust.result$mean.occupancy
m.mean <- Mclust(data = cov.df, modelNames = 'V', G = 4)
plot(m.mean, what = "classification",
     xlab = "Mean Agricultural Occupancy")


# check out dimension reduction
t <- MclustDR(m4)
#plot(t)


##### Summary Outputs #####
ggplot(clust.summary) + geom_bar(aes(x = ag.class.roll, y = mean), stat = "identity") + 
  geom_errorbar(aes(x = ag.class.roll, ymin = lower, ymax = upper), width = .2) + ylab("Mean Occupancy") + xlab("Tactic Cluster")

ggplot(clust.result, aes(x = mean.occupancy, y = ag.class.roll)) + geom_boxplot() + 
  xlab('Mean Occupancy') + ylab('Tactic Cluster')

##### Export #####
# new used df with cluster classification. 
# has ag.class.mean and ag.class.both classifications included

# cluster results - lumped over collar lifetime
{outfile <- paste0("./GMM/results/clustResult_GME_", Sys.Date(), ".csv")
  write.csv(clust.result, outfile)}


saveRDS(m1, './GMM/results/mSelect_2020-07-14.rds')


