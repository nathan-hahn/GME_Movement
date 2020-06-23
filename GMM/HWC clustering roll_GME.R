## Gaussian Mixture Model clustering by roll.max

library(ggplot2)
theme_set(theme_bw() + theme(panel.border = element_rect(), panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title=element_text(size=14), axis.text = element_text(size=12))) 
library(mclust)
library(lubridate)
library(dplyr)

source("GME_functions.R")

####Prep Data####
# Filter to individuals that meet temporal and spatial criteria

# movdat <- readRDS('./movdata/GMEcollars_002_used_2020-06-02.rds')
# movdat$subject_name <- ifelse(is.na(movdat$subject_name), movdat$collar_id, movdat$subject_name)

# dataframe for fitting
df <- movdat.filter %>%
  mutate(ag.used = if_else(lc.estes == 1, 1, 0),
         month = month(date)) #%>%
  #filter(minute(date) <= 10 | minute(date) >= 50)

#### Mean occupancy ####
ag.mean <- df %>%
  group_by(subject_name) %>%
  summarise(n = n(),
            mean.ag = mean(ag.used))

(ag.mean)

#### Occupancy by Month ####
# calculate occupancy by month
ag.month <- df %>%
  group_by(site, subject_name, month) %>%
  summarise(n.month = n(),
            month.occupancy = mean(ag.used)) 

# max monthly occupancy

max <- ag.month %>%
  group_by(subject_name) %>%
  slice(which.max(month.occupancy)) %>%
  ungroup() %>%
  mutate(max.occupancy = round(month.occupancy, 5)) %>% dplyr::select(-month.occupancy) #%>%
  #filter(n.month >= 5*24)

# bar plot of ag occupancy by month
p <- ggplot(ag.month, aes(month, month.occupancy)) + geom_bar(stat = "identity") + facet_wrap(.~subject_name) + 
  scale_x_discrete(name ="month", breaks=c(3,6,9))
p


##### Rolling Occupancy #####
# filter for rolling mean calculations - at least 4 months of data
roll.filter <- movdat.filter %>%
  group_by(subject_name) %>%
  tally() %>% 
  filter(n >= 2880) 

roll.df <- filter(df, subject_name %in% roll.filter$subject_name)
split <- split(roll.df, roll.df$subject_name)

# window sizes
window <- c(30*24, 90*24) # in hours. Split to before and after when align = center

roll.output <- NULL
for(i in 1:length(window)){
  # calculate moving window stats - see rollstats function
  roll.output[[i]] <- window_stats(df.list = split, window = window[[i]])
  # remove NAs from ag windows - created by window cutoff
  roll.output[[i]] <- filter(roll.output[[i]], is.na(ag.window) == FALSE) %>% droplevels()
}

# extract max window output. Test for window size
roll.max <- roll.output[[2]] %>%
  group_by(subject_name) %>%
  summarise(roll.max = max(mean)) 


##### Tabulate Occupancy Results #####
# Join the results from each occupancy method into a single table of individuals for model fitting
result <- full_join(max, ag.mean) %>% full_join(., roll.max)
(result)

result$diff <- result$roll.max - result$mean.ag

####Model Fitting####

# fit model with mean ag occupancy
m1 <- Mclust(result$mean.ag)
plot(m1, what = "classification",
     xlab = "mean ag occupancy",
     ylab = "cluster")

# fit model with max monthly ag occupancy
# 2 and 3 group results are almost identical
m2 <- Mclust(result$max.occupancy)
summary(m2, classification = TRUE)
plot(m2, what = "classification",
     xlab = "max monthly ag occupancy",
     ylab = "cluster")

m2.roll <- Mclust(result$roll.max)
plot(m2.roll)

# fit model with mean and max occupancy data
both.df <- cbind(result$mean.ag, result$max.occupancy)
m3 <- Mclust(both.df, G = 4)
plot(m3, what = "classification",
     xlab = "mean ag occupancy",
     ylab = "max monthly ag occupancy")

# fit model with mean and roll.max occupancy data - test window sizes
both.df <- cbind(result$mean.ag, result$max.occupancy)
m3.3 <- Mclust(both.df, G = 3)
plot(m3.3, what = "classification",
     xlab = "mean ag occupancy",
     ylab = "max monthly ag occupancy")


roll.df <- cbind(result$roll.max, result$mean.ag)
roll.df <- roll.df[complete.cases(roll.df), ]
m4 <- Mclust(roll.df, G = 4, modelNames = 'VVE')
plot(m4, what = "classification",
     xlab = "mean ag occupancy",
     ylab = "roll.max ag occupancy")

####Results####

# select top based on BIC values and biological realism (m1 and m2)
plot(m1, what = "BIC")
plot(m2, what = "BIC")
plot(m3, what = "BIC")
plot(m4, what = "BIC")

# add classification to results
ag.class.mean <- as.factor(m1$classification)
ag.class.both.3 <- as.factor(m3.3$classification)
ag.class.both.3 <- factor(ag.class.both.3, levels(ag.class.both.3)[c(1,3,2)], labels = c(1,2,3)) # FOR GME 002 ONLY
ag.class.both.4 <- as.factor(m3$classification)
ag.class.both.4 <- factor(ag.class.both.4, levels(ag.class.both.4)[c(1,3,4,2)], labels = c(1,2,3,4)) # FOR GME 002 ONLY
ag.class.roll <- as.factor(m4$classification)
ag.class.roll <- factor(ag.class.roll, levels(ag.class.roll)[c(1,3,2,4)], labels = c(1,2,3,4))
#ag.class.both <- factor(ag.class.both, levels(ag.class.both)[c(1,2,4,3)], labels = c(1,2,3,4)) # FOR GME  001 ONLY
#ag.class.both <- factor(ag.class.both,levels(ag.class.both)[c(1,4,2,3)], labels = c(1,2,3,4)) #FOR MEP ONLY

result$ag.class.mean <- ag.class.mean
result$ag.class.both.3 <- ag.class.both.3
result$ag.class.both.4 <- ag.class.both.4

roll.df <- result[complete.cases(result), ]
roll.df$ag.class.roll <- ag.class.roll
roll.df <- dplyr::select(roll.df, c(subject_name, ag.class.roll)) 
result <- full_join(result, roll.df)

# Summarize by individual
clust.result <- result %>%
  dplyr::select(-c(month)) 
(clust.result)

##### Model validation/inspecition #####

### bootstrap - m3 vs m1

boot.m3 <- MclustBootstrap(m3, type = "bs")
summary(boot.m3, what = "ci")
par(mfrow=c(4,3))
plot(boot.m3, what = "mean")

boot.m1 <- MclustBootstrap(m1, type = "bs")
summary(boot.m1, what = "ci")

### crossvalidation - m3 with 3 vs 4 groups

# G = 3
trainData <- cbind(clust.result$mean.ag, clust.result$max.occupancy)
trainClass <- clust.result$ag.class.both.3

models <- mclust.options()$emModelNames
tab <- matrix(NA, nrow = length(models), ncol = 5)
rownames(tab) <- models
colnames(tab) <- c("model","BIC", "CV_err", "se", "G")
for(i in seq(models)) {
  tab[i,1] <- models[i]
  mod <- MclustDA(trainData, trainClass,
                  modelType = "EDDA", modelNames = models[i])
  tab[i,2] <- mod$bic
  t <- cvMclustDA(mod, nfold = 10, verbose = FALSE)
  tab[i,3] <- t$error
  tab[i,4] <- t$se
  tab[i,5] <- "3"
}
tab.3 <- as.data.frame(tab)

# G = 4
trainData <- cbind(clust.result$mean.ag, clust.result$max.occupancy)
trainClass <- clust.result$ag.class.both.4

models <- mclust.options()$emModelNames
tab <- matrix(NA, nrow = length(models), ncol = 5)
rownames(tab) <- models
colnames(tab) <- c("model","BIC", "CV_err", "se", "G")
for(i in seq(models)) {
  tab[i,1] <- models[i]
  mod <- MclustDA(trainData, trainClass,
                  modelType = "EDDA", modelNames = models[i])
  tab[i,2] <- mod$bic
  t <- cvMclustDA(mod, nfold = 10, verbose = FALSE)
  tab[i,3] <- t$error
  tab[i,4] <- t$se
  tab[i,5] <- "4"
}
tab.4 <- as.data.frame(tab)

# roll, G = 4
trainData <- cbind(clust.result$mean.ag, clust.result$roll.max)
trainData <- trainData[complete.cases(trainData), ]
trainClass <- clust.result$ag.class.roll
trainClass <- trainClass[complete.cases(trainClass)]

models <- mclust.options()$emModelNames
tab <- matrix(NA, nrow = length(models), ncol = 5)
rownames(tab) <- models
colnames(tab) <- c("model","BIC", "CV_err", "se", "G")
for(i in seq(models)) {
  tab[i,1] <- models[i]
  mod <- MclustDA(trainData, trainClass,
                  modelType = "EDDA", modelNames = models[i])
  tab[i,2] <- mod$bic
  t <- cvMclustDA(mod, nfold = 10, verbose = FALSE)
  tab[i,3] <- t$error
  tab[i,4] <- t$se
  tab[i,5] <- "roll.4"
}

tab.roll <- as.data.frame(tab)

# create single table and filter to 5 largest BIC values
tab.cv <- rbind(tab.3, tab.4, tab.roll) %>%
  mutate(BIC = as.numeric(as.character(BIC)), CV_err = as.numeric(as.character(CV_err))) %>%
  top_n(n = 5, BIC) %>% arrange(desc(BIC))
tab.cv


##### Fit Top Model to Data #####
cov.df <- cbind(clust.result$mean.ag, clust.result$roll.max)
cov.df <- cov.df[complete.cases(cov.df), ]
mSelect <- Mclust(data = cov.df, modelNames = 'VVE', G = 4)
plot(mSelect,  what = "classification",
     xlab = "mean ag occupancy",
     ylab = "max monthly ag occupancy")



##### Fit Top Model to data subsets $$$$$
test.df <- cbind(test$mean.occupancy, test$max.occupancy)

pred <- predict.Mclust(mSelect, test.df)

plot(test.df, col = mclust.options("classPlotColors")[pred$classification], pch = 15, cex = 0.5)
points(cov.df, pch = mSelect$classification)

##### Summary Outputs #####
# Summarize by class
clust.summary <- result %>%
  group_by(ag.class.both.4) %>%
  summarise(n = length(ag.class.both.4),
            mean = mean(mean.occupancy),
            sd = sd(mean.occupancy),
            se = sd/sqrt(n),
            lower = mean - sd,
            upper = mean + sd)
(clust.summary)

clust.df <- result %>%
  dplyr::select(subject_name, ag.class.both.4) %>%
  inner_join(.,df, by = "subject_name") %>%
  dplyr::select(-month) %>%
  as.data.frame()


ggplot(clust.summary) + geom_bar(aes(x = ag.class.both.4, y = mean), stat = "identity") + 
  geom_errorbar(aes(x = ag.class.both.4, ymin = lower, ymax = upper), width = .2) + ylab("mean ag occupancy") + xlab("ag usage cluster")


##### Export #####
# new used df with cluster classification. 
# has ag.class.mean and ag.class.both classifications included
{outfile <- paste0("./GMM/GMEcollars_002_usedClust_", Sys.Date(), ".rds")
saveRDS(clust.df, outfile)}

{outfile <- paste0("./GMM/GMEcollars_002_usedClust_", Sys.Date(), ".csv")
  write.csv(clust.df, outfile)}

# cluster results
{outfile <- paste0("./GMM/clustResult_GME_", Sys.Date(), ".csv")
  write.csv(clust.result, outfile)}

# cluster summary (mean or both)
{outfile <- paste0("./GMM/GMEcollars_002_clustSummary_", Sys.Date(), ".csv")
  write.csv(clust.result, outfile)}

{outfile <- paste0("./GMM/GMEcollars_002_clustSummary_", Sys.Date(), ".csv")
  write.csv(clust.summary, outfile)}


