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
movdat.filter <- readRDS('./movdata/GMEcollars_002_usedFilter_2020-06-30.rds')

# dataframe for fitting
df <- movdat.filter %>%
  mutate(ag.used = if_else(lc.estes == 1, 1, 0),
         month = month(date)) # %>% group_by(subject_name) %>%
  # # filter to indiv with 1+ year of data
  # mutate(start.time = min(date), end.time = max(date), 
  #        difftime = difftime(end.time, start.time, units = 'days')) %>% filter(difftime > 365) %>%
  # ungroup()

t <- df %>% group_by(subject_name, start.time, end.time, difftime) %>% tally()
  




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
# filter for rolling mean calculations - at least 4 months of data
roll.filter <- df %>%
  group_by(subject_name) %>%
  tally() %>% 
  filter(n >= 2880) 

roll.df <- filter(df, subject_name %in% roll.filter$subject_name)
split <- split(roll.df, roll.df$subject_name)

# assign window sizes - check alignment of window
window <- c(90*24) # days*hrs

roll.output <- NULL
for(i in 1:length(window)){
  # calculate moving window stats - see rollstats function
  roll.output[[i]] <- window_stats(df.list = split, window = window[[i]], align = "center")
  # remove NAs from ag windows - created by window cutoff
  roll.output[[i]] <- filter(roll.output[[i]], is.na(ag.window) == FALSE) %>% droplevels()
}

# extract max window output. Test for window size
roll.max <- roll.output[[1]] %>%
  group_by(subject_name) %>%
  summarise(roll.max = max(mean)) 


##### Tabulate Occupancy Results #####
# Join the results from each occupancy method into a single table of individuals for model fitting
result <- full_join(max, ag.mean) %>% full_join(., roll.max) 
(result)

result$diff <- result$roll.max - result$mean.occupancy

####Model Fitting####

# fit model with mean ag occupancy

m1 <- Mclust(result$mean.occupancy, G = 3)
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
m2.roll <- Mclust(t, G = 3)
plot(m2.roll, what = "classification",
     xlab = "roll max",
     ylab = "cluster")

# fit model with mean and max occupancy data
both.df <- cbind(result$mean.occupancy, result$max.occupancy)
m3 <- Mclust(both.df)
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
m5 <- Mclust(diff.df, G = 3)
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
ag.class.both <- factor(ag.class.both, levels(ag.class.both)[c(1,3,2,4)], labels = c(1,2,3,4))

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

# Mean Ag, G = 3
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
  tab[i,5] <- "3"
  tab[i,6] <- "m1"
}
tab.m1 <- as.data.frame(tab)

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
  tab[i,5] <- "3"
  tab[i,6] <- "m4"
}

tab.m4 <- as.data.frame(tab)

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
# Summarize by class
clust.summary <- clust.result %>%
  group_by(ag.class.both) %>%
  summarise(n = length(mean.occupancy),
            mean = mean(mean.occupancy),
            sd = sd(mean.occupancy),
            se = sd/sqrt(n),
            lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
            upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se)
(clust.summary)

# clust.df <- clust.result %>%
#   dplyr::select(subject_name, ag.class.both.4) %>%
#   inner_join(.,df, by = "subject_name") %>%
#   dplyr::select(-month) %>%
#   as.data.frame()

ggplot(clust.summary) + geom_bar(aes(x = ag.class.roll, y = mean), stat = "identity") + 
  geom_errorbar(aes(x = ag.class.roll, ymin = lower, ymax = upper), width = .2) + ylab("Mean Occupancy") + xlab("Tactic Cluster")

ggplot(clust.result, aes(x = mean.occupancy, y = ag.class.roll)) + geom_boxplot() + 
  xlab('Mean Occupancy') + ylab('Tactic Cluster')

##### Export #####
# new used df with cluster classification. 
# has ag.class.mean and ag.class.both classifications included
{outfile <- paste0("./GMM/GMEcollars_002_usedClust_", Sys.Date(), ".rds")
saveRDS(clust.df, outfile)}

{outfile <- paste0("./GMM/GMEcollars_002_usedClust_", Sys.Date(), ".csv")
  write.csv(clust.df, outfile)}

# cluster results
{outfile <- paste0("./GMM/results/clustResult_GME_", Sys.Date(), ".csv")
  write.csv(clust.result, outfile)}

# cluster summary (mean or both)
{outfile <- paste0("./GMM/GMEcollars_002_clustSummary_", Sys.Date(), ".csv")
  write.csv(clust.result, outfile)}

{outfile <- paste0("./GMM/GMEcollars_002_clustSummary_", Sys.Date(), ".csv")
  write.csv(clust.summary, outfile)}


