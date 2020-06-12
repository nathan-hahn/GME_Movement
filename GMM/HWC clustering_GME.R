## Gaussian Mixture Model clustering for agricultural usage

## DONE: Make summary table by ag class
## TODO: Add startStop dates to the metadata file and attach to individual-level summary
## DONE: Make summary table by individual, with mean ag class rating, and add metadata (name and sex)
## DONE: Add cluster results as a covariate to movement df

library(ggplot2)
theme_set(theme_bw() + theme(panel.border = element_rect(), panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title=element_text(size=14), axis.text = element_text(size=12))) 
library(mclust)
library(lubridate)
library(dplyr)


####Calculate Ag occupancy####
movdat <- readRDS('./movdata/GMEcollars_002_used_2020-06-02.rds')
# mep adjustments
movdat <- movdat[-which(movdat$subject_name%in%c("Tunai","Bettye","Nancy", "Wilbur", 
                             "Julia", "Ndorre", "Nkoidila", "Santiyan", "Earhart",
                             "Rudisha", "Naeku", "Olkeri", "ST2010-1441")),]
# grumeti adjustments
movdat <- movdat[-which(movdat$subject_name%in%c("Nyanza")),]
movdat$subject_name <- ifelse(is.na(movdat$subject_name), movdat$collar_id, movdat$subject_name)
movdat <-subset(movdat, minute(movdat$date) <= 10 | minute(movdat$date) >= 50)

# dataframe for fitting
df <- movdat
# df <- movdat %>%
#   filter(site == "mep")

####Mean occupancy####
fixes <- df %>%
  group_by(subject_name) %>%
  tally()

ag.mean <- df %>%
  group_by(subject_name) %>%
  tally(lc.estes == 1) %>%
  droplevels() %>%
  ungroup() %>%
  mutate(mean.ag = n/fixes$n*100)

(ag.mean)

####Occupancy by Month####
# get total number of fixes for each month - summed across all fixes for each individual
# get number of fixes within ag for each month
# bind with total fixes per month
# calculate occupancy by month

df$month <- format(df$date, "%m")
fixes.m <- df %>%
  group_by(site, subject_name, month) %>%
  tally() 

ag.month <- df %>%
  group_by(site, subject_name, month) %>%
  tally(lc.estes == 1) %>% # tally fixes in ag
  droplevels() %>%
  ungroup() %>%
  rename(n.ag = n) %>%
  mutate(month = as.factor(month)) %>%
  mutate(n = fixes.m$n) %>%
  arrange(subject_name, month) %>%
  mutate(month.occupancy = n.ag/n)

(ag.month) # the month for each individual where ag occupancy is highest (across all)

# bar plot of ag occupancy by month
p <- ggplot(ag.month, aes(month, month.occupancy)) + geom_bar(stat = "identity") + facet_wrap(.~name) + 
  scale_x_discrete(name ="month", breaks=c(3,6,9))
p

# filter to each individuals max monthly ag occupancy 
result <- ag.month %>%
  #filter(site == "mep") %>% ## filter by site if needed! 
  group_by(subject_name) %>%
  slice(which.max(month.occupancy)) %>%
  ungroup() %>%
  mutate(
    max.occupancy = round(month.occupancy, 5), 
    mean.occupancy = round(ag.mean$mean.ag/100, 5), 
    month = as.numeric(month))
result$month.occupancy <- NULL
(result)
####Model Fitting####

# fit model with mean ag occupancy
m1 <- Mclust(result$mean.occupancy, G = 3)
plot(m1, what = "classification",
     xlab = "mean ag occupancy (*100)",
     ylab = "cluster")

# fit model with max monthly ag occupancy
# 2 and 3 group results are almost identical
m2 <- Mclust(result$max.occupancy)
summary(m2, classification = TRUE)
plot(m2, what = "classification",
     xlab = "max monthly ag occupancy (*100)",
     ylab = "cluster")

# fit model with mean and max occupancy data (2D)
# produces 9 clusters
both.df <- cbind(result$mean.occupancy*100, result$max.occupancy*100)
m3 <- Mclust(both.df)
plot(m3, what = "classification",
     xlab = "mean ag occupancy",
     ylab = "max monthly ag occupancy")


m3.3 <- Mclust(both.df, G = 3)


####Results####

# select top based on BIC values and biological realism (m1 and m2)
plot(m1, what = "BIC")
plot(m2, what = "BIC")
plot(m3, what = "BIC")

# add classification to results
ag.class.mean <- as.factor(m1$classification)
ag.class.both.3 <- as.factor(m3.3$classification)
ag.class.both.3 <- factor(ag.class.both.3, levels(ag.class.both.3)[c(1,3,2)], labels = c(1,2,3)) # FOR GME 002 ONLY
ag.class.both.4 <- as.factor(m3$classification)
ag.class.both.4 <- factor(ag.class.both.4, levels(ag.class.both.4)[c(1,4,3,2)], labels = c(1,2,3,4)) # FOR GME 002 ONLY
#ag.class.both <- factor(ag.class.both, levels(ag.class.both)[c(1,2,4,3)], labels = c(1,2,3,4)) # FOR GME  001 ONLY
#ag.class.both <- factor(ag.class.both,levels(ag.class.both)[c(1,4,2,3)], labels = c(1,2,3,4)) #FOR MEP ONLY

result$ag.class.mean <- ag.class.mean
result$ag.class.both.3 <- ag.class.both.3
result$ag.class.both.4 <- ag.class.both.4

# Summarize by individual
clust.result <- result %>%
  dplyr::select(-c(n.ag, n)) %>%
  dplyr::select(-month, max.occupancy) 
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

both.df <- cbind(clust.result$mean.occupancy*100, clust.result$max.occupancy*100)
class <- clust.result$ag.class.both
m3.da <- MclustDA(both.df, class, modelType = 'EDDA') # fit VVV model to the 4 groups
m3.cv <- cvMclustDA(m3.da, nfold = 5, verbose = FALSE) # cross validate 
m3.cv$error

# G = 3
trainData <- cbind(clust.result$mean.occupancy*100, clust.result$max.occupancy*100)
trainClass <- clust.result$ag.class.both.3

models <- mclust.options()$emModelNames
tab <- matrix(NA, nrow = length(models), ncol = 4)
rownames(tab) <- models
colnames(tab) <- c("model","BIC", "CV_err", "G")
for(i in seq(models)) {
  tab[i,1] <- models[i]
  mod <- MclustDA(trainData, trainClass,
                  modelType = "EDDA", modelNames = models[i])
  tab[i,2] <- mod$bic
  tab[i,3] <- cvMclustDA(mod, nfold = 5, verbose = FALSE)$error
  tab[i,4] <- "3"
}
tab.3 <- as.data.frame(tab)

# G = 4
trainData <- cbind(clust.result$mean.occupancy*100, clust.result$max.occupancy*100)
trainClass <- clust.result$ag.class.both.4

models <- mclust.options()$emModelNames
tab <- matrix(NA, nrow = length(models), ncol = 4)
rownames(tab) <- models
colnames(tab) <- c("model","BIC", "CV_err", "G")
for(i in seq(models)) {
  tab[i,1] <- models[i]
  mod <- MclustDA(trainData, trainClass,
                  modelType = "EDDA", modelNames = models[i])
  tab[i,2] <- mod$bic
  tab[i,3] <- cvMclustDA(mod, nfold = 5, verbose = FALSE)$error
  tab[i,4] <- "4"
}
tab.4 <- as.data.frame(tab)

# create single table and filter to 5 largest BIC values
tab.cv <- rbind(tab.3, tab.4) %>%
  mutate(BIC = as.numeric(as.character(BIC)), CV_err = as.numeric(as.character(CV_err))) %>%
  top_n(n = 5, BIC) %>% arrange(desc(BIC))
tab.cv
  
##### Summary Outputs #####
# Summarize by class
clust.summary <- result %>%
  group_by(ag.class.both.4) %>%
  summarise(n = length(ag.class.mean),
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


