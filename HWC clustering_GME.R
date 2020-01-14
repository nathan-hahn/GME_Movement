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
movdat <- readRDS('./movdata/GMEcollars_001_used_2019-12-18.rds')
movdat <- movdat[-which(movdat$name%in%c("Tunai","Bettye","Nancy", "Wilbur", 
                             "Julia", "Ndorre", "Nkoidila", "Santiyan", "Earhart",
                             "Rudisha", "Naeku", "Olkeri", "ST2010-1441")),]
movdat <- movdat[-which(movdat$name%in%c("Nyanza")),]

# dataframe for fitting
df <- movdat

####Mean occupancy####
fixes <- df %>%
  group_by(name) %>%
  tally()

ag.mean <- df %>%
  group_by(name) %>%
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
  group_by(site, name, month) %>%
  tally() 

ag.month <- df %>%
  group_by(site, name, month) %>%
  tally(lc.estes == 1) %>% # tally fixes in ag
  droplevels() %>%
  ungroup() %>%
  rename(n.ag = n) %>%
  mutate(month = as.factor(month)) %>%
  mutate(n = fixes.m$n) %>%
  arrange(name, month) %>%
  mutate(month.occupancy = n.ag/n)

(ag.month) # the month for each individual where ag occupancy is highest (across all)

# bar plot of ag occupancy by month
p <- ggplot(ag.month, aes(month, month.occupancy)) + geom_bar(stat = "identity") + facet_wrap(.~name) + 
  scale_x_discrete(name ="month", breaks=c(3,6,9))
p

# filter to each individuals max monthly ag occupancy 
result <- ag.month %>%
  #filter(site == "mep") %>% ## filter by site if needed! 
  group_by(name) %>%
  slice(which.max(month.occupancy)) %>%
  ungroup() %>%
  mutate(
    max.occupancy = round(month.occupancy, 5), 
    mean.occupancy = round(ag.mean$mean.ag/100, 5), 
    month = as.numeric(month))

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
both.df <- cbind(result$mean.occupancy, result$max.occupancy)
m3 <- Mclust(both.df, G = 4)
plot(m3, what = "classification",
     xlab = "mean ag occupancy (*100)",
     ylab = "max monthly ag occupancy *100")

####Results####

# select top based on BIC values and biological realism (m1 and m2)
plot(m1, what = "BIC")
plot(m2, what = "BIC")
plot(m3, what = "BIC")

# add classification to results
ag.class.mean <- as.factor(m1$classification)
ag.class.both <- as.factor(m3$classification)
ag.class.both <- factor(ag.class.both, levels(ag.class.both)[c(1,2,4,3)], labels = c(1,2,3,4)) # FOR GME ONLY
#ag.class.both <- factor(ag.class.both,levels(ag.class.both)[c(1,4,2,3)], labels = c(1,2,3,4)) #FOR MEP ONLY

result$ag.class.mean <- ag.class.mean
result$ag.class.both <- ag.class.both

# Summarize by individual
clust.result <- result %>%
  dplyr::select(-c(n.ag, n)) %>%
  dplyr::select(-month, -max.occupancy) 
(clust.result)

# Summarize by class
clust.summary <- result %>%
  group_by(ag.class.both) %>%
  summarise(n = length(ag.class.mean),
            mean = mean(mean.occupancy),
            sd = sd(mean.occupancy),
            se = sd/sqrt(n),
            lower = mean - se,
            upper = mean + se)
(clust.summary)

clust.df <- result %>%
  dplyr::select(name, ag.class.mean, ag.class.both) %>%
  inner_join(.,df, by = "name") %>%
  dplyr::select(-month) %>%
  as.data.frame()
clust.df$merge_id <- NULL


ggplot(clust.summary) + geom_bar(aes(x = ag.class.both, y = mean), stat = "identity") + 
  geom_errorbar(aes(x = ag.class.both, ymin = lower, ymax = upper), width = .2) + ylab("mean ag occupancy*100") + xlab("ag usage cluster")


####Export####
# new used df with cluster classification. 
# has ag.class.mean and ag.class.both classifications included
{outfile <- paste0("GMEcollars_001_usedClust_", Sys.Date(), ".rds")
saveRDS(clust.df, outfile)}

{outfile <- paste0("GMEcollars_001_usedClust_", Sys.Date(), ".csv")
write.csv(clust.df, outfile)}

# cluster results
{outfile <- paste0("./GMM/clustResult_GME_", Sys.Date(), ".csv")
  write.csv(clust.result, outfile)}

# cluster summary (mean or both)
{outfile <- paste0("GMEcollars_001_clustSummary_", Sys.Date(), ".csv")
  write.csv(clust.result, outfile)}

{outfile <- paste0("GMEcollars_001_clustSummary_", Sys.Date(), ".csv")
  write.csv(clust.summary, outfile)}







