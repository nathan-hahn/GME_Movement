##### GMM Predict #####

## Get top model and movement data
mSelect <- readRDS("./GMM/results/mSelect.rds")
#movdat.filter <- readRDS('./movdata/GMEcollars_002_usedFilter_2020-06-20.rds')

##### Create data subsets #####
# dataframe for fitting
df <- temp.filter %>%
  mutate(ag.used = if_else(lc.estes == 1, 1, 0),
         month = month(date)) 
# Mean Occ.
ag.mean <- df %>%
  group_by(subject_name) %>%
  summarise(n = n(),
            mean.ag = mean(ag.used)) 

(ag.mean)

# Month Occ
# calculate occupancy by month
ag.month <- df %>%
  group_by(site, subject_name, month) %>%
  summarise(n.month = n(),
            month.occupancy = mean(ag.used)) %>%
  filter(n.month >= 5*24)

##### Rolling Occupancy #####
# filter for rolling mean calculations - at least 4 months of data
roll.filter <- df %>%
  group_by(subject_name) %>%
  tally() %>% 
  filter(n >= 120*24) 

roll.df <- filter(df, subject_name %in% roll.filter$subject_name)
split <- split(roll.df, roll.df$subject_name)

# window sizes
window <- c(90*24) # in hours. Split to before and after when align = center

roll.output <- NULL
for(i in 1:length(window)){
  # calculate moving window stats - see rollstats function
  roll.output[[i]] <- window_stats(df.list = split, window = window[[i]])
  # remove NAs from ag windows - created by window cutoff
  roll.output[[i]] <- filter(roll.output[[i]], is.na(ag.window) == FALSE) %>% droplevels()
}

# extract max window output. Test for window size
roll.max <- roll.output[[1]] %>%
  group_by(subject_name) %>%
  summarise(roll.max = max(mean)) 


##### Tabulate Occupancy Results #####
# Join the results from each occupancy method into a single table of individuals for model fitting
t <- full_join(max, ag.mean) %>% full_join(., roll.max)
(t)

# TEST

# For 1D plot
test.df <- cbind(t$mean.ag)
test1D <- test.df[complete.cases(test.df), ]
pred1D <- predict.Mclust(m1, test1D)

# For 2D plot
test.df <- cbind(t$mean.ag, t$roll.max)
test2D <- test.df[complete.cases(test.df), ]
pred2D <- predict.Mclust(mSelect, test2D)


hist(pred1D$classification, breaks = 5)
mclust1Dplot(test1D, z = pred1D$z, classification = pred1D$classification)

hist(pred2D$classification, breaks = 5)
mclust2Dplot(test2D, z = pred2D$z, classification = pred2D$classification)


mclust2Dplot(test2D, z = pred2D$z, classification = pred2D$classification, what = 'uncertainty')

