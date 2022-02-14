#' #### Stage Tagging Using Ensemble Classification ####

#' The aim is to define staging behavior prior to crop raiding from the movement track and use these staging 
#' events to investigate behavioral and landscape factors that facilitate the use of staging for crop raiding. 
#' 
#' We developed an algorithm to distinguish staging prior to agricultural use based on three parameters: 1) The 
#' start and end time of the stage; 2) the size of the staging window; and 3) the percentage of encamped 
#' relocations inside the time window. We hypothesized that staging events should have a higher percentage of 
#' encamped relocations than normal because 24-hour activity budgets suggested that elephants were increasing 
#' their encamped and meandering movements during the day when in crop raiding bouts. Additionally, we restricted 
#' possible start and end times between 6am and 6pm because 24-hour activity budgets suggested that elephants were 
#' increasing their encamped movements during this time. We also restricted the staging definition to exclude 
#' any directed-walk movements during the staging event.
#' 
#' We tested all possible combinations of the three parameters and calculated the accuracy of each parameter set.
#' Accuracy was calculated by comparing the total number of events the algorithm identified to the number 
#' of staging events that preceded agricultural use by elephants. 
#' 
#' The following code evaluates the results of the algorithm and uses weighted majority voting to produce an ensemble
#' classification of staging events.

##### Set up results dataframe #####
library(tidyverse)
library(lubridate)

setwd('~/Dropbox/CSU/GME_Movement')

# import algorithm results matrix
result.matrix <- readRDS("~/Dropbox/CSU/GME_Movement/Staging Areas/nsd loop tests/NSDloop_w2_6to21_seq200_20211223.RDS")

# convert matrix to dataframe for viz
result.df <- as.data.frame(result.matrix)
result.df <- mutate_at(result.df, c('ag.window.ext', 'n.stage', 'n.stage.adj', 'nsd.seq', 'n.stage.ratio',
                                    'n.stage.err', 'n.stage.acc', 'win.siz'), .funs = as.numeric)
result.df$win.start <- as.numeric(unlist(lapply(strsplit(as.character(result.df$win.period), "\\-"), "[", 1)))
result.df$win.end <- as.numeric(unlist(lapply(strsplit(as.character(result.df$win.period), "\\-"), "[", 2)))

result.df <- result.df[-1,]

# remove repetitive sequences
result.df <- result.df %>%
  #filter out repeated sequences
  group_by(n.stage.err) %>%
  filter(nsd.seq == max(nsd.seq))

#### Plot Loop Results ####

plot.df <- filter(result.df, ag.window.ext == 1)

ggplot(plot.df, aes(x = nsd.seq, y = (1-n.stage.err), group = as.factor(win.start))) +
  geom_line(aes(color = as.factor(win.start))) +
  facet_wrap(.~win.siz) +
  xlab('Mean NSD in window') + ylab('accuracy') + labs(color = 'window start time') +
  ggtitle('NSD Algorithm accuracy stratified by window size (hours)')

#' Based on the apparent trend of higher staging accuracy in the middle of the day,
#' plot the accuracy as a function of window start (left) and end (right) time. 

## Plot accuracy by start and end of window
par(mfrow = c(1,2))
boxplot((1-plot.df$n.stage.err) ~ plot.df$win.start, xlab = 'window start time', ylab = 'accuracy')
boxplot((1-plot.df$n.stage.err) ~ plot.df$win.end, xlab = 'window end time', ylab = 'accuracy')

#' #### Tag Stage Events ####
#' Use ensemble learning approach to tag staging events based on the full set of 
#' candidate models. First, import the movement dataset and specify a variable
#' for 24-hour ag windows. These 1/0 events are then indexed as unique events.

##### Import dataset #####

# import dataset
#gme <- as.data.frame(data.table::fread('~/Dropbox/CSU/GME_Movement/Staging Areas/GMEcollars_004_HMMclassified_stage_20201012.csv'))
gme <- readRDS('./Staging Areas/nsd loop tests/tor.3.w2.filter.RDS')

## assign unique id to each ag window event
# define unique id for window/non-window phases
gme <- gme[!is.na(gme$ag.window.ext),]
eventIndex <- ifelse(gme$ag.window.ext == 1, gme$dayBurst, 0)
# assign event index
gme$ag.window.index <- eventIndex

#' ##### Ensemble Voting #####

#' Ensemble approach uses weighted voting to evaluate each GPS relocation and tag 
#' it as stage or non-stage. Each model gets a vote which is weighted by it's 
#' accuracy. Relocations with a vote total greater than the mean is tagged as a 
#' stage.

##### Ensemble Voting #####

# prep plot.df and create empty lists
plot.df <- arrange(plot.df, n.stage.err) # order by descending error
plot.df$n.stage.acc <- 1-plot.df$n.stage.err # accuracy = 1-error

# restrict candidate model set for speed and to exclude irrelevant models with low accuracy.
# 150 cutoff is based on the ggplot of accuracy by nsd sequence. Appears that models stabalize around 150 net displacement. Removing these avoids artificially inflating model voting scores
plot.df <- filter(plot.df, nsd.seq <= 150 & n.stage.acc >= 0.5)
dim(plot.df) # unique permutations

# create new dataframe
gme.stage <- gme

# dump dataframes
remove(gme)
remove(result.matrix)
remove(result.df)

# split up dataframes for memory
n <- 50
nr <- nrow(plot.df)
split <- split(plot.df, rep(1:n, length.out = nr, each = ceiling(nr/n)))

#eventWeight <- NULL # Accuracy of the specific model, to be used for weighting
#eventWeight <- matrix(data = NA, nrow = dim(gme)[1], ncol = dim(plot.df)[1])
#eventWeight.all <- matrix(data = NA, nrow = dim(gme)[1], ncol = dim(plot.df)[1])

# set sums storage df 
sums <- NULL
sums.all <- NULL
start <- Sys.time()
for (i in 1:length(split)){
    df <- split[[i]]
    eventWeight <- NULL
    eventWeight.all <- NULL
  
  for (j in 1:dim(split[[i]])[1]) { # voting by row of split df -- row corresponds with the parameter set and accuracy for each model
    pct.threshold = df$nsd.seq[j]
    win.size = df$win.siz[j]
    gme.stage$stage.period <- ifelse(hour(gme.stage$date) >= df$win.start[j] & hour(gme.stage$date) <= df$win.end[j], 1, 0)
    
    df2 <- gme.stage %>%
      filter(!is.na(ag.window.ext)) %>%
      group_by(dayBurst, stage.period, ag.window.ext) %>%
      # mean nsd
      mutate(mean.nsd = mean(r2n))
    
    # tag all staging events
    stage.event <- ifelse(df2$mean.nsd <= pct.threshold & df2$stage.period == 1, 1, 0)
    # tag ag staging events
    ag.stage.event <- ifelse(df2$mean.nsd <= pct.threshold & df2$stage.period == 1 & df2$ag.window.ext == 1, 1, 0)
    
    ## produce list of vectors. each vec is length of dataset with 0/1 corresponding to stage tagging for one parameter set
    
    # all events - use this to get ensemble accuracy
    eventWeight.all[[j]] <- ifelse(stage.event == 1, 1*(df$n.stage.acc[i]), 0)
    
    # only ag events - use this for analysis of results
    eventWeight[[j]] <- ifelse(ag.stage.event == 1, 1*(df$n.stage.acc[j]), 0)
    
    }
    
    # calculate rowSums for each j loop and store as list
    t <- as.data.frame(do.call(cbind, eventWeight))
    sums[[i]] <- rowSums(t) # vector of rowsums stored as list element. 
    t.all <- as.data.frame(do.call(cbind, eventWeight.all))
    sums.all[[i]] <- rowSums(t.all)
    
}  

end <- Sys.time()

end - start

#saveRDS(sums.all, 'Staging Areas/nsd loop tests/nsd_ensemble_sums.all_20211229.RDS')
sums.all <- readRDS('Staging Areas/nsd loop tests/nsd_ensemble_sums.all_20211229.RDS')

#' The resulting dataframe has a column for each model (xx) and each row corresponds
#' to a GPS relocation. If a stage is detected for a given model, the vote value 
#' is tied to the model's accuracy for weighting (1*model.accuracy). 

## weighted majority vote test ##
# create dataframe
df.weight <- as.data.frame(do.call(cbind, sums.all))

# drop
remove(eventWeight)

dim(df.weight)

# tally weighted votes - sum of sums
df.weight$stagesum <- rowSums(df.weight, na.rm = T)
# get mean vote value for majority voting
w.mu <- sum(df.weight$stagesum)/sum(!!df.weight$stagesum)
# tag ensemble stage events
df.weight$majority <- ifelse(df.weight$stagesum > w.mu, 1,0)

# relocs in and outside stage
table(df.weight$majority)

## Assign stage tags to relocs in the main dataset
gme.stage$vote <- df.weight$majority

##### Get Ensemble Accuracy #####
## Not Run

## index staging events
eventFlag <- ifelse(gme.stage$vote == 1, TRUE, FALSE)
eventIndex <- inverse.rle(within.list(rle(eventFlag), 
                                      values[values] <- seq_along(values[values])))
# assign a unique event index for each stage event
gme.stage$stage.event.index <- eventIndex

table <- gme.stage %>% group_by(ag.window.ext) %>%
  summarise(n.stage = length(unique(stage.event.index)), 
            n.stage.adj = n.stage/length(unique(dayBurst)) )

table <- table %>%
  mutate(n.stage.ratio = n.stage.adj[2]/n.stage.adj[1]) %>%
  mutate(n.stage.err = n.stage.adj[1]/sum(n.stage.adj)) %>%
  mutate(n.stage.acc = 1-n.stage.err)
table

# check omission error rate
gme.stage %>% ungroup() %>% filter(ag.window.ext == 1) %>% summarise(n.agdays = length(unique(ag.window.index)), 
                                                                     n.stage = length(unique(stage.event.index)),
                                                                     omission = 1 - n.stage/n.agdays )

##### Adjust for Ag Spatial Threshold #####

##### Adjust Using Spatial Threshold Filter #####

# Add dist2ag 
locs <- sf::st_as_sf(gme.stage, coords = c('x','y'), crs = 32736) %>%
  terra::vect()
dist2ag <- terra::rast('./spatial data/dist2ag_estes_32736_2019-11-21.tif')
dist2forest <- terra::rast('./spatial data/dist2forest_hansen_cover60_32736_30.tif')

gme.stage$dist2ag <- terra::extract(dist2ag, locs)[,2]
gme.stage$dist2forest <- terra::extract(dist2forest, locs)[,2]

# Define spatial threshold
#gme.stage$vote.thresh <- ifelse(gme.stage$dist2ag <= 3500, gme.stage$vote, 0)
gme.stage$vote.thresh <- ifelse(gme.stage$dist2forest <= 1200, gme.stage$vote, 0)

## index staging events
gme.stage$stage.event.index <- ifelse(gme.stage$vote.thresh == 1, gme.stage$dayBurst, 0)

table <- gme.stage %>% group_by(ag.window.ext) %>%
  summarise(n.stage = length(unique(stage.event.index)), 
            n.stage.adj = n.stage/length(unique(dayBurst)) )

table <- table %>%
  mutate(n.stage.ratio = n.stage.adj[2]/n.stage.adj[1]) %>%
  mutate(n.stage.err = n.stage.adj[1]/sum(n.stage.adj)) %>%
  mutate(n.stage.acc = 1-n.stage.err)
table

# spatial threshold - check omission error rate
gme.stage %>% ungroup() %>% filter(ag.window.ext == 1) %>% summarise(n.agdays = length(unique(ag.window.index)), 
                                                                     n.stage = length(unique(stage.event.index)),
                                                                     omission = 1-n.stage/n.agdays )


##### Create stage-tagged dataset #####

## apply ag filter to identify true staging
gme.stage$vote.ag <- ifelse(gme.stage$ag.window.ext == 1, gme.stage$vote, 0)
gme.stage$vote.nonag <- ifelse(gme.stage$ag.window.ext == 0, gme.stage$vote, 0)

## index ag staging events
eventFlag <- ifelse(gme.stage$vote.ag == 1, TRUE, FALSE)
eventIndex <- inverse.rle(within.list(rle(eventFlag), 
                                      values[values] <- seq_along(values[values])))
# assign a unique event index for each stage event
gme.stage$stage.event.index <- eventIndex

#' ##### Staging Area Data Summaries

#' Check the hourly distribution of staging relocations

# check distribution of stage relocs over the course of the day
t <- filter(gme.stage, vote == 1) %>%
  group_by(hour(date)) %>% tally()
ggplot(t, aes(as.factor(`hour(date)`), n)) + geom_bar(stat = 'identity') +
  xlab('hour of day') + ylab('n stage relocations')

#' The total number of stage events in the dataset is 5692.

## index staging events
eventFlag <- ifelse(gme.stage$vote == 1, TRUE, FALSE)
eventIndex <- inverse.rle(within.list(rle(eventFlag),
                                           values[values] <- seq_along(values[values])))
gme.stage$stage.event.index <- eventIndex

length(unique(gme.stage$stage.event.index))

#' Calculate distribution of staging frequency prior to crop raiding among individuals. 

# summarise overall % staging distribution
stage.summary <- gme.stage %>%
  #group_by(subject_name) %>%
  summarise(n.stage = length(unique(stage.event.index)),
            n.raid = length(unique(ag.window.index)),
            pct.stage = n.stage/n.raid) 

summary(stage.summary$pct.stage)

#' Summarise staging frequency distribution by individual and tactic. Use yearly tactics (n=101) 
#' Note the Rare group with 100% stage percentages have very small numbers of 
#' stages. See the end of the script for data summary.

# summarise % staging distribution by individual-year tactic
stage.summary <- gme.stage %>%
  group_by(subject_name, tactic.season, year.cuts) %>%
  summarise(n.stage = length(unique(stage.event.index)),
            n.raid = length(unique(ag.window.index)),
            pct.stage = n.stage/n.raid)

boxplot(stage.summary$pct.stage ~ stage.summary$tactic.season,
        main = 'distribution of staging frequency by individuals in the 4 tactics',
        xlab = 'Tactics (Rare to Habitual)', ylab = 'Pct. of Raids with a Stage')


#' ##### Plot the data on a map

##### Plot Staging Events #####
library(sf)
gme <- sf::st_read('./spatial data/GSE/GSE_2020.shp') 
stage.relocs <- filter(gme.stage, vote.ag == 1)

ggplot(data = gme) + geom_sf() + coord_sf(datum=st_crs(32736)) + 
  geom_point(data = stage.relocs, aes(x, y, color = 'red'), size = 0.2, alpha = 0.2) +
  labs(color = 'staging relocations')

#' Table of staging events and ag bouts for each individual-year
stage.summary <- gme.stage %>%
  group_by(subject_name, tactic.season, year.cuts) %>%
  summarise(n.stage = length(unique(stage.event.index)),
            n.raid = length(unique(ag.window.index)),
            pct.stage = n.stage/n.raid)
(stage.summary)


#### Save Staging-Tagged Data ####
write.csv(gme.stage,'./Staging Areas/movdata/GMEcollars_004_stageclassified_NSD_20211019.csv')


