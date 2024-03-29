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

setwd('~/Dropbox/CSU/GME_Movement') # for ssh

# import algorithm results matrix
result.matrix <- readRDS('./Staging Areas/persistence velocity loop tests/pvloop_6to21_seq0_agmask_20220606.RDS')

# convert matrix to dataframe for viz
result.df <- as.data.frame(result.matrix)
result.df <- mutate_at(result.df, c('ag.window.ext', 'n.stage', 'n.stage.adj', 'pv.seq', 'n.stage.ratio',
                                    'n.stage.err', 'win.siz'), .funs = as.numeric)
result.df$win.start <- as.numeric(unlist(lapply(strsplit(as.character(result.df$win.period), "\\-"), "[", 1)))
result.df$win.end <- as.numeric(unlist(lapply(strsplit(as.character(result.df$win.period), "\\-"), "[", 2)))

result.df <- result.df[-1,]

# remove repetitive sequences
result.df <- result.df %>%
  #filter out repeated sequences
  group_by(n.stage.err) %>%
  filter(pv.seq == max(pv.seq))

#### Plot Loop Results ####

plot.df <- filter(result.df, ag.window.ext == 1) # filter out repeats for non-ag windows
dim(plot.df) # number of permutations 

ggplot(plot.df, aes(x = pv.seq, y = (1-n.stage.err), group = as.factor(win.start))) +
  geom_line(aes(color = as.factor(win.start))) +
  #geom_point(aes(shape = as.factor(win.end))) +
  facet_wrap(.~win.siz) +
  xlab('mean pv in window') + ylab('accuracy') + labs(color = 'window start time') +
  ggtitle('Algorithm accuracy stratified by window size (hours)')

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

# ag mask update
gme <- as.data.frame(data.table::fread('~/Dropbox/CSU/GME_Movement/Staging Areas/pct loop tests/GMEcollars_004_HMMclassified_stage_20220103.csv'))

## unique id for each individual-day
gme$dayBurst <- paste(gme$id, as.Date(gme$date))

## 24-hour ag window identifies days where they used ag - dayburst helps group by day and individual 
gme <- gme %>%
  group_by(dayBurst) %>%
  mutate(ag.window.ext = ifelse(sum(ag.window.mask) >= 1, 1, 0)) # if ag.window identifies ag use phase, classify that day as ag use. Helps to pick up stages at the beginning of an ag use phase that might otherwise be split up.

## assign unique id to each ag window event
# define unique id for window/non-window phases
gme <- gme[!is.na(gme$ag.window.ext),]
eventIndex <- ifelse(gme$ag.window.ext == 1, gme$dayBurst, 0)
# assign event index
gme$ag.window.index <- eventIndex

## add persistence velocity ##
gme <- gme %>% group_by(id) %>%
  mutate(velocity = dist/(dt/60/60)) %>% # meters per hour
  mutate(pv = velocity * cos(rel.angle)) %>% # turning angle (not heading)
  filter(!is.na(pv))

tt <- gme[gme$ag.window.ext == 1,]
ggplot(data = gme.shp) + geom_sf() + coord_sf(datum=st_crs(32736)) +
  geom_point(data = tt, aes(x, y, color = 'red'), size = 0.2, alpha = 0.2) +
  labs(color = 'staging relocations')

#' ##### Ensemble Voting #####

#' Ensemble approach uses weighted voting to evaluate each GPS relocation and tag 
#' it as stage or non-stage. Each model gets a vote which is weighted by it's 
#' accuracy. Relocations with a vote total greater than the mean is tagged as a 
#' stage.

##### Ensemble Voting #####

# prep plot.df and create empty lists
plot.df <- arrange(plot.df, n.stage.err) # order by descending error
plot.df$n.stage.acc <- 1-plot.df$n.stage.err # accuracy = 1-error

# remove crappy algorithms (<50% accuracy)
plot.df <- filter(plot.df, n.stage.acc >= 0.5)

# split up dataframes for memory
n <- 50
nr <- nrow(plot.df)
split <- split(plot.df, rep(1:n, length.out = nr, each = ceiling(nr/n)))

# set sums storage df 
sums <- NULL
sums.all <- NULL
start <- Sys.time()
for (i in 1:length(split)){
  df <- split[[i]]
  eventWeight <- NULL
  eventWeight.all <- NULL
  
  for (j in 1:dim(split[[i]])[1]) { # voting by row of split df -- row corresponds with the parameter set and accuracy for each model
    pv.threshold = df$pv.seq[j]
    win.size = df$win.siz[j]
    gme$stage.period <- ifelse(hour(gme$date) >= df$win.start[j] & hour(gme$date) <= df$win.end[j], 1, 0)
    
    df2 <- gme %>%
      group_by(dayBurst, stage.period, ag.window.ext) %>%
      mutate(pv.mean = mean(pv, na.rm = T))
    
    # tag all staging events
    stage.event <- ifelse(df2$pv.mean <= pv.threshold & df2$stage.period == 1, 1, 0)
    stage.event <- ifelse(df2$stage.period == 1, stage.event, NA) # NA if a reloc is outside of the time window
    # tag ag staging events
    ag.stage.event <- ifelse(df2$pv <= pv.threshold & df2$stage.period == 1 & df2$ag.window.ext == 1, 1, 0)
    ag.stage.event <- ifelse(df2$stage.period == 1, ag.stage.event, NA) # NA if a reloc is outside of the time window
    
    ## produce list of vectors. each vec is length of dataset with 0/1 corresponding to stage tagging for one parameter set
    # all events - use this to get ensemble accuracy
    eventWeight.all[[j]] <- ifelse(stage.event == 1, 1*(plot.df$n.stage.acc[i]), stage.event)
    # only ag events - use this for analysis of results
    eventWeight[[j]] <- ifelse(ag.stage.event == 1, 1*(plot.df$n.stage.acc[i]), ag.stage.event)
   
  }
  
  # calculate rowSums for each j loop and store as list
  t <- as.data.frame(do.call(cbind, eventWeight))
  sums[[i]] <- rowSums(t, na.rm = T) # vector of rowsums stored as list element. 
  t.all <- as.data.frame(do.call(cbind, eventWeight.all))
  sums.all[[i]] <- rowSums(t.all, na.rm = T)
  
}  

end <- Sys.time()

end - start

saveRDS(sums, 'Staging Areas/persistence velocity loop tests/pv_ensemble_sums_mask_20220614.RDS')
saveRDS(sums.all, 'Staging Areas/persistence velocity loop tests/pv_ensemble_sums.all_mask_20220614.RDS')

# # for updated ag mask - removal of algorithms <50% accurate# sums.all <- readRDS('Staging Areas/step loop tests/step_ensemble_sums.all_mask_20220602.RDS')# #' The resulting dataframe has a column for each model (xx) and each row corresponds
#' to a GPS relocation. If a stage is detected for a given model, the vote value
#' is tied to the model's accuracy for weighting (1*model.accuracy).

sums.ag <- readRDS('Staging Areas/persistence velocity loop tests/pv_ensemble_sums_mask_20220614.RDS')
sums.all <- readRDS('Staging Areas/persistence velocity loop tests/pv_ensemble_sums.all_mask_20220614.RDS')


## weighted majority vote test ##
# create dataframe
#sums.all <- sums
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
gme$vote <- df.weight$majority

##### Get Ensemble Accuracy #####

## index staging events
gme$stage.event.index <- ifelse(gme$vote == 1, gme$dayBurst, 0)

# check stage summary
table <- gme %>% group_by(ag.window.ext) %>%
  summarise(n.stage = length(unique(stage.event.index)),
            n.stage.adj = n.stage/length(unique(dayBurst)) )

table <- table %>%
  mutate(n.stage.ratio = n.stage.adj[2]/n.stage.adj[1]) %>%
  mutate(n.stage.err = n.stage.adj[1]/sum(n.stage.adj)) %>%
  mutate(n.stage.acc = 1-n.stage.err)
table

# check omission error rate
gme %>% ungroup() %>% filter(ag.window.ext == 1) %>% summarise(n.agdays = length(unique(ag.window.index)),
                                                               n.stage = length(unique(stage.event.index)),
                                                               omission = 1 - n.stage/n.agdays )


##### Adjust Using Spatial Threshold Filter #####
gme$vote.thresh <- ifelse(gme$dist2ag <= 3500, gme$vote, 0)
gme$vote.thresh <- ifelse(gme$dist2forest <= 1200, gme$vote, 0)

## index staging events
gme$stage.event.index <- ifelse(gme$vote.thresh == 1, gme$dayBurst, 0)


table <- gme %>% group_by(ag.window.ext) %>%
  summarise(n.stage = length(unique(stage.event.index)),
            n.stage.adj = n.stage/length(unique(dayBurst)) )

table <- table %>%
  mutate(n.stage.ratio = n.stage.adj[2]/n.stage.adj[1]) %>%
  mutate(n.stage.err = n.stage.adj[1]/sum(n.stage.adj)) %>%
  mutate(n.stage.acc = 1-n.stage.err)
table

# spatial threshold - check omission error rate
gme %>% ungroup() %>% filter(ag.window.ext == 1) %>% summarise(n.agdays = length(unique(ag.window.index)),
                                                               n.stage = length(unique(stage.event.index)),
                                                               omission = 1 - n.stage/n.agdays )


##### Create stage-tagged dataset #####

# remove non-ag staging events
gme.stage <- gme
gme.stage$vote.thresh <- NULL #remove old column
gme.stage$vote.ag <- ifelse(gme$ag.window.ext == 1, gme$vote == 1, 0)
# index ag staging events
gme.stage$stage.event.index <- ifelse(gme.stage$vote.ag == 1, gme.stage$dayBurst, 0)
length(unique(gme.stage$stage.event.index))

# # results of ag-only ensemble model (does not include non-ag stages)
# sums.ag <- readRDS('Staging Areas/persistence velocity loop tests/pv_ensemble_sums_mask_20220614.RDS')
# df.weight <- as.data.frame(do.call(cbind, sums.ag))
# # tally weighted votes - sum of sums
# df.weight$stagesum <- rowSums(df.weight, na.rm = T)
# # get mean vote value for majority voting
# w.mu <- sum(df.weight$stagesum)/sum(!!df.weight$stagesum)
# # tag ensemble stage events
# df.weight$majority <- ifelse(df.weight$stagesum > w.mu, 1,0)
# 
# ## Assign stage tags to relocs in the main dataset
# gme.stage <- as.data.frame(cbind(gme, 'vote.ag' = df.weight$majority))
# sum(gme.stage$vote.ag)
# 
# ## index ag staging events with a unique number
# eventFlag <- ifelse(gme.stage$vote.ag == 1, TRUE, FALSE)
# eventIndex <- inverse.rle(within.list(rle(eventFlag),
#                                       values[values] <- seq_along(values[values])))
# # assign a unique event index for each stage event
# gme.stage$stage.event.index <- eventIndex


#' ##### Staging Area Data Summaries

#' Check the hourly distribution of staging relocations

# check distribution of stage relocs over the course of the day
t <- filter(gme.stage, vote.ag == 1) %>%
  group_by(hour(date)) %>% tally()
ggplot(t, aes(as.factor(`hour(date)`), n)) + geom_bar(stat = 'identity') +
  xlab('hour of day') + ylab('n stage relocations')

#' The total number of stage events in the dataset is 6969.

gme.stage %>% group_by(vote.ag) %>%
  summarise(n = length(unique(stage.event.index)))

#' Calculate distribution of staging frequency prior to crop raiding among individuals.

# summarise overall % staging distribution
stage.summary <- gme.stage %>%
  group_by(subject_name) %>%
  summarise(n.stage = length(unique(stage.event.index)),
            n.raid = length(unique(ag.window.index)),
            pv.stage = n.stage/n.raid)

summary(stage.summary$pv.stage)

#' Summarise staging frequency distribution by individual and tactic. Use yearly tactics (n=101)
#' Note the Rare group with 100% stage percentages have very small numbers of
#' stages. See the end of the script for data summary.

# summarise % staging distribution by individual-year tactic
stage.summary <- gme.stage %>%
  group_by(subject_name, tactic.season, year.cuts) %>%
  summarise(n.stage = length(unique(stage.event.index)),
            n.raid = length(unique(ag.window.index)),
            step.stage = n.stage/n.raid)

boxplot(stage.summary$step.stage ~ stage.summary$tactic.season,
        main = 'distribution of staging frequency by individuals in the 4 tactics',
        xlab = 'Tactics (Rare to Habitual)', ylab = 'Pct. of Raids with a Stage')


#' ##### Plot the data on a map

##### Plot Staging Events #####
library(sf)
gme.shp <- sf::st_read('./spatial data/GSE/GSE_2020.shp')
stage.relocs <- filter(gme.stage, vote.ag == 1) #filter(gme.stage, ag.window.mask == 1)

ggplot(data = gme.shp) + geom_sf() + coord_sf(datum=st_crs(32736)) +
  geom_point(data = stage.relocs, aes(x, y, color = 'red'), size = 0.2, alpha = 0.2) +
  labs(color = 'staging relocations')

#' Table of staging events and ag bouts for each individual-year
stage.summary <- gme.stage %>%
  group_by(subject_name, tactic.season, year.cuts) %>%
  summarise(n.stage = length(unique(stage.event.index)),
            n.raid = length(unique(ag.window.index)),
            step.stage = n.stage/n.raid)
(stage.summary)


#### Save Staging-Tagged Data ####
write.csv(gme.stage,'./Staging Areas/movdata/GMEcollars_004_stageclassified_pv_20220614.csv')


