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

setwd('~/Dropbox (Personal)/CSU/GME_Movement')

# import algorithm results matrix
result.matrix <- readRDS('stage_loop_result_pct T2.1 seq0 GME_004.RDS')

# convert matrix to dataframe for viz
result.df <- as.data.frame(result.matrix)
result.df <- mutate_at(result.df, c('ag.window.ext', 'n.stage', 'n.stage.adj', 'pct.seq', 'n.stage.ratio',
                                    'n.stage.err', 'win.siz'), .funs = as.numeric)
result.df$win.start <- as.numeric(unlist(lapply(strsplit(as.character(result.df$win.period), "\\-"), "[", 1)))
result.df$win.end <- as.numeric(unlist(lapply(strsplit(as.character(result.df$win.period), "\\-"), "[", 2)))

result.df <- result.df[-1,]

# remove repetitive sequences
result.df <- result.df %>%
  #filter out repeated sequences
  group_by(n.stage.err) %>%
  filter(pct.seq == max(pct.seq))

#### Plot Loop Results ####

plot.df <- filter(result.df, ag.window.ext == 1)

ggplot(plot.df, aes(x = pct.seq, y = (1-n.stage.err), group = as.factor(win.start))) +
  geom_line(aes(color = as.factor(win.start))) +
  #geom_point(aes(shape = as.factor(win.end))) +
  facet_wrap(.~win.siz) +
  xlab('pct encamped in window') + ylab('accuracy') + labs(color = 'window start time') +
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
gme <- as.data.frame(data.table::fread('~/Dropbox (Personal)/CSU/GME_Movement/Staging Areas/GMEcollars_004_HMMclassified_stage_20201012.csv'))

## unique id for each individual-day
gme$dayBurst <- paste(gme$id, as.Date(gme$date))

## 24-hour ag window identifies days where they used ag - dayburst helps group by day and individual 
gme <- gme %>%
  group_by(dayBurst) %>%
  mutate(ag.window.ext = if_else(sum(ag.window) >= 1, 1, 0)) # if ag.window identifies ag use phase, classify that day as ag use. Helps to pick up stages at the beginning of an ag use phase that might otherwise be split up.

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
eventFlag <- NULL # TRUE/FALSE whether a relocation is part of a stage
eventIndex <- NULL # Unique ID for stage event
eventWeight <- NULL # Accuracy of the specific model, to be used for weighting
eventWeight.all <- NULL

# create new dataframe
gme.stage <- gme

for (i in 1:dim(plot.df)[1]) { # voting by row of plot.df -- row corresponds with the parameter set and accuracy for each model
  pct.threshold = plot.df$pct.seq[i]
  win.size = plot.df$win.siz[i]
  gme.stage$stage.period <- ifelse(hour(gme.stage$date) >= plot.df$win.start[i] & hour(gme.stage$date) <= plot.df$win.end[i], 1, 0)
  
  t <- gme.stage %>%
    group_by(dayBurst, stage.period, ag.window.ext) %>%
    mutate(enc.day = sum(viterbi==1), meander.day = sum(viterbi==2), dw.day = sum(viterbi==3), n.day = n()) %>%
    mutate(pct = round((enc.day/win.size), 3))
  
  # tag all staging events
  stage.event <- ifelse(t$pct >= pct.threshold & t$stage.period == 1 & t$dw.day == 0, 1, 0)
  # tag ag staging events
  ag.stage.event <- ifelse(t$pct >= pct.threshold & t$stage.period == 1 & t$dw.day == 0 & t$ag.window.ext == 1, 1, 0)
  
  ## index staging events
  # eventFlag[[i]] <- ifelse(ag.stage.event == 1, TRUE, FALSE)
  # eventIndex[[i]] <- inverse.rle(within.list(rle(eventFlag[[i]]),
  #                                            values[values] <- seq_along(values[values])))
  eventWeight[[i]] <- ifelse(ag.stage.event == 1, 1*(plot.df$n.stage.acc[i]), 0)
  }

#' The resulting dataframe has a column for each model (418) and each row corresponds
#' to a GPS relocation. If a stage is detected for a given model, the vote value 
#' is tied to the model's accuracy for weighting (1*model.accuracy). 

## weighted majority vote test ##
# create dataframe
df.weight <- as.data.frame(do.call(cbind, eventWeight))
dim(df.weight)

# tally weighted votes
df.weight$stagesum <- rowSums(df.weight)
# get mean vote value for majority voting
w.mu <- sum(df.weight$stagesum)/sum(!!df.weight$stagesum)
# tag ensemble stage events
df.weight$majority <- ifelse(df.weight$stagesum > w.mu, 1,0)

# relocs in and outside stage
table(df.weight$majority)

## Assign stage tags to relocs in the main dataset
gme.stage$vote <- df.weight$majority

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
  group_by(subject_name) %>%
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
gme <- sf::st_read('~/Dropbox (Personal)/CSU/GME_Movement/spatial data/GSE/GSE_2020.shp') 
stage.relocs <- filter(gme.stage, vote == 1)

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
write.csv(gme.stage,'./Staging Areas/movdata/GMEcollars_004_stageclassified_20211019.csv')


