#### Ensemble Learning ####


library(tidyverse)
library(lubridate)

# import results matrix
result.matrix <- readRDS('stage_loop_result_pct T2.1 GME_004.RDS')

# convert matrix to dataframe for viz
result.df <- as.data.frame(result.matrix)
result.df <- mutate_at(result.df, c('ag.window.ext', 'n.stage', 'n.stage.adj', 'pct.seq', 'n.stage.ratio',
                                    'n.stage.err', 'win.siz'), .funs = as.numeric)
result.df <- result.df[-1,]

# remove repetitive sequences
result.df <- result.df %>%
  group_by(ag.window.ext) %>%
  filter(n.stage > 2) %>%
  mutate(n.stage.ratio = ifelse(n.stage.ratio == lag(n.stage.ratio), NA, n.stage.ratio)) %>%
  mutate(n.stage.err = ifelse(n.stage.err == lag(n.stage.err), NA, n.stage.err))

#### Plot Loop Results ####

plot.df <- filter(result.df, !is.na(n.stage.err) & ag.window.ext == 1)
plot.df$win.start <- as.numeric(unlist(lapply(strsplit(as.character(plot.df$win.period), "\\-"), "[", 1)))
plot.df$win.end <- as.numeric(unlist(lapply(strsplit(as.character(plot.df$win.period), "\\-"), "[", 2)))

ggplot(plot.df, aes(x = pct.seq, y = (1-n.stage.err), group = as.factor(win.start))) +
  geom_line(aes(color = as.factor(win.start))) +
  #geom_point(aes(color = win.period), position = 'jitter') +
  facet_wrap(.~win.siz)

## Plot accuracy by start and end of window
par(mfrow = c(1,2))
plot(as.numeric(plot.df$win.start), (1-plot.df$n.stage.err))
plot(as.numeric(plot.df$win.end), (1-plot.df$n.stage.err))

boxplot((1-plot.df$n.stage.err) ~ plot.df$win.start)
boxplot((1-plot.df$n.stage.err) ~ plot.df$win.end)

#### Tag Staging Events for Each Parameter Set ####

# import dataset
set.seed(1992)
gme <- as.data.frame(data.table::fread("HMM/movdata/GMEcollars_004_HMMclassified_20211012.csv"))
colnames(gme)[1] <- 'uid'
unique(gme$subject_name)

# filter out outlier steps 
split <- split(gme, gme$viterbi)
split[[1]] <- filter(split[[1]], dist <= median((split[[2]]$dist + runif(1, -200, 200)), na.rm = T))
split[[2]] <- filter(split[[2]], dist <= median((split[[3]]$dist + runif(1, -200, 200)), na.rm = T))
split[[3]] <- filter(split[[3]], dist >= median((split[[2]]$dist + runif(1, -200, 200)), na.rm = T))
gme <- do.call(rbind, split)
gme <- gme[with(gme, order(uid)), ]

## unique id for each individual-day
gme$dayBurst <- paste(gme$id, as.Date(gme$date))

## 24-hour ag window identifies days where they used ag - dayburst helps group by day and individual 
gme <- gme %>%
  group_by(dayBurst) %>%
  mutate(ag.window.ext = ifelse(sum(ag.window) >= 1, 1, 0)) # if ag.window identifies ag use phase, classify that day as ag use. Helps to pick up stages at the beginning of an ag use phase that might otherwise be split up.

## assign unique id to each ag window event
# define unique id for window/non-window phases
gme <- gme[!is.na(gme$ag.window),]
eventFlag <- ifelse(gme$ag.window == 1, TRUE, FALSE)
eventIndex <- inverse.rle(within.list(rle(eventFlag), 
                                      values[values] <- seq_along(values[values])))
# assign event index
gme$ag.window.index <- eventIndex

##### Tag staging events with optimal sequence #####
plot.df <- arrange(plot.df, n.stage.err)
plot.df$n.stage.acc <- 1-plot.df$n.stage.err
eventFlag <- NULL
eventIndex <- NULL
eventWeight <- NULL

gme.stage <- gme %>%
  filter(!is.na(ag.window.ext)) 

for (i in 1:3) {
  pct.threshold = plot.df$pct.seq[i]
  win.size = plot.df$win.siz[i]
  gme.stage$stage.period <- ifelse(hour(gme.stage$date) >= plot.df$win.start[i] & hour(gme.stage$date) <= plot.df$win.end[i], 1, 0)
  
  t <- gme.stage %>%
    group_by(dayBurst, stage.period, ag.window.ext) %>%
    mutate(enc.day = sum(viterbi==1), meander.day = sum(viterbi==2), dw.day = sum(viterbi==3), n.day = n()) %>%
    mutate(pct = round((enc.day/win.size), 3))
  
  # tag all staging events
  stage.event <- ifelse(t$pct > pct.threshold & t$stage.period == 1 & t$dw.day == 0, 1, 0)
  # tag ag staging events
  ag.stage.event <- ifelse(t$pct > pct.threshold & t$stage.period == 1 & t$dw.day == 0 & t$ag.window.ext == 1, 1, 0)
  
  ## index staging events
  eventFlag[[i]] <- ifelse(ag.stage.event == 1, TRUE, FALSE)
  eventIndex[[i]] <- inverse.rle(within.list(rle(eventFlag[[i]]),
                                             values[values] <- seq_along(values[values])))
  eventWeight[[i]] <- ifelse(ag.stage.event == 1, 1*(plot.df$n.stage.acc[i]), 0)
}

df.weight <- as.data.frame(do.call(cbind, eventWeight))
df.weight$stagesum <- rowSums(df.weight)
