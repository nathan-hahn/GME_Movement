#### Stage Test - E-M ratio ####

library(tidyverse)
library(lubridate)


# import state-classified data
gme <- as.data.frame(data.table::fread("HMM/movdata/GMEcollars_003_HMMclassified_20201206.csv"))
colnames(gme)[1] <- 'uid'
unique(gme$subject_name)

# filter out outlier steps 
split <- split(gme, gme$viterbi)
split[[1]] <- filter(split[[1]], dist <= median((split[[2]]$dist + runif(1, -200, 200)), na.rm = T))
split[[2]] <- filter(split[[2]], dist <= median((split[[3]]$dist + runif(1, -200, 200)), na.rm = T))
split[[3]] <- filter(split[[3]], dist >= median((split[[2]]$dist + runif(1, -200, 200)), na.rm = T))
gme <- do.call(rbind, split)
gme <- gme[with(gme, order(uid)), ]

## define the period for staging assessment (10am - 2pm, 5 hrs)
gme$stage.period <- ifelse(hour(gme$date) >= 10 & hour(gme$date) <= 14, 1, 0)
# unique id for each individual-day
gme$dayBurst <- paste(gme$id, as.Date(gme$date))

# 24-hour ag window identifies days where they used ag - dayburst helps group by day and individual 
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

#### Viz of E-C ratios ####
## Mean E-C ratio for each hour of the day for ag and non-ag days

# ratio with individual variation -- to get confidence intervals
t <- gme %>%
  filter(!is.na(ag.window.ext)) %>%
  mutate(hour = hour(date)) %>%
  group_by(subject_name, hour, ag.window.ext) %>%
  mutate(enc.day = sum(viterbi==1), meander.day = sum(viterbi==2), dw.day = sum(viterbi==3)) %>%
  mutate(ratio = enc.day/meander.day) %>%
  mutate(ratio = ifelse(ratio == Inf, enc.day, ratio))

# get confidence intervals
t <- t %>%
  group_by(hour, ag.window.ext) %>%
  summarise(mean.ratio = mean(ratio, na.rm = TRUE),
            sd = sd(ratio, na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd / sqrt(n),
         lwr.ci = mean.ratio - qt(1 - (0.05 / 2), n - 1) * se,
         upr.ci = mean.ratio + qt(1 - (0.05 / 2), n - 1) * se)

t$ag.window <- as.factor(t$ag.window.ext)
levels(t$ag.window) <- c('non-window', 'ag-window')

ggplot(t, aes(as.factor(hour), mean.ratio)) + geom_pointrange(
  aes(ymin = lwr.ci, ymax = upr.ci, color = ag.window), 
  position = position_dodge(0.5), size = .2) + 
  ggtitle("Encamped/Meandering Ratio by Hour and Ag Window")


#### Stage Loops ####

# stage tagging function (applied in k loop)
tag_stage <- function(df, ratio.seq){
  # tag all staging events
  df$stage.event <- ifelse(df$ratio > ratio.seq & df$stage.period == 1 & df$dw.day == 0, 1, 0)
  # tag ag staging events
  df$ag.stage.event <- ifelse(df$ratio > ratio.seq & df$stage.period == 1 & df$dw.day == 0 & df$ag.window.ext == 1, 1, 0)
  
  ## index staging events
  eventFlag <- ifelse(df$stage.event == 1, TRUE, FALSE)
  eventIndex <- inverse.rle(within.list(rle(eventFlag), 
                                        values[values] <- seq_along(values[values])))
  # assign a unique event index for each stage event
  df$stage.event.index <- eventIndex
  
  table <- df %>% group_by(ag.window.ext) %>%
    summarise(n.stage = length(unique(stage.event.index)),
              n.stage.adj = n.stage/length(unique(dayBurst)),
              ratio.seq = ratio.seq)
  
  table <- table %>%
    mutate(n.stage.ratio = n.stage.adj[2]/n.stage.adj[1]) %>%
    mutate(n.stage.err = n.stage.adj[1]/sum(n.stage.adj))
  
  return(table)
} 

hr.start <- c(8:16) # hour start times (used in i loop)
ratio.seq = seq(0.1, 5, 0.1) # sequence of possible E-M ratios (used in k loop)


# Run the Loop! 
system.time({

stage.table <- NULL # store accuracy results (used in k loop)
stage.table.list <- NULL

for(i in 2:9){ # window size loop
  
  # create matrix for possible start and end times 
  hr.end <- hr.start + (i-1) #subtract 1 from current index for inclusive window size
  hr.win <- subset(cbind(hr.start, hr.end), hr.end <= max(hr.start))
  
  for (j in 1:dim(hr.win)[1]) { # window period loop indexes the matrix row to define stage period
  
    for (k in 1:length(ratio.seq)){
      # define the period for staging assessment based on current window size and period
      gme$stage.period <- ifelse(hour(gme$date) >= hr.win[j,1] & hour(gme$date) <= hr.win[j,2], 1, 0)
      gme.stage <- gme %>%
        filter(!is.na(ag.window.ext)) %>%
        group_by(dayBurst, stage.period, ag.window.ext) %>%
        mutate(enc.day = sum(viterbi==1), meander.day = sum(viterbi==2), dw.day = sum(viterbi==3), n.day = n()) %>%
        mutate(ratio = enc.day/meander.day) %>%
        mutate(ratio = ifelse(ratio == Inf, enc.day, ratio)) %>%
        mutate(ratio = ifelse(is.na(ratio), 0, ratio)) # for periods with all directed-walk
      
      
      stage.table[[k]] <- tag_stage(gme.stage, ratio.seq[[k]])
      
      # add reference for window size and ratio
      stage.table[[k]]$win.siz <- i
      stage.table[[k]]$win.period <- paste(hr.win[j,1], hr.win[j,2], sep = '-')
      
    }
  }
  
  
  stage.table.list[[i]] <- stage.table
  
}
})


