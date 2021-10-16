#### Stage Test - E-M ratio ####

library(tidyverse)
library(lubridate)

set.seed(15)
#setwd('~/Dropbox/CSU/GME_Movement') # for ssh

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
#gme$stage.period <- ifelse(hour(gme$date) >= 10 & hour(gme$date) <= 14, 1, 0)
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

#hr.start <- c(8:16) # hour start times (used in i loop) (i in 1:8)
#ratio.seq = seq(0.1, 10, 0.1) # sequence of possible E-M ratios (used in k loop)

# get testing dataframe
#gme <- filter(gme, subject_name == "Fred")

hr.start <- c(6:17) # increase to all daylight hours (i in 1:11)
ratio.seq = seq(0.1, 13, 0.1) # sequence of possible E-M ratios (used in k loop)

#### Stage Test Loop ####
system.time({

# store each result in matrix
result.matrix <- matrix(ncol = 8)

for(i in 1:11){ # window size loop - index adds to hr.start 
  
  # create matrix for possible start and end times 
  hr.end <- hr.start + (i) #add i index to hr.start hour for inclusive window size (e.g. 6am fix represents before from 6-7am)
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
      
      
      stage.table <- tag_stage(gme.stage, ratio.seq[[k]])
      
      # add reference for window size and ratio
      stage.table$win.siz <- i+1
      stage.table$win.period <- paste(hr.win[j,1], hr.win[j,2], sep = '-')
      
      # bind each stage.table
      result.matrix <- rbind(result.matrix, as.matrix(stage.table))
      
    } # close k loop
  } # close j loop
} # close i loop
}) # close system.time (175 minutes) -- what is increase when adding 4 hours and 4 extra ratios

# Save/Read result
saveRDS(result.matrix, 'stage_loop_result i = 1:11, hrstart=6:17, seq = 13, seed15.RDS')
saveRDS(gme, 'stage_loop_result_ratio_df_15.RDS')
#result.matrix <- readRDS('stage_loop_result i = 1:11, hrstart=6:17, seq = 13.RDS')

# # convert matrix to dataframe for viz
# result.df <- as.data.frame(result.matrix)
# result.df <- mutate_at(result.df, c('ag.window.ext', 'n.stage', 'n.stage.adj', 'ratio.seq', 'n.stage.ratio',
#                                     'n.stage.err', 'win.siz'), .funs = as.numeric)
# result.df <- result.df[-1,]
# 
# # remove repetitive sequences
# result.df <- result.df %>%
#   group_by(ag.window.ext) %>%
#   filter(n.stage > 2) %>%
#   mutate(n.stage.ratio = ifelse(n.stage.ratio == lag(n.stage.ratio), NA, n.stage.ratio)) %>%
#   mutate(n.stage.err = ifelse(n.stage.err == lag(n.stage.err), NA, n.stage.err))
# 
# 
# 
# #### Plot Loop Results ####
# 
# plot.df <- filter(result.df, !is.na(n.stage.err) & ag.window.ext == 1)
# 
# ggplot(plot.df, aes(x = ratio.seq, y = (1-n.stage.err), group = win.period)) +
#   geom_line(aes(color = win.period)) +
#   #geom_point(aes(color = win.period), position = 'jitter') +
#   facet_wrap(.~win.siz)
# 
# 
# ### 3D plot
# library(plotly)
# 
# split <- split(plot.df, plot.df$win.siz)
# names(split) <- as.character(c(2:9))
# 
# 
# fig <- plot_ly(plot.df, x = ~ratio.seq, y = ~win.period, z = ~(1-n.stage.err),
#                marker = list(color = ~n.stage, colorscale = c('#FFE1A1', '#683531'), showscale = TRUE))
# fig <- fig %>% add_markers()
# fig <- fig %>% layout(scene = list(xaxis = list(title = 'E-C ratio'),
#                                    yaxis = list(title = 'Window Period', type = 'category'),
#                                    zaxis = list(title = 'Accuracy')),
#                       annotations = list(
#                         x = 1.04,
#                         y = 1.02,
#                         text = 'Number of Stage Events',
#                         xref = 'paper',
#                         yref = 'paper',
#                         showarrow = TRUE
#                       ))
# fig
# 
# fig %>% add_surface()
# 
# 
# 
# 
# #### Optimal Sequence Tagging ####
# 
# gme$stage.period <- ifelse(hour(gme$date) >= 10 & hour(gme$date) <= 14, 1, 0)
# 
# gme.stage <- gme %>%
#   filter(!is.na(ag.window.ext)) %>%
#   group_by(dayBurst, stage.period, ag.window.ext) %>%
#   mutate(enc.day = sum(viterbi==1), meander.day = sum(viterbi==2), dw.day = sum(viterbi==3), n.day = n()) %>%
#   mutate(ratio = enc.day/meander.day) %>%
#   mutate(ratio = ifelse(ratio == Inf, enc.day, ratio)) %>%
#   mutate(ratio = ifelse(is.na(ratio), 0, ratio)) # for periods with all directed-walk
# 
# ##### Tag staging events with optimal sequence #####
# ratio.threshold = 1
# 
# # tag all staging events
# gme.stage$stage.event <- ifelse(gme.stage$ratio > ratio.threshold & gme.stage$stage.period == 1 & gme.stage$dw.day == 0, 1, 0)
# # tag ag staging events
# gme.stage$ag.stage.event <- ifelse(gme.stage$ratio > ratio.threshold & gme.stage$stage.period == 1 & gme.stage$dw.day == 0 & gme.stage$ag.window.ext == 1, 1, 0)
# 
# ## index staging events
# eventFlag <- ifelse(gme.stage$stage.event == 1, TRUE, FALSE)
# eventIndex <- inverse.rle(within.list(rle(eventFlag),
#                                       values[values] <- seq_along(values[values])))
# # assign a unique event index for each stage event
# gme.stage$stage.event.index <- eventIndex
# 
# length(unique(gme.stage$stage.event.index))
# 
# 
# ## Plot
# plot(gme.stage[gme.stage$ag.stage.event == 1,]$x, gme.stage[gme.stage$ag.stage.event == 1,]$y)
# 
# 
# ## Ag Staging Density
# library(sp)
# library(adehabitatHR)
# library(raster)
# 
# # create dataframe with all ag stage event relocs
# ag.sp <- gme.stage[gme.stage$ag.stage.event == 1,]
# coordinates(ag.sp) <- ~x+y # create a SpatialPointsDataFrame
# proj4string(ag.sp) <- CRS("+init=epsg:32636")
# 
# # kde - eval h using LSCV
# kde <- kernelUD(ag.sp[,'ag.stage.event'], h = "LSCV",
#                 kern = 'bivnorm', grid = 1000) # reference
# 
# #par(mfrow = c(2,2))
# for(i in 1:length(kde)){
#   rast <- raster(as(kde[[i]], "SpatialPixelsDataFrame"))
#   plot(rast, main = paste(names(kde[i])))
# }
# 
# # get the UD
# ud <- getvolumeUD(kde)
# 
# # plot on map
# library(mapview)
# rast <- raster(as(ud$`1`,"SpatialPixelsDataFrame"))
# rast[rast>99] <- NA
# mapview(rast) + mapview(ag.sp, cex = 1.5, layer.name = 'stage.locs') #change legend: layer.name = 'X47801a'
# 
# 
# 
# 
# 
