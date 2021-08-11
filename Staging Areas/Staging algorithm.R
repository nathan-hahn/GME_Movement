#### Staging Area Algorithm ####
library(tidyvere)
library(sf)
library(adehabitatHR)
library(recurse)

## Steps
# 1. Identify sequence runs (rle)
# 2. Tag sequences for encamped >= 6 (0/1)
# 3. Tag sequences for directed-walk >= 2 (0/1)
# 4. Algorithm uses encamped and directed walk 0/1 data to tag 'staging' events


# import state-classified data
gme <- read.csv("HMM/movdata/GMEcollars_003_HMMclassified_20201206.csv")
unique(gme$subject_name)

# filter out outlier steps 
split <- split(gme, gme$viterbi)
split[[1]] <- filter(split[[1]], dist <= median((split[[2]]$dist + runif(1, -200, 200)), na.rm = T))
split[[2]] <- filter(split[[2]], dist <= median((split[[3]]$dist + runif(1, -200, 200)), na.rm = T))
split[[3]] <- filter(split[[3]], dist >= median((split[[2]]$dist + runif(1, -200, 200)), na.rm = T))
gme <- do.call(rbind, split)
gme <- gme[with(gme, order(X)), ]


# check step length distributions
boxplot(gme$dist ~ gme$viterbi, xlab = 'cohort:state', ylab = 'step length')



#### 1. Calculate sequences runs ####
# split by ele - collar (same as HMM fitting)
split <- split(gme, gme$id)

# Identifies sequences (rle) and their length and creates a new column with the length
state.seq <- function(df, state, state.name) {
  t <- rle(as.numeric(df$viterbi))# == state)
  seqCol <- rep(t$lengths, t$lengths)
  df$seqCol <- seqCol
  # specify seqCol by state name
  colnames(df)[which(names(df) == "seqCol")] <- paste(state.name)
  return(df)
}

# add sequence counts for encamped
split.seq <- lapply(split, state.seq, state = 1, state.name = 'seq')

#### 2. Tag windows of encamped and directed walk sequences ####

# tag sequences of a specific state and of a specific length. ifelse specifies relocs that are part of a sequence of interest by 0/1 values
seq.tag <- function(df, state, window, state.name, ref = 1) {
  t <- ifelse(df$viterbi == state & df$seq >= as.numeric(window), as.numeric(ref), 0)
  df$t <- t
  # specify seqCol by state name
  colnames(df)[which(names(df) == "t")] <- paste(state.name)
  return(df)
}

# define relocs associated with encamped and directed walk sequences
split.tag <- lapply(split.seq, seq.tag, state = 1, window = 4, state.name = 'enc.seq')
split.tag <- lapply(split.tag, seq.tag, state = 3, window = 1, state.name = 'dw.seq', ref = 1)

#### 3. Define Staging Events ####
# Identify sequences of x+ encamped, followed by y+ sequences of dw
stage_id <- function(df, enc.seq = enc.seq, dw.seq = dw.seq){
  require(tidyverse)
  # define unique id for encamped/non-encamped sequences
  session <- cumsum(c(TRUE,as.logical(diff(df$enc.seq))))
  # remove IDs for the non-encamped sequences
  df$session <- ifelse(df$enc.seq == 1, session, 0)
  
  # Tag the end of staging events by shifting dw.seq back one. The end of the staging event is identified as an overlap between 1 in 'enc.seq' and 1 in 'lead' vector
  lead <- lead(df$dw.seq)
  df$stage.end <- ifelse(df$enc.seq == 1 & lead == 1, 1, 0)
  
  # Define staging events
  df <- df %>%
    group_by(session) %>%
    mutate(stage = if_else(sum(stage.end) == 1, 1, 0))
  
  # Unique ID for stage/non-stage events 
  stage.session <- cumsum(c(TRUE,as.logical(diff(df$stage))))
  # Remove IDs for non-staging events 
  df$stage.session <- ifelse(df$stage == 1, stage.session, 0)
  
  # deselect columns
  df$seq <- NULL
  df$session <- NULL
  df$stage.end <- NULL
  
  return(df)
  
}

# Apply algorithm
split.stage <- lapply(split.tag, stage_id)

# check output
df <- split.stage$`Lempiris-e0db`
t <- as.data.frame(cbind(df$enc.seq, df$stage, df$stage.session))
View(t)
unique(df$stage.session)

##### 4. Test Encamped Sequence distributions within ag windows  #####

## define 24-hour ag window definitions 


## assign unique id to each ag window event
# define unique id for window/non-window phases
eventFlag <- ifelse(gme$ag.window == 1, TRUE, FALSE)
eventIndex <- inverse.rle(within.list(rle(eventFlag), 
                                          values[values] <- seq_along(values[values])))
# remove IDs for the non-window phases
gme$ag.window.index <- eventIndex
window <- gme[gme$ag.window == 1,]
nonwindow <- gme[gme$ag.window == 0,]

# check encamped sequence runs
split <- split(window, window$ag.window.index)
split <- lapply(split, state.seq, state.name = 'seq')
window <- do.call(rbind, split)

boxplot(window$seq ~ window$viterbi)


## Calculate staging numbers
# all staging events
t <- lapply(split.stage, function(x) length(unique(x$stage.session)))
stage.ind <- as.data.frame(do.call(rbind, t))
sum(stage.ind$V1)
# ag window
t <- lapply(split.stage, function(x) filter(x, ag.window == 1))
t <- lapply(t, function(x) length(unique(x$stage.session)))
stage.ind <- as.data.frame(do.call(rbind, t))
sum(stage.ind$V1)

## Plot staging output as all staging events and ag-phase staging
map <- do.call(rbind, split.stage)
map <- filter(map, stage == 1)
map$ag.window <- as.factor(map$ag.window)

library(sf)
library(mapview)

map <- st_as_sf(map, coords = c('x', 'y'), crs = 32736)
map <- st_transform(map, crs = 4326)
mapview(map, zcol = "ag.window", burst = TRUE, cex = .5)

# only look at stages when individual was in an ag window
map2 <- filter(map, ag.window == 1)
points(map2$x, map2$y, pch = 19, col = 'red', cex = 0.5)

#### 4. Build Curves ####

## Directed Walk Sequence Testing

all <- vector()
ag.only <- vector()
dw <- seq(1,6,1)

for(i in 1:length(dw)){
  # Test window sizes
  y <- lapply(split.seq, seq.tag, state = 1, window = 4, state.name = 'enc.seq')
  y <- lapply(y, seq.tag, state = 3, window = i, state.name = 'dw.seq', ref = 1)
  z <- lapply(y, stage_id)
  
  # all staging events
  t <- lapply(z, function(x) length(unique(x$stage.session)))
  stage.ind <- as.data.frame(do.call(rbind, t))
  all[i] <- sum(stage.ind$V1)
  
  # ag window
  t <- lapply(z, function(x) filter(x, ag.window == 1))
  t <- lapply(t, function(x) length(unique(x$stage.session)))
  stage.ind <- as.data.frame(do.call(rbind, t))
  ag.only[i] <- sum(stage.ind$V1)
}

# create dataframe of results
df.dw <- as.data.frame(cbind(dw, all, ag.only))
colnames(df.dw) <- c('dw.window', 'all', 'ag.only')
melt.dw <- melt(df.dw, id.vars = c('dw.window'))

ggplot(melt.dw, aes(y = value, color = variable)) + geom_point(aes(x = factor(dw.window))) + geom_line(aes(x = dw.window)) +
  xlab('directed walk sequence length (encamped = 4)') + ylab('n staging events tagged') 


## Encamped Sequence Testing

all <- vector()
ag.only <- vector()
ec <- seq(1,12,1)
  
  for(j in 1:length(ec)) {
    # Test window sizes
    y <- lapply(split.seq, seq.tag, state = 1, window = j, state.name = 'enc.seq')
    y <- lapply(y, seq.tag, state = 3, window = 1, state.name = 'dw.seq', ref = 1)
    z <- lapply(y, stage_id)
  
    # all staging events
    t <- lapply(z, function(x) length(unique(x$stage.session)))
    stage.ind <- as.data.frame(do.call(rbind, t))
    all[j] <- sum(stage.ind$V1)
    
    # ag window
    t <- lapply(z, function(x) filter(x, ag.window == 1))
    t <- lapply(t, function(x) length(unique(x$stage.session)))
    stage.ind <- as.data.frame(do.call(rbind, t))
    ag.only[j] <- sum(stage.ind$V1)
  
  }

# create dataframe of results
df.ec <- as.data.frame(cbind(ec, all, ag.only))
colnames(df.ec) <- c('ec.window', 'all', 'ag.only')
melt.ec <- melt(df.ec, id.vars = c('ec.window'))

## 
ggplot(melt.ec, aes(y = value, color = variable)) + geom_point(aes(x = factor(ec.window))) + geom_line(aes(x = ec.window)) +
  xlab('encamped sequence length (directed walk = 1)') + ylab('n staging events tagged')


#### 5. Re-run with optimal sequence lengths ####

## Tag stages
split.tag <- lapply(split.seq, seq.tag, state = 1, window = 3, state.name = 'enc.seq')
split.tag <- lapply(split.tag, seq.tag, state = 3, window = 1, state.name = 'dw.seq', ref = 1)
split.stage <- lapply(split.tag, stage_id)


## Count staging events
# all staging events
t <- lapply(split.stage, function(x) length(unique(x$stage.session)))
events.all <- as.data.frame(do.call(rbind, t))
sum(events.all$V1)

# ag window
t <- lapply(split.stage, function(x) filter(x, ag.window == 1))
t <- lapply(t, function(x) length(unique(x$stage.session)))
events.ag <- as.data.frame(do.call(rbind, t))
sum(events.ag$V1)

## Combine lists to single dataframe
bind <- do.call(rbind, split.stage) 
bind$stage.id <- paste(bind$id, bind$stage.session, sep = '-')

bind.filter <- filter(bind, stage.session != 0) %>% droplevels()
bind.ag <- filter(bind.filter, ag.window == 1)# %>%
  #group_by(stage.id) %>% filter(n() >= 5) # make sure staging groups are greater than 5 locs for MCP building

# Write gps locs for viewing in QGIS 
#write.csv(bind.filter, './Staging Areas/events_all.csv' )
#write.csv(bind.ag, './Staging Areas/events_ag.csv' )

#### 6. Staging Polygons ####
library(sp)
library(adehabitatHR)

## All staging areas
sp <-  filter(bind.ag, !is.na(x)) %>%
  dplyr::select(x, y, stage.id) 
coordinates(sp) <-~x+y

mcp.stage <- mcp(sp, percent = 100)
mcp.stage@data <- as.data.frame(bind.ag)


## Buffer and Union of staging hulls
# convert data frame to sf object
sf.object <- filter(bind.ag, !is.na(x)) %>%
  st_as_sf(coords = c("x", "y"))

# convert MCPs to sf polygon objects - tied to the stage id
hulls <- sf.object %>%
  group_by(stage.id) %>%
  summarise( geometry = st_combine( geometry ) ) %>%
  st_convex_hull() 

# Buffer and union of hulls
hull.union <- hulls %>%
  st_cast("POLYGON") %>%
  ungroup() %>%
  # buffer
  st_buffer(dist = 600) %>% # buffer by median step length of the dataset
  # union polygons
  st_union() %>%
  # set crs
  st_set_crs(st_crs(32736))
  

##### Plotting #####
library(mapview)
mapview(hull.union, col.regions = 'red')


# save shapefiles
st_write(hull.union, 'staging.area.poly_20210808.shp')


##### Recursion #####
library(recurse)


# TODO: Make id the ag.window unique id - need to make unique id's for the ag.windows.
df.recurse <- gme %>%
  filter(ag.window == 1, subject_name == 'Fitz') %>%
  dplyr::select(c('x', 'y', 't'='date', 'id'='subject_name')) %>%
  mutate(t = as.POSIXct(t), id = as.factor(id))


# test with ag windows - best
test <- getRecursions(df.recurse, radius = 350, threshold = 100, timeunits = 'hours')

par(mfrow = c(2,1))
plot(test, df.recurse, legendPos = c(13, -10), pch = 19, cex = 0.5)
hist(test$revisits, breaks = 25, main = "", xlab = "Revisits (radius = 350m)")

boxplot(test$revisitStats$timeInside ~ hour(test$revisitStats$entranceTime))

# add to spatial object for plotting
df.recurse$revisits <- test$revisits

sf <- st_as_sf(df.recurse, coords = c('x', 'y'), crs=32736)


library(RColorBrewer)
mapview(sf, zcol = "revisits", col.regions = rev(brewer.pal(11, 'RdYlBu')), cex = 1.3, stroke = 0)

# test with all data - compare recursions stats relative to the entire dataset
df.recurse <- gme %>%
  filter(subject_name == 'Ivy') %>%
  dplyr::select(c('x', 'y', 't'='date', 'id'='subject_name')) %>%
  mutate(t = as.POSIXct(t), id = as.factor(id))

test <- getRecursions(df.recurse, radius = 350, threshold = 100, timeunits = 'hours')

plot(test, df.recurse, legendPos = c(13, -10), pch = 19, cex = 0.5)

hist(test$revisits, breaks = 20, main = "", xlab = "Revisits (radius = 2)")

boxplot(test$revisitStats$timeInside ~ hour(test$revisitStats$entranceTime))

# test with recursions to stage locations
df.recurse <- gme %>%
  filter(ag.window == 1, subject_name == 'Ivy') %>%
  dplyr::select(c('x', 'y', 't'='date', 'id'='subject_name')) %>%
  mutate(t = as.POSIXct(t), id = as.factor(id))

locations <- bind %>%
  filter(stage.session > 0) %>% dplyr::select(c('x', 'y'))

test <- getRecursionsAtLocations(df.recurse, locations = as.data.frame(locations), radius = 350, threshold = 100)
plot(test, df.recurse, pch = 19, cex = 0.5)

