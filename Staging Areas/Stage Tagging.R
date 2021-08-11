##### Ag Window Tagging #####
library(dplyr)
library(lubridate)
library(sf)

# import state-classified data
gme <- read.csv("HMM/movdata/GMEcollars_003_HMMclassified_20201206.csv") 

gme <- gme %>%
  arrange(subject_name, date)
gme$X <- rownames(gme) # unique id for each relocation
unique(gme$subject_name)

# filter out outlier steps 
split <- split(gme, gme$viterbi)
split[[1]] <- filter(split[[1]], dist <= median((split[[2]]$dist + runif(1, -200, 200)), na.rm = T))
split[[2]] <- filter(split[[2]], dist <= median((split[[3]]$dist + runif(1, -200, 200)), na.rm = T))
split[[3]] <- filter(split[[3]], dist >= median((split[[2]]$dist + runif(1, -200, 200)), na.rm = T))
gme <- do.call(rbind, split)

# check step length distributions
boxplot(gme$dist ~ gme$viterbi, xlab = 'cohort:state', ylab = 'step length')

# Moving Window - Expanded
# roll.df <- gme
# split <- split(roll.df, roll.df$subject_name)

# window <- c(90*24) # in hours. Split to before and after when align = center
# 
# output <- NULL
# for(i in 1:length(window)){
#   # calculate moving window stats - see rollstats function
#   output[[i]] <- window_stats(df.list = split, window = window[[i]], align = 'center')
#   
# }



# Dataset for tagging - Ivy 2018
filter <- filter(gme, subject_name %in% c('Ivy') & year(date) == 2018)

head(filter)


# Unique ID for ag/non-ag events 
stage.session <- cumsum(c(TRUE,as.logical(diff(filter$ag.window))))
# Remove IDs for non-staging events 
filter$stage.session <- ifelse(filter$ag.window == 1, stage.session, 0)
# create list of ag.window trajectories
ag.filter <- filter(filter, ag.window == 1)


##### TRAJ - create ag.window trajectories #####

## full Ivy 2018 traj
ele.sf <- st_as_sf(filter, coords = c("x", "y"), crs = 32736) # for some reason doesn't work with ..36
split <- split(ele.sf, ele.sf$viterbi)
# create sldf of trajectories. Key: each step must be it's own line!! 
# Creates begin and end coords that are used to build individual steps
# https://stackoverflow.com/questions/20531066/convert-begin-and-end-coordinates-into-spatial-lines-in-r
coords <- ele.sf %>% st_coordinates()
begin.coord <- as.matrix(coords[1:nrow(coords)-1,])
end.coord <- as.matrix(coords[2:nrow(coords),])

l <- vector("list", nrow(begin.coord))
for (i in seq_along(l)) {
  l[[i]] <- Lines(list(Line(rbind(begin.coord[i, ], end.coord[i,]))), as.character(i))
}

sl <- SpatialLines(l, proj4string = CRS('+proj=utm +init=epsg:32736'))
ele.sldf <- SpatialLinesDataFrame(sl, data = ele.df)

## ag window traj

# split dataframe and create ID column
df.split <- split(ag.filter, ag.filter$stage.session)
df.split <- lapply(df.split, function(x) x %>% mutate(row.id = seq.int(nrow(.))) %>% drop_na(x,y))

# prep spatial data as a list
sf.split <- lapply(df.split, st_as_sf, coords = c("x", "y"), crs = 32736)
coords <- lapply(sf.split, st_coordinates)
begin.coord <- lapply(coords, function(x) as.matrix(x[1:nrow(x)-1,]))
end.coord <- lapply(coords, function(x) as.matrix(x[2:nrow(x),]))

# for loop creates lines in inner loop, and builds the trajectories in outter loop
sl.list <- vector("list", length(df.split))
sldf.list <- vector("list", length(df.split))

system.time({
  for(i in 1:length(df.split)){
    l <- vector("list", nrow(begin.coord[[i]])) # adaptive list for lines 
    #create lines
    for (j in seq_along(l)) {
      l[[j]] <- Lines(list(Line(rbind(begin.coord[[i]][j, ], end.coord[[i]][j,]))), as.character(j))
    }
    # create trjaectories
    sl.list[[i]] <- SpatialLines(l, proj4string = CRS('+proj=utm +init=epsg:32736'))
    sldf.list[[i]] <- SpatialLinesDataFrame(sl.list[[i]], data = df.split[[i]])
  }
})


traj.sf <- lapply(sldf.list, st_as_sf) # convert to sf 
names(traj.sf) <- names(df.split)

##### Write to Disk #####

## Points
# write full traj - Ivy 2018
st_write(st_as_sf(filter), 'Full_Points_Ivy2019', driver = 'ESRI Shapefile')

# by ag window
locs.sf <- lapply(df.split, st_as_sf, coords = c("x", "y"), crs = 32736)
sapply(names(locs.sf),
       function (x) st_write(locs.sf[[x]], paste(x, 'locs_ivy2018_agwindow', sep='_'), driver = 'ESRI Shapefile'))

## Traj
# by ag window
sapply(names(traj.sf), 
       function (x) st_write(traj.sf[[x]], paste(x, "traj_Ivy2018_agwindow", sep="_"), driver = "ESRI Shapefile"))

# write full trajectory - Ivy 2018
st_write(st_as_sf(ele.sldf), 'Full_Traj_Ivy2018', driver = 'ESRI Shapefile')

# traj dataframe to rds
saveRDS(traj.sf, "EleCollars_211119_trajectories_2019-11-30.rds")
