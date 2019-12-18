## HMM Data Vis Code
library(tidyverse)
library(lubridate)
library(sf)
library(sp)
library(mapview)
library(anipaths)

# load dataset as output from HMM.
# viterbi is the state: 1 - encamped, 2 - foraging, 3 - exploratory
output <- readRDS("./movdata/GMEcollars_001_used_2019-12-18.rds")



## Trajectory creation
# output is a dataframe with collar ID + locations + covariates + states
# GSME is protected area polygon
# requires package 'mapview'
# ele.sldf is the trajectory, ele.sf are the points
GSME <- st_read("./spatial data/GSE/GSE.shp")
ele.df <- filter(output, name == "Shorty") %>%
  mutate(row.id = seq.int(nrow(.))) %>%
  drop_na(x,y) %>%
  droplevels()
ele.sf <- st_as_sf(ele.df, coords = c("x", "y"), crs = 32737) # for some reason doesn't work with ..36

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

sl <- SpatialLines(l, proj4string = CRS('+proj=utm +init=epsg:32737'))
ele.sldf <- SpatialLinesDataFrame(sl, data = ele.df)

## Visualize the result by individual in mapview
# hacky way to plot. Using zcol not cooperating
mapview(ele.sf, cex = 0.1, alpha = 0, legend = TRUE) + # for legend, make points super small
  mapview(GSME, col.regions="darkgreen", color="darkgreen", alpha.regions=0.25) +
  mapview(ele.sldf, color = "grey", cex = 0.5) +
  mapview(ele.sf, color = "black", cex = 0.5)

## Path animation
library(anipaths)
source("~/Dropbox (Personal)/CSU/R/Functions/Traj_Function.R")

# Convert to latlong
output.latlong <- AddLatLong(output)

ele <- filter(output.latlong, name == "Shorty")
#ele <- filter(output.latlong, date > as.POSIXct("2018-01-01 00:00:00"))
delta.t <- "day"

# get background map
# https://console.cloud.google.com/google/maps-apis/overview
ggmap::register_google(key = "AIzaSyClJFcr4SDm0xmhdp5iuBpVE4Rqvo5vB1k")
background <- ggmap::get_googlemap(center = c(34.83504, -1.68927),
                                   zoom = 11,
                                   maptype = "satellite")

# animation function
animate_paths(paths = ele,
              delta.t = delta.t,
              coord = c("location.long", "location.lat"),
              Time.name = "date",
              #covariate = "viterbi", covariate.colors = (c("#E69F00", "#56B4E9", "#009E73")),
              ID.name = "ID",
              background = background,
              method = "mp4")



## Data Export as shapefile - takes a while to run

# LOCS
setwd("./movdata")

df <- output %>%
  drop_na(x,y) %>%
  droplevels()
locs.sf <- st_as_sf(df, coords = c("x", "y"), crs = 32737) # for some reason doesn't work with ..36
locs.split <- split(locs.sf, ele.sf$ID)

# write to disk
# full 
st_write(ele.sf, "EleCollars_211119_HMMclassified_2019-11-30", driver = "ESRI Shapefile")
# by individual
sapply(names(locs.split), 
       function (x) st_write(locs.split[[x]], paste(x, "2019-11-30", sep="_"), driver = "ESRI Shapefile"))
# rds
saveRDS(ele.sf, "EleCollars_211119_locs_2019-11-30.rds")

# TRAJ
# Have to split dataframe and build trajectories for each individual
# Final output converted to an sf object

# split dataframe and create ID column
df.split <- split(df, df$ID)
df.split <- lapply(df.split, function(x) x %>% mutate(state = as.numeric(viterbi), row.id = seq.int(nrow(.))))

# prep spatial data as a list
sf.split <- lapply(df.split, st_as_sf, coords = c("x", "y"), crs = 32737)
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
  sl.list[[i]] <- SpatialLines(l, proj4string = CRS('+proj=utm +init=epsg:32737'))
  sldf.list[[i]] <- SpatialLinesDataFrame(sl.list[[i]], data = df.split[[i]])
}
})

traj.sf <- lapply(sldf.list, st_as_sf) # convert to sf 
names(traj.sf) <- names(df.split)

# write trajectories to disk

# by individual
sapply(names(traj.sf), 
       function (x) st_write(traj.sf[[x]], paste(x, "traj_2019-11-30", sep="_"), driver = "ESRI Shapefile"))

# traj dataframe to rds
saveRDS(traj.sf, "EleCollars_211119_trajectories_2019-11-30.rds")

setwd('..')

head(output)
# pct. fixes in ag in state 2
tally(filter(output, viterbi == 1, lc.estes == 2))/tally(filter(output, lc.estes == 2))

# pct. fixes in ag in state 2
tally(filter(output, viterbi == 2, lc.estes == 2))/tally(filter(output, lc.estes == 2))

# pct. fixes in ag in state 3
tally(filter(output, viterbi == 3, lc.estes == 2))/tally(filter(output, lc.estes == 2))




