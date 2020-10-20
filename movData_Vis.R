## HMM Data Vis Code
library(tidyverse)
library(lubridate)
library(sf)
library(sp)
library(mapview)
library(ggmap)
library(anipaths)

#TODO: Paralelize trajectory building with foreach %dopar%
#TODO: Animate paths. Issue with google API key

# load dataset as output from HMM.
# viterbi is the state: 1 - encamped, 2 - foraging, 3 - exploratory

hmm.df <- readRDS("./HMM/model results/Full Pop July2020/m8_pop_Julyfinal.rds") # output from hmm model (no squared term)
output <- readRDS("./HMM/model results/Full Pop July2020/GMEcollars_002_population_original.rds") # dataset prior to data transforms

# add viterbi estimates to the original dataframe -- revert to non-standardized covs and steps
library(momentuHMM)
output$viterbi <- viterbi(hmm.df)

## Trajectory creation
# output is a dataframe with collar ID + locations + covariates + states
# GSME is protected area polygon
# requires package 'mapview'
# ele.sldf is the trajectory, ele.sf are the points
GSME <- st_read("./spatial data/GSE/GSE.shp")
GMF <- st_read("./spatial data/ag/ag_vec_change03_reclassMara.shp")
ele.df <- filter(output, subject_name == "Ivy") %>%
  mutate(row.id = seq.int(nrow(.))) %>%
  drop_na(x,y) %>%
  droplevels()
ele.sf <- st_as_sf(ele.df, coords = c("x", "y"), crs = 32736) # for some reason doesn't work with ..36
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

t1 <- split[[1]]
t2 <- split[[2]]
t3 <- split[[3]]

## Visualize the result by individual in mapview
# hacky way to plot. Using zcol not cooperating
mapview(GSME, col.regions="darkgreen", color="black", alpha.regions=0, cex = 3) +
  mapview(ele.sldf, color = "grey", cex = 0.1) +
  mapview(t3, color = "#009E73", cex = 1.5) +
  mapview(t2, color = "#56B4E9", cex = 1.5) +
  mapview(t1, color = "#E69F00", cex = 1.5)

pal <- c("#E69F00", "#56B4E9", "#009E73") # color pallet 
#mapview(ele.sf, zcol = "viterbi", col.region = pal, cex = 0.1, alpha = 0, legend = TRUE) + # for legend, make points super small
  mapview(GSME, col.regions="darkgreen", color="black", alpha.regions=0) +
  mapview(GMF, col.regions="yellow", color="yellow", alpha.regions=0.5) +
  mapview(ele.sldf, color = "darkgrey", cex = 0.1) + 
  mapview(filter(ele.sf, viterbi == 2), color = pal[2], cex = 1.5, alpha = 1, legend = FALSE) +
  mapview(filter(ele.sf, viterbi == 3), color = pal[3], cex = 1.5, alpha = 1, legend = FALSE) +
  mapview(filter(ele.sf, viterbi == 1), color = pal[1], cex = 1.5, alpha = 1, legend = FALSE) 


st_write(ele.sf, 'Ivy_sf.shp', driver = 'ESRI Shapefile')
st_write(st_as_sf(ele.sldf), 'Ivy_sldf', driver = 'ESRI Shapefile')

## Path animation
source("C:/Users/nhahn/Dropbox (Personal)/CSU/R/Functions/Traj_Function.R")

# Convert to latlong
output.latlong <- AddLatLong(output)

ele <- filter(output.latlong, name == "Shorty")
#ele <- filter(output.latlong, date > as.POSIXct("2018-01-01 00:00:00"))
delta.t <- "day"

# get background map
# https://console.cloud.google.com/google/maps-apis/overview
ggmap::register_google(key = "AIzaSyAyZnyz0E5teo9SbaKvyvoxkV1sAz-ty10", write = TRUE)
background <- ggmap::get_googlemap(center = c(34.83504, -1.68927),
                                   zoom = 11,
                                   maptype = "satellite")

# animation function
animate_paths(paths = ele,
              delta.t = delta.t,
              coord = c("location.long", "location.lat"),
              Time.name = "date",
              #covariate = "viterbi", covariate.colors = (c("#E69F00", "#56B4E9", "#009E73")),
              ID.name = "name",
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
df.split <- split(df, df$name)
df.split <- lapply(df.split, function(x) x %>% mutate(row.id = seq.int(nrow(.))))

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




