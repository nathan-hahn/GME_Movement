##### Stage KDE #####

library(tidyverse)
library(adehabitatHR)

# import state-classified data
gme <- read.csv("HMM/movdata/GMEcollars_003_HMMclassified_20201206.csv") 

gme <- gme %>%
  arrange(subject_name, date)
gme$X <- rownames(gme) # unique id for each relocation
unique(gme$subject_name)

ag.phase <- filter(gme, ag.window == 1)

study.area <- CRS('+proj=utm +init=epsg:32736')
coordinates(ag.phase) = c("x", "y")
proj4string(ag.phase) <- CRS("+init=epsg:32736")


## KDE for ag phase areas
kud <- kernelUD(ag.phase[,'subject_name'], grid = 250)

ud <- getvolumeUD(kud)

library(mapview)
library(raster)
# quickly plot the raster on the map
rast <- raster(as(ud$Fred,"SpatialPixelsDataFrame"))
rast[rast>30] <- NA
mapview(rast, alpha = 0.7)


kerneloverlaphr(kud, method = 'HR') # proportion of home range overlap

overlap.30 <- kerneloverlap(ag.phase[,'subject_name'], grid = 250, method = 'HR', percent = 30)
overlap.30


traj <- as.ltraj(xy = cbind(ag.phase$x, ag.phase$y), date = as.POSIXct(ag.phase$date), id = ag.phase$subject_name)
move <- move(traj)
# fit the dynamic brownian bridge using 500 meter raster and 60 minute time steps - takes ~10 hrs
system.time(
dBB.ele <- brownian.bridge.dyn(move, ext=1.2, raster=250, location.error=20, time.step = 480/15)
)


