
#### Ag Network Connectivity ####

# Set environment to EAT
Sys.setenv(TZ="Africa/Nairobi") 

# install packages
#install.packages('devtools') # this allows you to download packages from github
#devtools::install_github("BastilleRousseau/moveNT")
#install.packages('adehabitatLT')

# load packages
library(moveNT)
library(mapview)
library(adehabitatLT)
library(dplyr)

#' #2 - Applying network theory to movement data (package *moveNT*)
#' This training module presents a new approach to study animal movement 
#' based on network theory. Network (graph) theory is a popular analytical framework to 
#' characterize the structure and dynamics among discrete objects (nodes) and the connections 
#' among them (edges) and is particularly effective at identifying critical hubs and patterns 
#' of connectivity. For example, a network can be represented from the different airports
#' within a country (hubs) and the connections among them (edges). Network theory can be used to 
#' calculate network metrics that would help identify airports with higher connectivity or 
#' importance. Similarly, the identification of attributes related to connectivity is a critical
#' component of understanding animal movement, yet network theory has rarely been applied 
#' directly to animal relocation data. We develop an approach that allows the analysis of 
#' movement data using network theory by defining occupied pixels as nodes and connections 
#' among these pixels as edges. By identifying critical nodes, our approach provides a robust 
#' quantitative framework to identify local properties of space use. The approach is intuitive, 
#' and can be implemented across imperfectly sampled or large-scale data sets efficiently, 
#' providing a new valuable framework. The main steps of the approach are presented in Figure 2.
#' Basically, a grid is overlayed over the GPS locations and transitions among all pixels 
#' (estimated from the movement data) are tallied into and adjacency matrix. From this matrix, 
#' it is possible to calculate various network metrics. 

#' #A- Opening data and creation of trajectory object. 
#' We will load and open the elephant data and create a trajectory object using the 
#' `adeHabitatLT` package

#setwd("./Move-Course/Lab 5 Networks") # set working directory as needed

## Load the data
df <- as.data.frame(data.table::fread('./movdata/GMEcollars_004_stageclassified_20211019.csv'))
df <- df %>%
  dplyr::select(uid, subject_name, year.cuts, tactic.season, tactic.agg, x, y, date, dist, id, subject_sex, region, 
                subject_ageClass, ag.used, viterbi, ag.window, dayBurst, ag.window.ext, vote)

head(df)

# how many elephants are there?
unique(df$id)

# filter to ag windows
elephants <- df[df$ag.window == 1,]
#elephants <- elephants[elephants$fixType != 'drift',]

# what is the size of the dataset?
dim(elephants)

## Create trajectory using adehabitatLT package
library(adehabitatLT)
# format dates
date<-as.POSIXct(strptime(as.character(elephants$date),"%Y-%m-%d %H:%M:%S"), tz = 'Africa/Nairobi')
date2 <- as.POSIXct(round(date, units="hours"))
elephants$date <- date2

# check for duplicates
temp_id<-paste(elephants$id, elephants$date)
elephants<-elephants[!duplicated(temp_id),]

# create trajectory - traj cut has already been conducted
traj<-as.ltraj(xy=elephants[,6:7], date=elephants$date, id=elephants$id)
ref.dat<-strptime("2000-01-01 00:00","%Y-%m-%d %H:%M")
traj2<- sett0(traj, ref.dat, 1, units = "hour", tol=0.5)
traj3<-cutltraj(traj2, "dt > 3600*24",nextr = T) # will essentially create a burst for each ag.window -- important for interpolation step
traj4<-redisltraj(traj3, 3600, type="time")

#' #D- Calculating network metrics for multiple individuals

# create user defined grid - first we need to define the bounding area of the GPS data 
# I'm going to use the min and max of the x and y coordinates in the tracking data to define
# the bounding box
tt<-SpatialPoints(ld(traj4)[,1:2])
tt1<-apply(coordinates(tt), 2, min)
tt2<-apply(coordinates(tt), 2, max)

# create the raster 
# I am using the floor and ceiling commands to round the coordinates down/up
# We also define the grid resolution based on rounded median step length (250m)
median(ld(traj4)$dist, na.rm=T) # check median step length
grid.res = 250 # round down a bit for a standard grid size
ras250<-raster(xmn=floor(tt1[1]), ymn=floor(tt1[2]),xmx=ceiling(tt2[1]), ymx=ceiling(tt2[2]), res= grid.res) 

# create adjacency matrix
adj_patches2<-traj2adj(traj, res=250, grid=ras250) #Grid size based on median

#' #E- Looping over all individuals *loop*
#' The function *loop* is a wrapper of *traj2adj* and *adj2stack* applied to all individuals 
#' within a trajectory. The function will keep the same grid for all individuals. The user 
#' simply need to specify the trajectory object and the grid size. The loop function also adds 
#' additional movement properties regarding speed, absolute angle, and turning angle. We will 
#' use *loop* to apply the network approach to our XX elephants and display the results for the
#' second and third elephants. 

# apply the loop function
out1<-loop(traj, res = 250)

# name the list index
names(out1) <- unique(elephants$id)

# plot results
plot(out1[[2]]) #Plot the second elephant
plot(out1[[3]]) #Plot the third elephant


#' ##F- Mosaic individual
#' 
#' Even if the the function *loop* perform the analysis to every individuals, the outputs 
#' produced are at the individual-level. We can see this by looking at the list elements in the
#' object `out1` - there is one for each individual. The function *mosaic_network* can combine 
#' the different individual levels into a single raster representation. When multiple 
#' individuals overlap, *mosaic_network* applies a function (mean or max) to calculate a 
#' population-level value for that pixel. To use the function, the user needs to specify which 
#' variable to mosaic (using index), whether to scale the individual layers (recommended) and 
#' the function to apply. We recommend to use mean for degree and weight and max for the 
#' betweenness. 
#' 
#' * QUESTION: Once you have run the code below, try changing the function (mean, median, max) 
#' to see how it affects the outputs. In terms of the ecology, what is the difference between 
#' summarizing individual differences in landscape use by mean vs. max? 

mean_weight<-mosaic_network(out1, index=2, sc=T, fun=mean) #Perform mean weight (not-interpolated)
plot(mean_weight, main = 'mean weight')

mean_degree<-mosaic_network(out1, index=4, sc=T, fun=mean) #Perform mean weight (not-interpolated)
plot(mean_degree, main = 'mean degree')

max_between<-mosaic_network(out1, index=5, sc=T, fun=max) #Perform max weight (not-interpolated)
plot(max_between, main = 'max betweeness')

max_loop <- mosaic_network(out1, index=3, sc=T, fun=max)

# for any of the rasters, use this code to plot in mapview for closer inspection
r.plot <- max_loop
crs(r.plot) <- '+init=epsg:32736' 

mapview(r.plot, col.regions = RColorBrewer::brewer.pal(9, "YlOrRd"))


## Add staging areas on top
stage.sp <- sf::st_as_sf(elephants, coords = c('x', 'y'), crs = 32736)
#crs(stage.sp) <- '+init=epsg:32736'

mapview(r.plot, col.regions = RColorBrewer::brewer.pal(9, "YlOrRd")) + mapview(stage.sp)

#' ##F- Linear interpolation 
#' As can be seen in the last plot produced, one of the limitations of the current approach is 
#' that it creates gaps in areas where no locations are observed (only pixels with gps locations
#' in them have values). This can sometimes limit interpretability or the visual appeal of the
#' maps produced. To assist with this, we created a linear interpolation approach that can be
#' applied to the individual level network calculation (i.e. after *loop*). The interpolation
#' linearly interpolate each step (i.e. straight line) and assign the network metric of each 
#' starting location to the whole step. When multiples overlap in a pixel, a function is 
#' applied to summarize these steps (e.g. mean or max). This function will take an output from
#' *loop* and performed the interpolation for five metrics (weight, degree, betweenness, speed,
#' and turning angles). We recommend to take the mean for weight, degree, betweenness, and 
#' speed, the max for betweenness, and the dot-product for the turning angles (default).   

# The interpolation function takes awhile to run, but you can try it if you have time
out2<-interpolation(traj4, out1)

# As a shortcut, I provided the interpolated data. Read in RData file with the interpolation
out2 <- readRDS('ag_network_interpolation.RDS')

mean_mean_degree <- mosaic_network(out2, index=2, sc=T, fun=mean)
max_max_between <- mosaic_network(out2, index=3, sc=T, fun=max)
mean_mean_speed <- mosaic_network(out2, index=4, sc=T, fun=mean)
mean_dot_TA <- mosaic_network(out2, index=5, sc=T, fun=mean)

par(mfrow=c(2,2))
plot(mean_mean_degree, "Degree", main = "Degree")
plot(max_max_between, "Betweenness", main = "Betweeness")
plot(mean_mean_speed, "Speed", main = "Speed")
plot(mean_dot_TA, "Directionality", main = "Directionality")

# plot interactive connectivity map
r.plot <- max_max_between
crs(r.plot) <- '+init=epsg:32736' 
stage.sp <- sf::st_as_sf(elephants, coords = c('x', 'y'), crs = 32736)

mapview(r.plot, col.regions = RColorBrewer::brewer.pal(9, "YlOrRd")) #+ mapview(stage.sp, cex = 1.5, col.region = black) ## Add staging areas on top
  

#' These four layers are showing interpolated and mosaicked population-level network or 
#' movement properties. Again, these raster could be exported to be opened in ArcGIS using 
#' *writeRaster*. 
#' 
#' * QUESTION: What is the difference between the interpolated and non-interpolated rasters?
#' What are the pros and cons of each approach to displaying the data? Hint: Plot the 
#' interpolated raster with the GPS movement data.


library(sf)
library(ggplot2)
gme <- sf::st_read('~/Dropbox (Personal)/CSU/GME_Movement/spatial data/GSE/GSE_2020.shp') 
stage.relocs <- filter(df, vote == 1) %>%
  st_as_sf(coords = c('x','y'), crs = 32736)
ag.relocs <- elephants %>% 
  filter(vote == 0) %>%
  st_as_sf(coords = c('x','y'), crs = 32736) 

mapview::mapview(gme, col.region = , alpha = 0.2) + mapview::mapview(stage.relocs, cex = 1.5, col.region = 'red') 

library(tmap)

# define network metric to plot - betweeness - and give it a projection system
r.plot <- max_max_between
crs(r.plot) <- '+init=epsg:32736' 

stage.points <- tm_shape(gme) + tm_polygons(col = 'pa_status', palette = RColorBrewer::brewer.pal(3, 'Greens'), title = 'protected areas') +
  tm_shape(ag.relocs) + tm_dots(size = 0.01, alpha = 0.03, col = '#FEE0D2') + 
  tm_shape(stage.relocs) + tm_dots(size = 0.01, alpha = 0.02, col = '#DE2D26', title = 'staging relocs') + 
  # add manual legend for tm_dot layers
  tm_add_legend(title = 'staging relocs', type = 'symbol', labels = c('ag day relocs', 'staging relocs'), col = c('#FEE0D2', '#DE2D26')) +
  # add map title
  tm_layout(title = 'Staging Relocations')
stage.points


ag.betweeness <- tm_shape(gme) + tm_polygons(col = 'pa_status', palette = RColorBrewer::brewer.pal(3, 'Greens'), title = 'serengeti-mara') +
  tm_shape(r.plot, raster.downsample = FALSE) + tm_raster(palette = RColorBrewer::brewer.pal(9, 'YlOrRd'), title = 'max_max_betweeness') + 
  tm_scale_bar() + tm_layout(title = 'Network Betweeness')
ag.betweeness

t <- tmap_arrange(stage.points, ag.betweeness)
tmap_save(t, "Staging Panel.png", dpi = 300)
tmap_save(stage.points, "Staging Relocs Map.png", dpi = 300)


