#### rgee Covariate Extraction ####
#' Nathan Hahn
#' 
#' Covariates List (limited to earth engine products):
#' - MODIS NDVI - 16day
#' - Slope
#' - human footprint index (gHM)
#' - Pct forest cover (extrapolated from hansen forest loss/gain data)


library(tidyverse)
#remotes::install_github("r-spatial/rgee") # follow the directions the first time to install 
library(rgee)

ds <- as.data.frame(data.table::fread('./SSF/lab6_randsteps_20211112.csv'))
ds$uid <- ds$V1
head(ds)

# Prep dataframe for extraction -- using the step taken (t2) for cov extraction
trackingdata <- as.data.frame(cbind(ds$uid, ds$id, ds$x2_, ds$y2_))
trackingdata$timestamp <- as.POSIXct(ds$t2_) # add timestamp too!
# set column names for extraction. only necessary columns are an ID, x/y coordinates in any projection, and a time stamp.
#           Ideally should also have a uid or other unique identifer for each row for merging later
#           This code is designed to work with the following column names:
colnames(trackingdata) <- c('uid','id', 'x', 'y', 'timestamp')
trackingdata$x <- as.numeric(trackingdata$x)
trackingdata$y <- as.numeric(trackingdata$y)

# earth Engine stores time in UTC, so we will convert our timestamps to UTC 
lubridate::tz(trackingdata$timestamp) <- 'Africa/Nairobi' # set the timezone
attr(trackingdata$timestamp, "tzone") <- "UTC" # convert to UTC

# format date/time for geoJSON to work with Earth Engine
trackingdata$Date <- as.factor(trackingdata$timestamp)
trackingdata$Date <- sub(" ", "T", trackingdata$Date) #Put in a format that can be read by javascript

# get date range for subset
tmin <- min(trackingdata$timestamp)
tmax <- max(trackingdata$timestamp)
trackingdata$timestamp <- NULL # remove timestamp column

# Our data is in UTM coordinates, and earth engine works in lat-long WGS84. Lets convert it.
library(sp)
library(raster)

# create a spatial points dataframe with our current projection system
tt<-SpatialPoints(trackingdata[,3:4])
crs(tt) <- '+init=epsg:32736'
# transform to lat-long
tt <- spTransform(tt, CRS = "+proj=longlat +datum=WGS84")
tt <- as.data.frame(tt)
# add new coordinates to the tracking dataset
trackingdata$location.long <- tt$x
trackingdata$location.lat <- tt$y
# remove the old UTM coordinates
trackingdata$x <- NULL
trackingdata$y <- NULL


##### Initialize rgee #####
library(rgee)
rgee::ee_Initialize()


##### Define Functions #####
## These functions are needed for rgee extraction -- see vignette in Crego et al. 2021

#Function to add property with time in milliseconds
add_date<-function(feature) {
  date <- ee$Date(ee$String(feature$get("Date")))$millis()
  feature$set(list(date_millis=date))
}

#Join Image and Points based on a maxDifference Filter within a temporal window

#Set the filter
maxDiffFilter<-ee$Filter$maxDifference(
  difference=tempwin*24*60*60*1000, #days * hr * min * sec * milliseconds
  leftField= "date_millis", #Timestamp of the telemetry data
  rightField="system:time_start" #Image date
)

# Define the join. We implement the saveBest function for the join, which finds the image that best matches the filter (i.e., the image closest in time to the particular GPS fix location). 
saveBestJoin<-ee$Join$saveBest(
  matchKey="bestImage",
  measureKey="timeDiff"
)

#Function to add property with raster pixel value from the matched image
add_value<-function(feature){
  #Get the image selected by the join
  img1<-ee$Image(feature$get("bestImage"))$select(band)
  #Extract geometry from the feature
  point<-feature$geometry()
  #Get pixel value for each point at the desired spatial resolution (argument scale)
  pixel_value<-img1$sample(region=point, scale=250, tileScale = 16, dropNulls = F) 
  #Return the data containing pixel value and image date.
  feature$setMulti(list(PixelVal = pixel_value$first()$get(band), DateTimeImage = img1$get('system:index')))
}

# Function to remove image property from features
removeProperty<- function(feature) {
  #Get the properties of the data
  properties = feature$propertyNames()
  #Select all items except images
  selectProperties = properties$filter(ee$Filter$neq("item", "bestImage"))
  #Return selected features
  feature$select(selectProperties)
}

##### Extract Temporally-Match Covariates (NDVI) #####

## This feature in rgee package allows you to extract values for temporally varying covariates. For example matching timestamped GPS data with
##              rainfall, temperature, or other data. Here, we'll use 16-day NDVI from MODIS


### Load Image Collection - 16-day MODIS NDVI

# check data range in data
tmin
tmax

# give at least of month of buffer in case we are missing values from lots of cloudy days
start<-"2012-03-01"
end<-"2016-05-01"

# Select MODIS imagery 
# NOTE: The MODIS NDVI and EVI products are computed from atmospherically 
#       corrected bi-directional surface reflectances that have been masked 
#       for water, clouds, heavy aerosols, and cloud shadows.
imagecoll<-ee$ImageCollection('MODIS/006/MOD13Q1')$filterDate(start,end)
band <- "NDVI" #Name of the band to use. Also can select EVI from this product

#Set temporal window in days for filter. This will depend on the remote sensing data used.
tempwin <- 16 # 16 days for MODIS

## rgee needs data as an sf object. Convert tracking frame into an sf object
library(sf)
datasf <- st_as_sf(trackingdata, coords = c("location.long", "location.lat"), crs = 4326)

# Add a unique identifier to groups of data. rgee will process the data in batches because of limits to how much can be uploaded to Earth Engine at once
#           This code will make 1000 batches of 1000 rows each. 
datasf$uniq <- rep(1:1000, each=1000)[1:nrow(datasf)] #This is for up to 1 million points. To increase the max number of points, increase the value for max repetitions. To change the number of points to run per time, change the value in the argument each.

# Start the NDVI extraction and time it 
start_time <- Sys.time()
dataoutput <- data.frame()
for(x in unique(datasf$uniq)){
  data1 <- datasf %>% filter(uniq == x)
  # Send sf to GEE
  data <- sf_as_ee(data1)
  # Transform day into milliseconds
  data<-data$map(add_date)
  # Apply the join
  Data_match<-saveBestJoin$apply(data, imagecoll, maxDiffFilter)
  # Add pixel value to the data
  DataFinal<-Data_match$map(add_value)
  # Remove image property from the data
  DataFinal<-DataFinal$map(removeProperty)
  # Move GEE object into R
  temp<- ee_as_sf(DataFinal, via = 'getInfo')
  # Append
  dataoutput <- rbind(dataoutput, temp)
}
end_time <- Sys.time()

# time needed to run 100000 points -- 15 minutes
end_time - start_time

# update band name
names(dataoutput)[3] <- band
dataoutput

#### Extract Static Covariates ####
# We can also use it to extract static landscape covariates. Here, we'll use two global layers available in Earth Engine,
#         and another custom layer that has already been made in Earth Engine and saved as an asset online. As a note, 
#         you can make layers in R or other GIS software programs and then upload them to Earth Engine as an asset to use
#         this same workflow. This may be faster if you are working with very large ares (e.g. country, continent, or global scale)

## make a slope image -- using SRTM 
ee.drm <- ee$Image("USGS/SRTMGL1_003")
ee.slope <- ee$Terrain$slope(ee.drm)
# check
vis_params <- list(min=0, max=45, palette='white,black')
Map$addLayer(ee.slope, vis_params)

## make a human modification image -- using gHM layer
ee.gHM <- ee$ImageCollection("CSP/HM/GlobalHumanModification")
ee.gHM <- ee$Image(ee.gHM$first())
#check
vis_params <- list(min=0, max=1, palette='white,black')
Map$addLayer(ee.gHM, vis_params)

## add a layer from earth engine assets -- 
##          this will pull layers that I have already calculated and have stored in 
##          my Earth Engine account online
ee.forest.250 <- ee$Image('users/nhahnwa/Hansen_forest30_2019_focal_r250')
vis_params <- list(min=0, max=1, palette='white,black')
Map$addLayer(ee.forest.250, vis_params)

# Again, we split the data into 1000 chunks of 1000 rows. 
datasf$uniq <- rep(1:1000, each=1000)[1:nrow(datasf)] #This is for up to 1 million points. To increase the max number of points, increase the value for max repetitions. To change the number of points to run per time, change the value in the argument each.

# Start the extraction and time it 
start_time <- Sys.time()
dataoutput <- data.frame()
for(x in unique(datasf$uniq)){
  data1 <- datasf %>% filter(uniq == x)
  data.gHM <- ee_extract(x = ee.gHM, y = data1, sf = FALSE, scale = 1000)
  data.slope <- ee_extract(x = ee.slope, y = data1, sf = FALSE, scale = 30)
  # Append - gHM used to index uid, slope variable appended on
  temp <- inner_join(data.gHM, data.slope, by = c('uid', 'id', 'Date','uniq'), keep = F)
  dataoutput <- rbind(dataoutput, temp)
}
end_time <- Sys.time()

# time for 100000 points -- 12 minutes
end_time - start_time


#### Combine Temporal and Static Extraction ####

## As another example, we can combine extraction of all covariates at once. This makes the code more compact.
##          We will use the same layers as before but extract all covariates in the same loop.         

# Again, split the data into 1000 groups of 1000 rows each.
datasf$uniq <- rep(1:1000, each=1000)[1:nrow(datasf)] #This is for up to 1 million points. To increase the max number of points, increase the value for max repetitions. To change the number of points to run per time, change the value in the argument each.

# Start the extraction and time it 
start_time <- Sys.time()
dataoutput <- data.frame()
for(x in unique(datasf$uniq)){
  
  ## Filter to chunk by uniq index
  chunk <- datasf %>% filter(uniq == x)
  
  ## Perform static extraction
  data.gHM <- ee_extract(x = ee.gHM, y = chunk, sf = TRUE, scale = 1000)
  data.slope <- ee_extract(x = ee.slope, y = chunk, sf = TRUE, scale = 30)
  data.forest <- ee_extract(x = ee.forest.250, y = chunk, sf = TRUE, scale = 30)
  # Append back to chunk - gHM used to index uid, slope variable appended on
  slope <- data.slope$slope
  forest <- data.forest$forest
  chunk <- cbind(data.gHM, slope)
  
  ## Perform temporal extraction
  # Send sf to GEE
  data <- sf_as_ee(chunk)
  # Transform day into milliseconds
  data <- data$map(add_date)
  # Apply the join
  Data_match <- saveBestJoin$apply(data, imagecoll, maxDiffFilter)
  # Add pixel value to the data
  DataFinal <- Data_match$map(add_value)
  # Remove image property from the data
  DataFinal <- DataFinal$map(removeProperty)
  # Move GEE object into R
  temp <- ee_as_sf(DataFinal, via = 'getInfo')
  # Append
  dataoutput <- rbind(dataoutput, temp)
}
end_time <- Sys.time()

# time for 100000 points -- 35 minutes
end_time - start_time


# update band name for NDVI -- check indexing first!!!
dataoutput
names(dataoutput)[3] <- band
dataoutput

#### Save your results ####
#write.csv('')



