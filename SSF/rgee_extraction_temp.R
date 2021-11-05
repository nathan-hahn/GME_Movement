#### rgee Covariate Extraction ####
#' Nathan Hahn
#' 
#' Covariates List (limited to earth engine products):
#' - MODIS NDVI - 16day
#' - Slope
#' - human footprint index (gHM)


library(tidyverse)
remotes::install_github("r-spatial/rgee")


df <- as.data.frame(data.table::fread('./SSF/movdata_004_randsteps_2021-10-12.csv'))
head(df)
df$uid <- df$V1

# Explore with random subset (for speed)

ds <- df[sample(nrow(df), 1000000), ]

# Prep dataframe for extraction -- using the step taken (t2) for cov extraction
trackingdata <- as.data.frame(cbind(ds$uid, ds$id, ds$x2_, ds$y2_))

# rename to x and y for streamlined updating later
colnames(trackingdata) <- c('uid','id', 'x', 'y')
trackingdata$x <- as.numeric(trackingdata$x)
trackingdata$y <- as.numeric(trackingdata$y)
# use timestamp from second step, but rename to a generic col name
trackingdata$timestamp <- as.POSIXct(ds$t2_)
lubridate::tz(trackingdata$timestamp) <- 'Africa/Nairobi'
# convert to UTC
attr(trackingdata$timestamp, "tzone") <- "UTC"

# format date/time for geoJSON
trackingdata$Date <- as.factor(trackingdata$timestamp)
trackingdata$Date <- sub(" ", "T", trackingdata$Date) #Put in a format that can be read by javascript

# get date range for subset
tmin <- min(trackingdata$timestamp)
tmax <- max(trackingdata$timestamp)
trackingdata$timestamp <- NULL # remove timestamp column

# convert to lat-long for earth engine
library(sp)
library(raster)
tt<-SpatialPoints(trackingdata[,3:4])
crs(tt) <- '+init=epsg:32736'
tt <- spTransform(tt, CRS = "+proj=longlat +datum=WGS84")
tt <- as.data.frame(tt)

trackingdata$location.long <- tt$x
trackingdata$location.lat <- tt$y

# remove UTM
trackingdata$x <- NULL
trackingdata$y <- NULL

# convert tracking frame into sf object
library(sf)
datasf <- st_as_sf(trackingdata, coords = c("location.long", "location.lat"), crs = 4326)

##### Initialize rgee #####
library(rgee)
rgee::ee_Initialize()


##### Define Functions #####
#Function to add property with time in milliseconds
add_date<-function(feature) {
  date <- ee$Date(ee$String(feature$get("Date")))$millis()
  feature$set(list(date_millis=date))
}

#Join Image and Points based on a maxDifference Filter within a temporal window

#Set temporal window in days for filter. This will depend on the remote sensing data used.
tempwin <- 16 

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

#### Load Image Collection ####

## 16-day MODIS NDVI

# check range
tmin
tmax

# give some buffer
start<-"2011-08-01"
end<-"2021-04-01"

imagecoll<-ee$ImageCollection('MODIS/006/MOD13Q1')$filterDate(start,end)
band <- "NDVI" #Name of the band to use. You can change to EVI for instance when using MOD13Q1.

##### Extract NDVI Values #####
#datasf <- datasf[1:100000,]
datasf$uniq <- rep(1:1000, each=1000)[1:nrow(datasf)] #This is for up to 1 million points. To increase the max number of points, increase the value for max repetitions. To change the number of points to run per time, change the value in the argument each.

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

# time needed to run 100000 points
end_time - start_time

# update band name
names(dataoutput)[3] <- band
dataoutput

#### Extract Multiple Covariates ####

# make a slope image -- using SRTM 
drm <- ee$Image("USGS/SRTMGL1_003")
slope <- ee$Terrain$slope(drm)

# check
visualization_params <- list(min=0, max=45, palette='white,black')
Map$addLayer(slope, visualization_params)



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