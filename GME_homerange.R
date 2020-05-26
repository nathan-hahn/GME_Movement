##### GME Home Ranges #####

library(adehabitatHR)

# GME movement data - 2019
gme <- readRDS("./movdata/GMEcollars_002_clean_2020-05-21.rds")

# prep spdf
t <-  filter(gme, !is.na(x)) %>%
  dplyr::select(x, y, subject_name) %>%
  filter(subject_name == c("Ivy", "Fred")) 
coordinates(t) <-~x+y

# build MCP's 
mcp <- mcp(t[,1], percent = 100)
plot(mcp)

# build MCP's by subject_name
split <- split(t, t$subject_name)
mcp <- lapply(split, mcp)


##### cut raster stack by 
s <- readRDS("./spatial data/GME_rasterStack_20200522.rds")

# mcp.sf <- st_as_sf(t)
# st_crs(mcp.sf) = 32736

# clip rasters to homeranges
cr <- crop(s, extent(mcp), snap="out")                    
fr <- rasterize(mcp, cr)   
lr <- mask(x=cr, mask=fr)



mask.poly.raster <- function(poly, raster) {
  # crop raster to polygon extent
  cr <- crop(raster, extent(poly), snap = "out")
  # # rasterize the polygon for masking
  fr <- rasterize(poly, cr)
  # # mask the 
  lr <- mask(x = cr, mask = fr)
  return(lr)
}


system.time(mcp.rast <- lapply(mcp, mask.poly.raster, raster = s$change03_181_reclassMara_2019.11.22))




# landscape metrics
library(landscapemetrics)





