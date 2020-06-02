##### GME Home Ranges #####

library(adehabitatHR)

# GME movement data - 2019
gme <- readRDS("./movdata/GMEcollars_002_clean_2020-05-21.rds")

# prep spdf
t <-  filter(gme, !is.na(x)) %>%
  dplyr::select(x, y, subject_name) 
  #filter(subject_name %in% c("Tressa", "Fred", "Olchoda")) 
coordinates(t) <-~x+y

# build MCP's by subject_name
split <- split(t, t$subject_name)
mcp <- lapply(split, mcp)

##### cut raster stack by 
s <- readRDS("./spatial data/GME_rasterStack_20200522.rds")

mask.poly.raster <- function(poly, raster) {
  # crop raster to polygon extent
  cr <- crop(raster, extent(poly), snap = "out")
  # mask the raster by mcp 
  lr <- mask(x = cr, mask = poly)
  return(lr)
}

# mask landcover raster to homeranges
system.time(
  mcp.rast <- lapply(mcp, mask.poly.raster, raster = s$change03_181_reclassMara_2019.11.22)
)

saveRDS(mcp.rast, paste0("./homerange/GME_mcp_rasters_002", Sys.Date(),".rds"))


# landscape metrics
library(landscapemetrics)
library(landscapetools)

lsm_l_ed(mcp.rast$Olchoda)
lsm_l_ed(mcp.rast$Kimbizwa)

metrics <- lapply(mcp.rast, lsm_l_ed)
names(metrics) <- levels(as.factor(gme$subject_name))
t <- do.call("rbind", metrics)
t$subject_name <- levels(as.factor(gme$subject_name))

par(mfrow = c(2,2))
cuts = c(0,1,2,3,4)
pal = colorRampPalette(c('white', "yellow", "dark green", "dark grey", "grey"))
plot(mcp.rast$Alina, main = paste0("Sporadic: Alina, ed=", round(metrics$Alina$value, 2)),
     col = pal(5))
plot(mcp.rast$Fred, main = paste0("Mixed Seasonal: Fred, ed=", round(metrics$Fred$value,2)),
     col = pal(5))
plot(mcp.rast$Kimbizwa, main = paste0("Seasonal: Kimbizwa, ed=", round(metrics$Kimbizwa$value,2)),
     col = pal(5))
plot(mcp.rast$Olchoda, main = paste0("Habitual: Olchoda, ed=", round(metrics$Olchoda$value, 2)),
     col = pal(5))

