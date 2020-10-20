#' ## Home Range Analysis
#' Test for differences in landscape and ag structure between individual homeranges
#' Are there significant differences in the homerange landscape structure between the ag tactic groups?
#' 
#' 1. MCP homeranges
#' 2. Landscape Metrics
#' 3. ANOVA on Ag Tactics


library(dplyr, warn.conflicts = FALSE, quietly = TRUE)
library(adehabitatHR, warn.conflicts = FALSE, verbose = FALSE, quietly = TRUE)
library(raster, warn.conflicts = FALSE, quietly = TRUE)
library(landscapemetrics, warn.conflicts = FALSE, quietly = TRUE)
library(car, warn.conflicts = FALSE, quietly = TRUE)
library(ggplot2, warn.conflicts = FALSE, quietly = TRUE)

source("./GME_functions.R")

# GME movement data - 2019
gme <- readRDS("./GMM/GMEcollars_002_usedClust_2020-07-14.rds")
#stack <- readRDS("./spatial data/GME_rasterStack_20200602.rds")


#' ### Build homeranges
# prep spdf
t <-  filter(gme, !is.na(x)) %>%
  dplyr::select(x, y, subject_name, year.cuts) 
coordinates(t) <-~x+y

# build MCP's by subject_name and year.cut
#splt.by <- c('subject_name','year.cuts')
split <- split(t, list(t$subject_name, t$year.cuts), drop = TRUE)
mcp <- lapply(split, mcp)
split.area <- lapply(split, cbind, year.mcp = mcp$area)

split.area <- list()
for(i in 1:length(split)){
  split.area[[i]] <- cbind(split[[i]], newcol = mcp[[i]]$area)
  colnames(split.area[[i]]@data) <- c("subject_name", "year.cuts", "year.mcp.area")
}
mcp.areas <- do.call(rbind, split.area)
gme$year.mcp.area <- mcp.areas@data$year.mcp.area


##### cut raster stack by 
# mask landcover raster to homeranges
system.time(
  mcp.rast <- lapply(mcp, mask_poly_raster, raster = stack$change03_181_reclassMara_2019.11.22)
)

#' ### Landscape Metrics
##### Landscape Metrics #####

# calculate a metric for all individual homeranges
metrics <- lapply(mcp.rast, lsm_l_contag)

# add relevant metadata to each list element
t <- gme %>% group_by(subject_name, ag.class.both) %>% tally() 
l <- Map(cbind, metrics, subject_name = names(metrics), ag.class = t$ag.class.both)
# create metrics dataframe
metrics.df <- do.call("rbind", l)

#' Plots of homeranges with contagion metric results
#' Yellow = Ag, Green = Natural Habitat
par(mfrow = c(2,2))
cuts = c(0,1,2,3,4)
pal = colorRampPalette(c('white', "yellow", "dark green", "dark grey", "grey"))
plot(mcp.rast$Alina, main = paste0("Sporadic: Alina, cont=", round(metrics$Alina$value, 2)),
     col = pal(5))
plot(mcp.rast$Fred, main = paste0("Mixed Seasonal: Fred, cont=", round(metrics$Fred$value,2)),
     col = pal(5))
plot(mcp.rast$Kimbizwa, main = paste0("Seasonal: Kimbizwa, cont=", round(metrics$Kimbizwa$value,2)),
     col = pal(5))
plot(mcp.rast$Olchoda, main = paste0("Habitual: Olchoda, cont=", round(metrics$Fitz$value, 2)),
     col = pal(5))

#' ### ANOVA
#' Looked at 5 metrics. Many had issues meeting the variance assumption due to the small sample size of the 
#' Habitual group (n = 8). Contagion produced the best results, although it may be a poor metric for comparing 
#' landscapes which have differnt proportions of habitat (Wong et al. 2014)
#' 
#' Metrics recommended by Wong et al. 2014 for comparing across landscapes with unequal proportions of habitat 
#' did not produce significant results. 
#' 
#' The way the homerange is defined also likely influences these metrics. MCP is course and for some individuals
#' can include a much larger amount of ag than they are actually accessing (see above plots). Could look at KDE 
#' or LoCoH to improve on this. 

##### ANOVA #####
# Data Summary
summary <- metrics.df %>%
  group_by(ag.class) %>%
  summarise(n = n(),
            contagion.median = median(value),
            contagion.sd = sd(value))
(summary)


# ANOVA
# Compute the analysis of variance
res.aov <- aov(value ~ ag.class, data = metrics.df)
# Summary of the analysis
summary(res.aov)

# Tukey
tukey <- TukeyHSD(res.aov)
(tukey)

#' #### Check assumptions
#' 1. Homogeneity of variances
#' Small number of individuals in habitual ag class create some issues with variance assumption
#' Plot of residuals shows megaphone shape - not good! 
library(car)
leveneTest(value ~ ag.class, data = metrics.df)

plot(res.aov, 1)

#' 2. Check normality of metric distribution among the tactic groups
#' Looks ok
p <- ggplot(metrics.df, aes(x = value, fill = ag.class)) + geom_density(alpha = 0.6)
p
