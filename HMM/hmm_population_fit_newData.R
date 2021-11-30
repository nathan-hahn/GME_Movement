#'-----
#'title: HMM Model Fitting - Population Level
#'author: Nathan Hahn
#'-----

# This code is for implementation of HMM pop-level models across multiple machines. momentuHMM object prepped in the HMM model pre_pop.R file. Import here and proceed with model fitting. 

# set wd - for server run
setwd('~/Dropbox/CSU/GME_Movement') 

# libraries
library(momentuHMM)
library(forcats)

# Set environment to EAT
Sys.setenv(TZ="Africa/Nairobi") 

pop.log.veloc <- readRDS("./HMM/GMEcollars_004_population_logVeloc.rds")
pop.log.veloc$subject_name <- as.factor(pop.log.veloc$subject_name)

#' ---- Specify parameters ---- 

# label states
stateNames3 <- c("Encamped", "Area-Restricted Search", "Exploratory")
stateNames2 <- c("Encamped", "Exploratory")

# set distributions
distNorm = list(step = "norm", angle = "vm")

# ****top model - step****
par1 <- list(step = c(1, 4, 6, 2, 1, 1),
             angle = c(.2,.2,.2)) 

#' ---- Fit models ----
#' Each model runs, reports sys.time, and saves 

m.pop <- momentuHMM::fitHMM(data = pop.log.veloc, nbStates = 3, dist = distNorm,
                                Par0 = par1,
                                #retryFits = 3,
                                stateNames = stateNames3,
                                formula = ~ dist2agedge + dist2agedge^2 + dist2forest + dist2permwater + dist2seasonalwater + gHM + slope + (1:subject_name),
                                modelName = "dist2agedge + dist2agedge^2 + dist2forest + dist2permwater + dist2seasonalwater + gHM + slope")

saveRDS(m.pop, "m_pop_GME_004.rds")


