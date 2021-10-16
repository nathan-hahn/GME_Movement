#'-----
#'title: HMM Model Fitting - Population Level
#'author: Nathan Hahn
#'-----
#'
#'This code is for implementation of HMM pop-level models across multiple machines. momentuHMM object prepped in the HMM model pre_pop.R file. Import here and proceed with model fitting. 
#'
#'
#' libraries
library(momentuHMM)
library(forcats)

# Set environment to EAT
Sys.setenv(TZ="Africa/Nairobi") 


pop.log.veloc <- readRDS("./HMM/GMEcollars_004_population_logVeloc.rds")

#' ---- Specify parameters ---- 

# label states
stateNames3 <- c("encamped", "exploritory", "goal-oriented movement")
stateNames2 <- c("encamped", "exploratory")

# set distributions
# distGam = list(step = "gamma", angle = "vm") # zero-inflated gamma distributions
distNorm = list(step = "norm", angle = "vm")


# ****top model - step****
par1 <- list(step = c(1, 4, 6, 2, 1, 1),
             angle = c(.2,.2,.2)) 

# par2 <- list(step = c(1, 4, 6, 3, 2, 1),
#              angle = c(.1,.3,.7)) 


#' ---- Fit models ----
#' Each model runs, reports sys.time, and saves 


##### Ag Edge #####


#' ---- Specify parameters ---- 

# label states
stateNames3 <- c("Encamped", "Area-Restricted Search", "Exploratory")
stateNames2 <- c("Encamped", "Exploratory")

# set distributions
# distGam = list(step = "gamma", angle = "vm") # zero-inflated gamma distributions
distNorm = list(step = "norm", angle = "vm")


# ****top model - step****
par1 <- list(step = c(1, 4, 6, 2, 1, 1),
             angle = c(.2,.2,.2)) 

# par2 <- list(step = c(1, 4, 6, 3, 2, 1),
#              angle = c(.1,.3,.7)) 


#' ---- Fit models ----
#' Each model runs, reports sys.time, and saves 

m2.pop <- momentuHMM::fitHMM(data = pop.log.veloc, nbStates = 3, dist = distNorm,
                       Par0 = par1,
                       #retryFits = 4,
                       stateNames = stateNames3,
                       formula = ~ dist2agedge + I(dist2agedge^2) + dist2forest + dist2water + gHM + slope,
                       modelName = "dist2agedge + dist2agedge^2 + dist2forest + dist2water + gHM + slope")
  


saveRDS(m2.pop, "m2_pop.rds")

m3.pop <- momentuHMM::fitHMM(data = pop.log.veloc, nbStates = 3, dist = distNorm,
                               Par0 = par1,
                               #retryFits = 4,
                               stateNames = stateNames3,
                               formula = ~ dist2ag + I(dist2ag^2) + dist2forest + dist2permwater + dist2seasonalwater + gHM + slope,
                               modelName = "dist2ag + dist2ag^2 + dist2forest + dist2water + gHM + slope")



saveRDS(m2.pop, "m3_pop.rds")


m4.pop <- momentuHMM::fitHMM(data = pop.log.veloc, nbStates = 3, dist = distNorm,
                       Par0 = par1,
                       #retryFits = 4,
                       stateNames = stateNames3,
                       formula = ~ dist2agedge + dist2forest + dist2water + gHM + slope,
                       modelName = "dist2agedge + dist2forest + dist2water + gHM + slope")

saveRDS(m4.pop, "m4_pop.rds")


m5.pop <- momentuHMM::fitHMM(data = pop.log.veloc, nbStates = 3, dist = distNorm,
                               Par0 = par1,
                               #retryFits = 4,
                               stateNames = stateNames3,
                               formula = ~ dist2ag + dist2forest + dist2permwater + dist2seasonalwater + gHM + slope,
                               modelName = "dist2ag + dist2ag^2 + dist2permwater + dist2seasonalwater + gHM + slope")



saveRDS(m5.pop, "m5_pop.rds")


m6.pop <- momentuHMM::fitHMM(data = pop.log.veloc, nbStates = 3, dist = distNorm,
                       Par0 = par1,
                       #retryFits = 4,
                       stateNames = stateNames3,
                       formula = ~ dist2agedge + dist2forest + dist2permwater + dist2seasonalwater + gHM + slope,
                       modelName = "dist2agedge + dist2forest + dist2permwater + dist2seasonalwater + gHM + slope")

saveRDS(m6.pop, "m6_pop.rds")



 



m10.pop <- momentuHMM::fitHMM(data = pop.log.veloc, nbStates = 3, dist = distNorm,
                       Par0 = par1,
                       #retryFits = 4,
                       stateNames = stateNames3,
                       formula = ~ dist2agedge + dist2forest + dist2permwater + gHM + slope,
                       modelName = "dist2agedge +  dist2forest + dist2permwater + gHM + slope")

saveRDS(m10.pop, "m10_pop.rds")


m12.pop <- momentuHMM::fitHMM(data = pop.log.veloc, nbStates = 3, dist = distNorm,
                                Par0 = par1,
                                #retryFits = 4,
                                stateNames = stateNames3,
                                formula = ~ dist2agedge + dist2forest + dist2permwater + dist2seasonalwater + slope,
                                modelName = "dist2agedge +  dist2forest + dist2permwater + dist2seasonalwater + slope")

saveRDS(m12.pop, "m12_pop.rds")



m14.pop <- momentuHMM::fitHMM(data = pop.log.veloc, nbStates = 3, dist = distNorm,
                                Par0 = par1,
                                #retryFits = 4,
                                stateNames = stateNames3,
                                formula = ~ dist2agedge + dist2permwater + dist2seasonalwater + slope,
                                modelName = "dist2agedge +  dist2permwater + dist2seasonalwater + slope")

saveRDS(m14.pop, "m14_pop.rds")



m16.pop <- momentuHMM::fitHMM(data = pop.log.veloc, nbStates = 3, dist = distNorm,
                                Par0 = par1,
                                #retryFits = 4,
                                stateNames = stateNames3,
                                formula = ~ dist2agedge + dist2forest + slope,
                                modelName = "dist2agedge +  dist2forest + slope")

saveRDS(m16.pop, "m16_pop.rds")

