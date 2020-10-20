#'-----
#'title: HMM Model Fitting - Individual Level
#'author: Nathan Hahn
#'-----
#'
#'This code is for implementation of HMM individual-level models. 
#'momentuHMM object prepped seperately
#'
#'
#' libraries
library(momentuHMM)
library(forcats)

# Set environment to EAT
Sys.setenv(TZ="Africa/Nairobi") 

indiv.log.velocity <- readRDS("./HMM/GMEcollars_002_individual_logVeloc_train.rds")

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

library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)

 
library(foreach)
system.time({
  m1.indiv <- foreach::foreach(i = 1:3) %dopar%
    momentuHMM::fitHMM(data = indiv.log.velocity[[i]], nbStates = 3, dist = distNorm,
                       Par0 = par1,
                       #retryFits = 4,
                       stateNames = stateNames3,
                       formula = ~ dist2ag + I(dist2ag^2) + dist2forest + dist2water + gHM + slope,
                       modelName = "dist2ag + dist2ag^2 + dist2forest + dist2water + gHM + slope")

})

saveRDS(m1.indiv, "m1_individual.rds")



library(foreach)
system.time({
  m3.indiv <- foreach::foreach(i = 1:3) %dopar%
    momentuHMM::fitHMM(data = indiv.log.velocity[[i]], nbStates = 3, dist = distNorm,
                       Par0 = par1,
                       #retryFits = 4,
                       stateNames = stateNames3,
                       formula = ~ dist2ag + I(dist2ag^2) + dist2forest + dist2permwater + dist2seasonalwater + gHM + slope,
                       modelName = "dist2ag + dist2ag^2 + dist2forest + dist2permwater + dist2seasonalwater + gHM + slope")

})

saveRDS(m3.indiv, "m3_individual.rds")


library(foreach)
system.time({
  m5.indiv <- foreach::foreach(i = 1:3) %dopar%
    momentuHMM::fitHMM(data = indiv.log.velocity[[i]], nbStates = 3, dist = distNorm,
                       Par0 = par1,
                       #retryFits = 4,
                       stateNames = stateNames3,
                       formula = ~ dist2ag + I(dist2ag^2) + dist2forest + dist2permwater + dist2seasonalwater + gHM + slope,
                       modelName = "dist2ag + dist2forest + dist2permwater + dist2seasonalwater + gHM + slope")

})

saveRDS(m5.indiv, "m5_individual.rds")


library(foreach)
system.time({
  m7.indiv <- foreach::foreach(i = 1:3) %dopar%
    momentuHMM::fitHMM(data = indiv.log.velocity[[i]], nbStates = 3, dist = distNorm,
                       Par0 = par1,
                       #retryFits = 4,
                       stateNames = stateNames3,
                       formula = ~ dist2ag + I(dist2ag^2) + dist2forest + dist2permwater + gHM + slope,
                       modelName = "dist2ag + dist2ag^2 + dist2dist2forest + dist2permwater + gHM + slope")
  
})

saveRDS(m7.indiv, "m7_individual.rds")



library(foreach)
system.time({
  m9.indiv <- foreach::foreach(i = 1:3) %dopar%
    momentuHMM::fitHMM(data = indiv.log.velocity[[i]], nbStates = 3, dist = distNorm,
                       Par0 = par1,
                       #retryFits = 4,
                       stateNames = stateNames3,
                       formula = ~ dist2ag + dist2forest + dist2water + gHM + slope,
                       modelName = "dist2ag + dist2forest + dist2water + gHM + slope")

})

saveRDS(m9.indiv, "m9_individual.rds")


