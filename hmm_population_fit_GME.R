#'-----
#'title: HMM Model Fitting - Individual Level
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


ele.log.step <- readRDS("GMEcollars_001_logStep_train_TEST.rds")
split <- split(ele.log.step, ele.log.step$ID)

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


# null model - log step length
# system.time(
# mod1 <- fitHMM(data = ele.log.velocity, nbStates = 3, dist = distNorm,
#                Par0 = par1,
#                retryFits = 0,
#                stateName = stateNames3,
#                modelName = "dist2ag + dist2forest + gHM + (1:ID)")
# )
# 
# saveRDS(mod1, file = "mod1_L.rds")
# 
# system.time(
# mod2 <- fitHMM(data = ele.log.velocity, nbStates = 3, dist = distNorm,
#                Par0 = par1,
#                retryFits = 0,
#                stateName = stateNames3,
#                formula = ~ dist2ag + dist2forest + gHM + (1:name),
#                modelName = "dist2ag + dist2forest + gHM + pct.ag + pct.forest + (1:ID)")
# )
# 
# saveRDS(mod2, file = "mod2_L.rds")
# 
# system.time(
#   mod3 <- fitHMM(data = ele.log.velocity, nbStates = 3, dist = distNorm,
#                  Par0 = par1,
#                  retryFits = 0,
#                  stateName = stateNames3,
#                  formula = ~ dist2ag + I(dist2ag^2) + dist2forest + gHM + pct.ag + pct.forest + (1:name),
#                  modelName = "dist2ag + dist2ag^2 + dist2forest + gHM + pct.ag + pct.forest + (1:ID)")
# )
# 
# saveRDS(mod3, file = "mod3_L.rds")
# 

system.time(
  mod1 <- fitHMM(data = ele.log.step, nbStates = 3, dist = distNorm,
                 Par0 = par1,
                 retryFits = 3,
                 stateName = stateNames3,
                 modelName = "Null")
)

saveRDS(mod1, file = "mod1_null_log.rds")

system.time(
  mod2 <- fitHMM(data = ele.log.step, nbStates = 3, dist = distNorm,
                 Par0 = par1,
                 retryFits = 0,
                 stateName = stateNames3,
                 formula = ~ dist2ag + I(dist2ag^2) + dist2forest + gHM + slope + (1:Name),
                 modelName = "dist2ag + dist2ag^2 + dist2forest + gHM + slope + (1:ID)")
)

saveRDS(mod2, file = "mod2_full_log.rds")


system.time(
  mod3 <- fitHMM(data = ele.log.step, nbStates = 3, dist = distNorm,
                 Par0 = par1,
                 retryFits = 0,
                 stateName = stateNames3,
                 formula = ~ dist2ag + I(dist2ag^2) + dist2forest + gHM + slope,
                 modelName = "dist2ag + dist2ag^2 + dist2forest + gHM + slope")
)

saveRDS(mod3, file = "mod3_full_log.rds")


# system.time(
#   mod6 <- fitHMM(data = ele.log.velocity, nbStates = 3, dist = distNorm,
#                  Par0 = par1,
#                  retryFits = 0,
#                  stateName = stateNames3,
#                  formula = ~ dist2ag + I(dist2ag^2) + dist2forest + gHM + pct.forest + (1:name),
#                  modelName = "dist2ag + dist2ag^2 + dist2forest + gHM + pct.forest + (1:ID)")
# )
# 
# saveRDS(mod6, file = "mod6_L.rds")
