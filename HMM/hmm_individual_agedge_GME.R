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

# library(doParallel)
# cl <- makeCluster(3)
# registerDoParallel(cl)
# 
# 
# library(foreach)
# system.time({
#   m1.indiv <- foreach::foreach(i = 1:3) %dopar% 
#     momentuHMM::fitHMM(data = indiv.log.velocity[[i]], nbStates = 3, dist = distNorm,
#                        Par0 = par1,
#                        #retryFits = 4,
#                        stateNames = stateNames3,
#                        formula = ~ dist2ag + I(dist2ag^2) + dist2forest + dist2water + gHM,
#                        modelName = "dist2ag + dist2ag^2 + dist2forest + dist2water + gHM")
#   
# })
# 
# saveRDS(m1.indiv, "m1_individual.rds")
# 
# # get viterbi estimates and extract data frames
# temp <- list()
# 
# for (i in 1:length(m1.indiv)) {
#   m.indiv[[i]]$data$viterbi <- viterbi(m.indiv[[i]])
#   temp[[i]] <- m.indiv[[i]]$data
# }
# 
# result.mfit <- do.call(rbind, temp)
# 
# # save dataframes  
# saveRDS(result.mfit, "df_mfit_indiv.rds")



library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)


# library(foreach)
# system.time({
#   m2.indiv <- foreach::foreach(i = 1:3) %dopar% 
#     momentuHMM::fitHMM(data = indiv.log.velocity[[i]], nbStates = 3, dist = distNorm,
#                        Par0 = par1,
#                        #retryFits = 4,
#                        stateNames = stateNames3,
#                        formula = ~ dist2agedge + I(dist2agedge^2) + dist2forest + dist2water + gHM,
#                        modelName = "dist2agedge + dist2agedge^2 + dist2forest + dist2water + gHM")
#   
# })
# 
# saveRDS(m2.indiv, "m2_individual.rds")
# 
# 
# 
# library(foreach)
# system.time({
#   m4.indiv <- foreach::foreach(i = 1:3) %dopar% 
#     momentuHMM::fitHMM(data = indiv.log.velocity[[i]], nbStates = 3, dist = distNorm,
#                        Par0 = par1,
#                        #retryFits = 4,
#                        stateNames = stateNames3,
#                        formula = ~ dist2agedge + dist2forest + dist2water + gHM,
#                        modelName = "dist2agedge + dist2forest + dist2water + gHM")
#   
# })
# 
# saveRDS(m4.indiv, "m4_individual.rds")



library(foreach)
system.time({
  m6.indiv <- foreach::foreach(i = 1:3) %dopar% 
    momentuHMM::fitHMM(data = indiv.log.velocity[[i]], nbStates = 3, dist = distNorm,
                       Par0 = par1,
                       #retryFits = 4,
                       stateNames = stateNames3,
                       formula = ~ dist2agedge + dist2forest + dist2permwater + dist2seasonalwater + gHM,
                       modelName = "dist2agedge + dist2forest + dist2permwater + dist2seasonalwater + gHM")
  
})

saveRDS(m6.indiv, "m6_individual.rds")



library(foreach)
system.time({
  m8.indiv <- foreach::foreach(i = 1:3) %dopar% 
    momentuHMM::fitHMM(data = indiv.log.velocity[[i]], nbStates = 3, dist = distNorm,
                       Par0 = par1,
                       #retryFits = 4,
                       stateNames = stateNames3,
                       formula = ~ dist2agedge + I(dist2agedge^2) + dist2forest + dist2permwater + dist2seasonalwater + gHM,
                       modelName = "dist2agedge + dist2agedge^2 + dist2forest + dist2permwater + dist2seasonalwater + gHM")
  
})

saveRDS(m8.indiv, "m8_individual.rds")



library(foreach)
system.time({
  m10.indiv <- foreach::foreach(i = 1:3) %dopar% 
    momentuHMM::fitHMM(data = indiv.log.velocity[[i]], nbStates = 3, dist = distNorm,
                       Par0 = par1,
                       #retryFits = 4,
                       stateNames = stateNames3,
                       formula = ~ dist2agedge + I(dist2agedge^2) + dist2forest + dist2permwater + gHM,
                       modelName = "dist2agedge + dist2agedge^2 + dist2forest + dist2permwater + gHM")
  
})

saveRDS(m10.indiv, "m10_individual.rds")
