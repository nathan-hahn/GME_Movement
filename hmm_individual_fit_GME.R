#'-----
#'title: HMM Model Fitting
#'author: Nathan Hahn
#'-----
#'
#' libraries
library(momentuHMM)
library(plyr)
library(ggplot2)
library(dplyr)

source("~/Dropbox (Personal)/R/Functions/velocity.R")

# Set environment to EAT
Sys.setenv(TZ="Africa/Nairobi") 

#' ---- Prep data for model fitting ----
ele.log.step <- readRDS("./HMM/GMEcollars_001_logStep_train_TEST.rds")
split <- split(ele.log.step, ele.log.step$ID)
# log velocity - replaces values in the step variable
ele.log.velocity <- lapply(split, log.velocity)

#' ---- Specify parameters ---- 

# label states
stateNames3 <- c("encamped", "foraging", "exploritory")
stateNames2 <- c("encamped", "exploratory")

# set distributions
distGam = list(step = "gamma", angle = "vm") # zero-inflated gamma distributions
distNorm = list(step = "norm", angle = "vm")


# ****top model - step****
par1 <- list(step = c(1, 4, 6, 2, 1, 1),
                angle = c(.1,.3,.7)) 



#' ---- Fit model ----
#' Fit top model at the individual level
 

# for (i in 1:length(ele.log.velocity)) {
# mfits1 [[i]] <- fitHMM(data = ele.log.velocity[[i]], nbStates = 3, dist = distNorm,
#                        Par0 = par1,
#                        retryFits = 1,
#                        stateNames = stateNames3,
#                        formula = ~ dist2forest + I(dist2forest^2) + dist2ag + I(dist2ag^2) + pa + ToD)
# }

# use the foreach function
library(doParallel)
cl <- makeCluster(6)
registerDoParallel(cl)


library(foreach)
system.time({
  m.indiv <- foreach::foreach(i = 1:65) %dopar% 
    momentuHMM::fitHMM(data = ele.log.velocity[[i]], nbStates = 3, dist = distNorm,
                       Par0 = par1,
                       retryFits = 4,
                       stateNames = stateNames3,
                       formula = ~ dist2ag + I(dist2ag^2) + dist2forest + gHM + pct.forest + (1:name),
                       modelName = "dist2ag + dist2ag^2 + dist2forest + gHM + pct.forest + (1:ID)")
  
})

saveRDS(m.indiv, "m_individual.rds")

# get viterbi estimates and extract data frames
temp <- list()

for (i in 1:length(m.indiv)) {
  m.indiv[[i]]$data$viterbi <- viterbi(m.indiv[[i]])
  temp[[i]] <- m.indiv[[i]]$data
}

result.mfit <- do.call(rbind, temp)


# save dataframes  
saveRDS(result.mfit, "df_mfit_indiv.rds")

