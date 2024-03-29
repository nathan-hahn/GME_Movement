
library(tidyverse)
library(momentuHMM)
library(lubridate)

source("./GME_functions.R")

#' ---- Prep data for HMM population model fitting ----

# Add data with raster covariates. 
used.all <- readRDS("./movdata/GMEcollars_004_usedClust_2021-10-05.rds")
#used.all <- used.all[!is.na(used.all$x),]
used.all <- used.all %>%
  filter(fixType != "irregular") %>%
  mutate(date = ymd_hms(date, tz = 'Africa/Nairobi')) %>%
  arrange(burst, date)

#' ---- Create Test/Train Sets ----
#' Split the data into testing and training sets
#' Last n% of each burst for the testing, defined by 'cut'


# filter to bursts that are large enough for test/train
t <- used.all %>%
  group_by(burst, subject_name) %>% tally() %>%
  filter(n > 500)
  
  # mutate(train = n*.90, test = n*.10) %>%
  # filter(train > 500 & test > 500) %>% droplevels()

used.all <- filter(used.all, burst %in% t$burst) %>% droplevels()

# divide datasets by individual
# split <- split(used.all, used.all$burst)
# split <- c(split[1], split[15], split[6]) #TEMP -- only use 3 elephants
# 
# # withold data
# cut = .10
# train <- lapply(split, withold, cut = cut, type = "train")
# train <- do.call("rbind", train)
# test <- lapply(split, withold, cut = cut, type = "test")
# test <- do.call("rbind", test)
# 
# # save dataframe to be able to revert to unstandadized and non-log distances
# saveRDS(train, "./HMM/TEST_traindf_original.rds")
# saveRDS(test, "./HMM/TEST_testdf_original.rds")


#' ---- Standardize Covariates ----

##### Adjust ag edge #####
used.all$dist2agedge <- ifelse(used.all$lc.estes == 1, -(used.all$dist2agedge), used.all$dist2agedge)

# standardize covariates
cor.vars <- c("dist2ag", "dist2agedge", "dist2forest", "dist2water", "dist2permwater","dist2seasonalwater", "slope")
st.covs <- used.all %>%
  dplyr::select(cor.vars)%>%
  apply(., 2, function(x) (x - mean(x)) / sd(x)) %>%
  as.data.frame()

# normalize covariates - check
nm.covs <- used.all %>%
  dplyr::select(cor.vars) %>%
  apply(., 2, function(x) (x - min(x)) / (max(x) - min(x))) %>%
  as.data.frame()

# check dimensions
dim(st.covs)

# merge with reloc data
used.st <- used.all %>%
  dplyr::select(-cor.vars) %>%
  bind_cols(., st.covs) %>%
  rename(ID = burst) %>% #for HMM package
  dplyr::select(ID, x, y, date, all_of(cor.vars), gHM, subject_name) %>%
  droplevels()
# must be a dataframe to work!! 
used.st <- as.data.frame(used.st)


#data_ordered <- used.st[with(used.st, order(ID, date)),]

#' ---- create momentu objects ----
#' Apply the velocity function by group. Create list by burst, apply function to each data.frame, and then unlist.

# Create objects for the MomentuHMM package: step, log step, velocity, log velocity 
library(momentuHMM)
library(rlist)

ele.step <- prepData(data = used.st, coordNames = c("x", "y"), covNames = c(cor.vars, 'subject_name', 'gHM'))
split <- split(ele.step, ele.step$ID)

individual.velocity <- lapply(split, log_velocity)
population.velocity <- do.call("rbind", individual.velocity)

# save as rdata for us on secondary machines
saveRDS(individual.velocity, "./HMM/GMEcollars_004_individual_logVeloc.rds")
saveRDS(population.velocity, "./HMM/GMEcollars_004_population_logVeloc.rds")

# save original dataframe to attach after fitting
saveRDS(used.all, "./HMM/GMEcollars_004_population_original.rds")



