
library(tidyverse)
library(momentuHMM)

source("./GME_functions.R")

#' ---- Prep data for HMM population model fitting ----

# Add data with raster covariates. 
used.all <- readRDS("./GMM/GMEcollars_002_usedClust_2020-06-05.rds")
used.all <- used.all[!is.na(used.all$x),]
used.all <- used.all %>%
  filter(fixType != "irregular")

#' ---- Create Test/Train Sets ----
#' Split the data into testing and training sets
#' Last n% of each burst for the testing, defined by 'cut'


# filter to bursts that are large enough for test/train
t <- used.all %>%
  group_by(burst, subject_name, ag.class.both) %>% tally() %>%
  mutate(train = n*.90, test = n*.10) %>%
  filter(train > 500 & test > 500) %>% droplevels()

used.all <- filter(used.all, burst %in% t$burst) %>% droplevels()

# divide datasets by individual
split <- split(used.all, used.all$burst)
split <- c(split[1], split[2], split[6], split[14], split[23]) #TEMP -- only use Alina and Kegol

# withold data
cut = .10
train <- lapply(split, withold, cut = cut, type = "train")
train <- do.call("rbind", train)
test <- lapply(split, withold, cut = cut, type = "test")
test <- do.call("rbind", test)

# save dataframe to be able to revert to unstandadized and non-log distances
saveRDS(train, "./HMM/TEST_traindf_original.rds")


#' ---- Standardize Covariates ----
# standardize covariates
cor.vars <- c("dist2ag", "dist2forest", "dist2water", "slope")
st.covs <- train %>%
  dplyr::select(cor.vars)%>%
  apply(., 2, function(x) (x - mean(x)) / sd(x)) %>%
  as.data.frame()

# check dimensions
dim(st.covs)

# merge with reloc data
train.st <- train %>%
  dplyr::select(-cor.vars) %>%
  bind_cols(., st.covs) %>%
  rename(ID = burst) %>% #for HMM package
  droplevels()
# must be a dataframe to work!! 
train.st <- as.data.frame(train.st)


#' ---- create momentu objects ----
#' Apply the velocity function by group. Create list by burst, apply function to each data.frame, and then unlist.

# Create objects for the MomentuHMM package: step, log step, velocity, log velocity 
library(momentuHMM)

ele.step <- prepData(data = train.st, coordNames = c("x", "y"))
split <- split(ele.step, ele.step$ID)

individual.velocity <- lapply(split, log_velocity)
population.velocity <- do.call("rbind", individual.velocity)

# save as rdata for us on secondary machines
saveRDS(individual.velocity, "./HMM/GMEcollars_002_individual_logVeloc_train.rds")
saveRDS(population.velocity, "./HMM/GMEcollars_002_population_logVeloc_train.rds")
