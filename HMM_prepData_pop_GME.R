
library(dplyr)
library(momentuHMM)

source(velocity.R)

#' ---- Prep data for HMM population model fitting ----

# Add data with raster covariates. 
used.all <- readRDS("./movdata/GMEcollars_001_usedClust_2020-01-14.rds")
used.all <- used.all[!is.na(used.all$x),]

# standardize covariates
cor.vars <- c("dist2ag", "dist2forest", "dist2water", "slope")
st.covs <- used.all %>%
  dplyr::select(all_of(cor.vars))%>%
  apply(., 2, function(x) (x - mean(x)) / sd(x)) %>%
  as.data.frame()

# check dimensions
dim(st.covs)

# merge with reloc data
used.st <- used.all %>%
  dplyr::select(-cor.vars) %>%
  bind_cols(., st.covs) %>%
  rename(ID = burst) %>% #for HMM package
  droplevels()


# must be a dataframe to work!! 
used.st <- as.data.frame(used.st)
saveRDS(used.st, "./movdata/GMEcollars_001_usedClustST_2020-05-14.rds")

#' ---- create momentu objects ----
#' Apply the velocity function by group. Create list by burst, apply function to each data.frame, and then unlist.
#' Faster way to do this without losing momentuHMM object? 

# Create objects for the MomentuHMM package: step, log step, velocity, log velocity 
library(momentuHMM)

ele.step <- prepData(data = used.st, coordNames = c("x", "y"))
ele.log.step <- ele.step
ele.log.step$step <- log(ele.log.step$step + 0.001) # add constant for zero steps

# save as rdata for us on secondary machines
saveRDS(ele.log.step, "ele_logStep_mmt.rds")


