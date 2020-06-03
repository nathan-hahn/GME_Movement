
library(dplyr)
library(momentuHMM)

source("~/Dropbox (Personal)/R/Functions/velocity.R")

#' ---- Prep data for HMM population model fitting ----

# Add data with raster covariates. 
used.all <- readRDS("./GMM/GMEcollars_002_usedClust_2020-06-02.rds")
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
  dplyr::select(-all_of(cor.vars)) %>%
  bind_cols(., st.covs) %>%
  rename(ID = burst) %>% #for HMM package
  droplevels()


# must be a dataframe to work!! 
used.st <- as.data.frame(used.st)
saveRDS(used.st, "./HMM/GMEcollars_001_usedClustST_2020-05-14.rds")


#' ---- Create Test/Train Sets ----
#' Split the data into testing and training sets
#' Last n% of each burst for the testing, defined by 'cut'


# filter to bursts that are large enough for test/train
t <- used.st %>%
  group_by(ID, subject_name) %>% tally() %>%
  mutate(train = n*.90, test = n*.10) %>%
  filter(train > 500 & test > 500) %>% droplevels()

used.st <- filter(used.st, ID %in% t$ID) %>% droplevels()

# divide datasets by individual
split <- split(used.st, used.st$ID)
split <- c(split[1], split[6]) #TEMP -- only use Alina and Kegol

withold <- function(x, cut, type) {
  n <- nrow(x)*cut
  
  if (type == "train") {
    y <- head(x, (nrow(x)-n))
  }
  
  if (type == "test") {
    y <- tail(x, n)
  }
  
  return(y)
}

cut = .10
train <- lapply(split, withold, cut = cut, type = "train")
train <- do.call("rbind", train)
test <- lapply(split, withold, cut = cut, type = "test")
test <- do.call("rbind", test)

#' ---- create momentu objects ----
#' Apply the velocity function by group. Create list by burst, apply function to each data.frame, and then unlist.

# Create objects for the MomentuHMM package: step, log step, velocity, log velocity 
library(momentuHMM)

ele.step <- prepData(data = train, coordNames = c("x", "y"))
split <- split(ele.step, ele.step$ID)

individual.velocity <- lapply(split, log.velocity)
population.velocity <- do.call("rbind", individual.velocity)

# save as rdata for us on secondary machines
saveRDS(individual.velocity, "./HMM/GMEcollars_002_individual_logVeloc_train.rds")
saveRDS(population.velocity, "./HMM/GMEcollars_002_population_logVeloc_train.rds")
