#### HMM Model Assessments ####
##### Activity Moving Windows #####

library(tidyverse)
theme_set(  theme_bw() + # set theme with no legend of strip text
              theme(panel.grid.major = element_blank(),
                    strip.background = element_blank(),
                    panel.border = element_rect(colour = "black"),
                    strip.text = element_text(size = 12),
                    legend.text = element_text(size = 10),
                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title=element_text(size=14), axis.text = element_text(size=12))
)
library(lubridate)
library(xts)
library(zoo)
library(momentuHMM)

source('GME_functions.R')

##### Prep Data #####
hmm.df <- readRDS('m8_pop_003.rds') # output from hmm model (global)
output <- readRDS("./HMM/GMEcollars_003_population_original.rds") # dataset prior to data transforms

# add viterbi estimates to the original dataframe -- revert to non-standardized covs and steps
output$viterbi <- viterbi(hmm.df)
output$v.dist <- hmm.df$data$step
output$v.angle <- hmm.df$data$angle


# add stationary state props
t <- stationary(hmm.df)
df <- as.data.frame(t[[1]])
output <- cbind(output, df)


##### Rose Diagrams #####
library(CircStats)
plotRA <- function(x) {
  rose.diag(x[!is.na(x$rel.angle),]$rel.angle, bins=24, prop=1.8,
            main=paste0("Relative Angles: ", unique(x$state)))
}

{df <- output %>% mutate(state = recode_factor(viterbi,
                                               "1" = "Encamped", "2" = "R-A Search", "3" = "Exploratory"))
  state.split <- split(df, df$viterbi)
  lapply(state.split, plotRA)}

##### State classification scatter plot #####
# need to figure out filtering for bad state classifications
# 1, > 300; 2, > 700

split <- split(output, output$subject_name)
plot <- list()
for(i in 1:length(split)){
  
  # sample 500 random points
  sample <- split[[i]][sample(nrow(split[[i]]), 500), ]
  sample <- filter(sample, exp(v.dist) >300 && viterbi != 1)
  #sample <- filter(sample, exp(v.dist) > 800 && viterbi != 2)
  
  plot[[i]] <- ggplot(sample, aes(exp(v.dist), v.angle)) + geom_point(aes(color = factor(viterbi)))
}

# Simulation