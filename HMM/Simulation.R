##### Simulation #####

#' Window Size = 12 hours in ag
#' Window rules:
#' 1. Any point in ag creates a +/- 12 hour window around the point
#' 2. Overlapping windows condensed into a single window

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
library(DescTools)
library(adehabitatLT)

source('GME_functions.R')

##### Prep Data #####
hmm <- readRDS('m8_pop_003.rds') # output from hmm model (global)

hmm.df <- hmm$data
hmm.df$states <- viterbi(hmm)

##### Simulation #####
# simulate data with 10,000
sim <- simData(model = hmm, obsPerAnimal = 500000, states = TRUE)


## Model ability to reproduce overall movement process (step and turn distributions)
sim.angle <- sim$angle
sim.step <- sim$step

obs.angle <- hmm.df$angle
obs.step <- hmm.df$step

# angle densities
ggplot() + geom_density(aes(obs.angle), alpha = 0.6) + geom_density(aes(sim.angle), color = "red", alpha = 0.6) +
  xlab("turning angle")

# step densities
ggplot() + geom_density(aes(obs.step), alpha = 0.6) + geom_density(aes(sim.step), color = "red", alpha = 0.6) +
  xlab("log(step)")

## Model ability to reproduce movement process at state level

# slightly overestimates exploratory movement and underestimates encamped. This could be a result of simulated landscape covariates though
prop.table(table(sim$states))
prop.table(table(hmm.df$states))


# plots 
# combine state-level datasets

sim.step <- as.data.frame(cbind(sim$step, sim$states, "sim"))
colnames(sim.step) <- c('step', 'states', 'key')
sim.angle <- as.data.frame(cbind(sim$angle, sim$states, "sim"))
colnames(sim.angle) <- c('angle', 'states', 'key')

obs.step <- as.data.frame(cbind(hmm.df$step, hmm.df$states, "obs"))
colnames(obs.step) <- c('step', 'states', 'key')
obs.angle <- as.data.frame(cbind(hmm.df$angle, hmm.df$states, "obs"))
colnames(obs.angle) <- c('angle', 'states', 'key')

# dataframes for plotting
step <- rbind(obs.step, sim.step) %>%
  mutate(step = as.numeric(as.character(step))) %>%
  filter(!is.na(step)) %>%
  mutate(states = recode_factor(states, "1" = "encamped",
                                "2" = "meandering",
                                "3" = "directed walk"))

angle <- rbind(obs.angle, sim.angle) %>%
  mutate(angle = as.numeric(as.character(angle))) %>%
  filter(!is.na(angle)) %>%
  mutate(states = recode_factor(states, "1" = "encamped",
                                "2" = "meandering",
                                "3" = "directed walk"))


# plot steps
ggplot(data = step, aes(x = step, color = key)) + geom_density() + facet_wrap(.~states) +
  xlab("log(step)") + scale_color_manual(values = c("black", "red"))

# plot angles
ggplot(data = angle, aes(x = angle, color = key)) + geom_density() + facet_wrap(.~states) +
  xlab("turning angle") + scale_color_manual(values = c("black", "red"))


# K-S tests

split.angle <- split(angle, list(angle$key, angle$states))

ks.test(split.angle[[1]]$angle, split.angle[[2]]$angle)
chisq.test(t.test[[1]]$angle, split.angle[[2]]$angle)

ks.test(split.angle[[3]]$angle, split.angle[[4]]$angle)

ks.test(split.angle[[5]]$angle, split.angle[[6]]$angle)









