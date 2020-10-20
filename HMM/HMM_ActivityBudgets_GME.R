##### Activity Moving Windows #####

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

source('GME_functions.R')

##### Prep Data #####
hmm.df <- readRDS("./HMM/model results/Full Pop July2020/m8_pop_Julyfinal.rds") # output from hmm model (no squared term)
output <- readRDS("./HMM/model results/Full Pop July2020/GMEcollars_002_population_original.rds") # dataset prior to data transforms

# add viterbi estimates to the original dataframe -- revert to non-standardized covs and steps
output$viterbi <- viterbi(hmm.df)


# add stationary state props
t <- stationary(hmm.df)
df <- as.data.frame(t[[1]])
output <- cbind(output, df)


# create new varibles to differentiate ag and non-ag in non-protected areas
output$ag.used <- ifelse(output$lc.estes == 1, 1, 0)
output$pa <- factor(output$pa, levels = c(3:1))
# 4 is ag in non-protected areas
output$pa.2 <- ifelse(output$pa == 1 & output$ag.used == 1, 4, output$pa)

output <- output %>%
  dplyr::select(-c(mean, stdev, hi.95, lo.95, ag.window))

###### Identify Ag Use Windows #####

# Apply rollstat. Create time series and apply rollstat function
window <- 25 # in hours. Split to before and after when align = center
split <- split(output, output$burst) # careful of what you split by. May cuase problems splitting by name if there are large gaps between recollars. Let na.rm = FALSE in rollstat

system.time({
  ts <- lapply(split, as.ts, start = split[[i]]$date, frequency = 1) # frequency set with 30min fixes
  tx <- lapply(ts, as.xts)
  res <- lapply(seq_along(tx), function(i) rollapply(tx[[i]]$ag.used, width = window, align = "center", by.column = FALSE,
                                                     FUN = rollstat))
  merge <- lapply(seq_along(res), function(i) cbind(split[[i]], as.data.frame(res[[i]])))
})

output <- do.call("rbind", merge)

write.csv(output, "./HMM/model results/agWindow_temporary_df.csv")

##### Activity Budgets #####

output.plot <- mutate(output, viterbi = recode_factor(viterbi, 
                                                 "1" = "encamped",
                                                 "2" = "foraging",
                                                 "3" = "transit"),
                tactic.season = recode_factor(tactic.season,
                                              "1" = "Rare",
                                              "2" = "Sporadic",
                                              "3" = "Seasonal",
                                              "4" = "Habitual"),
                 pa.2 = recode_factor(pa.2, "1" = "protected", "2" = "limited use", 
                                      "3" = "not protected", "4" = "ag")) %>%
                filter(!(is.na(ag.window))
) 


# t <- ifelse(is.na(output.plot$ag.window) == TRUE, 0, output.plot$ag.window)
# output.plot$ag.window <- as.factor(t)

# check data before plotting -- are there enough points for a density plot?
t <- output.plot %>%
  group_by(burst, ag.window, viterbi) %>%
  tally()
(t)

#output.plot <- filter(output.plot, burst == "Shorty-3016.1" | burst == "ST2010-2767.1")

## Density plots v2
# budget by ag window
agWindow <- plot_budget(t = output.plot, facet = viterbi ~ ag.window, title = "GME: Activity budget by ag window")
agWindow
# budget by land use
# landUse <- plot_budget(t = output.plot, facet = viterbi ~ pa.2, title = "GME: Activity budget by land use type")
# landUse
# budget by ag tactic
tactic <- plot_budget(t = output.plot, facet = viterbi ~ tactic.season, title = "GME: Activity budget by land use type")
tactic





## Budget by tactic and ag window
t0 <- filter(output.plot, ag.window == 0)
t1 <- filter(output.plot, ag.window == 1) # rare individuals do not have enough points

ag.tactic <- ggplot() +
  facet_grid(viterbi~tactic.season) +
  geom_density(data = t0, aes(x = hour(date),
                   fill = factor(viterbi), colour = factor(ag.window)), alpha = 0.3, adjust = 1.5) +
  geom_density(data = t1, aes(x = hour(date),
                              fill = factor(viterbi), colour = factor(ag.window)), alpha = 0.3, adjust = 1.5) +
  
  # add sunrise/sunset
  geom_vline(xintercept=6
             ,color="dark grey", linetype="dashed", size=1) + 
  geom_vline(xintercept=18,
           color="dark grey", linetype="dashed", size=1) +

  # add colors
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73"), name = c("state")) +
  scale_colour_manual(values = c("black", 'red' ), name = c("ag window")) +
  xlab("hour (0-23)") + ggtitle('Activity Budget by Tactic and Ag Use')
ag.tactic

##### Overlap Tests #####

## Testing overlapping package
d1 <- output.plot %>%
  filter(viterbi == "transit" & ag.window == 0) %>%
  filter(tactic.season == 'Seasonal')
#d1 <- density(as.numeric(hour(d1$date))) #check

d0 <- output.plot %>%
  filter(viterbi == "transit" & ag.window == 0) %>%
  filter(tactic.season == 'Habitual')
#d0 <- density(as.numeric(hour(d0$date))) #check

d <- list(as.numeric(hour(d0$date)), as.numeric(hour(d1$date)))

library(overlapping)
t <- overlap(d, plot = T)
boot <- boot.overlap(d, B = 1000)
boot$OVboot_stats

q <- quantile(boot$OVboot_dist, probs = c(.05, .95))
  

#### function - working
overlap.test <- function(df, comp0, comp1, boot.it = 1000) {
  require(tidyverse)
  require(lubridate)
  require(overlapping)
  d0 <- dplyr::right_join(df, as.data.frame(comp0))
  d1 <- dplyr::right_join(df, as.data.frame(comp1))
  d <- list(as.numeric(hour(d0$date)), as.numeric(hour(d1$date)))
  
  # overlap and bootstrap test
  t <- overlap(d, plot = T)
  boot <- boot.overlap(d, B = boot.it)
  boot$OVboot_stats
  
  return(list(d = df, bootstrap = boot))

}


## Tactic comps loop

x <- c(1:4)
tactic.comp <- combn(x, 2)
comp <- list()
OV_results <- list()
# viterbi loop
for(i in 1:3){
  # tactic comp loop
  for(j in 1:6){
    #create comps
    comp0 <- list(viterbi = as.integer(i), tactic.season = tactic.comp[,j][1])
    comp1 <- list(viterbi = as.integer(i), tactic.season = tactic.comp[,j][2])
    
    # run and store overlap output
    t <- overlap.test(df = output, comp0 = comp0, comp1 = comp1)
    name <- paste(as.integer(i), tactic.comp[,j][1], tactic.comp[,j][2], sep='-')
    q <- quantile(t$bootstrap$OVboot_dist, probs = c(0.05, 0.95))
    
    # store state and comp name
    comp[[name]] <- t$bootstrap$OVboot_stats
    comp[[name]]$lwr <- as.numeric(q[1])
    comp[[name]]$upr <- as.numeric(q[2])
    comp[[name]]$state <- as.integer(i)
    comp[[name]]$comp <- paste(tactic.comp[,j][1], tactic.comp[,j][2], sep='-')
    
    # store OV bootstrap results
    OV_results[[j]] <- t$bootstrap$OVboot_dist
  }
}

# comp output should be 3*6 - every 6 interations is a viterbi group
state.comps <- do.call(rbind, comp)


## Ag Window comps look
agWindow <- c(0,1)
x <- c(1:4)
tactic.comp <- combn(x, 2)
comp <- list()
#comp <- vector(mode = "list", length = 2*3*6)
# state loop
for(i in 1:2){
  # viterbi loop
  for(j in 1:3){
    # tactic comp loop
    for(k in 1:6){
      #create comps
      comp0 <- list(ag.window = agWindow[i], viterbi = as.integer(j), tactic.season = tactic.comp[,k][1])
      comp1 <- list(ag.window = agWindow[i], viterbi = as.integer(j), tactic.season = tactic.comp[,k][2])
      
      # run and store overlap output
      t <- overlap.test(df = output, comp0 = comp0, comp1 = comp1)
      name <- paste(as.integer(j), tactic.comp[,k][1], tactic.comp[,k][2], agWindow[i], sep='-')
      q <- quantile(t$bootstrap$OVboot_dist, probs = c(0.05, 0.95))
      
      # store state and comp name
      comp[[name]] <- t$bootstrap$OVboot_stats
      comp[[name]]$lwr <- as.numeric(q[1])
      comp[[name]]$upr <- as.numeric(q[2])
      comp[[name]]$AgWindow <- as.integer(agWindow[i])
      comp[[name]]$state <- as.integer(j)
      comp[[name]]$comp <- paste(tactic.comp[,k][1], tactic.comp[,k][2], sep='-')
      
    }
  }
}

ag.comp <- do.call(rbind, comp)


## Vizualise results ##

# prep ag comps
ag.comp$compType <- "ag"
ag.comp <- mutate(ag.comp, AgWindow = if_else(AgWindow == 1, "Ag Window", "Non-Ag Window"))

# prep baseline comps
state.comps$compType <- "base"
state.comps$AgWindow <- "Baseline"

# bind all comps into long form
comps <- rbind(state.comps, ag.comp) 
comps$AgWindow <- factor(comps$AgWindow,levels(factor(comps$AgWindow))[c(1,3,2)])
comps$state.names <- as.factor(comps$state)
levels(comps$state.names) <- c("Encamped", "R-A Search", "Exploratory")

library(DescTools)
# get tactic-level mean overlaps
comp.sum <- comps %>%
  #mutate(AgWindow.2 = ifelse(AgWindow == "Baseline", "Baseline", "Ag.Weighted")) %>%
  #group_by(AgWindow.2, state.names) %>%
  group_by(AgWindow, state.names) %>%
  summarise(n = n(),
            grand.mean.OV = mean(estOV),
            lwr.ci = MeanCI(estOV, method = "boot", type = "norm", R=1000)[2],
            upr.ci = MeanCI(estOV, method = "boot", type = "norm", R=1000)[3])

comp.sum


# plot comps
ggplot(data = comp.sum, aes(x = state.names, y = grand.mean.OV)) +
  geom_pointrange(aes(ymin = lwr.ci, ymax = upr.ci, color = as.factor(AgWindow)), position = position_dodge(0.4))

comps.2 <- comps %>%
  mutate(AgWindow.2 = ifelse(AgWindow == "Baseline", "Baseline", "Window"))
ggplot(data = comps.2, aes(x = state.names, y = estOV)) + geom_boxplot(aes(color = as.factor(AgWindow.2)))

# plot
ggplot(data = comps, aes(x = comp, y = estOV)) + 
  geom_pointrange(aes(ymin = lwr, ymax = upr, color = as.factor(AgWindow)), position = position_dodge(0.4)) + 
  facet_wrap(~state) + xlab("Tactic Comparison") + ylab("Overlap Estimate") +
  labs(color = "Comparison Type") +
  scale_color_manual(values = c("#F8766D", "#00BFC4", "grey50")) + theme_bw()


ggplot(data = filter(comps, AgWindow != "Non-Ag Window"), aes(x = comp, y = estOV)) + 
  geom_pointrange(aes(ymin = lwr, ymax = upr, color = as.factor(AgWindow)), position = position_dodge(0.4)) + 
  facet_wrap(~state.names) + xlab("Tactic Comparison") + ylab("Overlap Estimate") +
  labs(color = "Comparison Type") +
  scale_color_manual(values = c("#F8766D", "grey50")) + theme_bw()



##### Stationary States #####
t <- stationary(hmm.df, covs = data.frame(dist2agedge = 0, dist2forest = hmm.df$rawCovs$dist2forest, dist2permwater = 0, dist2seasonalwater = 0, gHM = mean(hmm.df$rawCovs$gHM), slope = 0))
plot(output$dist2forest, t[[1]][,1], col = "orange", ylim = c(0,1))
points(output$dist2forest, t[[1]][,2], col = "blue")
points(output$dist2forest, t[[1]][,3], col = "green")


plot(output$dist2forest, output$Encamped, pch = 19)

plot(output$dist2agedge, output$Encamped, pch = 19)


ggplot(output, aes(dist2agedge, Encamped)) + geom_point() + geom_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1)
ggplot(output, aes(dist2agedge, `Area-Restricted Search`)) + geom_point() + geom_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1)
ggplot(output, aes(dist2agedge, Exploratory)) + geom_point() + geom_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1)

##### Activity Budget %%s#####
t <- output.plot %>%
  group_by(viterbi) %>%
  filter(dist2forest == 0) %>%
  summarise(n = n())
t$prop <- t$n/sum(t$n)
t


t <- output.plot %>%
  group_by(viterbi) %>%
  filter(dist2ag == 0) %>%
  summarise(n = n())
t$prop <- t$n/sum(t$n)
t

## Ag Window Stats
t <- output.plot %>%
  group_by(ag.window) %>%
  summarise(n = n(),
            median = median(dist2agedge),
            iqr.upper = median + IQR(dist2agedge),
            iqr.lower = median - IQR(dist2agedge),
            min = min(dist2agedge),
            max = max(dist2agedge))

t

median(output$dist2agedge)
IQR(output$dist2agedge)
