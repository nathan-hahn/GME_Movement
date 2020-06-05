##### Activity Moving Windows #####

#' Window Size = 12 hours in ag
#' Window rules:
#' 1. Any point in ag creates a +/- 12 hour window around the point
#' 2. Overlapping windows condensed into a single window

library(tidyverse)
library(xts)
library(zoo)


hmm.df <- readRDS("./HMM/TEST_m1_indiv_df.rds")
output <- readRDS("./HMM/TEST_traindf_original.rds")

output$viterbi <- hmm.df$viterbi
output$ag.used <- ifelse(output$lc.estes == 1, 1, 0)
output$pa <- factor(output$pa, levels = c(3:1))
output$pa.2 <- ifelse(output$pa == 1 & output$ag.used == 1, 4, output$pa)


###### Identify Ag Use Windows #####

# function to return stats from roll apply
rollstat <- function(x, na.rm = TRUE) {
  # x     = numeric vector
  # na.rm = boolean, whether or not to remove NA's
  
  m  <- mean(x, na.rm = na.rm)
  s  <- sd(x, na.rm = na.rm)
  hi <- m + 2*s
  lo <- m - 2*s
  max <- max(x) # all points within the window will be labeled 1 - signals a "raiding period" for activity budgets/plotting

  ret <- c(mean = m, stdev = s, hi.95 = hi, lo.95 = lo, ag.window = max) 
  return(ret)
}

# Apply rollstat. Create time series and apply rollstat function
window <- 25 # in hours. Split to before and after when align = center
split <- split(output, output$burst) # careful of what you split by. May cuase problems splitting by name if there are large gaps between recollars. Let na.rm = FALSE in rollstat

system.time({
  ts <- lapply(split, as.ts, start = split[[i]]$date, frequency = 1) # frequency set with 30min fixes
  tx <- lapply(ts, as.xts)
  res <- lapply(seq_along(tx), function(i) rollapply(tx[[i]]$ag.used, width = window, align = "center", by.column = FALSE,
                                                     FUN = rollstat, na.rm = FALSE))
  merge <- lapply(seq_along(res), function(i) cbind(split[[i]], as.data.frame(res[[i]])))
})

output <- do.call("rbind", merge)

##### Activity Budgets #####

t <- output %>%
  group_by(burst, pa, lc.estes) %>%
  filter(lc.estes == 1) %>%
  tally()

# output$ag <- ifelse(output$lc.estes == 1, 1, 0) # calculate ag.used (0,1)
# output$p <- ifelse(output$pa == "protected", 1, 0) # calculate ag.used (0,1)
# output$np <- ifelse(output$pa == "not protected", 1, 0) # calculate ag.used (0,1)
# output$lu <- ifelse(output$pa == "limited use", 1, 0) # calculate ag.used (0,1)

output <- mutate(output, viterbi = recode_factor(viterbi, 
                                                 "1" = "encamped",
                                                 "2" = "foraging",
                                                 "3" = "transit"),
                 pa.2 = recode_factor(pa.2, "1" = "protected", "2" = "limited use", 
                                      "3" = "not protected", "4" = "ag")
) 


# fix NAs for the ag window
t <- ifelse(is.na(output$ag.window) == TRUE, 0, output$ag.window)
output$ag.window <- as.factor(t)


## Density plots v2
theme_set(  theme_bw() + # set theme with no legend of strip text
              theme(panel.grid.major = element_blank(),
                    strip.background = element_blank(),
                    panel.border = element_rect(colour = "black"),
                    strip.text = element_text(size = 12),
                    legend.text = element_text(size = 10),
                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title=element_text(size=14), axis.text = element_text(size=12))
)


plot_budget <- function(t=t, facet, title){
  ggplot(t, aes(hour(date))) +
    facet_grid(facet) +
    geom_density(aes(x = hour(date),
                     fill = factor(viterbi), colour = factor(viterbi)), alpha = 0.3, adjust = 1.5)
  
  # add day/night lines - code for shading below
  geom_vline(aes(xintercept=6),
             color="blue", linetype="dashed", size=1) +
    geom_vline(aes(xintercept=18),
               color="blue", linetype="dashed", size=1) +
    
    # add colors
    scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73"), name = c("state")) +
    scale_colour_manual(values = c("#E69F00", "#56B4E9", "#009E73"), name = c("state")) +
    xlab("hour (0-23)") + ggtitle(title)
}

# budget by ag window
plot_budget(t = output, facet = viterbi ~ as.factor(max), title = "ag.class.mean")

# budget by land use
plot_budget(t = output, facet = viterbi ~ pa.2, title = "GME: Activity budget by land use type")







