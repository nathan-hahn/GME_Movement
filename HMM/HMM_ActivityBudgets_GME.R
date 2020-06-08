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

##### Prep Data #####
hmm.df <- readRDS("./HMM/model results/TEST_m2_indiv_df.rds") # output from hmm model (no squared term)
output <- readRDS("./HMM/model results/TEST_traindf_original.rds") # dataset prior to data transforms

# add viterbi estimates to the original dataframe -- revert to non-standardized covs and steps
output$viterbi <- hmm.df$viterbi
# create new varibles to differentiate ag and non-ag in non-protected areas
output$ag.used <- ifelse(output$lc.estes == 1, 1, 0)
output$pa <- factor(output$pa, levels = c(3:1))
output$pa.2 <- ifelse(output$pa == 1 & output$ag.used == 1, 4, output$pa)

###### Identify Ag Use Windows #####

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

write.csv(output, "./HMM/model results/agWindow_temporary_df.csv")

##### Activity Budgets #####

output.plot <- mutate(output, viterbi = recode_factor(viterbi, 
                                                 "1" = "encamped",
                                                 "2" = "foraging",
                                                 "3" = "transit"),
                 pa.2 = recode_factor(pa.2, "1" = "protected", "2" = "limited use", 
                                      "3" = "not protected", "4" = "ag")
) 


t <- ifelse(is.na(output.plot$ag.window) == TRUE, 0, output.plot$ag.window)
output.plot$ag.window <- as.factor(t)

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
landUse <- plot_budget(t = output.plot, facet = viterbi ~ pa.2, title = "GME: Activity budget by land use type")
landUse






