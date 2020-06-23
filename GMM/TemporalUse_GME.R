##### Ag Temporal Use #####
library(tidyverse)
source("GME_functions.R")

##### Data Prep #####
movdat <- readRDS('./movdata/GMEcollars_002_usedFilter_2020-06-20.rds')
movdat$ag.used <- ifelse(movdat$lc.estes == 1, 1, 0) 

##### Rolling Stats #####
split <- split(movdat, movdat$subject_name)

window <- c(1*24, 7*24, 30*24, 90*24) # in hours. Split to before and after when align = center

output.plot <- NULL
for(i in 1:length(window)){
  # calculate moving window stats - see rollstats function
  output.plot[[i]] <- window_stats(df.list = split, window = window[[i]])
  # remove NAs from ag windows - created by window cutoff
  output.plot[[i]] <- filter(output.plot[[i]], is.na(ag.window) == FALSE) %>% droplevels()
}

##### Visualize #####
res.90 <- output.plot[[3]] 
res.90$ag.used <- ifelse(res.90$ag.used == 1, 0.75, NA) # adjust ag.used values for plotting for visibility
p90 <- ggplot(res.90, aes(x = date, color = subject_name)) + #, color = name
  geom_point(aes(y = ag.used), color = "grey40", size = .2, alpha = 0.5) +
  geom_ribbon(aes(ymin = lo.95, ymax = hi.95), size = .5, alpha = 0.4) +
  geom_line(aes(y = mean), size = 1, alpha = 0.7) +
  facet_wrap(~ subject_name) +
  labs(title = "Ag Usage: Volatility and Trend",
       subtitle = "60-Day Moving Average with 95% CI Bands") +
  theme(legend.position="none")

p90


