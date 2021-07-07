library(tidyverse)
library(ggplot2)
theme_set

##### Data Prep #####
movdat <- readRDS('./GMM/movdata/GMEcollars_003_usedClust_2020-10-30.rds')
movdat$ag.used <- ifelse(movdat$lc.estes == 1, 1, 0) 


## All eles
res.90 <- movdat
res.90$ag.used <- ifelse(res.90$ag.used == 1, 0.5, NA) # adjust ag.used values for plotting for visibility



p90 <- ggplot(res.90, aes(x = date, color = factor(subject_name))) + 
  #geom_point(aes(y = ag.used), color = "grey40", size = .2, alpha = 0.5) +
  #geom_ribbon(aes(ymin = lo.95, ymax = hi.95), size = .5, alpha = 0.4) +
  geom_line(aes(y = mean), size = .5) +
  facet_wrap(region ~ tactic.agg) +
  labs(title = "Ag Use Trend",
       subtitle = "90-Day Moving Average") +
  xlab("Date") + ylab("90-Day Mean Occupancy") +
  theme(legend.position="none")
p90 + ylim(0, 0.5)



## Sample
res.90 <- movdat %>%
  filter(subject_name %in% c('Bonchugu', 'Maddy', 'Hangzhou', 'Pepper', 'Mytene',
                             'Jerahapembe', 'Matobo', 'Marima',
                             'Mchuri', 'Chuma', 'Imara', 'Ivy', 'Fred', 'Naibosho',
                             'Kimbizwa', 'Lowana', 'Lucy', 'Olchoda', 'Bobo')) %>%
  droplevels()
res.90$ag.used <- ifelse(res.90$ag.used == 1, 0.5, NA) # adjust ag.used values for plotting for visibility
res.90$tactic.agg <- as.factor(res.90$tactic.agg)
levels(res.90$tactic.agg) <- c('Rare', 'Sporadic', 'Seasonal', 'Habitual')
res.90$tactic.season <- as.factor(res.90$tactic.season)
levels(res.90$tactic.season) <- c('Rare', 'Sporadic', 'Seasonal', 'Habitual')

p90 <- ggplot(res.90, aes(x = date, color = region, shape = subject_name)) + 
  geom_line(aes(y = mean), size = 0.5) +
  facet_wrap(~tactic.season) +
  xlab("Date") + ylab("90-Day Mean Occupancy") +
  theme(legend.position="none")
p90 + ylim(0, 0.5) + theme_bw() + theme(legend.position='none')

ggsave('Fig S1 temporal use plot.tiff', dpi = 300)

