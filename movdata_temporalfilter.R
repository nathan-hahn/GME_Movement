#### Used Filter #####

library(dplyr)
library(lubridate)

movdat <- readRDS('./movdata/GMEcollars_003_used_2020-10-30.rds')

gr <- movdat %>%
  filter(site == "gr") %>%
  filter(!(subject_name %in% c("Nyanza"))) %>%
  filter(minute(date) <= 10 | minute(date) >= 50) %>%
  droplevels()

mep <- movdat %>%
  filter(site == "mep") %>%
  filter(!(subject_name %in% c("Heritage", "Tunai","Bettye","Nancy", "Wilbur", 
                               "Julia", "Ndorre", "Nkoidila", "Santiyan", "Earhart",
                               "Rudisha", "Naeku", "Olkeri", "ST2010-1441", "Rudisha"))) 

# bind and limit up to 2019 data
downsampled <- rbind(gr, mep) %>% 
  filter(year(date) <= 2019) %>% droplevels()
max(downsampled$date)


saveRDS(downsampled, paste0("./movdata/GMEcollars_003_usedFilter_", Sys.Date(),".rds"))

# general mov stats using filtered data
downsampled %>% mutate(ag.used = ifelse(lc.estes == 1, 1, 0)) %>%
  summarise(mean(ag.used))

t <- downsampled %>% group_by(subject_name, subject_sex, subject_ageClass) %>% tally()
table(t$subject_sex)

