#### Used Filter #####

library(dplyr)
library(lubridate)

movdat <- readRDS('./movdata/GMEcollars_004_used_2021-10-05.rds')

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

# bind and limit up to filterdate. 

filterdate <- as.POSIXct('2021-04-01 00:00:00') # for gme_004, filter date is marked as latest ag cut date (april 1 2021)

downsampled <- rbind(gr, mep) %>% 
  filter(date < filterdate) %>% droplevels()
max(downsampled$date)


saveRDS(downsampled, paste0("./movdata/GMEcollars_004_usedFilter_", Sys.Date(),".rds"))



