#### Used Filter #####

library(dplyr)
library(lubridate)

movdat <- readRDS('./movdata/GMEcollars_002_used_2020-06-30.rds')
movdat$subject_name <- ifelse(is.na(movdat$subject_name), movdat$collar_id, movdat$subject_name)

gr <- movdat %>%
  filter(site == "gr") %>%
  filter(!(subject_name %in% c("Nyanza"))) %>%
  filter(minute(date) <= 10 | minute(date) >= 50) %>%
  droplevels()

mep <- movdat %>%
  filter(site == "mep") %>%
  filter(!(subject_name %in% c("Tunai","Bettye","Nancy", "Wilbur", 
                               "Julia", "Ndorre", "Nkoidila", "Santiyan", "Earhart",
                               "Rudisha", "Naeku", "Olkeri", "ST2010-1441"))) 

downsampled <- rbind(gr, mep) 

saveRDS(downsampled, "./movdata/GMEcollars_002_usedFilter_2020-06-30.rds")
