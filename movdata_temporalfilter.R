#### Used Filter #####

library(dplyr)
library(lubridate)

movdat <- readRDS('./movdata/GMEcollars_002_used_2020-07-14.rds')

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

downsampled <- rbind(gr, mep) 

saveRDS(downsampled, paste0("./movdata/GMEcollars_002_usedFilter_", Sys.Date(),".rds"))
