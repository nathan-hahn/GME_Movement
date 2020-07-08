#### Used Filter #####

library(dplyr)
library(lubridate)

movdat <- readRDS('./movdata/GMEcollars_002_used_2020-06-30.rds')
movdat$subject_name <- ifelse(is.na(movdat$subject_name), movdat$collar_id, movdat$subject_name)

gr <- movdat %>%
  filter(site == "gr") %>%
  # mutate(collar_id = tolower(collar_id)) %>%
  # mutate(subject_name = collar_id) %>% ## TEMPORARY
  filter(!(subject_name %in% c("Nyanza"))) %>% #"st2010-3047"
  filter(minute(date) <= 10 | minute(date) >= 50) %>%
  droplevels()

mep <- movdat %>%
  filter(site == "mep") %>%
  filter(!(subject_name %in% c("Tunai","Bettye","Nancy", "Wilbur", 
                               "Julia", "Ndorre", "Nkoidila", "Santiyan", "Earhart",
                               "Rudisha", "Naeku", "Olkeri", "ST2010-1441"))) 


downsampled <- rbind(gr, mep) 

fixes <- downsampled %>%
  group_by(subject_name) %>%
  tally() %>%
  filter(n >= 120*24)

movdat.filter <- filter(downsampled, subject_name %in% fixes$subject_name)

saveRDS(movdat.filter, "./movdata/GMEcollars_002_usedFilter_2020-06-30.rds")


### TEMPORARY

# gr <- movdat %>%
#   filter(site == "gr") %>%
#   mutate(collar_id = tolower(collar_id)) %>%
#   mutate(subject_name = collar_id) %>% ## TEMPORARY
#   filter(!(subject_name %in% c("st2010-3047"))) %>% #"st2010-3042"
#   filter(minute(date) <= 10 | minute(date) >= 50) %>%
#   droplevels()
# 
# mep <- movdat %>%
#   filter(site == "mep") %>%
#   filter(!(subject_name %in% c("Tunai","Bettye","Nancy", "Wilbur", 
#                                "Julia", "Ndorre", "Nkoidila", "Santiyan", "Earhart",
#                                "Rudisha", "Naeku", "Olkeri", "ST2010-1441"))) 
# 
# 
# downsampled <- rbind(gr, mep) 
# 
# fixes <- downsampled %>%
#   group_by(subject_name) %>%
#   tally() %>%
#   filter(n >= 120*24)
# 
# temp.filter <- filter(downsampled, subject_name %in% fixes$subject_name)
