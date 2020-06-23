#### Used Filter #####

movdat <- readRDS('./movdata/GMEcollars_002_used_2020-06-20.rds')
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

fixes <- bind %>%
  group_by(subject_name) %>%
  tally() %>%
  filter(n >= 120*24)

movdat.filter <- filter(downsampled, subject_name %in% fixes$subject_name)

saveRDS(movdat.filter, "./movdata/GMEcollars_002_usedFilter_2020-06-20.rds")
