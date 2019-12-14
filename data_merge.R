##Merge datasets for GME 
# Ingest cleaned datasets
# MEP
# Grumeti

library(tidyverse)

mep <- readRDS("MEPcollars_001_clean_2019-06-02.rds")
gr <- readRDS("GRcollars_001_clean_2019-11-21.rds")

metaGR <- read.csv("EleCollars_211119_metadata.csv")

# prep datasets
gr <- gr %>%
  inner_join(., metaGR, by = "id") %>%
  mutate(collar_type = "Henrik GL200", fixType = "regular", site = "gr") %>%
  rename(name = Name, sex = Sex)
gr$ele.xy.TEMPERATURE <- NULL


mep <- mep %>%
  mutate(site = "mep") %>%
  rename(date = Fixtime) 
mep$collar_id <- NULL

# merge
output <- rbind(mep, gr)

# save
outfile <- paste0("GMEcollars_001_clean_",Sys.Date(),".csv" )
saveRDS(output, outfile)
write.csv(output, outfile)



