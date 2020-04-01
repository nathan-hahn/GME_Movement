##Merge datasets for GME 
# Ingest cleaned datasets
# MEP
# Grumeti

library(tidyverse)


##### LOAD DATA #####
mep <- readRDS("./movdata/MEPcollars_001_clean_2019-06-02.rds")
gr <- readRDS("./movdata/GRcollars_001_clean_2019-11-21.rds")
metaGR <- read.csv("./movdata/EleCollars_211119_metadata.csv")


##### MERGE #####
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

# filter out 
# mep <- mep[-which(mep$name%in%c("Tunai","Bettye","Nancy", "Wilbur", 
#                                 "Julia", "Ndorre", "Nkoidila", "Santiyan", "Earhart",
#                                 "Rudisha", "Naeku", "Olkeri", "ST2010-1441")),]

# merge
output <- rbind(mep, gr)


##### Tracking Summary #####
# summary table function
GME_summary <- function(data) {
  require(dplyr)
  require(lubridate)
  
  # 1. get summary stats by chronofile
  sum <- data %>%
    group_by(id, name, sex, site, fixType) %>%
    summarise(dataStart = min(date),
              dataStop = max(date),
              daysTracked = dataStop - dataStart,
              n_fixes = n(),
              fix_rate = round((1 - sum(is.na(x)/n_fixes))*100,2)) 

  # 2. format summary table
  table <- as.data.frame(sum) %>%
    mutate(dataStart = as_date(dataStart),
           dataStop = as_date(dataStop),
           daysTracked = round(as.numeric(daysTracked)))
  
  levels(table$sex)[levels(table$sex)=="M"] <- "Male"
  levels(table$sex)[levels(table$sex)=="F"] <- "Female"
  
  # 3. Return summary table
  return(table)
  
}

# summary table
tracking.summary <- GME_summary(data = output)


##### SAVE #####
outfile <- paste0("./movdata/GMEcollars_001_clean.csv" )
saveRDS(output, outfile)
write.csv(output, outfile)

outfile <- paste0("./movdata/GMEsummary_001.csv" )
write.csv(tracking.summary, outfile)





