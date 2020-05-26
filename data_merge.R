##Merge datasets for GME 
# Ingest cleaned datasets
# MEP
# Grumeti

library(tidyverse)


##### LOAD DATA #####
mep <- readRDS("./movdata/MEPcollars_20200516_clean_2020-05-21.rds")
metaMEP <- read.csv("./movdata/MEPcollars_20200516_metadata.csv")
metaMEP$X <- NULL
gr <- readRDS("./movdata/GRcollars_20191231_clean_2020-05-21.rds")
metaGR <- read.csv("./movdata/GRCollars_211119_metadata.csv")


##### MERGE #####
# prep datasets
mep <- mep %>%
  mutate(ele.xy.TEMPERATURE = NA) %>%
  inner_join(., metaMEP, by = "id") %>%
  mutate(site = "mep") %>%
  rename(date = Fixtime) 
mep$n <- NULL
mep.col.names <- colnames(mep)


gr <- gr %>%
  left_join(., metaGR, by = "id") %>%
  rename(subject_name = Name, subject_sex = Sex) %>%
  mutate(collar_id = id, 
    id = paste(subject_name, substr(collar_id, nchar(collar_id)-3, nchar(collar_id)), sep = "-"), 
    region = "Serengeti", fixType = "regular", site = "gr") %>%
  mutate(subject_sex  = ifelse(subject_sex == "F", "female", "male")) %>%
  dplyr::select(mep.col.names)


# filter out 
# mep <- filter(mep, region )

# merge
output <- bind_rows(mep, gr) %>%
  filter(year(date) < 2020)


##### Tracking Summary #####
# summary table function
GME_summary <- function(data) {
  require(dplyr)
  require(lubridate)
  
  # 1. get summary stats by id
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
outfile <- paste0("./movdata/GMEcollars_002_clean_", Sys.Date(), ".rds" )
#write.csv(output, outfile)
saveRDS(output, outfile)

# outfile <- paste0("./movdata/GMEsummary_001.csv" )
# write.csv(tracking.summary, outfile)





