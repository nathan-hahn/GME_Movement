##Merge datasets for GME 
# Ingest cleaned datasets
# MEP
# Grumeti


library(tidyverse)
library(lubridate)


##### LOAD DATA #####
mep <- readRDS("./movdata/mep/MEPcollars_20211001_clean_2021-10-02.rds")
metaMEP <- read.csv("./movdata/mep/MEPcollars_20211001_metadata.csv") 
metaMEP$X <- NULL
gr <- as.data.frame(data.table::fread("./movdata/grumeti/GRcollars_20210625_clean_2021-09-16.csv"))
metaGR <- read.csv("./movdata/grumeti/GRcollars_20191231_metadata.csv")


##### MERGE #####
# prep datasets
mep <- mep %>%
  mutate(ele.xy.TEMPERATURE = NA) %>%
  inner_join(., metaMEP, by = "id") %>%
  mutate(site = "mep",
         Fixtime = ymd_hms(Fixtime)) %>%
  rename(date = Fixtime)
mep$n <- NULL
mep$region_raw <- NULL
mep.col.names <- colnames(mep)


gr <- gr %>%
  mutate(id = tolower(id)) %>%
  left_join(., metaGR, by = "id") %>%
  rename(subject_name = Name, subject_sex = Sex, subject_age = Age) %>%
  mutate(source_id = as.character(id), 
    id = paste(subject_name, substr(source_id, nchar(source_id)-3, nchar(source_id)), sep = "-"), 
    region = "serengeti", fixType = "regular", site = "gr",
    date = ymd_hms(date)) %>%
  dplyr::select(mep.col.names)

# merge
output <- bind_rows(mep, gr) %>%
  #filter(year(date) < 2020) %>%
  filter(region %in% c("serengeti", "masai mara")) %>%
  mutate(subject_sex = if_else(subject_sex %in% c("male", "M"), "male", "female"))

# create a uid value for GME dataset
output$uid <- c(1:(dim(output)[1]))


##### Tracking Summary #####
# summary table function
GME_summary <- function(data) {
  require(dplyr)
  require(lubridate)
  
  # 1. get summary stats by id
  sum <- data %>%
    group_by(subject_name, id, subject_sex, subject_ageClass, subject_dob, site, fixType) %>%
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
  
  # 3. Return summary table
  return(table)
  
}

# summary table
tracking.summary <- GME_summary(data = output)
View(tracking.summary)

##### SAVE #####
outfile <- paste0("./movdata/GMEcollars_004_clean_", Sys.Date(), ".rds" )
#write.csv(output, outfile)
saveRDS(output, outfile)

outfile <- paste0("./movdata/GMEsummary_004.csv" )
write.csv(tracking.summary, outfile)





