## CLEAR THE ENVIRONMENT FIRST ==================================

rm(list = ls(all.names = TRUE))


## START TIME ===================================================

start <- Sys.time()


## LOAD PACKAGES ================================================

library(dbplyr)         # for database data manipulation
library(dplyr)          # standard database manipulation
library(data.table)     # for fread() - super fast load-in of data
library(RSQLite)        # SQL database connection


## MAKE A LIST OF FILES ===============

sst.files <- list.files('/Volumes/home/GL-Seasonal-Environmental-Change/data/lake-ontario-surface-water-temperature/raw',
                        full.names = TRUE, recursive = TRUE, pattern = ".csv")[7579:14328]


sst.lo <- do.call(rbind, lapply(sst.files, function(i) {
  ## Print current place in loop
  print(i)
  
  ## Load data and join season dates
  sst.data <- fread(file = i) %>% 
    filter(lat == 44.05, lon == -76.17) %>% 
    mutate(sst = ifelse(sst < 0, 0, sst))
}))

write.csv(sst.lo, "LakeOntario_ChaumontBay.csv", row.names = FALSE)
