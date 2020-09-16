## ===========================================================
## Clear the environment first
## ===========================================================
rm(list = ls(all.names=TRUE))


## ===========================================================
## Load packages
## ===========================================================
library(dplyr)
library(magrittr)
library(readxl)
library(ggplot2)


## ===========================================================
## Load data
## ===========================================================
climate.temp.lo.1 <- read_excel("data/HOBO/2019-Artedi-Temperature-Embryo.xlsx", sheet = "LO-Trt2.0") %>% 
  mutate(date = as.POSIXct(paste0(substr(datetime, 1, 4), "-", substr(datetime, 6, 7), "-", substr(datetime, 9, 10)))) %>% 
  select(population, temperature, date, temp_c_mean)

climate.temp.lo.2 <- read_excel("data/HOBO/2020-Artedi-Temperature-Embryo.xlsx", sheet = "LO-Trt4.5") %>% 
  mutate(date = as.POSIXct(paste0(substr(datetime, 1, 4), "-", substr(datetime, 6, 7), "-", substr(datetime, 9, 10)))) %>% 
  select(population, temperature, date, temp_c_mean)

## Combine all treatments into one dataframe
climate.temp.all <- bind_rows(climate.temp.lo.1, climate.temp.lo.2)


## ===========================================================
## Calculate accumulative degree days
## ===========================================================
climate.temp.summary <- climate.temp.all %>% group_by(population, temperature, date) %>% 
  summarize(mean.daily.temp = mean(temp_c_mean)) %>% ungroup() %>% 
  group_by(population, temperature) %>% 
  mutate(ADD = cumsum(mean.daily.temp)) %>% ungroup()

write.csv(climate.temp.summary, "data/2019-Artedi-ADD.csv", row.names = FALSE)


