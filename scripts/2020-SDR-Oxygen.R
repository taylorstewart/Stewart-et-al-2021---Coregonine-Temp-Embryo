## ===========================================================
## Clear the environment first
## ===========================================================
rm(list = ls(all.names = TRUE))


## ===========================================================
## Load packages
## ===========================================================
library(tidyverse)
library(readxl)
library(lubridate)
library(magrittr)


## ===========================================================
## Load data
## ===========================================================
sdr.files <- list.files('data/SDR_Oxygen', full.names = TRUE, pattern = "xlsx")

data.oxygen <- do.call(rbind, lapply(sdr.files, function(i) {
  print(i)
  num <- which(sdr.files == i)
    
  tmp <- read_excel(i, sheet = "Oxygen Orginal", skip = 12)[c(1, 3:26)] %>% 
    mutate(num = num)
  
  n <- ceiling(nrow(tmp) * 0.10)
  
  tmp %<>% 
    ## remove first 60 minutes when temperature may be unstable
    slice(-1:-n) %>% 
    ## remove last 60 minutes when temperature may be unstable
    slice(1:(n()-n)) %>%
    ## format datetime
    rename(datetime = "Date/Time") %>% 
    mutate(datetime = gsub("[.]", "-", datetime),
           datetime = as.POSIXct(datetime, format = "%d-%m-%y %H:%M:%S", tz = "America/New_York"),
           datetime = round_date(datetime, unit = "minute"),
           ## add population from file name
           population = substr(i, 24, 25),
           ## change numerical population to long name
           population = gsub("LS", "Superior", population),
           population = gsub("LO", "Ontario", population),
           ## add temperature from file name
           temperature = substr(i, 30, 30),
           ## change numerical treatment to temperature
           temperature = as.numeric(gsub("4", "9.0", temperature)),
           temperature = as.numeric(gsub("3", "7.0", temperature)),
           temperature = as.numeric(gsub("2", "4.5", temperature)),
           temperature = as.numeric(gsub("1", "2.0", temperature)))
}))

## load survival data
data.survival <- read_excel("data/2020-Artedi-SDR-Oxygen.xlsx", sheet = "OxoDishes") %>% 
  mutate(survival = factor(survival))


## ===========================================================
## Data Manipulation - 30 Minutes
## ===========================================================
## Create new datetime to summarize by every 30 minutes
data.oxygen.30 <- data.oxygen %>% mutate(date = as.Date(datetime, format = "%Y/%m/%d", tz = "America/New_York"),
                                         hour = hour(datetime),
                                         minute = minute(datetime),
                                         minute.30 = ifelse(minute <= 30, 0, 30),
                                         datetime = as.POSIXct(paste0(date, " ", hour, ":", minute.30), format = "%Y-%m-%d %H:%M")) %>% 
  select(datetime, population, temperature, num, 2:25)

## gather well columns into tidy format
data.oxygen.30 %<>% pivot_longer(5:28, "well")

## Summarize by every 30 minutes
data.oxygen.30 %<>% group_by(datetime, population, temperature, num, well) %>% 
  summarize(value = mean(value))


## ===========================================================
## Data Manipulation - 1 Minute
## ===========================================================
## gather well columns into tidy format
data.oxygen.1 <- data.oxygen %>% pivot_longer(2:25, "well")


## ===========================================================
## 
## ===========================================================
data <- left_join(data.oxygen.30, data.survival) %>% 
  mutate(group = paste0(num, "-", well)) %>% 
  filter(population == "Ontario", temperature == 7)


## ===========================================================
## Visualization
## ===========================================================
ggplot(data, aes(x = datetime, y = value)) +
  geom_line(aes(group = num, color = survival)) +
  scale_color_manual(values = c("#ca0020", "#0571b0")) +
  labs(x = "Date", y = "Oxygen [cO2 [mg/L]]", title = paste0(data$temperature, " ", data$population)) +
  theme_bw() +
  facet_wrap(~well, scales = "free_x", ncol = 6)




ggplot(data, aes(x = datetime, y = value)) +
  geom_line(aes(group = group, color = well)) +
  theme_bw() +
  facet_wrap(population~temperature, scales = "free_x", ncol = 4)
